import multiprocessing
import re
# import numpy

try:
    from classMiner import Miner, ExpMiner, testIni, DummyLog, RCollection

except ModuleNotFoundError:
    from .classMiner import Miner, ExpMiner, testIni, DummyLog, RCollection

import pdb

class MultiprocMiner(Miner):

    def full_run(self, cust_params={}):
        self.rcollect.resetFinal()
        self.count = 0
        self.reinitOnNew = False
        self.workers = {}
        self.rqueue = multiprocessing.Queue()
        self.pstopqueue = multiprocessing.Queue()

        self.logger.printL(1, "Start mining", "status", self.getLogId())  # todo ID

        # progress initialized after listing pairs
        self.logger.clockTic(self.getLogId(), "pairs")
        self.initializeRedescriptions()

        self.logger.clockTic(self.getLogId(), "full run")
        self.initializeExpansions()

        self.keepWatchDispatch()
        self.logger.clockTac(self.getLogId(), "full run", "%s" % self.questionLive())
        if not self.questionLive():
            self.logger.printL(1, "Interrupted!", "status", self.getLogId())
        else:
            self.logger.printL(1, "Done.", "status", self.getLogId())
        self.logger.sendCompleted(self.getLogId())
        self.logger.logAllTracksEnd(self.rcollect, "F", self.getLogId())
        return self.rcollect

    def initializePairs(self, ids=None):
        self.pairs = 0
        # Loading pairs from file if filename provided
        loaded, done = self.initial_candidates.loadFromFile(self.data)
        if not loaded or done is not None:

            self.logger.printL(1, "Searching for initial pairs", "status", self.getLogId())
            explore_list = self.getInitExploreList(ids, done)
            self.logger.initProgressFull(self.constraints, explore_list, self.rcollect.getNbAvailableCols(), 1, self.getLogId())
            self.total_pairs = len(explore_list)

            min_bsize = 25
            if self.pe_balance == 0:
                # Finish all pairs before exhausting,
                # split in fixed sized batches, not sorted by cost
                K = self.max_processes
                # self.total_pairs / (self.max_processes-1)
                batch_size = max(self.total_pairs // (5*K), min_bsize)
                ## print("Batch size=", batch_size)
                pointer = 0
            else:
                # Sort from easiest to most expensive
                tot_cost = sum([c[-1] for c in explore_list])
                explore_list.sort(key=lambda x: x[-1], reverse=True)
                thres_cost = (tot_cost/self.max_processes) * (1-self.pe_balance/10.)

                cost = 0
                off = 0
                while off < len(explore_list)-1 and cost < thres_cost:
                    off += 1
                    cost += explore_list[-off][-1]
                if len(self.initial_candidates) > off:
                    off = 0
                ## print("OFF", off)
                K = self.max_processes-1
                batch_size = max((self.total_pairs-off) // K, min_bsize)

                # Launch last worker
                if K*batch_size < len(explore_list):
                    self.logger.printL(3, "Initial pairs process [%s] with %d batch" %
                                       (K, len(explore_list[K*batch_size:])), "status", self.getLogId())
                    # print("Init PairsProcess ", K, self.getId(), len(explore_list[K*batch_size:]))
                    self.workers[K] = PairsProcess(K, self.getId(), explore_list[K*batch_size:], self.data, self.charbon, self.constraints, self.rqueue)
                pointer = -1

            self.initial_candidates.setExploreList(explore_list, pointer, batch_size, done)
            # Launch other workers
            for k in range(K):
                ll = self.initial_candidates.getExploreNextBatch(pointer=k)
                if len(ll) > 0:
                    self.logger.printL(3, "Initial pairs process [%s] with %d batch" %
                                       (k, len(ll)), "status", self.getLogId())
                    # print("Init PairsProcess ", k, self.getId(), len(ll))
                    self.workers[k] = PairsProcess(k, self.getId(), ll, self.data, self.charbon, self.constraints, self.rqueue)
                else:
                    pointer = -1

            self.pairWorkers = len(self.workers)
            if pointer == 0:
                pointer = self.pairWorkers
            self.initial_candidates.setExplorePointer(pointer)
            self.initial_candidates.setTracking(False)
            self.reinitOnNew = True
        else:
            self.logger.initProgressFull(self.constraints, None, self.rcollect.getNbAvailableCols(), 1, self.getLogId())
            self.logger.printL(1, self.initial_candidates.msgLoaded(), "log", self.getLogId())

        return self.initial_candidates

    def initializeExpansions(self):
        self.reinitOnNew = False
        for k in set(range(self.max_processes)).difference(self.workers.keys()):
            initial_red = self.initial_candidates.getNextRed(self.data, testIni)

            if self.questionLive():
                if initial_red is None:
                    self.reinitOnNew = True
                else:
                    self.count += 1

                    self.logger.printL(1, "Expansion %d" % self.count, "log", self.getLogId())
                    self.logger.printL(3, "Expand process [%s]" % k, "status", self.getLogId())
                    self.logger.clockTic(self.getLogId(), "expansion_%d-%d" % (self.count, k))
                    ## print("Init ExpandProcess ", k, self.count)
                    self.workers[k] = ExpandProcess(k, self.getId(), self.count, self.data,
                                                    self.charbon, self.constraints,
                                                    self.rqueue, [initial_red],
                                                    rcollect=self.shareRCollect(), logger=self.shareLogger())

    def handlePairResult(self, m):
        pairs, idL, idR, pload = m["pairs"], m["idL"], m["idR"], m["pload"]
        self.pairs += 1

        self.logger.updateProgress({"rcount": 0, "pair": self.pairs, "pload": pload})
        if self.pairs % 100 == 0:
            self.logger.printL(3, "Searching pair %d/%d (%i <=> %i)" %
                               (self.pairs, self.total_pairs, idL, idR), "status", self.getLogId())
            self.logger.updateProgress(level=3, id=self.getLogId())
        elif self.pairs % 10 == 0:
            self.logger.printL(7, "Searching pair %d/%d (%i <=> %i)" %
                               (self.pairs, self.total_pairs, idL, idR), "status", self.getLogId())
            self.logger.updateProgress(level=7, id=self.getLogId())
        else:
            self.logger.printL(10, "Searching pair %d/%d (%i <=> %i)" %
                               (self.pairs, self.total_pairs, idL, idR), "status", self.getLogId())

        added = 0
        for cand in pairs:
            if self.initial_candidates.add(cand) is not None:
                added += 1
                self.logger.printL(6, str(cand), "log", self.getLogId())

        self.initial_candidates.addExploredPair((idL, idR))
        if added > 0 and self.reinitOnNew:
            self.initializeExpansions()

    def handleExpandResult(self, m):
        self.rcollect.addItems(m["out"].getItems())

        if m["out"].getLen("P") > 0:
            self.rcollect.selected("final", ids="F", new_ids=m["out"].getIidsList("P"), trg_lid="F", constraints=self.constraints)
        self.rcollect.importTracks(m["out"].getTracks(), m["id"])

        if not self.constraints.getCstr("amnesic"):
            self.rcollect.updateSeen(m["out"].getItems("P"))
        self.logger.clockTac(self.getLogId(), "expansion_%d-%d" % (m["count"], m["id"]), "%s" % self.questionLive())
        self.logger.logResults(self.rcollect, "F", self.getLogId())
        self.logger.updateProgress({"rcount": m["count"]}, 1, self.getLogId())

    def leftOverPairs(self):
        # print("LEFTOVER", self.workers)
        return self.initial_candidates.exhausted() and len(self.workers) > 0 and len([(wi, ww) for (wi, ww) in self.workers.items() if ww.isExpand()]) == 0

    def keepWatchDispatch(self):
        while len(self.workers) > 0 and self.questionLive():
            m = self.rqueue.get()
            if m["what"] == "done":
                del self.workers[m["id"]]

                if self.initial_candidates.getExplorePointer() >= 0:
                    ll = self.initial_candidates.getExploreNextBatch()
                    if len(ll) > 0:
                        ## print("Init Additional PairsProcess ", m["id"], self.getId(), self.explore_pairs["pointer"], len(ll))
                        self.workers[m["id"]] = PairsProcess(m["id"], self.getId(), ll, self.data, self.charbon, self.constraints, self.rqueue)
                        self.initial_candidates.incrementExplorePointer()
                    else:
                        self.initial_candidates.setExplorePointer(-1)

                if self.initial_candidates.getExplorePointer() == -1:
                    self.pairWorkers -= 1
                    if self.pe_balance > 0:
                        self.initializeExpansions()
                    # else:
                    #     print("Waiting for all pairs to complete")

                if self.pairWorkers == 0:
                    self.logger.updateProgress({"rcount": 0}, 1, self.getLogId())
                    self.logger.clockTac(self.getLogId(), "pairs")
                    self.logger.printL(1, self.initial_candidates.msgFound(), "log", self.getLogId())
                    self.logger.updateProgress(level=1, id=self.getLogId())
                    self.initial_candidates.setExploredDone()
                    self.initial_candidates.saveToFile()
                    if self.pe_balance == 0:
                        self.initializeExpansions()
                        ## print("All pairs complete, launching expansion")

            elif m["what"] == "pairs":
                self.handlePairResult(m)

            elif m["what"] == "expand":
                self.handleExpandResult(m)
                del self.workers[m["id"]]
                self.initializeExpansions()

            if self.leftOverPairs():
                self.logger.printL(1, self.initial_candidates.msgFound(), "log", self.getLogId())
                break


class PairsProcess(multiprocessing.Process):
    task = "pairs"

    def __init__(self, sid, ppid, explore_list, data, charbon, constraints, rqueue):
        multiprocessing.Process.__init__(self)
        self.daemon = True
        self.id = sid
        self.ppid = ppid
        self.explore_list = explore_list

        self.data = data
        self.charbon = charbon.copy()
        self.charbon.setStore(data.nbRows(), constraints.getSSetts(), constraints, "P")
        self.queue = rqueue

        self.start()

    def getId(self):
        return self.id

    def getParentId(self):
        return self.ppid

    def getLogId(self):
        return "%s-%s" % (self.getParentId(), self.getId())

    def isExpand(self):  # expand or init?
        return False

    def getTask(self):
        return self.task

    def run(self):
        # print("(X) PairsProcess Run", self.id, self.ppid)
        for pairs, (idL, idR, dtlsL, dtlsR, pload) in enumerate(self.explore_list):
            pairs = self.charbon.computePair(self.data.col(0, idL), self.data.col(1, idR), self.data.getColsC(), self.data)
            self.queue.put({"id": self.getId(), "what": self.getTask(), "pairs": list(pairs), "idL": idL, "idR": idR, "pload": pload})
        self.queue.put({"id": self.getId(), "what": "done"})


class ExpandProcess(multiprocessing.Process, ExpMiner):

    def __init__(self, sid, ppid, count, data, charbon, constraints, rqueue,
                 nextge, rcollect=None, logger=None, question_live=None, task="expand"):
        multiprocessing.Process.__init__(self)
        self.task = task
        self.daemon = True
        self.id = sid
        self.ppid = ppid
        self.count = count

        self.data = data
        self.charbon = charbon.copy()
        self.constraints = constraints
        self.queue = rqueue

        self.up_souvenirs = False
        self.question_live = question_live
        if logger is None:
            self.logger = DummyLog()
        else:
            self.logger = logger

        self.nextge = nextge
        if rcollect is None:
            self.rcollect = RCollection()
        else:
            self.rcollect = rcollect
        self.start()

    def getId(self):
        return self.id

    def getParentId(self):
        return self.ppid

    def getLogId(self):
        return "%s-%s" % (self.getParentId(), self.getId())

    def getTask(self):
        return self.task

    def isExpand(self):  # expand or init?
        return True

    def run(self):
        if self.getTask() == "improve":
            rcollect = self.improveRedescriptions(self.nextge, rcollect=self.rcollect)
        else:
            rcollect = self.expandRedescriptions(self.nextge, rcollect=self.rcollect)
        self.queue.put({"id": self.getId(), "what": self.getTask(), "out": rcollect, "count": self.count})
