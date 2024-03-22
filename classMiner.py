import re
import datetime
# import numpy

try:
    from toolLog import Log
    from classContent import BatchCollection
    from classRedescription import Redescription
    from classCandidates import initCands
    from classSouvenirs import Souvenirs
    from classConstraints import Constraints

    from classCharbonXFIM import CharbonXFIM
    from classCharbonGMiss import CharbonGMiss
    from classCharbonGStd import CharbonGStd
    from classCharbonTAlt import CharbonTCW, CharbonTSprit, CharbonTSplit
    from classCharbonTLayer import CharbonTLayer

except ModuleNotFoundError:
    from .toolLog import Log
    from .classContent import BatchCollection
    from .classRedescription import Redescription
    from .classCandidates import initCands
    from .classSouvenirs import Souvenirs
    from .classConstraints import Constraints

    from .classCharbonXFIM import CharbonXFIM
    from .classCharbonGMiss import CharbonGMiss
    from .classCharbonGStd import CharbonGStd
    from .classCharbonTAlt import CharbonTCW, CharbonTSprit, CharbonTSplit
    from .classCharbonTLayer import CharbonTLayer

import pdb

XAUST_CLASSES = {"fim": CharbonXFIM}
XAUST_DEF = CharbonXFIM
TREE_CLASSES = {"layeredtrees": CharbonTLayer,
                "cartwheel": CharbonTCW,
                "splittrees": CharbonTSplit,
                "sprit": CharbonTSprit}
TREE_DEF = CharbonTLayer

CHARBON_MISS_FORCE = False
# CHARBON_MISS_FORCE = True

# PAIR_LOADS = [[1,2,3],
#               [2,4,6],
#               [3,6,10]]


def testIni(pair):
    if pair is None:
        return False
    return True


class DummyLog(object):
    verbosity = 0

    def printL(self, level, message, type_message="*", source=None):
        pass

    def updateProgress(self, details, level=-1, id=None):
        pass

    def logResults(self, rcollect, lid, pid):
        pass
    
    def logAllTracksEnd(self, rcollect, lid, pid):
        pass

class RCollection(BatchCollection, Souvenirs):
    name_class = "R"

    def __init__(self, nAvailableMo=None, nAmnesic=False, tracks_switch="OFF"):        
        Souvenirs.__init__(self, nAvailableMo, nAmnesic)
        BatchCollection.__init__(self, tracks_switch=tracks_switch)
        # prepare buffer list
        self.newList(iid="F")
        self.newList(iid="P")
        self.newList(iid="W")

    def toShare(self):
        return RCollection(self.copyAvailableCols(), self.isAmnesic(), tracks_switch=self.getTracksSwitch())

    def __str__(self):
        svrs = ""
        if self.isActive():
            svrs = ", " + Souvenirs.__str__(self)
        return "RCollection: %d items in W=%d/P=%d/F=%d, %d tracks%s" % (self.nbItems(), self.getLen("W"), self.getLen("P"), self.getLen("F"), self.nbTracks(), svrs)

    def resetWorking(self, reds=[]):
        self.clearList("W")
        self.addItems(reds, "W")

    def resetPartial(self, reds=[]):
        self.clearList("P")
        self.resetWorking(reds)

    def resetFinal(self, reds=[]):
        self.clearList("F")
        self.resetPartial(reds)
        

class ExpMiner(object):
    def __init__(self, pid, count, data, charbon, constraints, logger=None, question_live=None, up_souvenirs=True):
        self.charbon = charbon
        self.data = data
        self.count = count
        self.pid = pid
        self.constraints = constraints
        self.up_souvenirs = up_souvenirs
        self.question_live = question_live
        if logger is None:
            self.logger = DummyLog()
        else:
            self.logger = logger

    def getId(self):
        return self.pid

    def getLogId(self):
        return "%s" % self.getId()

    def questionLive(self):
        if self.question_live is None:
            return True
        return self.question_live()

    def expandRedescriptions(self, nextge, rcollect=None):
        if len(nextge) > 0 and self.questionLive():
            if rcollect is None:
                rcollect = RCollection()
            rcollect.resetPartial(nextge)
            if self.charbon.isTreeBased():
                return self.expandRedescriptionsTree(nextge, rcollect)
            else:
                return self.expandRedescriptionsGreedy(nextge, rcollect)
        return rcollect

    def expandRedescriptionsTree(self, nextge, rcollect):
        for redi, red in enumerate(nextge):
            self.logger.printL(2, "Expansion %s.%s\t%s" % (self.count, redi, red), "log", self.getLogId())
            kids = self.charbon.computeExpandTree(-1, self.data, red)
            self.logger.printL(2, "Expansion %s.%s returned %s redescriptions" % (self.count, redi, len(kids)), "status", self.getLogId())
            kid_ids = []
            for kid in kids:
                kid.normalizeSides()
                self.logger.printL(2, kid, "log", self.getLogId())
                rcollect.addItem(kid, "W")
                kid_ids.append(kid.getUid())
            if len(kid_ids) > 0:
                track = {"do": "expand-tree", "trg": kid_ids, "src": [red.getUid()], "out": "W"}
                rcollect.addTrack(track)

            # self.logger.updateProgress({"rcount": self.count, "redi": redi}, 4, self.getLogId())
            self.logger.printL(4, "Candidate %s.%d.%d grown" % (self.count, len(red), redi), "status", self.getLogId())
        self.logger.printL(4, "Generation %s.%d expanded" % (self.count, len(red)), "status", self.getLogId())

        # extra selection step for tree-based to check not repeating itself
        sel = rcollect.selected("partial", ids="W", constraints=self.constraints)
        Xsel = rcollect.selected("tree_rectangles", ids="F", new_ids=sel, new_only=True, trg_lid="P", constraints=self.constraints)
        self.logger.logResults(rcollect, "P", self.getLogId())
        return rcollect

    def expandRedescriptionsGreedy(self, nextge, rcollect):
        first_round = True
        max_var = [self.constraints.getCstr("max_var", side=0), self.constraints.getCstr("max_var", side=1)]
        for r in nextge:
            r.initAvailable(rcollect, self.data, max_var)
        self.charbon.setStore(self.data.nbRows(), self.constraints.getSSetts(), self.constraints)

        while len(nextge) > 0:
            kids = set()
            redi = 0

            while redi < len(nextge):
                red = nextge[redi]
                # if first_round:
                # To know whether some of its extensions were found already, is it necessary?
                # print(">>>", red.lAvailableCols[0], red.lAvailableCols[1])
                red.updateAvailable(rcollect)

                if red.hasAvailableCols():
                    self.logger.printL(2, "Expansion %s.%s\t%s" % (self.count, redi, red), "log", self.getLogId())
                    exts, basis_red, modr = red.prepareExtElems(self.data, self.data.isSingleD(), souvenirs=rcollect)
                    if modr != 0:
                        if modr == 1:
                            rcollect.addItem(basis_red, "W")
                            track = {"do": "prepare-greedy", "trg": [basis_red.getUid()], "src": [red.getUid()], "out": "W"}
                            rcollect.addTrack(track)
                        red.cutOffAvailables()
                    self.charbon.clearStore(basis_red)

                    # WARNING DANGEROUS few extensions for DEBUG!
                    for (side, v, r) in exts:
                        if not self.questionLive():
                            nextge = []
                        else:
                            self.charbon.computeExpand(side, self.data.col(side, v), r, self.data.getColsC())

                    if self.logger.verbosity >= 4:
                        self.logger.printL(4, str(self.charbon.getStore()), "log", self.getLogId())

                    kids = self.charbon.getStore().getFoundReds(self.data)
                    self.logger.printL(2, "Expansion %s.%s returned %s redescriptions" % (self.count, redi, len(kids)), "status", self.getLogId())
                    kid_ids = []
                    for kid in kids:
                        rcollect.addItem(kid, "W")
                        kid_ids.append(kid.getUid())
                    if len(kid_ids) > 0:
                        track = {"do": "expand-greedy", "trg": kid_ids, "src": [basis_red.getUid()], "out": "W"}
                        rcollect.addTrack(track)

                    # SOUVENIRS
                    if self.up_souvenirs:
                        rcollect.updateSeen(kids)

                    # parent has been used remove availables
                    basis_red.cutOffAvailables()
                # self.logger.updateProgress({"rcount": self.count, "generation": len(red), "cand": redi}, 4, self.getLogId())
                self.logger.printL(4, "Candidate %s.%d.%d expanded" % (self.count, len(red), redi), "status", self.getLogId())
                redi += 1

            first_round = False
            self.logger.printL(4, "Generation %s.%d expanded" % (self.count, len(red)), "status", self.getLogId())
            nextge_keys = rcollect.selected("nextge", ids="W", constraints=self.constraints, mute_tracking=True)
            if self.constraints.getCstr("ext_once"):
                nextge = []
            else:
                nextge = [rcollect.getItem(i) for i in nextge_keys]
            for iid in rcollect.getIids():
                if iid not in nextge_keys:
                    rcollect.getItem(iid).cutOffAvailables()

        # Do before selecting next gen to allow tuning the beam
        # ask to update results
        rcollect.selected("partial", ids="W", trg_lid="P", constraints=self.constraints)
        self.logger.logResults(rcollect, "P", self.getLogId())
        return rcollect

    def improveRedescriptions(self, nextge, rcollect=None):
        if len(nextge) == 0 or not self.questionLive() or self.charbon.isTreeBased():  # improvement is only possible greedily
            return rcollect

        if rcollect is None:
            rcollect = RCollection()
        rcollect.resetPartial()
        self.charbon.setStore(self.data.nbRows(), self.constraints.getSSetts(), self.constraints)

        nbS = 0
        for redi, red in enumerate(nextge):
            if self.questionLive():
                rcollect.resetWorking()
                non_anons = [red]
                if red.containsAnon():
                    self.logger.printL(2, "The redescription contains anonymous literal(s)\t%s" % red, "log", self.getLogId())
                    pids = self.setAnonRedescription(red, rcollect, self.charbon, trg_lid="W")
                    non_anons = [rcollect.getItem(iid) for iid in pids]

                for xi, r in enumerate(non_anons):
                    self.logger.printL(2, "Improvement attempt %d.%d\t%s" % (redi, xi, r), "status", self.getLogId())
                    self.improveRedescription(r, rcollect, self.charbon, trg_lid="W")

                rcollect.addItem(red, "P")
                if rcollect.getLen("W") > 0:
                    sids = rcollect.selected("partial", ids="W", trg_lid="P", constraints=self.constraints)
                    self.logger.printL(2, "Improvement found in round %d [%d/%d]" % (redi, len(sids), rcollect.getLen("W")), "status", self.getLogId())
                    if self.logger.verbosity >= 4:
                        for rr in rcollect.getItems("W"):
                            self.logger.printL(4, "%3s %s" % ("*"*(rr.getUid() in sids), rr), "log", self.getLogId())

                else:
                    self.logger.printL(2, "No improvement found in round %d" % redi, "status", self.getLogId())

                self.logger.logResults(rcollect, "P", self.getLogId())
        return rcollect

    def improveRedescription(self, red, rcollect, charbon, skip={}, round_id=0, trg_lid=None, track_rounds=[], try_shorten=True):
        # print("-----------------", round_id, track_rounds)
        # print("IMPRV -(0:0)-\t%s" % red)
        bests = []
        best_score = red.score()
        got_better = False
        for side in [0, 1]:
            queries = [red.query(0), red.query(1)]
            org_q = red.query(side)
            for (ls, q) in org_q.minusOneLiteral():
                if (side, ls, round_id) in skip:
                    # no need to try improving the same element as previous round
                    continue

                xps = self.improveRedescriptionElement(charbon, queries, side, org_q, ls)
                for cand in xps:
                    if cand["red"].score() > best_score:
                        got_better = True
                        bests = [cand]
                        best_score = cand["red"].score()
                    elif got_better and cand["red"].score() == best_score:
                        bests.append(cand)

                if len(q) > 0 and try_shorten:
                    queries[side] = q
                    cand_red = Redescription.fromQueriesPair(queries, self.data)
                    if cand_red.score() > red.score():
                        bests.append({"red": cand_red, "side": side, "ls": ls, "lit": None})

        kid_ids = []
        desc_ids = []
        for ci, best in enumerate(bests):
            # if round_id < 50:
            #     print("%s|_ (%d:%d:%d)%s\t%s" % (" "*round_id, ci, best["lit"] is None, best["red"].getUid(), " "*(50-round_id), best["red"]))
            # else:
            #     print("[...]")
            kid_ids.append(best["red"].getUid())
        if len(kid_ids) > 0:
            track = {"do": "improve", "trg": kid_ids, "src": [red.getUid()], "out": trg_lid}
            rcollect.addTrack(track)

        for ci, best in enumerate(bests):
            rcollect.addItem(best["red"], trg_lid)

            skip_next = {(best["side"], best["ls"], round_id+1): best["red"].score()}
            desc_ids.extend(self.improveRedescription(best["red"], rcollect, charbon, skip_next, round_id+1,
                                                      trg_lid=trg_lid, track_rounds=[ci]+track_rounds, try_shorten=try_shorten))

        return kid_ids + desc_ids

    def setAnonRedescription(self, red, rcollect, charbon, first_round=True, trg_lid=None):
        kid_ids = []
        desc_ids = []

        bests = []
        best_score = 0.
        qs, ls = red.minusAnon()
        if (len(qs[0]) * len(qs[1]) > 0) and (len(ls[0]) + len(ls[1]) > 0):
            noan_red = Redescription.fromQueriesPair(qs, self.data)
            best_score = noan_red.score()

            if first_round:
                kid_ids.append(noan_red.getUid())
                rcollect.addItem(noan_red, trg_lid)

        for side in [0, 1]:
            queries = [red.query(0), red.query(1)]
            org_q = red.query(side)
            for (q, org_ls, lit, ls) in org_q.minusAnonButOne():

                xps = self.improveRedescriptionElement(charbon, queries, side, q, ls)
                for cand in xps:
                    cand["org_ls"] = org_ls
                    if cand["red"].score() > best_score:
                        bests = [cand]
                        best_score = cand["red"].score()
                    elif len(bests) > 0 and cand["red"].score() == best_score:
                        bests.append(cand)

        for ci, best in enumerate(bests):
            kid_ids.append(best["red"].getUid())
        if len(kid_ids) > 0:
            track = {"do": "set-anonymous", "trg": kid_ids, "src": [red.getUid()], "out": trg_lid}
            rcollect.addTrack(track)

        for ci, best in enumerate(bests):
            rcollect.addItem(best["red"], trg_lid)
            queries = [red.query(0), red.query(1)]
            cand_q = queries[best["side"]].copy()
            cand_q.setBukElemAt(best["lit"], best["org_ls"])
            queries[best["side"]] = cand_q
            cand_red = Redescription.fromQueriesPair(queries, self.data)
            desc_ids.extend(self.setAnonRedescription(cand_red, rcollect, charbon, first_round=False, trg_lid=trg_lid))
        return kid_ids+desc_ids

    def improveRedescriptionElement(self, charbon, queries, side, org_q, ls):
        xps = []

        lit = org_q.getBukElemAt(ls)
        op, qc, qd = org_q.partsOCD(ls)

        ext_op = False  # default use conjunction
        offsets = (0, 0)
        if len(qc) == 0:
            ext_op = True  # use disjunction
            queries[side] = qd
            red_basis = Redescription.fromQueriesPair(queries, self.data)
            supports = red_basis.supports().copy()
        else:
            queries[side] = qc
            red_basis = Redescription.fromQueriesPair(queries, self.data)
            supports = red_basis.supports().copy()
            if len(qd) > 0:
                supp_mask, miss_mask = qd.recompute(side, self.data)
                mlparts = supports.moveInterAllOut(side, supp_mask)
                num = supports.getSSetts().sumPartsId(side, supports.getSSetts().IDS_num[True], mlparts)
                den = supports.getSSetts().sumPartsId(side, supports.getSSetts().IDS_den[True], mlparts)
                offsets = (num, den)

        # #### SEARCH FOR IMPROVEMENT
        col = self.data.col(side, lit.colId())
        tmp_cands = charbon.getCandidatesImprov(side, col, red_basis, ext_op, supports, offsets)
        for cci, cand in enumerate(tmp_cands):
            if cci > 0:
                print("Multi")
                pdb.set_trace()
            cand_q = org_q.copy()
            cand_q.setBukElemAt(cand.getLit(), ls)
            queries[side] = cand_q
            cand_red = Redescription.fromQueriesPair(queries, self.data)
            xps.append({"red": cand_red, "side": side, "ls": ls, "lit": cand.getLit()})
        return xps


class Miner(object):

    # INITIALIZATION
    ##################
    def __init__(self, data, params, logger=None, mid=None, qin=None, cust_params={}, filenames={}):
        self.count = "-"
        self.qin = qin
        self.org_data = None
        self.up_souvenirs = True
        self.want_to_live = True
        if mid is not None:
            self.id = mid
        else:
            self.id = 1
        self.data = data

        self.max_processes = params["nb_processes"]
        self.pe_balance = params["pe_balance"]

        # SETTING UP DATA
        row_ids = None
        if "area" in cust_params:
            inw, outw = cust_params.get("in_weight", 1), cust_params.get("out_weight", 1)
            if "in_weight" in params:
                inw = params["in_weight"]
            if "out_weight" in params:
                outw = params["out_weight"]
            weights = dict([(r, outw) for r in range(self.data.nbRows())])
            for old in cust_params["area"]:
                weights[old] = inw
            cust_params["weights"] = weights

        keep_rows = None
        if self.data.hasSelectedRows() or self.data.hasLT():
            keep_rows = self.data.getVizRows({"rset_id": "learn"})
            if "weights" not in cust_params:
                row_ids = dict([(v, [k]) for (k, v) in enumerate(keep_rows)])

        if "weights" in cust_params:
            row_ids = {}
            off = 0
            for (old, mul) in cust_params["weights"].items():
                if keep_rows is None or old in keep_rows:
                    row_ids[old] = [off+r for r in range(mul)]
                    off += mul

        if row_ids is not None:
            self.org_data = self.data
            self.data = self.data.subset(row_ids)

        if logger is not None:
            self.logger = logger
        else:
            self.logger = Log()
        self.constraints = Constraints(params, self.data, filenames=filenames)

        self.charbon = self.initCharbon()
        self.rcollect = RCollection(self.data.usableIds(self.constraints.getCstr("min_itm_c"), self.constraints.getCstr("min_itm_c")),
                                    self.constraints.getCstr("amnesic"), tracks_switch=self.constraints.getCstr("tracks_switch"))
        self.logger.printL(1, "Miner set up (%s)" % self.charbon.getAlgoName(), "log", self.getLogId())
        self.logger.printL(1, "\t%s" % self.data.getInfo(), "log", self.getLogId())
        self.logger.printL(1, "\t%s" % self.rcollect, "log", self.getLogId())

    def getId(self):
        return self.id

    def getLogId(self):
        return "%s" % self.getId()

    def shareRCollect(self):
        return self.rcollect.toShare()

    def shareLogger(self):
        return None
        # return self.logger
        # if not self.logger.usesOutMethods():
        #     return self.logger
        # return None

    def initCharbon(self):
        if self.constraints.getCstr("mining_algo") in TREE_CLASSES:
            if self.data.hasMissing():
                self.logger.printL(1, "THE DATA CONTAINS MISSING VALUES, FALLING BACK ON GREEDY REREMI", "log", self.getLogId())
                return CharbonGMiss(self.constraints, logger=self.shareLogger())
            else:
                return TREE_CLASSES.get(self.constraints.getCstr("mining_algo"), TREE_DEF)(self.constraints, logger=self.shareLogger())
        elif self.constraints.getCstr("mining_algo") in XAUST_CLASSES:
            return XAUST_CLASSES.get(self.constraints.getCstr("mining_algo"), XAUST_DEF)(self.constraints, logger=self.shareLogger())

        else:  # INIT GREEDY CHARBON
            if self.constraints.getCstr("add_condition"):
                if not self.data.isConditional() and self.data.isGeospatial():
                    self.data.prepareGeoCond()

            if CHARBON_MISS_FORCE or self.data.hasMissing():
                return CharbonGMiss(self.constraints, logger=self.shareLogger())
            else:
                return CharbonGStd(self.constraints, logger=self.shareLogger())

    def kill(self):
        self.want_to_live = False

    def questionLive(self):
        if self.want_to_live and self.qin is not None:
            try:
                piece_result = self.qin.get_nowait()
                if piece_result["type_message"] == "progress" and piece_result["message"] == "stop":
                    self.want_to_live = False
            except:
                pass
        return self.want_to_live

# RUN FUNCTIONS
################################

    def filter_run(self, redescs):
        return BatchCollection(redescs).selectedItems(self.constraints.getActionList("redundant"), constraints=self.constraints)

    def part_run(self, cust_params):
        if "reds" in cust_params:
            reds = cust_params["reds"]
        elif "red" in cust_params:
            reds = [cust_params["red"]]
        else:
            reds = []
        if self.org_data is not None:
            for red in reds:
                red.recompute(self.data)

        rcollect = self.shareRCollect()
        if "side" in cust_params:
            rcollect.cutOffSide(1-cust_params["side"])

        self.count = "C"

        self.logger.initProgressPart(self.constraints, reds, rcollect.getNbAvailableCols(), 1, self.getLogId())
        self.logger.clockTic(self.getLogId(), "part run")

        if cust_params.get("task") == "improve":
            self.logger.printL(1, "Improving...", "status", self.getLogId())  # todo ID
            rcollect = self.improveRedescriptions(reds, rcollect)
        else:
            self.logger.printL(1, "Expanding...", "status", self.getLogId())  # todo ID
            rcollect = self.expandRedescriptions(reds, rcollect)

        self.logger.clockTac(self.getLogId(), "part run", "%s" % self.questionLive())
        if not self.questionLive():
            self.logger.printL(1, "Interrupted!", "status", self.getLogId())
        else:
            self.logger.printL(1, "Done.", "status", self.getLogId())
        self.logger.sendCompleted(self.getLogId())
        # print(rcollect)
        return rcollect

    def full_run(self, cust_params={}):
        self.rcollect.resetFinal()
        self.count = 0
        self.logger.printL(1, "Start mining", "status", self.getLogId())

        # progress initialized after listing pairs
        self.logger.clockTic(self.getLogId(), "pairs")
        self.initializeRedescriptions()
        self.logger.clockTac(self.getLogId(), "pairs")

        if not self.constraints.getCstr("only_pairs"):
            self.logger.clockTic(self.getLogId(), "full run")
            self.doExpansions(cust_params)
            self.logger.clockTac(self.getLogId(), "full run", "%s" % self.questionLive())
        else:
            initial_red = self.initial_candidates.getNextRed(self.data, testIni)
            while initial_red is not None and self.questionLive():
                self.rcollect.addItem(initial_red, "W")
                initial_red = self.initial_candidates.getNextRed(self.data, testIni)
            self.rcollect.selected("final", ids="W", trg_lid="F", constraints=self.constraints)
            
        if not self.questionLive():
            self.logger.printL(1, "Interrupted!", "status", self.getLogId())
        else:
            self.logger.printL(1, "Done.", "status", self.getLogId())
                
        self.logger.sendCompleted(self.getLogId())
        self.logger.logAllTracksEnd(self.rcollect, "F", self.getLogId())
        return self.rcollect

    def doExpansions(self, cust_params={}):
        if self.charbon.isIterative():
            self.doExpansionsIterative(cust_params)
        else:
            self.doExpansionsGlobal(cust_params)

    def doExpansionsGlobal(self, cust_params={}):
        initial_terms = self.initial_candidates.getFoundCands()
        if len(initial_terms) > 0 and self.questionLive():

            self.logger.clockTic(self.getLogId(), "fim_exp")
            self.logger.printL(1, "Global expansion", "status", self.getLogId())

            kids = self.charbon.computeExpansions(self.data, initial_terms, logger=self.logger)
            self.logger.printL(2, "Global expansion returned %s redescriptions" % len(kids), "status", self.getLogId())
            kid_ids = []
            for kid in kids:
                self.logger.printL(2, kid, "log", self.getLogId())
                self.rcollect.addItem(kid, "W")
                kid_ids.append(kid.getUid())
            if len(kid_ids) > 0:
                track = {"do": "fim-expand", "trg": kid_ids, "src": [], "out": "W"}
                self.rcollect.addTrack(track)
            self.logger.printL(4, "FIM expansion done", "status", self.getLogId())

            # ################################
            # self.rcollect.selected(self.constraints.getActionList("partial", action_substitutions=[("cut", None)]), ids="W", trg_lid="P", constraints=self.constraints)
            # self.logger.logResults(self.rcollect, "P", self.getLogId())

            # # self.logger.clockTic(self.getLogId(), "select")
            # if self.rcollect.getLen("P") > 0:
            #     self.rcollect.selected("final", ids="F", new_ids="P", trg_lid="F", constraints=self.constraints)

            # self.logger.clockTac(self.getLogId(), "fim_exp", "%s" % self.questionLive())
            # self.logger.logResults(self.rcollect, "F", self.getLogId())
            # ################################

            if self.data.isSymmetric() and self.constraints.getCstr("min_fin_acc") == -1:
                self.rcollect.selected("fim", ids="W", trg_lid="F", constraints=self.constraints)
            else:
                self.rcollect.selected("fimg", ids="W", trg_lid="F", constraints=self.constraints)
            self.logger.clockTac(self.getLogId(), "fim_exp", "%s" % self.questionLive())
            self.logger.logResults(self.rcollect, "F", self.getLogId())

    def doExpansionsIterative(self, cust_params={}):
        nb_round = 0
        initial_red = self.initial_candidates.getNextRed(self.data, testIni)
        while initial_red is not None and self.questionLive():
            self.count += 1

            self.logger.clockTic(self.getLogId(), "expansion_%d-%d" % (self.count, 0))
            self.expandRedescriptions([initial_red], self.rcollect)
            self.logger.updateProgress({"rcount": self.count}, 1, self.getLogId())

            # self.logger.clockTic(self.getLogId(), "select")
            if self.rcollect.getLen("P") > 0:
                self.rcollect.selected("final", ids="F", new_ids="P", trg_lid="F", constraints=self.constraints)

            # self.final["results"] = self.final["batch"].selected("final", constraints=self.constraints)
            # if (self.final["results"] != ttt):
            #     pdb.set_trace()
            #     print("Not same")

            # DEBUG self.final["results"] = range(len(self.final["batch"]))
            # self.logger.clockTac(self.getLogId(), "select")

            self.logger.clockTac(self.getLogId(), "expansion_%d-%d" % (self.count, 0), "%s" % self.questionLive())
            self.logger.logResults(self.rcollect, "F", self.getLogId())
            initial_red = self.initial_candidates.getNextRed(self.data, testIni)


# HIGH LEVEL CALLING FUNCTIONS
################################

####################################################
# INITIAL PAIRS
####################################################


    def initializeRedescriptions(self, ids=None):
        if self.charbon.withInitTerms():
            self.initial_candidates = initCands("T", self.data, self.constraints)
            self.initializeTerms(ids)
        else:
            self.initial_candidates = initCands("P", self.data, self.constraints, save_filename=self.constraints.getCstr("pairs_store"))
            self.initializePairs(ids)

    def initializeTerms(self, ids=None):
        self.logger.printL(1, "Searching for initial terms", "status", self.getLogId())
        self.logger.initProgressFull(self.constraints, None, self.rcollect.getNbAvailableCols(), 1, self.getLogId())

        if ids is None:
            ids = self.data.usableIds(self.constraints.getCstr("min_itm_c"), self.constraints.getCstr("min_itm_c"))

        sides = [0, 1]
        if self.constraints.getSSetts().isDataSymmetric():
            # if self.data.isSingleD() and set(ids[0]) == set(ids[1]):
            sides = [0]
        for side in sides:
            for idl in ids[side]:
                cands = self.charbon.computeInitTerms(side, self.data.col(side, idl))
                self.logger.printL(4, "Generated %d initial terms from variable %d %d" % (len(cands), side, idl), "status", self.getLogId())
                for cand in cands:
                    self.logger.printL(6, str(cand), "log", self.getLogId())
                    self.initial_candidates.add(cand)
        self.logger.printL(1, self.initial_candidates.msgFound(), "log", self.getLogId())
        self.initial_candidates.setExploredDone()
        # self.logger.sendCompleted(self.getLogId())

    def initializePairs(self, ids=None):
        ########################################
        ### collecting pairs traces, initialize
        timings = {}
        candidate_details = None
        if len(self.constraints.getCstr("pairs_traces")) != 0:
            candidate_details = {"unaccurate_pairs": [], "stats": {}}
        ########################################
        # SELECTION USING FOLDS
        # folds = numpy.array(self.data.col(0,-1).getVector())
        # counts_folds = 1.*numpy.bincount(folds)
        # nb_folds = len(counts_folds)
        # Loading pairs from file if filename provided
        loaded, done = self.initial_candidates.loadFromFile(self.data)
        if not loaded or done is not None:

            self.logger.printL(1, "Searching for initial pairs", "status", self.getLogId())
            self.logger.clockTic(self.getLogId(), "pairs_candidates")
            explore_list = self.getInitExploreList(ids, done)
            if candidate_details is not None and "stats" in candidate_details:
                candidate_details["stats"]["nb_explored_pairs"] = len(explore_list)
            timings["pairs_candidates"] = self.logger.clockTac(self.getLogId(), "pairs_candidates")
            self.logger.clockTic(self.getLogId(), "pairs_eval")
            self.logger.initProgressFull(self.constraints, explore_list, self.rcollect.getNbAvailableCols(), 1, self.getLogId())

            self.initial_candidates.setExploreList(explore_list, done=done)
            self.charbon.setStore(self.initial_candidates)
            total_pairs = len(explore_list)
            for pairs, (idL, idR, dtlsL, dtlsR, pload) in enumerate(explore_list):
                if not self.questionLive():
                    self.initial_candidates.saveToFile()
                    return

                self.logger.updateProgress({"rcount": self.count, "pair": pairs, "pload": pload})
                if pairs % 100 == 0:
                    self.logger.printL(3, "Searching pair %d/%d (%i <=> %i)" %
                                       (pairs, total_pairs, idL, idR), "status", self.getLogId())
                    self.logger.updateProgress(level=3, id=self.getLogId())
                elif pairs % 10 == 0:
                    self.logger.printL(7, "Searching pair %d/%d (%i <=> %i)" %
                                       (pairs, total_pairs, idL, idR), "status", self.getLogId())
                    self.logger.updateProgress(level=7, id=self.getLogId())
                else:
                    self.logger.printL(10, "Searching pair %d/%d (%i <=> %i)" %
                                       (pairs, total_pairs, idL, idR), "status", self.getLogId())

                ### COMPUTE PAIRS
                cands = self.charbon.computePair(self.data.col(0, idL), self.data.col(1, idR), self.data.getColsC(), self.data)
                if candidate_details is not None and len(list(cands)) == 0:
                    candidate_details["unaccurate_pairs"].append((idL, idR))
                for cand in cands:
                    self.logger.printL(6, str(cand), "log", self.getLogId())
                    # ########
                    # ######## Filter pair candidates on folds distribution
                    # rr = Redescription.fromInitialPair((literalsL[i], literalsR[i]), self.data)
                    # bcount = numpy.bincount(folds[list(rr.getSuppI())], minlength=nb_folds)
                    # if len(numpy.where(bcount > 0)[0]) > 1:
                    #     bb = bcount/counts_folds
                    #     # bpr = bcount/float(numpy.sum(bcount))
                    #     # entropS = -numpy.sum(numpy.log(bpr)*bpr)
                    #     bpr = bb/numpy.max(bb)
                    #     score = numpy.sum(bpr)
                    #     # entropM = -numpy.sum(numpy.log(bpr)*bpr)
                    #     if score > 1.5:
                    #         self.logger.printL(6, "\tfolds count: %f, %d, %s" % (score, numpy.sum(bcount), bcount), "log", self.getLogId())
                    #         ####
                    #         # self.initial_candidates.add(literalsL[i], literalsR[i], {"score": scores[i], 0: idL, 1: idR})
                    #         self.initial_candidates.add(literalsL[i], literalsR[i], {"score": score, 0: idL, 1: idR})
                    # # if pairs % 50 == 0 and pairs > 0:
                    # #     exit()
                self.initial_candidates.addExploredPair((idL, idR))

            self.logger.printL(1, self.initial_candidates.msgFound(), "log", self.getLogId())
            self.logger.updateProgress(level=1, id=self.getLogId())

            # Saving pairs to file if filename provided
            self.initial_candidates.setExploredDone()
            self.initial_candidates.saveToFile()
            timings["pairs_eval"] = self.logger.clockTac(self.getLogId(), "pairs_eval")
        else:
            self.logger.initProgressFull(self.constraints, None, self.rcollect.getNbAvailableCols(), 1, self.getLogId())
            self.logger.printL(1, self.initial_candidates.msgLoaded(), "log", self.getLogId())
        ########################################
        ### collecting pairs traces, logging
        if len(self.constraints.getCstr("pairs_traces")) != 0:            
            all_pairs = []
            if candidate_details is not None:
                if "stats" in candidate_details:
                    candidate_details["stats"]["nb_unaccurate_pairs"] = len(candidate_details["unaccurate_pairs"])
                    candidate_details["stats"]["nb_accurate_pairs"] = len(self.initial_candidates)

                for k, c in enumerate(candidate_details["unaccurate_pairs"]):
                    if type(c) is tuple:
                        if self.data.isTypeId(self.data.col(0, c[0]).typeId(), "Boolean") and self.data.isTypeId(self.data.col(1, c[1]).typeId(), "Boolean"):
                            cand, candid = self.charbon.makePairCand(self.data.col(0, c[0]), self.data.col(1, c[1]), store=False, no_const=True)
                            all_pairs.append((cand.getAcc(), cand.getLit(0), cand.getLit(1), "+", k, 0, 0))
                        else:
                            all_pairs.append((-1, c[0], c[1], "+", k, 0, 0))
                    else:
                        all_pairs.append((c.getAcc(), c.getLit(0), c.getLit(1), "+", k, 0, 0))
            for k, cand in self.initial_candidates.items():
                all_pairs.append((cand.getAcc(), cand.getLit(0), cand.getLit(1), "-", k, cand.ssizes.lpart(2), cand.pValRed()))
            all_pairs.sort()

            #### LOGGING RESULTS AND TIME
            import sys
            if self.constraints.getCstr("pairs_traces") == "-":
                f = sys.stdout
            else:
                f = open(self.constraints.getCstr("pairs_traces"), "w")

            timings["total"] = sum(timings.values(), datetime.timedelta())
            sum(timings.values(), datetime.timedelta())
            f.write("### "+"\t".join(["%s %s" % (k,v) for (k,v) in timings.items()])+"\n")
            f.write("### "+"\t".join(["%s %s" % (k,v) for (k,v) in candidate_details.get("stats", {}).items()])+"\n")

            for c in all_pairs:
                f.write("%s%s\t%s\t%s\t%f\t%s\t%s\n" % (c[3], c[4], c[1], c[2], c[0], c[5], c[6]))
                
            if self.constraints.getCstr("pairs_traces") == "-":
                f.close()
        # exit()  # STOP AFTER PAIRS
        return self.initial_candidates

    def getInitExploreList(self, ids, done=set()):
        explore_list = []
        if ids is None:
            ids = self.data.usableIds(self.constraints.getCstr("min_itm_c"), self.constraints.getCstr("min_itm_c"))

        # ### WARNING DANGEROUS few pairs for DEBUG!
        # if self.data.nbCols(0) > 100:
        # ids = [[13], [0]]
        # ids = [[1], [9]]
        # ids = [[2], [29]]
        # ids = [[0], [75]]
        for idL in ids[0]:
            for idR in ids[1]:
                if not self.data.arePairTypesIn(idL, idR, tset=self.constraints.getCstr("inits_types_exclude")) and self.data.areGroupCompat(idL, idR) and \
                        (not self.data.isSingleD() or idR > idL or idR not in ids[0] or idL not in ids[1]):
                    if done is None or (idL, idR) not in done:
                        explore_list.append((idL, idR, None, None, self.getPairLoad(idL, idR)))
                    else:
                        self.logger.printL(3, "Loaded pair (%i <=> %i)" % (idL, idR), "status", self.getLogId())
        return explore_list

    def getPairLoad(self, idL, idR):
        # pdb.set_trace()
        # print(idL, idR, eval("0b"+"".join(["%d" %(((idL+idR)%10)%i ==0) for i in [8,4,2]])))
        ## + ((idL + idR)%10)/1000.0
        return max(1, self.data.col(0, idL).getNbValues() * self.data.col(1, idR).getNbValues()/50) + 1./(1+((idL + idR) % 10))
        # return max(1, self.data.col(0, idL).getNbValues()* self.data.col(1, idR).getNbValues()/50)
        # return PAIR_LOADS[self.data.col(0, idL).type_id-1][self.data.col(1, idR).type_id-1]


####################################################
# REDS EXPANSIONS
####################################################

    def expandRedescriptions(self, nextge, rcollect=None):
        return ExpMiner(self.getId(), self.count, self.data, self.charbon, self.constraints, self.logger,
                        question_live=self.questionLive, up_souvenirs=self.up_souvenirs).expandRedescriptions(nextge, rcollect)

    def improveRedescriptions(self, nextge, rcollect=None):
        return ExpMiner(self.getId(), self.count, self.data, self.charbon, self.constraints, self.logger,
                        question_live=self.questionLive, up_souvenirs=self.up_souvenirs).improveRedescriptions(nextge, rcollect)
