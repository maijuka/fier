# Copyright 2024 Maiju Karjalainen and Esther Galbrun

# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at

#     http://www.apache.org/licenses/LICENSE-2.0

# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

import re
import datetime
from sympy import nextprime
import itertools
import numpy
import operator
import heapq 

try:
    from classCol import NumColM
    from classRedescription import Redescription
    from classMiner import Miner, testIni

except ModuleNotFoundError:
    from .classCol import NumColM
    from .classRedescription import Redescription
    from .classMiner import Miner, testIni

import pdb


class MinerLSH(Miner):

    # INITIALIZATION
    ##################
    def __init__(self, data, params, logger=None, mid=None, qin=None, cust_params={}, filenames={}):
        Miner.__init__(self, data, params, logger, mid, qin, cust_params, filenames)

        ########## LSH ADDITION BEGINS
        seed = self.constraints.getCstr("lsh_random_seed")
        if seed <= 0:
            rng = numpy.random.default_rng()
        else:
            rng = numpy.random.Generator(numpy.random.MT19937(numpy.random.SeedSequence(seed)))

        exx = self.constraints.getSSetts().Exx
        exo = self.constraints.getSSetts().Exo
        eox = self.constraints.getSSetts().Eox
        eoo = self.constraints.getSSetts().Eoo

        self.nb_bands = self.constraints.getCstr("lsh_nb_bands_ext")
        self.band_size = self.constraints.getCstr("lsh_nb_rows_ext")
        # if self.constraints.getCstr("lsh_nb_rows") == 0:
        #     self.band_size = max(1, int(numpy.log(1/self.nb_bands)/numpy.log(self.constraints.getCstr("init_minscore"))))
        # else:
        #     self.band_size = self.constraints.getCstr("lsh_nb_rows")
        max_id = self.data.nbRows()                
        next_prime = nextprime(max_id)
        self.hash_matrix = []
        for i in range(self.nb_bands):
            coeffA = rng.integers(0,max_id, self.band_size)
            coeffB = rng.integers(0,max_id, self.band_size)
            self.hash_matrix.append((numpy.outer(numpy.arange(max_id), coeffA) + numpy.tile(coeffB, (max_id, 1))) % next_prime)

        self.lsh_method = self.constraints.getCstr("lsh_ext_method")
        if self.lsh_method == "replace":
            self.lshext = {False:{exx:1,exo:0,eox:0,eoo:0},True:{exx:1,exo:1,eox:1,eoo:0}}
            self.hash_indices = rng.integers(0,self.data.nbRows(),(self.nb_bands,self.band_size))
        elif self.lsh_method == "drop":
            # and,left and,right or,left or,right
            self.lshext = {(False,0):{exx:1,exo:0},(False,1):{exx:1,eox:0},(True,0):{eox:1,eoo:0},(True,1):{exo:1,eoo:0}}

        self.drop = {(False,0):{exx:0,exo:0,eox:1,eoo:1},(False,1):{exx:0,eox:0,exo:1,eoo:1}}
        ########## LSH ADDITION ENDS

    def doExpansions(self, cust_params={}):
        buckets = self.makeLSHBuckets()
        if self.constraints.getCstr("method_lsh_ext") == "priorityqueue":
            self.expandRedescriptionsPQLSH(buckets)
        elif self.constraints.getCstr("method_lsh_ext") == "greedy":
            red = self.initial_candidates.getNextRed(self.data, testIni)
            while red is not None and self.questionLive():
                self.count += 1
                self.logger.clockTic(self.getLogId(), "expansion_%d-%d" % (self.count, 0))
                self.rcollect.resetPartial([red])
                self.expandRedescriptionsGreedyLSH([red], self.rcollect, buckets)
                self.logger.updateProgress({"rcount": self.count}, 1, self.getLogId())
                if self.rcollect.getLen("P") > 0:
                    self.rcollect.selected("final", ids="F", new_ids="P", trg_lid="F", constraints=self.constraints)

                self.logger.clockTac(self.getLogId(), "expansion_%d-%d" % (self.count, 0), "%s" % self.questionLive())
                self.logger.logResults(self.rcollect, "F", self.getLogId())
                red = self.initial_candidates.getNextRed(self.data, testIni)
    

    def expandRedescriptionsPQLSH(self, buckets):
        reds = []
        red = self.initial_candidates.getNextRed(self.data, testIni)
        max_var = [self.constraints.getCstr("max_var", side=0), self.constraints.getCstr("max_var", side=1)]
        while red is not None and self.questionLive():
            self.rcollect.addItem(red, "P")
            reds.append((-1*red.getAcc(),red))
            red.initAvailable(self.rcollect, self.data, max_var)
            red = self.initial_candidates.getNextRed(self.data, testIni)
        # reds.sort(key=operator.itemgetter(0))
        heapq.heapify(reds)
        self.charbon.setStore(self.data.nbRows(), self.constraints.getSSetts(), self.constraints)
        i = 0
        while reds and (i < self.constraints.getCstr("lsh_nb_total_exts")):
            # red = reds.pop()[1]
            red = heapq.heappop(reds)[1]
            red.updateAvailable(self.rcollect)
            if red.hasAvailableCols():
                exts = self.findLSHExts(red,buckets)
                self.logger.printL(2, "Expanding %s" % (str(red)), "status", self.getLogId())
                self.logger.printL(2, "Found %s candidate extensions" % (len(exts)), "status", self.getLogId())
                best = None
                self.charbon.clearStore(red)
                for ((side,col),(r,op)) in exts:
                    if (r.length(side)==self.constraints.getCstr("max_var_s%i" % (side))):
                        continue
                    if r.isAvailableCol(side,col.getId(),self.data.isSingleD()):
                        self.charbon.computeExpand(side, self.data.col(side, col.getId()), r, self.data.getColsC(), force_ops=[op])
                kids = self.charbon.getStore().getFoundReds(self.data)
                for kid in kids:
                    if best is None or (kid.getAcc()>best.getAcc()):
                        best = kid
                    self.rcollect.addItem(kid, "W")
                    self.logger.printL(2, "Extension: \n %s" % (str(kid)), "status", self.getLogId())
                if self.up_souvenirs:
                    self.rcollect.updateSeen(kids)
                if best is not None:
                    # best.initAvailable(self.rcollect, self.data, max_var)
                    self.rcollect.addItem(best, "P")
                    if (not self.constraints.getCstr("ext_once")) and ((best.length(0)<self.constraints.getCstr("max_var_s0")) or (best.length(1)<self.constraints.getCstr("max_var_s1"))):
                        heapq.heappush(reds,(-1*best.getAcc(),best))
                        # reds.append((best.getAcc(),best))
                        # reds.sort(key=operator.itemgetter(0))
                    i += 1
        self.logger.clockTic(self.getLogId(), "final_selection") 
        self.rcollect.selected("final", ids="P", trg_lid="F", constraints=self.constraints)
        self.logger.clockTac(self.getLogId(), "final_selection")

    def expandRedescriptionsGreedyLSH(self, nextge, rcollect, buckets): 
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
                red.updateAvailable(rcollect)

                if red.hasAvailableCols():
                    self.logger.printL(2, "Expansion %s.%s\t%s" % (self.count, redi, red), "log", self.getLogId())
                    exts = self.findLSHExts(red,buckets)
                    self.charbon.clearStore(red)
                    for ((side,col),(r,op)) in exts:
                        if r.isAvailableCol(side,col.getId(),self.data.isSingleD()):
                            self.charbon.computeExpand(side, self.data.col(side, col.getId()), r, self.data.getColsC(), force_ops=[op])
                    if self.logger.verbosity >= 4:
                        self.logger.printL(4, str(self.charbon.getStore()), "log", self.getLogId())

                    kids = self.charbon.getStore().getFoundReds(self.data)
                    self.logger.printL(2, "Expansion %s.%s returned %s redescriptions" % (self.count, redi, len(kids)), "status", self.getLogId())
                    kid_ids = []
                    for kid in kids:
                        rcollect.addItem(kid, "W")
                        kid_ids.append(kid.getUid())
                    if len(kid_ids) > 0:
                        track = {"do": "expand-greedy", "trg": kid_ids, "src": [red.getUid()], "out": "W"}
                        rcollect.addTrack(track)

                    # SOUVENIRS
                    if self.up_souvenirs:
                        rcollect.updateSeen(kids)

                    # parent has been used remove availables
                    red.cutOffAvailables()
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

    def findLSHExts(self,red,buckets):
        exts = set()
        for bi in range(self.nb_bands):
            bucket = buckets[bi]
            if self.constraints.getCstr("lsh_ext_jaccard"):
                hm = self.hash_matrix[bi]
            else:
                hash_inds = self.hash_indices[bi]
            for op,v in self.lshext.items():
                if self.constraints.getCstr("lsh_ext_jaccard"):
                    redsupp = [i for i,x in enumerate(red.supports().getVectorABCD()) if v.get(x) == 1]
                    sig_red = tuple(numpy.min(hm[redsupp,:],axis=0))
                else:
                    supps = [red.supports().getVectorABCD()[i] for i in hash_inds]
                    sig_red = tuple([v.get(item,item) for item in supps]) 
                if sig_red in bucket:
                    exts.update(itertools.product(bucket[sig_red],[(red,op)]))
        return exts

    def makeLSHBuckets(self):
        buckets = {}
        for bi in range(self.nb_bands):
            buckets[bi] = {} 
            for side in [0,1]:
                for col in self.data.colsSide(side):
                    if self.data.isTypeId(col.typeId(), "Numerical"):
                        method = self.constraints.getCstr("method_buckets")
                        nbbuk = self.constraints.getCstr("nb_ext_buckets")
                        for i in range(self.constraints.getCstr("lsh_bucket_mult")):
                            buks = col.buckets(method, {"nb_buckets": nbbuk})
                            nbB = len(buks[1])
                            for b in range(nbB):
                                if buks[2] == b:
                                    supp = col.modeSupp()
                                else:
                                    supp = buks[0][b]
                                sig = self.computeSignature(supp,bi)
                                if sig in buckets[bi]:
                                    buckets[bi][sig].append((side,col))
                                else: 
                                    buckets[bi][sig]= [(side,col)]
                            nbbuk *= 2
                    elif self.data.isTypeId(col.typeId(), "Boolean"):
                        if col.supp():
                            sig = self.computeSignature(col.supp(),bi)
                        else:
                            continue
                        if sig in buckets[bi]:
                            buckets[bi][sig].append((side,col))
                        else: 
                            buckets[bi][sig]= [(side,col)]
                    elif self.data.isTypeId(col.typeId(), "Categorical"):
                        for ci, cname in enumerate(col.cats()):
                            supp = col.getSCat(ci)
                            sig = self.computeSignature(supp,bi)
                            if sig in buckets[bi]:
                                buckets[bi][sig].append((side,col))
                            else: 
                                buckets[bi][sig]= [(side,col)]
        return buckets
    
    def computeSignature(self,supp,bi):
        if self.constraints.getCstr("lsh_ext_jaccard"):
            sig = tuple(numpy.min(self.hash_matrix[bi][list(supp), :], axis=0))
        else:
            sig = tuple([int(i in supp) for i in self.hash_indices[bi]])
        return sig

#############################################

# HIGH LEVEL CALLING FUNCTIONS
################################

####################################################
# INITIAL PAIRS
####################################################


        
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

            self.logger.printL(1, "Searching for initial pairs with LSH", "status", self.getLogId())
            self.logger.clockTic(self.getLogId(), "pairs_candidates")
            explore_list = self.getInitExploreListLSH(ids, candidate_details)
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

                ## Evaluate the candidate
                cand, candid = self.charbon.makePairCand(self.data.col(0, idL), self.data.col(1, idR), dtlsL, dtlsR) #, no_const=True)
                if candid is not None:
                    cands = [cand]
                else:
                    cands = []
                    if candidate_details is not None:
                        if cand is None:
                            candidate_details["unaccurate_pairs"].append((idL, idR))
                        else:
                            candidate_details["unaccurate_pairs"].append(cand)                            

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
        ########################################
        return self.initial_candidates

    def getInitExploreListLSH(self, ids=None, candidate_details=None):
        if ids is None:
            ids = self.data.usableIds(self.constraints.getCstr("min_itm_c"), self.constraints.getCstr("min_itm_c"))
            
        ### random sampling needs to be done with a random number generator, allowing to set the seed, so that repeatable experiments can be performed
        seed = self.constraints.getCstr("lsh_random_seed")
        if seed <= 0:
            rng = numpy.random.default_rng()
        else:
            rng = numpy.random.Generator(numpy.random.MT19937(numpy.random.SeedSequence(seed)))

        nb_bands = self.constraints.getCstr("lsh_nb_bands")
        if self.constraints.getCstr("lsh_nb_rows") == 0:
            band_size = max(1, int(numpy.log(1/nb_bands)/numpy.log(self.constraints.getCstr("init_minscore"))))
        else:
            band_size = self.constraints.getCstr("lsh_nb_rows")
        max_id = self.data.nbRows()                
        next_prime = nextprime(max_id)
        
        mout = self.constraints.getCstrSimple("min_itm_out")
        mscore = self.constraints.getCstr("init_minscore")
        
        if self.constraints.getCstr("lsh_min_itm_out") == 0 and mscore != 0:
            lmout = ((mscore*mout)-(mscore*max_id)+mout+max_id)/(2*mout)
        else:
            lmout = self.constraints.getCstr("lsh_min_itm_out")

        candidate_pairs = set()        
        for bi in range(nb_bands):
            coeffA = rng.integers(0,max_id, band_size)
            coeffB = rng.integers(0,max_id, band_size)

            hash_matrix = (numpy.outer(numpy.arange(max_id), coeffA) + numpy.tile(coeffB, (max_id, 1))) % next_prime
            
            col_buckets = {}
            col_signatures_0 = {}
            for side in [0, 1]:
                for cid in ids[side]:
                    if self.data.isTypeId(self.data.col(side, cid).typeId(), "Boolean") and len(self.data.col(side, cid).supp()) >= band_size: 
                        if side == 1 and self.data.isSingleD() and cid in col_signatures_0:
                            signature = col_signatures_0[cid]
                        else:                        
                            signature = tuple(numpy.min(hash_matrix[list(self.data.col(side, cid).supp()), :], axis=0)) 
                        if side == 0:
                            col_signatures_0[cid] = signature
                            if signature in col_buckets:
                                col_buckets[signature][0].append((cid,()))
                            else:
                                col_buckets[signature] = [[(cid,())], []]
                        elif signature in col_buckets:
                            col_buckets[signature][1].append((cid,()))
                    elif self.data.isTypeId(self.data.col(side, cid).typeId(), "Categorical"):
                        signatures = numpy.zeros((len(self.data.col(side, cid).cats()),band_size))
                        for ci, cname in enumerate(self.data.col(side, cid).cats()):
                            supp = self.data.col(side, cid).getSCat(ci)
                            sig = tuple(numpy.min(hash_matrix[list(supp), :], axis=0))
                            signatures[ci,:] = sig
                            if side == 0:
                                if sig in col_buckets:
                                    col_buckets[sig][0].append((cid, (ci)))
                                else:
                                    col_buckets[sig] = [[(cid, (ci))], []]
                            elif sig in col_buckets:
                                col_buckets[sig][1].append((cid, (ci)))
                        if self.constraints.getCstr("multi_cats"):
                            if self.constraints.getCstr("lsh_comb_length") > 0:
                                clen = min(self.constraints.getCstr("lsh_comb_length"), len(self.data.col(side, cid).cats())-1)
                            else:
                                clen = len(self.data.col(side, cid).cats())-1
                            if clen >= 2:
                                for sub in range(2,clen+1):
                                    for comb in itertools.combinations(numpy.arange(0,len(self.data.col(side, cid).cats()),1), sub):
                                        min_sig = tuple(numpy.min(signatures[comb,:] , axis=0).astype(int))
                                        if side == 0:
                                            if min_sig in col_buckets:
                                                col_buckets[min_sig][0].append((cid, comb))
                                            else:
                                                col_buckets[min_sig] = [[(cid, comb)], []]
                                        elif min_sig in col_buckets:
                                            col_buckets[min_sig][1].append((cid, comb))
                    elif self.data.isTypeId(self.data.col(side, cid).typeId(), "Numerical"):
                        method = self.constraints.getCstr("method_buckets") # bellman, collapsed etc.
                        if method == "bellman":
                            buckets = self.data.col(side, cid).buckets(method, {"nb_buckets": self.constraints.getCstr("nb_buckets"), "bellman_criterion": self.constraints.getCstr("bellman_criterion")})
                        else:
                            buckets = self.data.col(side, cid).buckets(method, {"nb_buckets": self.constraints.getCstr("nb_buckets")})
                        bUp = NumColM.buk_ind_maxes(buckets)
                        nbB = len(buckets[1])
                        # save the signature for each separate bucket
                        signatures = numpy.zeros((nbB,band_size))
                        for bi in range(nbB):
                            if buckets[2] == bi:
                                n_sig = tuple(numpy.min(hash_matrix[list(self.data.col(side, cid).modeSupp()), :], axis=0))
                            else:
                                n_sig = tuple(numpy.min(hash_matrix[list(buckets[0][bi]), :], axis=0))
                            signatures[bi,:] = n_sig
                        largest_bius = {} # to track subintervals with same signature
                        for bil in range (nbB-1):
                            biu = nbB-1
                            supp = sum(len(buckets[0][i]) for i in range(bil,biu+1))
                            while supp > max_id - (self.constraints.getCstrSimple("min_itm_out")*lmout) and biu > 0:
                                supp -= len(buckets[0][biu])
                                biu -= 1
                            while supp >= self.constraints.getCstrSimple("min_itm_in") and biu >= 0:
                                min_sig = tuple(numpy.min(signatures[bil:biu+1,:],axis=0).astype(int)) # take minimum hash values from the separate buckets
                                if biu>largest_bius.get(min_sig,-1):
                                    if side == 0:
                                        if min_sig in col_buckets:
                                            col_buckets[min_sig][0].append((cid, (bil,biu,method)))
                                        else:
                                            col_buckets[min_sig] = [[(cid, (bil,biu,method))], []]
                                    elif min_sig in col_buckets: # signature from right side is added only if the same signature was already found on the left side
                                        col_buckets[min_sig][1].append((cid, (bil,biu,method)))  
                                    largest_bius[min_sig] = biu
                                biu -= 1
                                if biu >= 0:
                                    supp -= len(buckets[0][biu+1])
                  
            for signature, (cidsL, cidsR) in col_buckets.items():
                candidate_pairs.update(itertools.product(cidsL, cidsR))

        if candidate_details is not None and "stats" in candidate_details:
            candidate_details["stats"]["nb_candidate_pairs"] = len(candidate_pairs)
        
        found_dict = {}
        for ((idL,cL),(idR,cR)) in candidate_pairs:
            if (idL,idR) in found_dict:
                found_dict[(idL,idR)].append((cL,cR))
            else:
                found_dict[(idL,idR)] = [(cL,cR)]
        
        found_pairs = []
        for (idL,idR) in list(found_dict): 
            if self.data.arePairTypesIn(idL, idR, tset=self.constraints.getCstr("inits_types_exclude")) or not self.data.areGroupCompat(idL, idR) or \
            (self.data.isSingleD() and not idR > idL and idR in ids[0] and idL in ids[1]):
                del found_dict[(idL,idR)]
            elif self.data.isTypeId(self.data.col(0, idL).typeId(), "Numerical") and self.data.isTypeId(self.data.col(1, idR).typeId(), "Numerical") and \
            self.constraints.getCstr("lsh_filtering") and len(found_dict[(idL,idR)]) > 1:
                pairs = numpy.zeros((len(found_dict[(idL,idR)]),6)).astype(int)
                for ip, ((bilL,biuL,mL),(bilR,biuR,mR)) in enumerate(found_dict[(idL,idR)]):
                    pairs[ip,:5] = [ip,bilL,biuL,bilR,biuR] 
                leftsorted = pairs[pairs[:,2].argsort()][::-1]
                leftsorted = leftsorted[leftsorted[:,1].argsort(kind="mergesort")]
                rightsorted = pairs[pairs[:,4].argsort()][::-1]
                rightsorted = rightsorted[rightsorted[:,3].argsort(kind="mergesort")]
                maxbiuL = -1
                maxbiuR = -1
                subs = {}
                for jp in range(len(found_dict[(idL,idR)])):
                    if leftsorted[jp,2] <= maxbiuL:
                        subs[leftsorted[jp,0]] = set(leftsorted[:jp][leftsorted[:jp,2]>=leftsorted[jp,2]][:,0])
                    else:
                        maxbiuL = leftsorted[jp,2]      
                for jp in range(len(found_dict[(idL,idR)])):
                    if rightsorted[jp,4] <= maxbiuR:
                        if bool(subs.get(rightsorted[jp,0],set()).intersection(set(rightsorted[:jp][rightsorted[:jp,4]>=rightsorted[jp,4]][:,0]))):
                            pairs[int(rightsorted[jp,0]),5] += 1
                    else:
                        maxbiuR = rightsorted[jp,4]
                for pair in pairs[pairs[:,5]==0]:
                    # is method always same? I guess it should be
                    found_pairs.append((idL,idR, {"val_range":(pair[1],pair[2],mL)}, {"val_range":(pair[3],pair[4],mR)},1))
            else:
                for (cL,cR) in found_dict[(idL,idR)]:
                    found_pairs.append((idL, idR, {"val_range":cL}, {"val_range":cR}, 1))
 
        # found_pairs = [(idL, idR, {"val_range":cL}, {"val_range":cR}, 1) for ((idL,cL), (idR,cR)) in candidate_pairs \
        #                 if not self.data.arePairTypesIn(idL, idR, tset=self.constraints.getCstr("inits_types_exclude")) and self.data.areGroupCompat(idL, idR) and \
        #                 (not self.data.isSingleD() or idR > idL or idR not in ids[0] or idL not in ids[1])]     
               
        return found_pairs

