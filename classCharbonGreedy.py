import numpy
import itertools
try:
    from classCol import ColM
    from classSParts import SParts, cmp_lower, cmp_greater, cmp_leq, cmp_geq, tool_ratio
    from classQuery import *
    from classRedescription import Redescription
    from classConstraints import Constraints
    from classCandidates import XStore, PairCand, ExtsBatch, PairsBatch, ExtensionWarning
    from classCharbon import Charbon
except ModuleNotFoundError:
    from .classCol import ColM
    from .classSParts import SParts, cmp_lower, cmp_greater, cmp_leq, cmp_geq, tool_ratio
    from .classQuery import *
    from .classRedescription import Redescription
    from .classConstraints import Constraints
    from .classCandidates import XStore, PairCand, ExtsBatch, PairsBatch, ExtensionWarning
    from .classCharbon import Charbon

import pdb


class CharbonGreedy(Charbon):
    # cannot be used directly, need complentary methods specific to Standard or Missing values subclasses

    name = "G"
    FMT_FLG = "{:04b}"
    FLG_IN = 0b1000  # broke min supp in constraint
    FLG_OUT = 0b0100  # broke min supp out constraint
    FLG_CONT = 0b0010  # broke min cont constraint
    STG_CMP = 1  # dropped upon comparison
    STG_OK = 0  # selected upon comparison

    @classmethod
    def is_raised_flag(tcl, v, flag):
        return (v & flag) > 0

    @classmethod
    def disp_status(tcl, status):
        return "".join([tcl.is_raised_flag(status, f)*l for f, l in [(tcl.FLG_CONT, "C"), (tcl.FLG_IN, "I"), (tcl.FLG_OUT, "O")]])
        # return tcl.FMT_FLG.format(status)

    def __init__(self, constraints, logger=None):
        Charbon.__init__(self, constraints, logger)
        self.setOffsets()
        self.store = None
        self.temp_store = None

    def isTreeBased(self):
        return False

    def handlesMiss(self):
        return False

    def withInitTerms(self):
        return False

    ###########################
    # Storage of results
    def startStoreDivert(self):
        self.temp_store = []

    def stopStoreDivert(self):
        tt = self.temp_store
        self.temp_store = None
        return tt

    def getStore(self):
        return self.store

    def addPairToStore(self, *pair_pieces):
        if self.temp_store is not None:
            self.temp_store.append(pair_pieces)
        elif self.store is not None:
            self.store.addPair(*pair_pieces)

    def addExtToStore(self, *cand_pieces):
        if self.temp_store is not None:
            self.temp_store.append(cand_pieces)
        elif self.store is not None:
            self.store.addExt(*cand_pieces)

    def clearStore(self, current=None):
        self.store.clear(current)

    def setStore(self, N=0, constraints=None, ssetts=None, current="E"):  # "E/P" for extension or pair...
        if isinstance(N, XStore):  # Already existing store
            self.store = N
        elif current is None:
            self.store = ExtsBatch(N, constraints, ssetts)
        elif type(current) is str:
            if current == "P":  # empty redescription, eg. starting from anon lits -> actually yields pairs
                self.store = PairsBatch(N, constraints, ssetts)
            else:
                self.store = ExtsBatch(N, constraints, ssetts)
        elif isinstance(current, Redescription):
            if len(current) == 0:  # empty redescription, eg. starting from anon lits -> actually yields pairs
                self.store = PairsBatch(N, constraints, ssetts, current)
            else:
                self.store = ExtsBatch(N, constraints, ssetts, current)
        else:
            raise Warning("Implementation hanging...")
            ## self.store = ExtensionsCombBatch(N, constraints, current)

    ###################
    # OFFSETS
    def getOffsets(self):
        return self.offsets

    def setOffsets(self, offsets=(0, 0)):
        self.offsets = offsets

    def offsetsNonZero(self):
        return (self.offsets[0]+self.offsets[1]) != 0

    ##################################################
    # TOOLS METHODS
    def unconstrained(self, no_const):
        return no_const or self.offsetsNonZero()

    def ratio(self, num, den):
        return tool_ratio(num, den)

    def offset_ratio(self, num, den):
        return tool_ratio(num+self.offsets[0], den+self.offsets[1])

    def updateACT(self, best, lit, side, op, neg, clp_tmi, is_cond=False):
        tmp_adv, tmp_clp, stt = self.getAC(side, op, neg, clp_tmi, is_cond)
        if tmp_adv is None:
            return best, stt
        elif cmp_lower(best[0], tmp_adv):
            return (tmp_adv, tmp_clp, [side, op, neg, lit]), self.STG_OK
        else:
            return best, self.STG_CMP

    def updateACTList(self, typs, best, lit, side, op, neg, clp_tmi, pre_multi=True):
        tmp_adv, tmp_clp, stt = self.getAC(side, op, neg, clp_tmi)
        if tmp_adv is None:
            return stt
        return self.insertBest(typs, best, tmp_adv, tmp_clp, lit, side, op, neg, pre_multi)

    def insertBest(self, typs, best, adv, clp, lit, side, op, neg, pre_multi=False):
        if self.constraints.getCstr("inits_productivity") == "extra":
            best.append((adv, clp, [side, op, neg, lit]))
            return self.STG_OK
            
        inserted = False
        i = 0
        while i < len(best):
            if cmp_greater(best[i][0], adv):
                # Best already contains conflicting of better quality
                if self.conflictPair(typs, best[i][-1][-1], lit, pre_multi):
                    return self.STG_CMP
            else:
                if not inserted:
                    best.insert(i, (adv, clp, [side, op, neg, lit]))
                    inserted = True
                # Best contains conflicting of lesser quality than inserted, remove
                elif self.conflictPair(typs, lit, best[i][-1][-1], pre_multi):
                    best.pop(i)
                    i -= 1
            i += 1
        if not inserted:
            best.append((adv, clp, [side, op, neg, lit]))
        return self.STG_OK

    def updateBests(self, typs, besti, pre_multi=False):
        if self.constraints.getCstr("inits_productivity") == "extra":
            besti.sort(reverse=True)
            return besti
        elif self.constraints.getCstr("inits_productivity") == "high":
            return besti
        elif self.constraints.getCstr("inits_productivity") == "low":
            if len(besti) > 1:
                del besti[1:]
            return besti
        # print("BEST BEFORE", " ".join(["%s:%s" % (i, len(best[i])) for i in range(len(best))]))
        ia = 1
        while ia < len(besti):
            ib = ia-1
            while ib >= 0:
                # Best contains conflicting of lesser quality than inserted, remove
                if self.conflictPair(typs, besti[ib][-1][-1], besti[ia][-1][-1], pre_multi):
                    besti.pop(ia)
                    ia -= 1
                    ib = 0
                ib -= 1
            ia += 1
        # print("BEST AFTER", " ".join(["%s:%s" % (i, len(best[i])) for i in range(len(best))]))
        return besti

    def conflictPair(self, typs, litA, litB, pre_multi=False):
        #     r = self.conflictPairX(typs, litA, litB, pre_multi)
        #     print("PAIR", r, "\t", litA, litB, "\t", self.constraints.getCstr("inits_productivity"), typs, pre_multi,self.constraints.getCstr("multi_cats"))
        #     return r
        # def conflictPairX(self, typs, litA, litB, pre_multi=False):
        prod = self.constraints.getCstr("inits_productivity")
        if prod == "extra":
            return False
        if typs == 22:
            if prod == "high" or (pre_multi and self.constraints.getCstr("multi_cats")):
                return (litA[0] == litB[0]) and (litA[1] == litB[1])
            elif prod == "medium":
                return (litA[0] == litB[0]) or (litA[1] == litB[1])
        elif typs == 23:
            if prod == "high" or (pre_multi and self.constraints.getCstr("multi_cats")):
                return (litA[0] == litB[0]) and (litA[1] <= litB[2] and litB[1] <= litA[2])
            elif prod == "medium":
                return (litA[0] == litB[0]) or (litA[1] <= litB[2] and litB[1] <= litA[2])
        elif typs == 33:
            if prod == "high":
                return (litA[0] <= litB[1] and litB[0] <= litA[1]) and (litA[2] <= litB[3] and litB[2] <= litA[3])
            elif prod == "medium":
                return (litA[0] <= litB[1] and litB[0] <= litA[1]) or (litA[2] <= litB[3] and litB[2] <= litA[3])
        return True

    ##################################################
    # CONDITIONAL

    def isCond(self, currentRStatus=0):
        return self.constraints.isStatusCond(currentRStatus)

    def getConditionCand(self, colsC, cond_sparts):
        cis = range(len(colsC))
        prev = None
        best = ([], 0)
        while len(cis) > 0:  # if there are multiple condition attributes, iteratively and greedily try to find a combination
            current = []
            for ci in cis:
                # collect the best candidate for each condition attribute, if fit only returns one
                self.startStoreDivert()
                clp_tm = self.prepareCountsTM(cond_sparts, colsC[ci])
                self.findCover(1, colsC[ci], clp_tm, cond_sparts, Constraints.getStatusCond())
                cands = self.stopStoreDivert()
                if len(cands) == 1:
                    cand = self.store.mkCondCand(*cands[0])
                    if cmp_lower(best[1], cand.getAcc()):
                        best = ([len(current)], cand.getAcc())
                    elif best[1] == cand.getAcc():
                        best[0].append(len(current))
                    current.append(cand)
            if len(best[0]) == 0:
                cis = []
            else:
                basis = (None, None, cond_sparts.nbRows(), 0.)
                for cc in best[0]:  # among the best candidates, select top w.r.t. VarRed
                    cand = current[cc]
                    supp = colsC[cand.getCid()].suppLiteral(cand.getLit())
                    if cand.getVarRed() > basis[-1] or (cand.getVarRed() == basis[-1] and len(supp) < basis[-2]):
                        basis = (cc, supp, len(supp), cand.getVarRed())
                cis = [ci for cii, ci in enumerate(cis) if basis[0] != cii]
                keep_cand, keep_supp = current[basis[0]], basis[1]
                if prev is None:
                    keep_cand.setLit([keep_cand.getLit()])
                else:
                    keep_cand.setLit([keep_cand.getLit()]+prev.getLit())
                prev = keep_cand
                cond_sparts.update(1, False, keep_supp)
                lparts = cond_sparts.lparts()
                best = ([], prev.getAcc())
        return prev

    def getCCandSupp(self, colsC, cond_cand):
        lits = cond_cand.getLit()
        if type(lits) is Literal:
            lits = [lits]
        return set.intersection(*[colsC[lit.colId()].suppLiteral(lit) for lit in lits])

    ##################################################
    # COVER METHODS

    def computeExpand(self, side, col, red, colsC=None, data=None, force_ops=None):
        if isinstance(red, ColM):
            (colL, colR) = (col, red)
            if side == 1:
                (colL, colR) = (col, red)
            return self.computePair(colL, colR, colsC, data)
        elif isinstance(red, Redescription):
            return self.getCandidates(side, col, red, colsC, data, force_ops)
        return []

    def getCandidates(self, side, col, red, colsC=None, data=None, force_ops=None):
        self.store.toNextRound()
        supports = red.supports()
        currentRStatus = Constraints.getStatusRed(red, side, force_ops)
        method_string = 'self.getCandidates%i' % col.typeId()
        try:
            method_compute = eval(method_string)
        except AttributeError:
            raise Exception('Oups No candidates method for this type of data (%i)!' % col.typeId())
        method_compute(side, col, supports, currentRStatus)

        if self.constraints.getCstr("debug_checks"):  # DEBUG
            for c in self.store.currentRoundIter():
                self.checkCountsExt(supports, col, c)
                self.checkRedExt(data, c)

        # compute additional condition
        if colsC is not None and self.constraints.getCstr("add_condition"):
            for c in self.store.currentRoundIter():
                ss = supports.copy()
                supp = col.suppLiteral(c.getLit())
                ss.update(side, c.getOp(), supp)
                cond_sparts = SParts(self.constraints.getSSetts(), ss.nbRows(), [ss.suppI(), ss.suppU()])
                cond_cand = self.getConditionCand(colsC, cond_sparts)
                if cond_cand is not None:
                    c.setCondition(cond_cand)
                    if self.constraints.getCstr("debug_checks"):  # DEBUG
                        # cond_sparts is modified when getting the condition
                        csp = SParts(self.constraints.getSSetts(), ss.nbRows(), [ss.suppI(), ss.suppU()])
                        self.checkCountsCond(colsC, csp, cond_cand)
                        self.checkRedExt(data, c)

        return self.store.currentRoundIter()

    def getCandidatesImprov(self, side, col, red, op, supports, offsets):
        self.store.toNextRound()
        self.setOffsets(offsets)
        currentRStatus = Constraints.getStatusRed(red, side, [op])
        method_string = 'self.getCandidates%i' % col.typeId()
        try:
            method_compute = eval(method_string)
        except AttributeError:
            raise Exception('Oups No candidates method for this type of data (%i)!' % col.typeId())
        method_compute(side, col, supports, currentRStatus)
        self.setOffsets()
        return self.store.currentRoundIter()

    def getCandidates1(self, side, col, supports, currentRStatus=0):
        clp_tm = self.prepareCountsTM(supports, col)
        lin = supports.lpartsInterX(col.supp())

        for op in self.constraints.getCstr("allw_ops", side=side, currentRStatus=currentRStatus):
            for neg in self.constraints.getCstr("allw_negs", side=side, type_id=col.typeId(), currentRStatus=currentRStatus):
                adv, clp, stt = self.getAC(side, op, neg, clp_tm+[lin], self.isCond(currentRStatus))
                if adv is not None:
                    self.addExtToStore(False, adv, clp, (side, op, neg, Literal(neg, BoolTerm(col.getId()))))

    def getCandidates2(self, side, col, supports, currentRStatus=0):
        self.getCandidatesNonBool(side, col, supports, currentRStatus)

    def getCandidates3(self, side, col, supports, currentRStatus=0):
        self.getCandidatesNonBool(side, col, supports, currentRStatus)

    def getCandidatesNonBool(self, side, col, supports, currentRStatus=0):
        clp_tm = self.prepareCountsTM(supports, col)
        self.findCover(side, col, clp_tm, supports, currentRStatus)

    def findCover(self, side, col, clp_tm, supports, currentRStatus=0):
        method_string = 'self.findCover%i' % col.typeId()
        try:
            method_compute = eval(method_string)
        except AttributeError:
            raise Exception('Oups No covering method for this type of data (%i)!' % col.typeId())
        method_compute(side, col, clp_tm, supports, currentRStatus)

    def findCover1(self, side, col, clp_tm, supports, currentRStatus=0):
        lin = supports.lpartsInterX(col.supp())
        if self.constraints.getCstr("neg_query_init", side=side, currentRStatus=currentRStatus):
            nclp_tm = [self.constraints.getSSetts().negateParts(1-side, p) for p in clp_tm]
            nlin = self.constraints.getSSetts().negateParts(1-side, lin)

        for op in self.constraints.getCstr("allw_ops", side=side, currentRStatus=currentRStatus):
            for neg in self.constraints.getCstr("allw_negs", side=side, type_id=col.typeId(), currentRStatus=currentRStatus):
                adv, clp, stt = self.getAC(side, op, neg, clp_tm+[lin], self.isCond(currentRStatus))
                if adv is not None:
                    self.addExtToStore(False, adv, clp, (side, op, neg, Literal(neg, BoolTerm(col.getId()))))

                # to negate the other side when looking for initial pairs
                if self.constraints.getCstr("neg_query_init", side=side, currentRStatus=currentRStatus):
                    adv, clp, stt = self.getAC(side, op, neg, nclp_tm+[nlin], self.isCond(currentRStatus))
                    if adv is not None:
                        self.addExtToStore(True, adv, clp, (side, op, neg, Literal(neg, BoolTerm(col.getId()))))

    def findCover2(self, side, col, clp_tm, supports, currentRStatus=0):
        allw_neg = True in self.constraints.getCstr("allw_negs", side=side, type_id=col.typeId(), currentRStatus=currentRStatus)
        negs = self.constraints.getCstr("allw_negs", side=side, type_id=col.typeId(), currentRStatus=currentRStatus)
        if self.constraints.getCstr("multi_cats"):  # redundant negation
            negs = [False]

        if self.constraints.getCstr("neg_query_init", side=side, currentRStatus=currentRStatus):
            nclp_tm = [self.constraints.getSSetts().negateParts(1-side, p) for p in clp_tm]

        for op in self.constraints.getCstr("allw_ops", side=side, currentRStatus=currentRStatus):
            for neg in negs:
                best = (None, None, None)
                bestNeg = (None, None, None)
                collect_goods = []
                collect_goodsNeg = []
                for (cat, supp) in col.iter_cats():
                    lin = supports.lpartsInterX(supp)
                    tmp_adv, tmp_clp, stt = self.getAC(side, op, neg, clp_tm+[lin], self.isCond(currentRStatus))
                    if tmp_adv is not None:
                        if cmp_lower(best[0], tmp_adv):
                            best = (tmp_adv, tmp_clp, [side, op, neg, Literal(neg, CatTerm(col.getId(), cat))])
                        collect_goods.append((tmp_adv, tmp_clp, cat))

                    # to negate the other side when looking for initial pairs
                    if self.constraints.getCstr("neg_query_init", side=side, currentRStatus=currentRStatus):
                        nlin = self.constraints.getSSetts().negateParts(1-side, lin)
                        tmp_adv, tmp_clp, stt = self.getAC(side, op, neg, nclp_tm+[nlin], self.isCond(currentRStatus))
                        if tmp_adv is not None:
                            if cmp_lower(bestNeg[0], tmp_adv):
                                bestNeg = (tmp_adv, tmp_clp, [side, op, neg, Literal(neg, CatTerm(col.getId(), cat))])
                            collect_goodsNeg.append((tmp_adv, tmp_clp, cat))

                if best[0] is not None:
                    bb = self.combCats(best, allw_neg, side, op, neg, col, collect_goods, currentRStatus=currentRStatus)
                    self.addExtToStore(False, *bb)
                if bestNeg[0] is not None:
                    bb = self.combCats(bestNeg, allw_neg, side, op, neg, col, collect_goodsNeg, currentRStatus=currentRStatus)
                    self.addExtToStore(True, *bb)

    def findCover3(self, side, col, clp_tm, supports, currentRStatus=0):
        counts = None  # so counts can be reused for ops, rather than recomputed
        for op in self.constraints.getCstr("allw_ops", side=side, currentRStatus=currentRStatus):

            if self.constraints.isStatusCond(currentRStatus) or self.inSuppBounds(side, op, clp_tm[0]):  # DOABLE
                cutpoints, counts = col.findCutpoints(self.constraints.getSSetts(), side, supports, op, counts)
                self.findCoverCutpoints(side, op, col, cutpoints, counts, clp_tm, currentRStatus, nego=False)

            # to negate the other side when looking for initial pairs
            if self.constraints.getCstr("neg_query_init", side=side, currentRStatus=currentRStatus):
                nclp_tm = [self.constraints.getSSetts().negateParts(1-side, p) for p in clp_tm]

                if self.inSuppBounds(side, op, clp_tm[0]):  # DOABLE (op is True, OR, for initial pairs)
                    ncutpoints, ncounts = col.findCutpoints(self.constraints.getSSetts(), side, supports.negate(1-side), op)
                    self.findCoverCutpoints(side, op, col, ncutpoints, ncounts, nclp_tm, currentRStatus, nego=True)

    def findCoverCutpoints(self, side, op, col, cutpoints, counts, clp_tm, currentRStatus=0, nego=None):
        if len(cutpoints) < self.constraints.getCstr("max_cutpoints"):
            self.findCoverFullSearch(side, op, col, cutpoints, counts, clp_tm, currentRStatus, nego)
        else:
            if (False in self.constraints.getCstr("allw_negs", side=side, type_id=col.typeId(), currentRStatus=currentRStatus)):
                self.findPositiveCover(side, op, col, cutpoints, counts, clp_tm, currentRStatus, nego)
            if (True in self.constraints.getCstr("allw_negs", side=side, type_id=col.typeId(), currentRStatus=currentRStatus)):
                self.findNegativeCover(side, op, col, cutpoints, counts, clp_tm, currentRStatus, nego)

    def findCoverFullSearch(self, side, op, col, cutpoints, counts, clp_tm, currentRStatus=0, nego=None):
        bests = {False: (None, None, None), True: (None, None, None)}
        for cutp_si, cutp_s in enumerate(cutpoints):
            for cutp_e in cutpoints[cutp_si:]:
                lin = numpy.sum(counts[cutp_s[0]:cutp_e[1]+1], axis=0)
                for neg in self.constraints.getCstr("allw_negs", side=side, type_id=col.typeId(), currentRStatus=currentRStatus):
                    bests[neg], stt = self.updateACT(bests[neg], (cutp_s[0], cutp_e[1]), side, op, neg, clp_tm+[lin], self.isCond(currentRStatus))

        for neg in self.constraints.getCstr("allw_negs", side=side, type_id=col.typeId(), currentRStatus=currentRStatus):
            if bests[neg][0] is not None:
                cutp = bests[neg][-1][-1]
                bests[neg][-1][-1] = col.getLiteralCutpoint(neg, bests[neg][-1][-1])
                if self.constraints.getCstr("debug_checks"):  # DEBUG
                    self.checkCountsCutpoint(counts, cutp, col, bests[neg][-1][-1], bests[neg][1], comp=neg)
                if bests[neg][-1][-1] is not None:
                    self.addExtToStore(nego, *bests[neg])

    def findPositiveCover(self, side, op, col, cutpoints, counts, clp_tm, currentRStatus=0, nego=None):
        is_cond = self.isCond(currentRStatus)

        lin_f = numpy.zeros(counts.shape[1], dtype=int)
        bound_ids_f = []
        best_f = (None, None, None)
        lin_b = numpy.zeros(counts.shape[1], dtype=int)
        bound_ids_b = []
        best_b = (None, None, None)
        for i, cutp_f in enumerate(cutpoints[:-1]):
            # FORWARD
            # cutp_f = cutpoints[i]
            lin_cutpf = numpy.sum(counts[cutp_f[0]:cutp_f[1]+1], axis=0)
            if i > 0 and cmp_lower(self.advAcc(side, op, False, clp_tm+[lin_cutpf], is_cond), self.advRatioVar(side, op, lin_f, is_cond)):
                lin_f += numpy.sum(counts[bound_ids_f[1]+1:cutp_f[1]+1], axis=0)
                bound_ids_f[1] = cutp_f[1]
            else:
                lin_f = lin_cutpf
                bound_ids_f = [cutp_f[0], cutp_f[1]]
            best_f, stt = self.updateACT(best_f, tuple(bound_ids_f), side, op,
                                         False, clp_tm+[list(lin_f)], is_cond)
            # print("FWD", i, best_f[0], best_f[1][0] if best_f[1] is not None else best_f[1], bound_ids_f, lin_f)

            # BACKWARD
            cutp_b = cutpoints[-(i+1)]
            lin_cutpb = numpy.sum(counts[cutp_b[0]:cutp_b[1]+1], axis=0)
            if i > 0 and cmp_lower(self.advAcc(side, op, False, clp_tm+[lin_cutpb], is_cond),
                                   self.advRatioVar(side, op, lin_b, is_cond)):
                lin_b += numpy.sum(counts[cutp_b[0]:bound_ids_b[0]], axis=0)
                bound_ids_b[0] = cutp_b[0]
            else:
                lin_b = lin_cutpb
                bound_ids_b = [cutp_b[0], cutp_b[1]]
            best_b, stt = self.updateACT(best_b, tuple(bound_ids_b), side, op,
                                         False, clp_tm+[list(lin_b)], is_cond)
            # print("BWD", i, best_b[0], best_b[1][0] if best_b[1] is not None else best_b[1], bound_ids_b, lin_b)

        # print("END", best_f, best_b)
        if best_b[0] is not None and best_f[0] is not None:
            # Attempt to intersect the two candidates if they overlap
            bests = [best_b, best_f]
            if cmp_greater(best_b[-1][-1][0], best_f[-1][-1][0]) and \
                    cmp_greater(best_b[-1][-1][1], best_f[-1][-1][1]) and cmp_leq(best_b[-1][-1][0], best_f[-1][-1][1]):
                bounds_m = (best_b[-1][-1][0], best_f[-1][-1][1])
                lin_m = list(numpy.sum(counts[bounds_m[0]:bounds_m[1]+1], axis=0))
                tmp_adv_m, tmp_clp_m, stt = self.getAC(side, op, False, clp_tm+[lin_m], is_cond)
                if tmp_adv_m is not None:
                    bests.append((tmp_adv_m, tmp_clp_m, [side, op, False, bounds_m]))

            bests.sort()
            best = bests[-1]

        elif best_f[0] is not None:
            best = best_f
        else:
            best = best_b

        if best[0] is not None:
            cutp = best[-1][-1]
            best[-1][-1] = col.getLiteralCutpoint(False, best[-1][-1])
            if self.constraints.getCstr("debug_checks"):  # DEBUG
                self.checkCountsCutpoint(counts, cutp, col, best[-1][-1], best[1])
            if best[-1][-1] is not None:
                self.addExtToStore(nego, *best)

    def findNegativeCover(self, side, op, col, cutpoints, counts, clp_tm, currentRStatus=0, nego=None):
        # negation is accounted for in counts in, so setting neg to False
        # should be returned to True for the produced candidates, to indicate they were generated with negation,
        # regardless of whether it was flipped when creating the literal, so duplicates with positive cover can be detected
        is_cond = self.isCond(currentRStatus)

        lin_f = numpy.zeros(counts.shape[1], dtype=int)
        bests_f = [(self.advAcc(side, op, False, clp_tm+[lin_f], is_cond), -1, lin_f)]
        lin_b = numpy.zeros(counts.shape[1], dtype=int)
        bests_b = [(self.advAcc(side, op, False, clp_tm+[lin_b], is_cond), -1, lin_b)]

        for i, cutp_f in enumerate(cutpoints):
            # FORWARD
            # find best cut points for intervals from bottom
            # cutp_f = cutpoints[i]
            if i > 0:
                lin_f += numpy.sum(counts[cutpoints[i-1][0]:cutp_f[0]], axis=0)
            else:
                lin_f = numpy.sum(counts[:cutp_f[0]], axis=0)
            if cmp_greater(self.advRatioVar(side, op, lin_f, is_cond), bests_f[-1][0]):
                lin_f += bests_f[-1][2]
                bests_f.append((self.advAcc(side, op, False, clp_tm+[lin_f], is_cond), cutp_f[0], lin_f))
                lin_f = numpy.zeros(counts.shape[1], dtype=int)

            # BACKWARD
            # find best cut points for intervals from top
            cutp_b = cutpoints[-(i+1)]
            if i > 0:
                lin_b += numpy.sum(counts[cutp_b[1]+1:cutpoints[-i][1]+1], axis=0)
            else:
                lin_b = numpy.sum(counts[cutp_b[1]+1:], axis=0)
            if cmp_greater(self.advRatioVar(side, op, lin_b, is_cond), bests_b[-1][0]):
                lin_b += bests_b[-1][2]
                bests_b.append((self.advAcc(side, op, False, clp_tm+[lin_b], is_cond), cutp_b[1], lin_b))
                lin_b = numpy.zeros(counts.shape[1], dtype=int)

        bests_f.pop(0)
        bests_b.pop(0)
        bests_b.reverse()

        # pair,
        # a forward, i.e. lower bound, cut point should be paired
        # with any backward, i.e. upper bound, cut point that appears before the next forward cut point
        # ## Pairing too close bounds might not leave enough values out!
        ixf, ixb = (-1, 0)
        best_t = (None, None, None)
        stt_b, stt_f, stt_p = (None, None, None)
        while ixf < len(bests_f) and ixb < len(bests_b):
            if ixf+1 == len(bests_f) or bests_b[ixb][1] < bests_f[ixf+1][1]:
                best_t, stt_b = self.updateACT(best_t, (None, bests_b[ixb][1]), side, op,
                                               False, clp_tm+[bests_b[ixb][2]], is_cond)

                if ixf > -1 and not self.is_raised_flag(stt_b, self.FLG_OUT):
                    best_t, stt_p = self.updateACT(best_t, (bests_f[ixf][1], bests_b[ixb][1]), side, op,
                                                   False, clp_tm+[bests_b[ixb][2] + bests_f[ixf][2]], is_cond)
                    # print("<<< Pair %s %s\tf:%s b:%s p:%s\t%s %s" % (bests_f[ixf][1], bests_b[ixb][1],
                    #                                                  self.disp_status(stt_f), self.disp_status(stt_b), self.disp_status(stt_p),
                    #                                                  self.is_raised_flag(stt_b, self.FLG_OUT), self.is_raised_flag(stt_p, self.FLG_OUT)))

                    # if self.is_raised_flag(stt_b, self.FLG_OUT) and not self.is_raised_flag(stt_p, self.FLG_OUT):
                    #     pdb.set_trace()
                    # store_status, store_best = stt_p, best_t
                    # off_f = ixf
                    # while off_f > 0 and self.is_raised_flag(stt_p, self.FLG_OUT) and not self.is_raised_flag(stt_p, self.FLG_IN):
                    #     # if out is too small, but in is not
                    #     # try pairing with previous fwds to increase out
                    #     off_f -= 1
                    #     best_t, stt_p = self.updateACT(best_t, (bests_f[off_f][1], bests_b[ixb][1]), side, op,
                    #                                       False, clp_tm+[bests_b[ixb][2] + bests_f[off_f][2]], is_cond)
                    # print(">> Pair %s %s\tf:%s b:%s p:%s\t%s %s" % (bests_f[off_f][1], bests_b[ixb][1],
                    #                                                 self.disp_status(stt_f), self.disp_status(stt_b), self.disp_status(stt_p),
                    #                                                 self.is_raised_flag(stt_b, self.FLG_OUT), self.is_raised_flag(stt_p, self.FLG_OUT)))

                    # if stt_p == 0 and store_status > 0:
                    #     print("--------------------")
                    #     print(">>>", store_status, store_best)
                    #     print("<<<", best_t)
                ixb += 1
            else:
                ixf += 1
                best_t, stt_f = self.updateACT(best_t, (bests_f[ixf][1], None), side, op,
                                               False, clp_tm+[bests_f[ixf][2]], is_cond)

        if best_t[0] is not None:
            cutp = best_t[-1][-1]
            best_t[-1][-1] = col.getLiteralCutpoint(True, best_t[-1][-1])
            if self.constraints.getCstr("debug_checks"):  # DEBUG
                self.checkCountsCutpoint(counts, cutp, col, best_t[-1][-1], best_t[1], comp=True)
            if best_t[-1][-1] is not None:
                best_t[-1][2] = True
                self.addExtToStore(nego, *best_t)

    ##################################################
    # PAIRS METHODS
    def makePairCand(self, colL, colR, dtlsL={}, dtlsR={}, colsC=None, dtlsC=None, store=True, no_const=False):
        c, candid = (None, None)
        litL, litR = (self.makeLiteral(0, colL, dtlsL), self.makeLiteral(1, colR, dtlsR))
        if litL is not None and litR is not None:
            supports = self.mkSPartsLits(colL, colR, litL)
            clp_tm = self.prepareCountsTM(supports, colR)
            lin = supports.lpartsInterX(colR.suppLiteral(litR))
            adv, clp, stt = self.getAC(1, True, False, clp_tm+[lin], no_const=no_const)
            if adv is not None:
                c = PairCand(self.constraints.getSSetts(), adv, clp, litL, litR)
                if self.constraints.getCstr("debug_checks"):  # DEBUG
                    if data is not None:
                        self.checkRedExt(data, c)        

                if colsC is not None and dtlsC is not None and self.constraints.getCstr("add_condition"):
                    raise Warning("Addition of condition on initial pair in CharboinGreedy makePair not implemented!!")
                    # # prepare additional condition from dtlsC
                    # cond_cand = ... to implement
                    # c.setCondition(cond_cand)
                    # if self.constraints.getCstr("debug_checks"):  # DEBUG
                    #     # cond_sparts is modified when getting the condition
                    #     rsparts = self.mkSPartsLits(colL, colR, c.getLit(0), c.getLit(1))
                    #     csp = SParts(self.constraints.getSSetts(), rsparts.nbRows(), [rsparts.suppI(), rsparts.suppU()])
                    #     self.checkCountsCond(colsC, csp, cond_cand)
                    #     self.checkRedExt(data, c)

                if store and stt == 0:
                    candid = self.store.addPairCand(c)
        return c, candid
        
    def computePair(self, colL, colR, colsC=None, data=None):
        self.store.toNextRound()
        min_type = min(colL.typeId(), colR.typeId())
        max_type = max(colL.typeId(), colR.typeId())
        method_string = 'self.do%i%i' % (min_type, max_type)
        try:
            method_compute = eval(method_string)
        except AttributeError:
            raise Exception('Oups this combination does not exist (%i %i)!' % (min_type, max_type))
        if colL.typeId() == min_type:
            method_compute(colL, colR, 1)
        else:
            method_compute(colL, colR, 0)

        if self.constraints.getCstr("debug_checks"):  # DEBUG
            if data is not None:
                for c in self.store.currentRoundIter():
                    self.checkRedExt(data, c)

        self.store.removeNegDuplicates()
        if colsC is not None and self.constraints.getCstr("add_condition"):
            for c in self.store.currentRoundIter():
                # compute additional condition
                rsparts = self.mkSPartsLits(colL, colR, c.getLit(0), c.getLit(1))
                cond_sparts = SParts(self.constraints.getSSetts(), rsparts.nbRows(), [rsparts.suppI(), rsparts.suppU()])
                cond_cand = self.getConditionCand(colsC, cond_sparts)
                if cond_cand is not None:
                    c.setCondition(cond_cand)
                    if self.constraints.getCstr("debug_checks"):  # DEBUG
                        # cond_sparts is modified when getting the condition
                        csp = SParts(self.constraints.getSSetts(), rsparts.nbRows(), [rsparts.suppI(), rsparts.suppU()])
                        self.checkCountsCond(colsC, csp, cond_cand)
                        self.checkRedExt(data, c)

        return self.store.currentRoundIter()

    def fit(self, col, supports, side, termX):
        clp_tm = self.prepareCountsTM(supports, col)
        currentRStatus = Constraints.getStatusPair(col, side, termX)
        self.store.setCurrentElement((termX, 1-side))
        self.findCover(side, col, clp_tm, supports, currentRStatus=currentRStatus)

    def do11(self, colL, colR, side):
        self.doBoolStar(colL, colR, side)

    def do12(self, colL, colR, side):
        self.doBoolStar(colL, colR, side)

    def do13(self, colL, colR, side):
        self.doBoolStar(colL, colR, side)

    def do22(self, colL, colR, side):
        self.subdo22Full(colL, colR, side)

    def do23(self, colL, colR, side):
        self.subdo23Full(colL, colR, side)

    def do33(self, colL, colR, side):
        if len(colL.interNonMode(colR.nonModeSupp())) >= self.constraints.getCstr("min_itm_in") or \
            (self.constraints.getCstr("support_separate_sides") and colL.lenNonMode() >= self.constraints.getCstr("min_itm_in") and colR.lenNonMode() >= self.constraints.getCstr("min_itm_in")):
            self.subdo33Alts(colL, colR, side)

    ##################################################
    # METHODS COMBINING CATEGORIES
    def matrixLParts(self, cmat, N, csF, csE, last_missing=False):
        if last_missing:
            lparts, lmiss, lin = (self.constraints.getSSetts().makeLParts(), self.constraints.getSSetts().makeLParts(), self.constraints.getSSetts().makeLParts())
            
            lparts[self.constraints.getSSetts().Exo] = cmat[csF, :].sum()
            lparts[self.constraints.getSSetts().Eoo] = N-lparts[self.constraints.getSSetts().Exo]
            
            lin[self.constraints.getSSetts().Exo] = cmat[csF, :][:, csE].sum()
            lin[self.constraints.getSSetts().Eoo] = cmat[:, csE].sum()-lin[self.constraints.getSSetts().Exo]
            
            if len(lparts) >= self.constraints.getSSetts().Emo:
                lparts[self.constraints.getSSetts().Emo] = cmat[-1, :].sum()
                lparts[self.constraints.getSSetts().Eoo] -= lparts[self.constraints.getSSetts().Emo]

                lin[self.constraints.getSSetts().Emo] = cmat[-1, csE].sum()
                lin[self.constraints.getSSetts().Eoo] -= lin[self.constraints.getSSetts().Emo]

                lmiss[self.constraints.getSSetts().Exo] = cmat[csF, -1].sum()
                lmiss[self.constraints.getSSetts().Emo] = cmat[-1, -1]
                lmiss[self.constraints.getSSetts().Eoo] = cmat[:, -1].sum()-lmiss[self.constraints.getSSetts().Exo]-lmiss[self.constraints.getSSetts().Emo]
            return [lparts, lmiss, lin]
        else:
            totF = cmat[csF, :].sum()
            interFE = cmat[csF, :][:, csE].sum()
            return [[totF, 0, 0, N-totF], [interFE, 0, 0, cmat[:, csE].sum() - interFE]]
    
    def additionsLParts(self, prevc, nextc, side=1, neg=False, other_side=False):
        # prevc and nextc are numpy array containing clp in rows lparts,[ lmiss,] lin
        c = numpy.array(nextc, dtype=int)
        if other_side:
            if side == 1:  # Exo(union) = Exo(A) + Exo(B); Eoo(union) = Eoo(A) - Exo(B); Emo(union) = Emo(A)
                p_from, p_to = (self.constraints.getSSetts().Eoo, self.constraints.getSSetts().Exo)
            else:
                p_from, p_to = (self.constraints.getSSetts().Eoo, self.constraints.getSSetts().Eox)
            if neg:
                p_from, p_to = (p_to, p_from)
            transv = prevc[:, p_to]
            c[:, p_from] -= transv
            c[:, p_to] += transv
        else:
            if neg:
                disj, comp = (-2, -1)  # adding disjoint lout, lparts[ and lmiss] are the same, adjusting lin
            else:
                disj, comp = (-1, -2)  # adding disjoint lin, lparts[ and lmiss] are the same, adjusting lout
            c[disj, :] += prevc[disj, :]
            if c.shape[0] == 4:
                c[comp, :] = c[0, :] - c[1, :] - c[disj, :]
            else:
                c[comp, :] = c[0, :] - c[disj, :]  # no lmiss

        return c

    def combCats(self, best, allw_neg, side, op, neg, col, collected, other_side=False, other_neg=False, currentRStatus=0):
        if not self.constraints.getCstr("multi_cats"):
            return best
        collected.sort(key=lambda x: (x[0], x[-1]))
        nextc = collected.pop()
        best_adv = nextc[0]
        cum_counts = numpy.array(nextc[1], dtype=int)
        best_cat = [nextc[-1]]
        # print("---", nextc[-1], "\t",  best_adv)
        while len(collected) > 0:
            nextc = collected.pop()
            # part counts have already been negated earlier where relevant,
            # this is accounted for in addition and neg is set to False in getAC,
            ccum_counts = self.additionsLParts(cum_counts, nextc[1], side, neg, other_side)
            tmp_adv, tmp_clp, stt = self.getAC(side, op, False, ccum_counts, self.isCond(currentRStatus), filled=True)
            # tmp_acc = self.advAcc(side, op, False, ccum_counts, self.isCond(currentRStatus), filled=True)
            # print(best_cat, nextc[-1], "\t",  tmp_adv, self.disp_status(stt))
            if tmp_adv is not None and cmp_lower(best_adv, tmp_adv):
                best_adv = tmp_adv
                cum_counts = ccum_counts
                best_cat.append(nextc[-1])
        if len(best_cat) > 1:  # otherwise best did not change
            # best_adv, best_clp, stt = self.getAC(side, op, False, cum_counts, self.isCond(currentRStatus), filled=True)
            if best_adv is not None:
                if col is None:
                    best = (best_adv, cum_counts, [side, op, neg, set(best_cat)])
                else:
                    lit = col.getLiteralCat(neg, best_cat, allw_neg)
                    if lit is not None:
                        best = (best_adv, cum_counts, [side, op, lit.isNeg(), lit])
        return best
        
    def combCatsPair(self, besti, colF, colE, nF, nE, allw_neg, cmat=None):
        if self.constraints.getCstr("exhaustive_comb_cats") and not allw_neg and cmat is not None:
            self.combCatsPairExhaustive(besti, colF, colE, cmat)
            return

        map_cat_rs = [[{}, {}], [{}, {}]]
        # collect together all candidates that have the same category
        for b in besti:
            for ss in [0, 1]:
                cat = b[-1][-1][ss]
                if cat not in map_cat_rs[0][ss]:
                    map_cat_rs[0][ss][cat] = []
                map_cat_rs[0][ss][cat].append((b[0], b[1], b[-1][-1][1-ss]))
                # adv, clp, cat

        # for ssx in [0,1]:
        #     print("SIDE", ssx)
        #     for z in map_cat_rs[0][ssx].keys():
        #         print("cat", z)
        #         for x in map_cat_rs[0][ssx][z]:
        #             print(x[-1], x[0][:2])
        # ss is the common side
        for ri in range(len(map_cat_rs)):
            for ss, nS, nX in [(0, nF, nE), (1, nE, nF)]:
                cats_loc = [None, None]
                kk = map_cat_rs[ri][ss].keys()
                for k in kk:
                    if len(map_cat_rs[ri][ss][k]) > 1:
                        bb = self.combCats(None, allw_neg, 1, True, nX, None, map_cat_rs[ri][ss][k], other_side=(ss == 1), other_neg=nS)
                        if bb is not None:
                            cats_loc = {ss: set(k) if type(k) is tuple else k, 1-ss: bb[-1][-1]}
                            if self.constraints.getCstr("debug_checks"):  # DEBUG
                                tF = colF.getLiteralCat(nF, cats_loc[0], allw_neg)
                                tE = colE.getLiteralCat(nE, cats_loc[1], allw_neg)
                                self.checkCountsPair(1, False, colF, colE, tF, tE, bb[1])
                            stt = self.insertBest(22, besti, bb[0], bb[1], (cats_loc[0], cats_loc[1]), 1, True, False, pre_multi=False)

                            # ADD TO map_cat of other side for further combination
                            if ri+1 < len(map_cat_rs):  # and stt == 0?
                                cat = tuple(sorted(bb[-1][-1]))
                                if cat not in map_cat_rs[ri+1][1-ss]:
                                    map_cat_rs[ri+1][1-ss][cat] = []
                                map_cat_rs[ri+1][1-ss][cat].append((bb[0], bb[1], k))

        # for ssx in [0,1]:
        #     print("SIDE", ssx)
        #     for z in map_cat_rs[1][ssx].keys():
        #         print("cat", z)
        #         for x in map_cat_rs[1][ssx][z]:
        #             print(x[-1], x[0][:2])

        
    def combCatsPairExhaustive(self, besti, colF, colE, cmat):
        for (ncF, ncE) in itertools.product(range(1, colF.nbCats()), range(1, colE.nbCats())):
            if ncF > 1 or ncE > 1:
                for csF in itertools.combinations(range(colF.nbCats()), ncF):                                        
                    for csE in itertools.combinations(range(colF.nbCats()), ncE):
                        ccum_counts = self.matrixLParts(cmat, colF.nbRows(), csF, csE, last_missing=self.handlesMiss())
                        tmp_adv, tmp_clp, stt = self.getAC(1, True, False, ccum_counts)
                        if tmp_adv is not None:
                            if self.constraints.getCstr("debug_checks"):  # DEBUG
                                tF = colF.getLiteralCat(False, csF, False)
                                tE = colE.getLiteralCat(False, csE, False)
                                self.checkCountsPair(1, False, colF, colE, tF, tE, tmp_clp)
                            stt = self.insertBest(22, besti, tmp_adv, tmp_clp, (csF, csE), 1, True, False, pre_multi=False)

    def combCatsNum(self, besti, colF, colE, nF, nE, allw_neg, side, buckets, bUp):
        map_cat = {}
        # collect together all candidates that have the same numerical interval
        for b in besti:
            buk = b[-1][-1][1:]
            if buk not in map_cat:
                map_cat[buk] = []
            map_cat[buk].append((b[0], b[1], b[-1][-1][0]))
            # adv, clp, cat

        kk = map_cat.keys()
        for k in kk:
            if len(map_cat[k]) > 1:
                bb = self.combCats(None, allw_neg, side, True, nF, None, map_cat[k], other_side=True, other_neg=nE)
                if bb is not None:
                    if self.constraints.getCstr("debug_checks"):  # DEBUG
                        tE = colE.getLiteralBuk(nE, buckets[1], k, buckets[bUp])
                        tF = colF.getLiteralCat(nF, bb[-1][-1], allw_neg)
                        self.checkCountsPair(side, False, colF, colE, tF, tE, bb[1])
                    stt = self.insertBest(23, besti, bb[0], bb[1], (bb[-1][-1], )+k, 1, True, False, pre_multi=False)

    ##################################################
    # TOOLS METHODS

    ##################################################
    # DOUBLE CHECKS FUNCTIONS

    def resultCheck(self, rv, msg):
        if not rv:
            print(msg)
            raise ExtensionWarning(msg)
        return rv

    def checkCountsCutpoint(self, counts, cutp, col, lit, clp, comp=False):
        if cutp[0] is None:
            cc = numpy.sum(counts[:cutp[1]+1], axis=0)
        elif cutp[1] is None:
            cc = numpy.sum(counts[cutp[0]:], axis=0)
        else:
            cc = numpy.sum(counts[cutp[0]:cutp[1]+1], axis=0)
        new_sums = str(list(cc))
        if comp:
            org_c = clp[-2]
        else:
            org_c = clp[-1]
        org_sums = str(list(org_c))
        org_in = sum(clp[-1])
        if lit is None:
            if org_in > 0:  # can't guess whether it's negated or not...
                new_in = col.nbRows() - col.lMiss()
            else:
                new_in = 0
        else:
            new_in = len(col.suppLiteral(lit))

        msg = "--- checkCountsCutpoint %s" % lit
        msg += "\nsums %s %s%s" % (org_sums, new_sums, (org_sums != new_sums)*" !!!")
        msg += "\nin %s %s%s" % (org_in, new_in, (org_in != new_in)*" !!!")
        rv = (org_sums == new_sums) and (org_in == new_in)
        return self.resultCheck(rv, msg)

    def checkCountsPair(self, side, neg, colF, colE, tF, tE, clp):
        part_ids = self.constraints.getSSetts().getInitPartIds(side)
        org_ltot = str([clp[0][pid] for pid in part_ids])
        if neg:
            org_lin = str([clp[-2][pid] for pid in part_ids])
        else:
            org_lin = str([clp[-1][pid] for pid in part_ids])

        if tE is None:
            if sum(clp[-1][:3]) > 0:  # can't guess whether it's negated or not...
                suppE = colE.rows() - colE.miss()
            else:
                suppE = set()
        else:
            suppE = colE.suppLiteral(tE)
        if tF is None:
            if sum(clp[-1][:3]) > 0:  # can't guess whether it's negated or not...
                suppF = colF.rows() - colF.miss()
            else:
                suppF = set()
        else:
            suppF = colF.suppLiteral(tF)
        if len(part_ids) == 2:
            new_ltot = str([len(suppF), colF.N - len(suppF)])
            new_lin = str([len(suppE & suppF), len(suppE) - len(suppE & suppF)])
            new_lmiss = str([0, 0])
            ok_miss = True
            msg_miss = ""
        else:
            missE = colE.miss()
            missF = colF.miss()
            new_ltot = str([len(suppF), colF.N - len(suppF) - len(missF), len(missF)])
            new_lin = str([len(suppE & suppF), len(suppE) - len(suppE & missF) - len(suppE & suppF), len(suppE & missF)])

            if len(clp) == 4:
                org_lmiss = str([clp[1][pid] for pid in part_ids])
                new_lmiss = str([len(missE & suppF), len(missE) - len(missE & missF) - len(missE & suppF), len(missE & missF)])
                ok_miss = (org_lmiss == new_lmiss)
                msg_miss = "\nlmiss %s %s%s" % (org_lmiss, new_lmiss, (org_lmiss != new_lmiss)*" !!!")
            else:
                ok_miss = True
                msg_miss = ""

        rv = (org_ltot == new_ltot) and (org_lin == new_lin) and ok_miss
        msg = "--- checkCountsPair %s %s %s" % (side, tF, tE)
        msg += "\nltot %s %s%s" % (org_ltot, new_ltot, (org_ltot != new_ltot)*" !!!")
        msg += "\nlin %s %s%s" % (org_lin, new_lin, (org_lin != new_lin)*" !!!")
        msg += msg_miss
        return self.resultCheck(rv, msg)

    def checkCountsCond(self, colsC, cond_sparts, c):
        supp = self.getCCandSupp(colsC, c)
        lparts = cond_sparts.lparts()
        lin = cond_sparts.lpartsInterX(supp)
        new_ltot = str(lparts)
        new_lin = str(lin)

        clp = c.getClp()
        org_ltot = str([cc for cc in clp[0]])
        org_lin = str([cc for cc in clp[-1]])

        rv = (org_ltot == new_ltot) and (org_lin == new_lin)
        msg = "--- checkCountsCond %s" % " and ".join([str(l) for l in c.getLit()])
        msg += "\nltot %s %s%s" % (org_ltot, new_ltot, (org_ltot != new_ltot)*" !!!")
        msg += "\nlin %s %s%s" % (org_lin, new_lin, (org_lin != new_lin)*" !!!")
        return self.resultCheck(rv, msg)

    def checkCountsExt(self, supports, col, c):
        supp = col.suppLiteral(c.getLit())
        lparts = supports.lparts()
        lin = supports.lpartsInterX(supp)

        if len(c.getClp()) == 4:
            lmiss = supports.lpartsInterX(set(col.miss()))
            lout = [lparts[i]-lmiss[i]-lin[i] for i in range(len(lparts))]

            org_clp = str([list(c.getClp()[0]), list(c.getClp()[1]), list(c.getClp()[-1])])
            new_clp = str([lparts, lmiss, lin])
        else:
            lout = [lparts[i]-lin[i] for i in range(len(lparts))]

            org_clp = str([list(c.getClp()[0]), list(c.getClp()[-1])])
            new_clp = str([lparts, lin])

        rv = org_clp == new_clp
        msg = "--- checkCountsExt %s" % self.store.dispCand(c)
        # c.disp(self.store.ssetts, self.store.N,
        #            self.store.c_vals["base_acc"], self.store.c_vals["prs"], self.store.constraints["score_coeffs"])
        msg += "\n<<< clp: %s" % org_clp
        msg += "\n>>> clp: %s" % new_clp
        return self.resultCheck(rv, msg)

    def checkRedExt(self, data, c):
        if data is not None:
            red = self.store.mkRedFromCand(c, data)
            chk = c.checkRed(red)
            rv = chk is None
            msg = "--- checkRedExt %s" % self.store.dispCand(c)
            msg += "\n%s" % red
            msg += "\n%s" % chk
            return self.resultCheck(rv, msg)
