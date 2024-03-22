import numpy
try:
    from classCol import NumColM
    from classData import Data
    from classConstraints import Constraints
    from classCharbonGreedy import CharbonGreedy
    from classSParts import SParts
    from classQuery import *
except ModuleNotFoundError:
    from .classCol import NumColM
    from .classData import Data
    from .classConstraints import Constraints
    from .classCharbonGreedy import CharbonGreedy
    from .classSParts import SParts
    from .classQuery import *

import pdb


class CharbonGMiss(CharbonGreedy):

    name = "GreedyMiss"

    def handlesMiss(self):
        return True

    ##################################################
    # TOOLS METHODS
    def prepareCountsTM(self, supports, col):
        lparts = supports.lparts()
        lmiss = supports.lpartsInterX(col.miss())
        return [lparts, lmiss]

    def mkSPartsLits(self, colL, colR, litL=None, litR=None):
        return SParts(self.constraints.getSSetts(), colL.nbRows(),
                      [colL.suppLiteral(litL) if litL is not None else set(), colR.suppLiteral(litR) if litR is not None else set(),
                       colL.miss() if litL is not None else set(), colR.miss() if litR is not None else set()])

    @classmethod
    def fillClp(tcl, clp_tmi, neg=False):  # ads lout and flips with lin if negated
        lparts, lmiss, lin = clp_tmi
        lout = [lparts[i] - lmiss[i] - lin[i] for i in range(len(lparts))]
        if neg:
            return (lparts, lmiss, lin, lout)
        else:
            return (lparts, lmiss, lout, lin)

    def inSuppBounds(self, side, op, lparts):
        ssetts = self.constraints.getSSetts()
        return ssetts.sumPartsId(side, ssetts.IDS_varnum[op] + ssetts.IDS_fixnum[op], lparts) >= self.constraints.getCstr("min_itm_in") and \
            ssetts.sumPartsId(side, ssetts.IDS_cont[op], lparts) >= self.constraints.getCstr("min_itm_c") and \
            ssetts.sumPartsId(side, ssetts.IDS_out[op], lparts) >= self.constraints.getCstr("min_itm_out")

    def advRatioVar(self, side, op, lin_f, is_cond=False):
        ssetts = self.constraints.getSSetts()
        if is_cond:
            num = ssetts.sumPartsId(side, ssetts.IDS_varnum[False], lin_f)
            den = ssetts.sumPartsId(side, ssetts.IDS_varden[False], lin_f)
            return self.ratio(num, den+num)
        num = ssetts.sumPartsId(side, ssetts.IDS_varnum[op], lin_f)
        den = ssetts.sumPartsId(side, ssetts.IDS_varden[op], lin_f)
        return self.ratio(num, den)

    def advAcc(self, side, op, neg, clp_tmi, is_cond=False, filled=False):
        ssetts = self.constraints.getSSetts()
        clp = clp_tmi if filled else self.fillClp(clp_tmi, neg)

        if is_cond:
            num = ssetts.sumPartsIdInOut(side, neg, ssetts.IDS_varnum[False], clp)
            den = num + ssetts.sumPartsIdInOut(side, neg, ssetts.IDS_varden[False], clp)
        else:
            num = ssetts.sumPartsIdInOut(side, neg, ssetts.IDS_varnum[op] + ssetts.IDS_fixnum[op], clp)
            den = ssetts.sumPartsIdInOut(side, neg, ssetts.IDS_varden[op] + ssetts.IDS_fixden[op], clp)
        return self.offset_ratio(num, den)

    def getAC(self, side, op, neg, clp_tmi, is_cond=False, no_const=False, filled=False):
        ssetts = self.constraints.getSSetts()
        clp = clp_tmi if filled else self.fillClp(clp_tmi, neg)

        if is_cond:
            cont = ssetts.sumPartsIdInOut(side, neg, ssetts.IDS_varnum[False], clp)
            if self.unconstrained(no_const) or (cont >= self.constraints.getCstr("min_itm_in")):
                fix_num, fix_den = (0, 0)
                var_num = ssetts.sumPartsIdInOut(side, neg, ssetts.IDS_varnum[False], clp)
                var_den = var_num + ssetts.sumPartsIdInOut(side, neg, ssetts.IDS_varden[False], clp)
                acc = self.offset_ratio(var_num, var_den)
                return (acc, var_num, var_den, cont, fix_num, fix_den), clp, 0
            return None, None, self.FLG_IN | self.FLG_CONT

        cstr_status = 0
        var_num = ssetts.sumPartsIdInOut(side, False, ssetts.IDS_varnum[op], clp)
        fix_num = ssetts.sumPartsIdInOut(side, False, ssetts.IDS_fixnum[op], clp)
        if var_num+fix_num < self.constraints.getCstr("min_itm_in"):
            cstr_status |= self.FLG_IN

        sout = ssetts.sumPartsIdInOut(side, False, ssetts.IDS_out[op], clp)
        if sout < self.constraints.getCstr("min_itm_out"):
            cstr_status |= self.FLG_OUT

        cont = ssetts.sumPartsIdInOut(side, False, ssetts.IDS_cont[op], clp)
        if cont < self.constraints.getCstr("min_itm_c"):
            cstr_status |= self.FLG_CONT
        # print(cstr_status, var_num+fix_num, sout, cont, clp_tmi[0], clp_tmi[-1])
        if self.unconstrained(no_const) or self.constraints.getCstr("support_separate_sides") or cstr_status == 0:
            var_den = ssetts.sumPartsIdInOut(side, False, ssetts.IDS_varden[op], clp)
            fix_den = ssetts.sumPartsIdInOut(side, False, ssetts.IDS_fixden[op], clp)
            acc = self.offset_ratio(var_num + fix_num, var_den + fix_den)
            # print("PIECES", sout, var_num, var_den, cont, fix_num, fix_den)
            return (acc, var_num, var_den, cont, fix_num, fix_den), clp, cstr_status
        return None, None, cstr_status

    ##################################################
    # SEARCH METHODS
    def doBoolStar(self, colL, colR, side):
        if side == 1:
            (supports, fixTerm, extCol) = (SParts(self.constraints.getSSetts(), colL.nbRows(),
                                                  [colL.supp(), set(), colL.miss(), set()]), BoolTerm(colL.getId()), colR)
        else:
            (supports, fixTerm, extCol) = (SParts(self.constraints.getSSetts(), colL.nbRows(),
                                                  [set(), colR.supp(), set(), colR.miss()]), BoolTerm(colR.getId()), colL)
        self.fit(extCol, supports, side, fixTerm)

    def subdo33Alts(self, colL, colR, side):
        # print("###### subdo33Alts",  colL.getUid(), colR.getUid())
        org_side = side
        best = []
        interMat = []
        max_prodbuckets = self.constraints.getCstr("max_prodbuckets")
        if self.constraints.getCstr("method_buckets") == "naive":
            bucketsL = colL.buckets()
            bucketsR = colR.buckets()
        elif self.constraints.getCstr("method_buckets") == "tails":
            tails_params = {"lower_tail_agg": self.constraints.getCstr("lower_tail_agg"),
                            "upper_tail_agg": self.constraints.getCstr("upper_tail_agg")}
            bucketsL = colL.buckets("tails", tails_params)
            bucketsR = colR.buckets("tails", tails_params)
        else:
            max_prodbuckets = 0
            bparams = {"nb_buckets": self.constraints.getCstr("nb_buckets")}
            if self.constraints.getCstr("method_buckets") == "bellman":
                bparams["bellman_criterion"] = self.constraints.getCstr("bellman_criterion")
            bucketsL = colL.buckets(self.constraints.getCstr("method_buckets"), bparams)
            bucketsR = colR.buckets(self.constraints.getCstr("method_buckets"), bparams)
        if len(bucketsL[0]) > len(bucketsR[0]):
            bucketsF = bucketsR
            colF = colR
            bucketsE = bucketsL
            colE = colL
            side = 1-side
        else:
            bucketsF = bucketsL
            colF = colL
            bucketsE = bucketsR
            colE = colR

        # DOABLE?
        nbb = max_prodbuckets / float(len(bucketsF[1]))
        # print("--- Nb buckets: %i x %i, max buckets = %s, max agg = %s" %
        #       (len(bucketsF[1]), len(bucketsE[1]), max_prodbuckets, self.constraints.getCstr("max_agg")))
        # print("--- nbb=%s\tnb B E/nb=%s" % (nbb, len(bucketsE[1])/nbb))
        if max_prodbuckets > 0 and len(bucketsE[1]) > nbb:
            ### TOO MANY BUCKETS !

            if len(bucketsE[1])/nbb < self.constraints.getCstr("max_agg"):
                # collapsing buckets on the largest side is enough to get within the reasonable size
                # print("subdo33Alts A --- Collapsing just E")
                bucketsE = colE.buckets("collapsed", {"max_agg": self.constraints.getCstr("max_agg"), "nbb": nbb})

            else:
                # collapsing buckets on the largest side is NOT enough to get within the reasonable size
                # print("subdo33Alts B --- As categories?")
                bucketsE = None

                # try cats
                exclL = NumColM.buk_excl_bi(bucketsL)
                exclR = NumColM.buk_excl_bi(bucketsR)
                bbs = [dict([("c:%d" % bi, es) for (bi, es) in enumerate(bucketsL[0])
                             if (len(es) > self.constraints.getCstrSimple("min_itm_in") and
                                 colL.nbRows() - len(es) > self.constraints.getCstrSimple("min_itm_out") and (exclL is None or bi != exclL))]),
                       dict([("c:%d" % bi, es) for (bi, es) in enumerate(bucketsR[0])
                             if (len(es) > self.constraints.getCstrSimple("min_itm_in") and
                                 colR.nbRows() - len(es) > self.constraints.getCstrSimple("min_itm_out") and (exclR is None or bi != exclR))])]

                # if len(bbs[0]) > 0 and ( len(bbs[1]) == 0 or len(bbs[0])/float(len(bucketsL[0])) < len(bbs[1])/float(len(bucketsR[0]))):
                nbes = [max(sum([len(v) for (k, v) in bbs[s].items()]), .5) for s in [0, 1]]
                sideN = None
                # Decide which side to make categorical, based on the number and cardinalities of the categories that would result
                # will get through in subdo23? ((len(buckets[1]) * len(colF.cats()) <= max_prodbuckets))
                if len(bbs[0]) > 0 and (len(bbs[1]) == 0 or nbes[0]/len(bbs[0]) > nbes[1]/len(bbs[1])):
                    sideN, ccN, ccC, ccCN, bucksC = (1, colR, Data.getColClassForName("Categorical")(bbs[0], colL.nbRows(), colL.miss()), colL, bucketsL)
                elif len(bbs[1]) > 0:
                    sideN, ccN, ccC, ccCN, bucksC = (0, colL, Data.getColClassForName("Categorical")(bbs[1], colR.nbRows(), colR.miss()), colR, bucketsR)

                if sideN is not None:
                    # working with variable as categories is doable
                    # print("Trying cats...", len(bucketsL[0]), len(bucketsR[0]), len(bbs[0]), len(bbs[1]))
                    self.startStoreDivert()
                    self.subdo23Full(ccC, ccN, 1, try_comb=False)
                    cands = self.stopStoreDivert()
                    for cand in cands:
                        # the column that was turned categorical
                        ltc = cand[2]
                        c = ltc.getTerm().getCat()
                        if type(c) is set and len(c) > 0:
                            c = sorted(c)[0]
                        c = int(c[2:])
                        valLow = bucksC[1][c]
                        valUp = bucksC[NumColM.buk_ind_maxes(bucksC)][c]
                        cand = list(cand)
                        cand[2] = Literal(ltc.isNeg(), NumTerm(ccCN.getId(), valLow, valUp))
                        cand[-1] = sideN
                        if self.constraints.getCstr("debug_checks"):  # DEBUG
                            # the column that was kept numerical
                            self.checkCountsPair(1, False, ccCN, ccN, cand[2], cand[3], cand[1])
                        self.addPairToStore(*cand)
                    return
                else:
                    # working with variable as categories is NOT doable
                    # the only remaining solution is aggressive collapse of buckets on both sides
                    # print("subdo33Alts C --- Collapsing both")
                    nbb = numpy.sqrt(max_prodbuckets)
                    bucketsE = colE.buckets("collapsed", {"max_agg": self.constraints.getCstr("max_agg"), "nbb": nbb})
                    bucketsF = colF.buckets("collapsed", {"max_agg": self.constraints.getCstr("max_agg"), "nbb": nbb})
                    # print("Last resort solution... Collapsing both E and F", nbb, len(bucketsL[0]), len(bucketsR[0]))

        # print("buckets lengths\t(0,%d) %d\t(1,%d) %d\tcollapsed %d -- product %d" % (colL.getId(), len(bucketsL[1]), colR.getId(), len(bucketsR[1]), len(bucketsE[1]), len(bucketsF[1]) * len(bucketsE[1])))
        if bucketsE is not None and ( max_prodbuckets == 0 or (len(bucketsF[1]) * len(bucketsE[1])) < max_prodbuckets ):
            # print("subdo33Alts D --- Trying buckets...", len(bucketsF[0]), len(bucketsE[0]))
            partsMubB = colF.lMiss()
            missMubB = len(colF.miss() & colE.miss())
            totInt = colE.nbRows() - colF.lMiss() - colE.lMiss() + missMubB
            # margE = [len(intE) for intE in bucketsE[0]]

            lmissFinE = [len(colF.miss() & bukE) for bukE in bucketsE[0]]
            lmissEinF = [len(colE.miss() & bukF) for bukF in bucketsF[0]]
            margF = [len(bucketsF[0][i]) - lmissEinF[i] for i in range(len(bucketsF[0]))]

            for bukF in bucketsF[0]:
                interMat.append([len(bukF & bukE) for bukE in bucketsE[0]])

            if bucketsF[2] is not None:
                lFmodeEmiss = len(colF.interMode(colE.miss()))
                margF[bucketsF[2]] += colF.lenMode() - lFmodeEmiss
                lmissEinF[bucketsF[2]] += lFmodeEmiss
                for bukEId in range(len(bucketsE[0])):
                    interMat[bucketsF[2]][bukEId] += len(colF.interMode(bucketsE[0][bukEId]))

            if bucketsE[2] is not None:
                # margE[bucketsE[2]] += colE.lenMode()
                lmissFinE[bucketsE[2]] += len(colE.interMode(colF.miss()))
                for bukFId in range(len(bucketsF[0])):
                    interMat[bukFId][bucketsE[2]] += len(colE.interMode(bucketsF[0][bukFId]))

            if bucketsF[2] is not None and bucketsE[2] is not None:
                interMat[bucketsF[2]][bucketsE[2]] += len(colE.interMode(colF.modeSupp()))

            totMissE = len(colE.miss())
            totMissEinF = sum(lmissEinF)

            # ### check marginals
            # totF = 0
            # for iF in range(len(bucketsF[0])):
            #     sF = sum(interMat[iF])
            #     if sF != margF[iF]:
            #         raise Error('Error in computing the marginals (1)')
            #     totF += sF

            # totE = 0
            # for iE in range(len(bucketsE[0])):
            #     sE = sum([intF[iE] for intF in interMat])
            #     if sE != margE[iE]:
            #         raise Error('Error in computing the marginals (2)')
            #     totE += sE

            # if totE != totF or totE != colE.nbRows():
            #     raise Error('Error in computing the marginals (3)')

            exclF = NumColM.buk_excl_bi(bucketsF)
            exclE = NumColM.buk_excl_bi(bucketsE)
            belowF = 0
            lowF = 0
            while lowF < len(interMat) and totInt - belowF >= self.constraints.getCstrSimple("min_itm_in"):
                aboveF = 0
                upF = len(interMat)-1
                if self.constraints.getCstrSimple("max_combined_buckets") > 0 and (lowF + self.constraints.getCstrSimple("max_combined_buckets")-1) < upF:
                    upF = lowF + self.constraints.getCstrSimple("max_combined_buckets") - 1
                    aboveF = numpy.sum(margF[upF+1:], dtype=int)

                if exclF is not None:
                    if lowF == exclF:  # basically, skip this value
                        upF = lowF-1
                    elif lowF < exclF:
                        upF = exclF-1
                        aboveF = numpy.sum(margF[upF+1:], dtype=int)

                while upF >= lowF and totInt - belowF - aboveF >= self.constraints.getCstrSimple("min_itm_in"):
                    if belowF + aboveF >= self.constraints.getCstrSimple("min_itm_out"):
                        EinF = [sum([interMat[iF][iE] for iF in range(lowF, upF+1)]) for iE in range(len(interMat[lowF]))]
                        EoutF = [sum([interMat[iF][iE] for iF in list(range(0, lowF))+list(range(upF+1, len(interMat)))]) for iE in range(len(interMat[lowF]))]
                        lmissE = sum(lmissEinF[lowF:upF+1])
                        # totEinF = sum(EinF)

                        clp_tm = [self.constraints.getSSetts().makeLParts([totInt - aboveF - belowF + lmissE,
                                                                           aboveF + belowF + totMissEinF - lmissE,
                                                                           partsMubB], side=side),
                                  self.constraints.getSSetts().makeLParts([lmissE,
                                                                           totMissEinF - lmissE,
                                                                           missMubB], side=side)]
                        belowEF = 0
                        outBelowEF = 0
                        lowE = 0
                        while lowE < len(interMat[lowF]) and ( totInt - belowF - aboveF - belowEF >= self.constraints.getCstrSimple("min_itm_in") or \
                                                               ( self.constraints.getCstr("support_separate_sides") and totInt - outBelowEF - belowEF >= self.constraints.getCstrSimple("min_itm_in") )):
                            aboveEF = 0
                            outAboveEF = 0
                            upE = len(interMat[lowF])-1

                            if self.constraints.getCstrSimple("max_combined_buckets") > 0 and (lowE + self.constraints.getCstrSimple("max_combined_buckets")-1) < upE:
                                upE = lowE + self.constraints.getCstrSimple("max_combined_buckets") - 1
                                aboveEF = numpy.sum(EinF[upE+1:], dtype=int)
                                outAboveEF = numpy.sum(EoutF[upE+1:], dtype=int)


                            if exclE is not None:
                                if lowE == exclE:  # basically, skip this value
                                    upE = lowE-1
                                elif lowE < exclE:
                                    upE = exclE-1
                                    aboveEF = numpy.sum(EinF[upE+1:], dtype=int)
                                    outAboveEF = numpy.sum(EoutF[upE+1:], dtype=int)

                            while upE >= lowE and ( totInt - belowF - aboveF - belowEF - aboveEF >= self.constraints.getCstrSimple("min_itm_in") or \
                                                     ( self.constraints.getCstr("support_separate_sides") and totInt - outBelowEF - belowEF - outAboveEF - aboveEF >= self.constraints.getCstrSimple("min_itm_in") )):
                                if not self.constraints.getCstr("support_separate_sides") or (outBelowEF + belowEF + outAboveEF + aboveEF >= self.constraints.getCstrSimple("min_itm_out")):
                                    lmissF = sum(lmissFinE[lowE:upE+1])
                                    lin = self.constraints.getSSetts().makeLParts([totInt - belowF - aboveF - belowEF - aboveEF,
                                                                               belowF + aboveF - outAboveEF - outBelowEF,
                                                                               lmissF], side=side)

                                    stt = self.updateACTList(33, best, (lowF, upF, lowE, upE), side, True, False, clp_tm+[lin])
                                aboveEF += EinF[upE]
                                outAboveEF += EoutF[upE]
                                upE -= 1
                            belowEF += EinF[lowE]
                            outBelowEF += EoutF[lowE]
                            lowE += 1
                    aboveF += margF[upF]
                    upF -= 1
                belowF += margF[lowF]
                lowF += 1

        bUpE = NumColM.buk_ind_maxes(bucketsE)  # in case of collapsed bucket the threshold is different
        bUpF = NumColM.buk_ind_maxes(bucketsF)  # in case of collapsed bucket the threshold is different
        self.updateBests(33, best)
        for b in best:
            tF = colF.getLiteralBuk(False, bucketsF[1], b[-1][-1][0:2], bucketsF[bUpF])
            tE = colE.getLiteralBuk(False, bucketsE[1], b[-1][-1][2:], bucketsE[bUpE])
            if self.constraints.getCstr("debug_checks"):  # DEBUG
                self.checkCountsPair(side, False, colF, colE, tF, tE, b[1])
            if tF is not None and tE is not None:
                self.addPairToStore(b[0], b[1], tF, tE, False, False, side)

    def subdo22Full(self, colL, colR, side):
        configs = [(False, False), (False, True), (True, False), (True, True)]

        allw_neg = True
        if True not in self.constraints.getCstr("neg_query", side=side, type_id=2):
            configs = configs[:1]
            allw_neg = False
        best = [[] for c in configs]

        # print("--------------------------------------")
        # print("\t".join(["", ""]+[catR for catR in colR.cats()]))
        # print("\t".join(["", ""]+["%d" % colR.lsuppCat(catR) for catR in colR.cats()]))
        # for catL in colL.cats():
        #     print("\t".join([catL, "%d" % colL.lsuppCat(catL)]+["%d" % len(colL.suppCat(catL).intersection(colR.suppCat(catR))) for catR in colR.cats()]))
        # print("--------------------------------------")

        cmat = numpy.zeros((colL.nbCats()+1, colR.nbCats()+1), dtype=int)
        for ciL, catL in enumerate(colL.cats()):
            supports = SParts(self.constraints.getSSetts(), colL.nbRows(),
                              [colL.suppCat(catL), set(), colL.miss(), set()])
            clp_tm = self.prepareCountsTM(supports, colR)
            cmat[ciL, -1] = -clp_tm[0][0]
            for ciR, catR in enumerate(colR.cats()):
                lin = supports.lpartsInterX(colR.suppCat(catR))
                cmat[ciL, ciR] = lin[self.constraints.getSSetts().Exo]
                if ciL > 0 and len(lin) >= self.constraints.getSSetts().Emo:
                    cmat[-1, ciR] = lin[self.constraints.getSSetts().Emo]
                # self.constraints.getSSetts().Exo
                for i, (nL, nR) in enumerate(configs):
                    if nL:
                        tmp_clp_tm = [self.constraints.getSSetts().negateParts(0, p) for p in clp_tm]
                        tmp_lin = self.constraints.getSSetts().negateParts(0, lin)
                    else:
                        tmp_clp_tm = clp_tm
                        tmp_lin = lin

                    stt = self.updateACTList(22, best[i], (catL, catR), side, True, nR, tmp_clp_tm+[tmp_lin], pre_multi=self.constraints.getCstr("multi_cats"))

        # compute the actual counts for the last column
        # (size of intersection of missing values for colR with categories of colL)
        cmat[:-1,-1] = -(cmat[:-1,-1]+cmat[:-1,:-1].sum(axis=1))
        cmat[-1,-1] = colL.nbRows() - cmat.sum()
                    
        for i, (nL, nR) in enumerate(configs):
            if self.constraints.getCstr("multi_cats"):
                self.combCatsPair(best[i], colL, colR, nL, nR, allw_neg, cmat=cmat)
            self.updateBests(22, best[i])

            for b in best[i]:
                tL = colL.getLiteralCat(nL, b[-1][-1][0], allw_neg)
                tR = colR.getLiteralCat(nR, b[-1][-1][1], allw_neg)
                # tL = Literal(nL, CatTerm(colL.getId(), b[-1][-1][0]))
                # tR = Literal(nR, CatTerm(colR.getId(), b[-1][-1][1]))
                if self.constraints.getCstr("debug_checks"):  # DEBUG
                    self.checkCountsPair(side, False, colL, colR, tL, tR, b[1])
                self.addPairToStore(b[0], b[1], tL, tR, nL, nR, side)

    def subdo23Full(self, colL, colR, side, try_comb=True):
        multi_cats = try_comb and self.constraints.getCstr("multi_cats")
        if side == 0:
            (colF, colE) = (colR, colL)
        else:
            (colF, colE) = (colL, colR)

        configs = [(False, False), (False, True), (True, False), (True, True)]
        allw_neg = True
        if True not in self.constraints.getCstr("neg_query", side=side, type_id=3):
            configs = configs[:1]
            allw_neg = False
        best = [[] for c in configs]

        buckets = colE.buckets()
        nbb = self.constraints.getCstr("max_prodbuckets") / len(colF.cats())
        if self.constraints.getCstr("max_prodbuckets") > 0 and len(buckets[1]) > nbb:  # self.constraints.getCstr("max_sidebuckets"):
            buckets = colE.buckets("collapsed", {"max_agg": self.constraints.getCstr("max_agg"), "nbb": nbb})

        # TODO DOABLE
        if buckets is not None and (self.constraints.getCstr("max_prodbuckets") == 0 or (len(buckets[1]) * len(colF.cats()) <= self.constraints.getCstr("max_prodbuckets"))):
            partsMubB = len(colF.miss())
            missMubB = len(colF.miss() & colE.miss())

            missMat = [len(colF.miss() & buk) for buk in buckets[0]]

            marg = [len(buckets[0][i]) - missMat[i] for i in range(len(buckets[0]))]
            if buckets[2] is not None:
                lEmodeFmiss = len(colE.interMode(colF.miss()))
                marg[buckets[2]] += colE.lenMode() - lEmodeFmiss
                missMat[buckets[2]] += lEmodeFmiss

            totMiss = sum(missMat)

            for cat in colF.cats():
                lFcEm = len(colF.suppCat(cat) & colE.miss())

                clp_tm = [self.constraints.getSSetts().makeLParts([colF.lsuppCat(cat),
                                                                   colF.nbRows() - colF.lsuppCat(cat) - partsMubB,
                                                                   partsMubB], side=side),
                          self.constraints.getSSetts().makeLParts([lFcEm,
                                                                   colE.lMiss() - lFcEm - missMubB,
                                                                   missMubB], side=side)]

                interMat = [len(colF.suppCat(cat) & buk) for buk in buckets[0]]
                if buckets[2] is not None:
                    interMat[buckets[2]] += len(colE.interMode(colF.suppCat(cat)))

                totIn = sum(interMat)
                below = 0
                missBelow = 0
                low = 0
                while low < len(interMat) and \
                        (totIn - below >= self.constraints.getCstr("min_itm_in") or
                         totIn - below >= self.constraints.getCstr("min_itm_out")):
                    above = 0
                    missAbove = 0
                    up = len(interMat)-1
                    while up >= low and \
                            (totIn - below - above >= self.constraints.getCstr("min_itm_in") or
                             totIn - below - above >= self.constraints.getCstr("min_itm_out")):
                        lin = self.constraints.getSSetts().makeLParts([totIn - below - above,
                                                                       sum(marg[low:up+1]) - totIn + below + above,
                                                                       totMiss - missBelow - missAbove], side=side)

                        for i, (nF, nE) in enumerate(configs):
                            if nF:
                                tmp_clp_tm = [self.constraints.getSSetts().negatePartsV(1-side, p) for p in clp_tm]
                                tmp_lin = self.constraints.getSSetts().negatePartsV(1-side, lin)
                            else:
                                tmp_clp_tm = clp_tm
                                tmp_lin = lin

                            stt = self.updateACTList(23, best[i], (cat, low, up), side, True, nE, tmp_clp_tm+[tmp_lin], pre_multi=multi_cats)

                        above += interMat[up]
                        missAbove += missMat[up]
                        up -= 1
                    below += interMat[low]
                    missBelow += missMat[low]
                    low += 1

        bUp = NumColM.buk_ind_maxes(buckets)
        for i, (nF, nE) in enumerate(configs):
            if multi_cats:
                self.combCatsNum(best[i], colF, colE, nF, nE, allw_neg, side, buckets, bUp)
            self.updateBests(23, best[i])

            for b in best[i]:
                tE = colE.getLiteralBuk(nE, buckets[1], b[-1][-1][1:], buckets[bUp])
                tF = colF.getLiteralCat(nF, b[-1][-1][0])
                if self.constraints.getCstr("debug_checks"):  # DEBUG
                    self.checkCountsPair(side, False, colF, colE, tF, tE, b[1])
                if tE is not None and tF is not None:
                    self.addPairToStore(b[0], b[1], tF, tE, nF, nE, side)
                    
