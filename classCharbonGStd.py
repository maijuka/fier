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


class CharbonGStd(CharbonGreedy):

    name = "GreedyStd"

    def handlesMiss(self):
        return False

    ##################################################
    # TOOLS METHODS
    def prepareCountsTM(self, supports, col):
        lparts = supports.lparts()
        return [lparts]

    def mkSPartsLits(self, colL, colR, litL=None, litR=None):
        return SParts(self.constraints.getSSetts(), colL.nbRows(),
                      [colL.suppLiteral(litL) if litL is not None else set(), colR.suppLiteral(litR) if litR is not None else set()])

    @classmethod
    def fillClp(tcl, clp_ti, neg=False):  # ads lout and flips with lin if negated
        if neg:
            return (clp_ti[0], clp_ti[1], [clp_ti[0][0] - clp_ti[1][0], clp_ti[0][1] - clp_ti[1][1], clp_ti[0][2] - clp_ti[1][2], clp_ti[0][3] - clp_ti[1][3]])
        else:
            return (clp_ti[0], [clp_ti[0][0] - clp_ti[1][0], clp_ti[0][1] - clp_ti[1][1], clp_ti[0][2] - clp_ti[1][2], clp_ti[0][3] - clp_ti[1][3]], clp_ti[1])
    # @classmethod
    # def fillClp(tcl, clp_ti, neg=False):  # ads lout and flips with lin if negated
    #     lparts, lin = clp_ti
    #     lout = [lparts[i] - lin[i] for i in range(len(lparts))]
    #     if neg:
    #         return (lparts, lin, lout)
    #     else:
    #         return (lparts, lout, lin)

    def inSuppBounds(self, side, op, lparts):
        if side == 1:
            (Eox, Exo, Exx, Eoo) = (0, 1, 2, 3)  # range(4)
        else:
            (Exo, Eox, Exx, Eoo) = (0, 1, 2, 3)  # range(4)
        if op:  # OR disjunction
            return (lparts[Exx] + lparts[Eox] >= self.constraints.getCstrSimple("min_itm_in")) and \
                   (lparts[Eox] >= self.constraints.getCstrSimple("min_itm_c")) and \
                   (lparts[Eoo] >= self.constraints.getCstrSimple("min_itm_out"))
        else:  # AND conjunction
            return (lparts[Exx] >= self.constraints.getCstrSimple("min_itm_in")) and \
                   (lparts[Exo] >= self.constraints.getCstrSimple("min_itm_c")) and \
                   (lparts[Eoo] + lparts[Exo] >= self.constraints.getCstrSimple("min_itm_out"))

    def clpToPieces(self, side, op, neg, lparts, lin, is_cond=False):
        if side == 1:
            (Eox, Exo, Exx, Eoo) = (0, 1, 2, 3)  # range(4)
        else:
            (Exo, Eox, Exx, Eoo) = (0, 1, 2, 3)  # range(4)

        if is_cond:  # Condition
            if neg:
                return 0, lparts[Exx] - lin[Exx], 0, lparts[Exo] - lin[Exo], -1, lparts[Exx] - lin[Exx]
            else:
                return 0, lin[Exx], 0, lin[Exo], -1, lin[Exx]
        elif op:  # OR disjunction
            if neg:
                return lparts[Exx], lparts[Eox] - lin[Eox], lparts[Exx] + lparts[Exo] + lparts[Eox], lparts[Eoo] - lin[Eoo], lin[Eoo], lparts[Eox] - lin[Eox]
            else:
                return lparts[Exx], lin[Eox], lparts[Exx]+lparts[Exo]+lparts[Eox], lin[Eoo], lparts[Eoo] - lin[Eoo], lin[Eox]
        else:  # AND conjunction
            if neg:
                return 0, lparts[Exx] - lin[Exx], lparts[Eox] + lparts[Exx], lparts[Exo] - lin[Exo], lin[Exo] + lparts[Eoo], lin[Exo]
            else:
                return 0, lin[Exx], lparts[Eox] + lparts[Exx], lin[Exo], lparts[Exo] - lin[Exo] + lparts[Eoo], lparts[Exo] - lin[Exo]

    def clpToPiecesLng(self, side, op, neg, lparts, lin, is_cond=False):
        if side == 1:
            (Eox, Exo, Exx, Eoo) = (0, 1, 2, 3)  # range(4)
        else:
            (Exo, Eox, Exx, Eoo) = (0, 1, 2, 3)  # range(4)

        if is_cond:  # Condition
            if neg:
                fix_num = 0
                var_num = lparts[Exx] - lin[Exx]
                fix_den = 0
                var_den = lparts[Exo] - lin[Exo]
                sout = -1
                cont = lparts[Exx] - lin[Exx]
            else:
                fix_num = 0
                var_num = lin[Exx]
                fix_den = 0
                var_den = lin[Exo]
                sout = -1
                cont = lin[Exx]
        elif op:  # OR disjunction
            if neg:
                fix_num = lparts[Exx]
                var_num = lparts[Eox] - lin[Eox]
                fix_den = lparts[Exx] + lparts[Exo] + lparts[Eox]
                var_den = lparts[Eoo] - lin[Eoo]
                sout = lin[Eoo]
                cont = lparts[Eox] - lin[Eox]
            else:
                fix_num = lparts[Exx]
                var_num = lin[Eox]
                fix_den = lparts[Exx]+lparts[Exo]+lparts[Eox]
                var_den = lin[Eoo]
                sout = lparts[Eoo] - lin[Eoo]  # lout[Eoo]
                cont = lin[Eox]
        else:  # AND conjunction
            if neg:
                fix_num = 0
                var_num = lparts[Exx] - lin[Exx]
                fix_den = lparts[Eox] + lparts[Exx]
                var_den = lparts[Exo] - lin[Exo]
                sout = lin[Exo] + lparts[Eoo]
                cont = lin[Exo]
            else:
                fix_num = 0
                var_num = lin[Exx]
                fix_den = lparts[Eox] + lparts[Exx]
                var_den = lin[Exo]
                sout = lparts[Exo] - lin[Exo] + lparts[Eoo]  # lout[Exo]+lparts[Eoo]
                cont = lparts[Exo] - lin[Exo]  # lout[Exo]
        return fix_num, var_num, fix_den, var_den, sout, cont

    def advRatioVar(self, side, op, lin_f, is_cond=False):
        if side == 1:
            (Eox, Exo, Exx, Eoo) = (0, 1, 2, 3)  # range(4)
        else:
            (Exo, Eox, Exx, Eoo) = (0, 1, 2, 3)  # range(4)

        if is_cond:
            num = lin_f[Exx]
            den = lin_f[Exo]
            return self.ratio(num, den+num)
        if op:  # OR disjunction
            num = lin_f[Eox]
            den = lin_f[Eoo]
        else:  # AND conjunction
            num = lin_f[Exx]
            den = lin_f[Exo]
        return self.ratio(num, den)

    def advAcc(self, side, op, neg, clp_ti, is_cond=False, filled=False):
        (fix_num, var_num, fix_den, var_den, sout, cont) = self.clpToPieces(side, op, neg, clp_ti[0], clp_ti[-1], is_cond)
        return self.offset_ratio(fix_num + var_num, fix_den + var_den)

    def getAC(self, side, op, neg, clp_ti, is_cond=False, no_const=False, filled=False):
        # try:
        (fix_num, var_num, fix_den, var_den, sout, cont) = self.clpToPieces(side, op, neg, clp_ti[0], clp_ti[-1], is_cond)
        # except IndexError:
        #     print(clp_ti[0], clp_ti[-1])
        #     pdb.set_trace()
        #     (fix_num, var_num, fix_den, var_den, sout, cont) = self.clpToPieces(side, op, neg, clp_ti[0], clp_ti[-1], is_cond)

        if is_cond:
            if self.unconstrained(no_const) or cont >= self.constraints.getCstrSimple("min_itm_in"):
                acc = self.offset_ratio(var_num, var_den)
                clp = clp_ti if filled else self.fillClp(clp_ti, neg)
                return (acc, var_num, var_den, cont, fix_num, fix_den), clp, 0
            return None, None, self.FLG_IN | self.FLG_CONT

        cstr_status = 0
        if not self.constraints.getCstr("support_separate_sides"):
            if var_num+fix_num < self.constraints.getCstrSimple("min_itm_in"):
                cstr_status |= self.FLG_IN
            if sout < self.constraints.getCstrSimple("min_itm_out"):
                cstr_status |= self.FLG_OUT
        if cont < self.constraints.getCstrSimple("min_itm_c"):
            cstr_status |= self.FLG_CONT
        # print(cstr_status, var_num+fix_num, sout, cont, clp_ti[0], clp_ti[-1])
        if self.unconstrained(no_const) or self.constraints.getCstr("support_separate_sides") or cstr_status == 0:
            acc = self.offset_ratio(var_num + fix_num, var_den + fix_den)
            clp = clp_ti if filled else self.fillClp(clp_ti, neg)
            return (acc, var_num, var_den, cont, fix_num, fix_den), clp, cstr_status
        return None, None, cstr_status

    ##################################################
    # SEARCH METHODS
    def doBoolStar(self, colL, colR, side):
        if side == 1:
            (supports, fixTerm, extCol) = (SParts(self.constraints.getSSetts(), colL.nbRows(),
                                                  [colL.supp(), set()]), BoolTerm(colL.getId()), colR)
        else:
            (supports, fixTerm, extCol) = (SParts(self.constraints.getSSetts(), colL.nbRows(),
                                                  [set(), colR.supp()]), BoolTerm(colR.getId()), colL)
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
            totInt = colE.nbRows()
            # margE = [len(intE) for intE in bucketsE[0]]

            margF = [len(bucketsF[0][i]) for i in range(len(bucketsF[0]))]

            for bukF in bucketsF[0]:
                interMat.append([len(bukF & bukE) for bukE in bucketsE[0]])

            if bucketsF[2] is not None:
                margF[bucketsF[2]] += colF.lenMode()
                for bukEId in range(len(bucketsE[0])):
                    interMat[bucketsF[2]][bukEId] += len(colF.interMode(bucketsE[0][bukEId]))

            if bucketsE[2] is not None:
                # margE[bucketsE[2]] += colE.lenMode()
                for bukFId in range(len(bucketsF[0])):
                    interMat[bukFId][bucketsE[2]] += len(colE.interMode(bucketsF[0][bukFId]))

            if bucketsF[2] is not None and bucketsE[2] is not None:
                interMat[bucketsF[2]][bucketsE[2]] += len(colE.interMode(colF.modeSupp()))

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
                        # totEinF = sum(EinF)
                        if side == 0:
                            clp_tm = [[0, totInt - aboveF - belowF, 0, aboveF + belowF]]
                        else:
                            clp_tm = [[totInt - aboveF - belowF, 0, 0, aboveF + belowF]]
                        # clp_tm = [self.constraints.getSSetts().makeLParts([totInt - aboveF - belowF,
                        #                                                    aboveF + belowF], side=side)]

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
                                    # lin = self.constraints.getSSetts().makeLParts([totInt - belowF - aboveF - belowEF - aboveEF,
                                    #                                                belowF + aboveF - outAboveEF - outBelowEF], side=side)
                                    if side == 0:
                                        lin = [0, totInt - belowF - aboveF - belowEF - aboveEF, 0, belowF + aboveF - outAboveEF - outBelowEF]
                                    else:
                                        lin = [totInt - belowF - aboveF - belowEF - aboveEF, 0, 0, belowF + aboveF - outAboveEF - outBelowEF]

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

        cmat = numpy.zeros((colL.nbCats(), colR.nbCats()), dtype=int)
        tot = colL.nbRows()
        for ciL, catL in enumerate(colL.cats()):
            totL = len(colL.suppCat(catL))
            clp_tm = [[totL, 0, 0, colL.nbRows()-totL]]

            for ciR, catR in enumerate(colR.cats()):
                totR = len(colR.suppCat(catR))
                interLR = len(colL.suppCat(catL) & colR.suppCat(catR))
                cmat[ciL, ciR] = interLR
                lin = [interLR, 0, 0, totR - interLR]
                for i, (nL, nR) in enumerate(configs):
                    if nL:
                        tmp_clp_tm = [self.constraints.getSSetts().negateParts(0, p) for p in clp_tm]
                        tmp_lin = self.constraints.getSSetts().negateParts(0, lin)
                    else:
                        tmp_clp_tm = clp_tm
                        tmp_lin = lin

                    stt = self.updateACTList(22, best[i], (catL, catR), side, True, nR, tmp_clp_tm+[tmp_lin], pre_multi=self.constraints.getCstr("multi_cats"))

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

            marg = [len(buk) for buk in buckets[0]]
            if buckets[2] is not None:
                marg[buckets[2]] += colE.lenMode()

            for cat in colF.cats():
                clp_tm = [self.constraints.getSSetts().makeLParts([colF.lsuppCat(cat),
                                                                   colF.nbRows() - colF.lsuppCat(cat)], side=side)]

                interMat = [len(colF.suppCat(cat) & buk) for buk in buckets[0]]
                if buckets[2] is not None:
                    interMat[buckets[2]] += len(colE.interMode(colF.suppCat(cat)))

                totIn = sum(interMat)
                below = 0
                low = 0
                while low < len(interMat) and \
                        (totIn - below >= self.constraints.getCstrSimple("min_itm_in") or
                         totIn - below >= self.constraints.getCstrSimple("min_itm_out")):
                    above = 0
                    up = len(interMat)-1
                    while up >= low and \
                            (totIn - below - above >= self.constraints.getCstrSimple("min_itm_in") or
                             totIn - below - above >= self.constraints.getCstrSimple("min_itm_out")):
                        lin = self.constraints.getSSetts().makeLParts([totIn - below - above,
                                                                       sum(marg[low:up+1]) - totIn + below + above], side=side)

                        for i, (nF, nE) in enumerate(configs):
                            if nF:
                                tmp_clp_tm = [self.constraints.getSSetts().negateParts(1-side, p) for p in clp_tm]
                                tmp_lin = self.constraints.getSSetts().negatePartsV(1-side, lin)
                            else:
                                tmp_clp_tm = clp_tm
                                tmp_lin = lin

                            stt = self.updateACTList(23, best[i], (cat, low, up), side, True, nE, tmp_clp_tm+[tmp_lin], pre_multi=multi_cats)

                        above += interMat[up]
                        up -= 1
                    below += interMat[low]
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
