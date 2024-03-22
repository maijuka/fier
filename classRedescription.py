import re
import string
import numpy
try:
    from classCol import ColM
    from classQuery import *
    from classSParts import SParts, tool_pValSupp, tool_pValOver, tool_ratio, SSetts
    from classProps import WithEVals, RedProps, mapSuppNames, ACTIVE_RSET_ID, HAND_SIDE
except ModuleNotFoundError:
    from .classCol import ColM
    from .classQuery import *
    from .classSParts import SParts, tool_pValSupp, tool_pValOver, tool_ratio, SSetts
    from .classProps import WithEVals, RedProps, mapSuppNames, ACTIVE_RSET_ID, HAND_SIDE

import pdb


class Redescription(WithEVals):
    diff_score = Query.diff_length + 1

    # PROPS WHAT
    info_what_dets = {}
    # info_what_dets = {"queryLHS": "self.prepareQueryLHS",
    #                   "queryRHS": "self.prepareQueryRHS",
    #                   "queryCOND": "self.prepareQueryCOND"}
    info_what = {"uid": "self.getUid()", "rid": "self.getShortId()",
                 "hasAvC": "self.hasAvailableCols()", "nbAvC": "self.totAvailableCols()", "diffLengthQs": "self.diffLengthQs()",
                 "containsAnon": "self.containsAnon()", "isTreeCompatible": "self.isTreeCompatible()",
                 "isBasis": "self.isBasis()", "status_enabled": "self.getStatus()",
                 "track": "self.getTrackStr()", "trackInit": "self.getTrackInitStr()", "trackExt": "self.getTrackExtStr()"}
    Pwhat_match = "(" + "|".join(["extra"]+list(info_what.keys())+list(info_what_dets.keys())) + ")"

    pair_map_side = {"L": 0, "R": 1}
    ## below, "O" is for the default way to compute the overlap, "overlap" is heritage
    pair_map_how = {"areaU": "union", "areaI": "inter", "areaR": "ratio", "areaF": "fractmin", "overlap": "ratio", "areaO": "ratio",
                    "colsU": ("cols", "union"), "colsI": ("cols", "inter"), "colsR": ("cols", "ratio"), "colsF": ("cols", "fractmin"), "colsO": ("cols", "fractmin"),
                    "rowsU": ("rows", "union"), "rowsI": ("rows", "inter"), "rowsR": ("rows", "ratio"), "rowsF": ("rows", "fractmin"), "rowsO": ("cols", "fractmin")}
    info_pair_what_dets = {"oneSideIdentical": "self.oneSideIdentical", "bothSidesIdentical": "self.bothSidesIdentical",
                           "equivalent": "self.equivalent", "superceding": "self.superceding", "compare": "self.compare",
                           "sameRectangles": "self.sameRectangles",
                           "overlapAreaTotal": "self.overlapAreaTotal", "overlapAreaMax": "self.overlapAreaMax",
                           "overlapRows": "self.overlapRows", "overlapCols": "self.overlapCols"}
    info_pair_what_side = {"interAreaSide": "self.interArea", "unionAreaSide": "self.unionArea", "overlapAreaSide": "self.overlapAreaSide"}
    Pwhat_pair_match = "(" + "|".join(list(pair_map_how.keys())+list(info_pair_what_side.keys())+list(info_pair_what_dets.keys())) + ")"

    class_letter = "r"
    RP = None
    @classmethod
    def setupRP(tcl, fields_fns=None):
        elems_typs = [(Query.class_letter, Query), (SParts.class_letter, SParts), (Redescription.class_letter, Redescription)]
        RedProps.setupProps(Query, Redescription, elems_typs)
        tcl.RP = RedProps(fields_fns)
    # getRP(tcl, rp=None) is defined in WithEVals class

    def __init__(self, nqueryL=None, nqueryR=None, nsupps=None, nN=-1, nPrs=[-1, -1], ssetts=None, iid=None):
        WithEVals.__init__(self, iid)
        if nqueryL is None:
            nqueryL = Query()
        if nqueryR is None:
            nqueryR = Query()
        self.queries = [nqueryL, nqueryR]
        if nsupps is not None:
            if isinstance(nsupps, SParts):
                self.sParts = nsupps
            else:
                self.sParts = SParts(ssetts, nN, nsupps, nPrs)
            self.dict_supp_info = None
        else:
            self.sParts = None
            self.dict_supp_info = {}
        self.lAvailableCols = [None, None]
        self.extras["status"] = 1
        self.extras["track"] = []
        self.condition = None

        
    @classmethod
    def fromInitialPair(tcl, initialPair, data, litC=None):
        if initialPair[0] is None and initialPair[1] is None:
            return None
        supps_miss = [set(), set(), set(), set()]
        queries = [None, None]
        for side in [0, 1]:
            suppS, missS = (set(), set())
            if type(initialPair[side]) is Query:
                queries[side] = initialPair[side]
                suppS, missS = initialPair[side].recompute(side, data)
            else:
                queries[side] = Query()
                if type(initialPair[side]) is Literal:
                    queries[side].extend(None, initialPair[side])
                    suppS, missS = data.literalSuppMiss(side, initialPair[side])
            supps_miss[side] = suppS
            supps_miss[side+2] = missS
        r = tcl(queries[0], queries[1], supps_miss, data.nbRows(), [len(supps_miss[0])/float(data.nbRows()), len(supps_miss[1])/float(data.nbRows())], data.getSSetts())
        r.initTrack()

        qC = None
        if litC is not None:
            if type(litC) is list:
                # if len(litC) > 1: pdb.set_trace()
                qC = Query(OR=False, buk=litC)
            else:
                qC = Query(buk=[litC])
        elif len(initialPair) > 2 and type(initialPair[-1]) is Query:
            qC = initialPair[-1]

        if qC is not None:
            supp_cond, miss_cond = qC.recompute(-1, data)
            r.setCondition(qC, supp_cond)
        WithEVals.recompute(r, data)
        return r

    @classmethod
    def fromQueriesPair(tcl, queries, data, copyQ=True, iid=None):
        if copyQ:
            r = tcl(queries[0].copy(), queries[1].copy(), iid=iid)
        else:
            r = tcl(queries[0], queries[1], iid=iid)
        r.recompute(data)
        r.initTrack()
        
        if len(queries) > 2 and queries[2] is not None:
            if copyQ:
                qC = queries[2].copy()
            else:
                qC = queries[2]
            supp_cond, miss_cond = qC.recompute(-1, data)
            r.setCondition(qC, supp_cond)
        return r

    @classmethod
    def fromCol(tcl, col, data):
        queries = [Query(), Query()]
        if col.getSide() == -1:
            queries.append(Query())
        if isinstance(col, ColM):
            queries[col.getSide()].extend(-1, Literal(False, col.getAnonTerm()))
        r = tcl.fromQueriesPair(queries, data, copyQ=False)
        # r = tcl(queries[0], queries[1])
        # r.recompute(data)
        r.initTrack()
        return r

    def isBasis(self):
        return (len(self.query(0)) == 0 and self.query(1).isBasis()) or \
               (len(self.query(1)) == 0 and self.query(0).isBasis())

    def isXpr(self):
        return (len(self.query(0)) == 0 and self.query(1).isXpr()) or \
               (len(self.query(1)) == 0 and self.query(0).isXpr())

    def getXprTerm(self):
        if self.isXpr():
            return self.query(0).getXprTerm() if self.query(0).isXpr() else self.query(1).getXprTerm()

    def containsAnon(self, side=None):
        if side is None:
            return self.query(0).containsAnon() or self.query(1).containsAnon()
        else:
            return self.query(side).containsAnon()

    def minusAnon(self, side=None):
        if side is None:
            q0noan, l0an = self.query(0).minusAnon()
            q1noan, l1an = self.query(1).minusAnon()
            return [q0noan, q1noan], [l0an, l1an]
        else:
            return self.query(side).minusAnon()

    def minusAnonRed(self, data):
        r, modr = self, 0
        if len(self) == 0:
            return r, [[], []], modr
        qsNoan, lsAnon = self.minusAnon()
        if (len(lsAnon[0]) + len(lsAnon[1])) > 0:  # contains anon, standard list of one-step extensions
            if (len(qsNoan[0]) + len(qsNoan[1])) == 0:  # only anon, need to start by mining pairs...
                modr = -1
            else:
                modr = 1
            r = Redescription.fromQueriesPair(qsNoan, data)
        return r, lsAnon, modr

    def typeId(self):
        if self.isBasis():
            return self.query(0).typeId() or self.query(1).typeId()
        return None

    def dropSupport(self):
        if self.sParts is not None:
            self.dict_supp_info.toDict()
            self.sParts = None

    # def __hash__(self):
    #      return int(hash(self.queries[0])+ hash(self.queries[1])*100*self.score())

    def __len__(self):
        return self.length(0) + self.length(1)
    
    def __lt__(self, other):
        return self.getAcc() < other.getAcc()

    def usesOr(self, side=None):
        if side is not None:
            return self.queries[side].usesOr()
        return self.queries[0].usesOr() or self.queries[1].usesOr()

    def isTreeCompatible(self, side=None):
        if side is not None:
            return self.queries[side].isTreeCompatible()
        return self.queries[0].isTreeCompatible() and self.queries[1].isTreeCompatible()

    def supp(self, side, rset_det=None):
        if rset_det is None:
            return self.supports().supp(side)
        set_parts = self.getRSetParts(rset_det)
        if set_parts is not None:
            return set_parts.supp(side)

    def miss(self, side, rset_det=None):
        if rset_det is None:
            return self.supports().miss(side)
        set_parts = self.getRSetParts(rset_det)
        if set_parts is not None:
            return set_parts.miss(side)

    def score(self, rset_det=None):
        return self.getAcc(rset_det)

    def supports(self):
        return self.sParts

    def nbRows(self):
        return self.supports().nbRows()

    def rows(self):
        return set(range(self.nbRows()))

    def partsAll(self):
        return self.supports().sParts

    def partsFour(self):
        return [self.supports().suppA(), self.supports().suppB(), self.supports().suppI(), self.supports().suppO()]

    def partsThree(self):
        return [self.supports().suppA(), self.supports().suppB(), self.supports().suppI()]

    def partsNoMiss(self):
        return self.supports().sParts[:4]

    def query(self, side=None):
        if side == -1:
            return self.getQueryC()
        return self.queries[side]

    def getQueries(self):
        return self.queries

    def getQueryC(self):
        if self.condition is not None:
            return self.condition.get("q")

    def getSupportsC(self):
        if self.condition is not None:
            return self.condition.get("sparts")

    def getSuppC(self):
        if self.condition is not None:
            return self.condition.get("supp")

    def hasCondition(self):
        return self.condition is not None

    def setCondition(self, qC=None, supp_cond=None):
        self.condition = None
        if qC is not None:
            if supp_cond is None:
                sparts = None
            else:
                sparts = self.supports().copy()
                sparts.update(0, False, supp_cond)
                sparts.update(1, False, supp_cond)
            self.condition = {"q": qC, "supp": supp_cond, "sparts": sparts}

    def probas(self):
        return self.supports().probas()

    def probasME(self, dbPrs, epsilon=0):
        return [self.queries[side].probaME(dbPrs, side, epsilon) for side in [0, 1]]

    def surpriseME(self, dbPrs, epsilon=0):
        # return [-numpy.sum(numpy.log(numpy.absolute(SParts.suppVect(self.supports().nbRows(), self.supports().suppSide(side), 0) - self.queries[side].probaME(dbPrs, side)))) for side in [0,1]]
        return -numpy.sum(numpy.log(numpy.absolute(SParts.suppVect(self.supports().nbRows(), self.supports().suppI(), 0) - self.queries[0].probaME(dbPrs, 0)*self.queries[1].probaME(dbPrs, 1))))

    def exME(self, dbPrs, epsilon=0):
        prs = [self.queries[side].probaME(dbPrs, side, epsilon) for side in [0, 1]]
        surprises = []
        tmp = [i for i in self.supports().suppI() if prs[0][i]*prs[1][i] == 0]
        surprises.append(-numpy.sum(numpy.log([prs[0][i]*prs[1][i] for i in self.supports().suppI()])))
        surprises.extend([-numpy.sum(numpy.log([prs[side][i] for i in self.supports().suppSide(side)])) for side in [0, 1]])

        return surprises + [len(tmp) > 0]

        N = self.supports().nbRows()
        margsPr = [numpy.sum([prs[side][i] for i in self.supports().suppSide(side)]) for side in [0, 1]]
        pvals = [tool_pValOver(self.supports().lenI(), N, int(margsPr[0]), int(margsPr[1])), tool_pValSupp(N, self.supports().lenI(), margsPr[0]*margsPr[1]/N**2)]
        return surprises, pvals

    def length(self, side, ex_anon=False):
        return self.queries[side].length(ex_anon)

    def diffLengthQs(self):
        return abs(self.length(0) - self.length(1))

    def prepareExtElems(self, data=None, single_dataset=False, souvenirs=None):
        r, lsAnon, modr = self.minusAnonRed(data)
        if modr != 0:  # contains anon, use to restrict available vars
            for side, lits in enumerate(lsAnon):
                org_available = self.lAvailableCols[side]
                if len(lits) > 0 or modr == 1:
                    still_available = [l[1].colId() for l in lits]
                    org_available = set()
                else:
                    still_available = None
                r.restrictAvailable(side, org_available, still_available)
        exts = []
        if modr == -1:  # only anon, need to start by mining pairs...
            for vLHS in r.lAvailableCols[0] or []:
                for vRHS in r.lAvailableCols[1] or []:
                    if (not single_dataset or vLHS != vRHS) and (not data.hasGroups(0) or not data.hasGroups(1) or data.areGroupCompat(vLHS, vRHS, 0, 1)):
                        exts.append((0, vLHS, data.col(1, vRHS)))
        elif len(r) > 0:
            for side in [0, 1]:
                if self.lAvailableCols[side] is not None and self.length(1-side) != 0:  # if other side is empty, that should be extended first
                    df = self.queries[1-side].invCols() if single_dataset else set()  # and self.lAvailableCols[1-side] is not None
                    # exts.extend([(side, v, r) for v in r.availableColsSide(side, data, single_dataset, souvenirs=souvenirs)])
                    exts.extend([(side, v, r) for v in r.lAvailableCols[side].difference(df)])
        return exts, r, modr

    def isAvailableCol(self, side, cid, single_dataset=False):
        if self.lAvailableCols[side] is not None and self.length(1-side) != 0:
            df = self.queries[1-side].invCols() if single_dataset else set()  # and self.lAvailableCols[1-side] is not None
            # exts.extend([(side, v, r) for v in r.availableColsSide(side, data, single_dataset, souvenirs=souvenirs)])
            return cid in self.lAvailableCols[side].difference(df)
        
    def hasAvailableCols(self):
        return any([(ss is not None and len(ss) > 0) for ss in self.lAvailableCols])

    def totAvailableCols(self):
        return sum([len(ss) if ss is not None else 0 for ss in self.lAvailableCols])

    def getNbAvC(self, details={}):
        if details.get("format") == "str":
            return " + ".join(["?" if self.lAvailableCols[side] is None else "%i" % len(self.lAvailableCols[side]) for side in self.getAvailableS()])
        else:
            return [0 if self.lAvailableCols[side] is None else len(self.lAvailableCols[side]) for side in self.getAvailableS()]

    def getAvailableS(self):
        return range(len(self.lAvailableCols))

    def copyAvailableCols(self, side=None):
        if side is None:
            return [self.copyAvailableCols(side) for side in self.getAvailableS()]
        else:
            if self.lAvailableCols[side] is None:
                return None
            else:
                return set(self.lAvailableCols[side])

    def initAvailable(self, souvenirs, data, max_var=None):
        invs = (self.invColsSide(0), self.invColsSide(1))
        for side in [0, 1]:
            if self.lAvailableCols[side] is None:
                self.lAvailableCols[side] = souvenirs.copyAvailableCols(side)
                for ss in [0, 1]:
                    data.upColsCompat(self.lAvailableCols[side], ss, invs[ss], side != ss)
        self.setFull(max_var)

    def updateAvailable(self, souvenirs):
        for side in [0, 1]:
            if self.lAvailableCols[side] is not None:
                self.lAvailableCols[side].difference_update(souvenirs.extOneStep(self, side))

                    

    def restrictAvailable(self, side=None, org_available=None, still_available=None, not_available=None):
        if side is None:
            for side in self.getAvailableS():
                ss = [subs[side] if subs is not None else None for subs in [org_available, still_available, not_available]]
                self.restrictAvailable(side, *ss)
        else:
            if org_available is not None:
                if still_available is not None:
                    if len(org_available) == 0:
                        self.lAvailableCols[side] = set(still_available)
                    else:
                        self.lAvailableCols[side] = org_available.intersection(still_available)
                else:
                    self.lAvailableCols[side] = set(org_available)
                if not_available is not None and self.lAvailableCols[side] is not None:
                    self.lAvailableCols[side].difference_update(not_available)

    def cutOffAvailables(self, side=None):
        if side is None:
            self.lAvailableCols = [None, None]
        else:
            self.lAvailableCols[side] = None

    def wasCutOffAvailables(self, side):
        return self.lAvailableCols[side] is None

    def update(self, data=None, side=-1, opBool=None, literal=None, suppX=None, missX=None):
        if side == -1:
            self.cutOffAvailables()
        else:
            op = Op(opBool)
            self.queries[side].extend(op, literal)
            self.supports().update(side, op.isOr(), suppX, missX)
            self.dict_supp_info = None
            if data is None:
                if self.lAvailableCols[side] is not None:
                    self.lAvailableCols[side].remove(literal.colId())
            else:
                for ss in [0, 1]:
                    data.upColsCompat(self.lAvailableCols[ss], side, [literal.colId()], side != ss)
            self.appendTrack(self.makeTrackTuple(side, op, literal.colId(), literal))

    def setFull(self, max_var=None):
        if max_var is not None:
            for side in [0, 1]:
                if self.length(side, ex_anon=True) >= max_var[side]:
                    self.cutOffAvailables(side)

    def mkRed(self, data, side=-1, op=None, literal=None, suppX=None, missX=None):
        red = self.copy(with_cond=False)  # condition is not inherited from parent red, adjusted separately
        red.update(data, side, op, literal, suppX, missX)
        return red

    def copy(self, iid=None, with_cond=True):
        r = Redescription(self.queries[0].copy(), self.queries[1].copy(),
                          self.supports().supparts(), self.supports().nbRows(), self.probas(), self.supports().getSSetts(), iid=iid)
        if self.hasCondition() and with_cond:
            r.setCondition(self.getQueryC().copy(), self.getSuppC().copy())

        r.lAvailableCols = self.copyAvailableCols()
        r.extras = self.copyExtras()
        r.restricted_sets = {}
        for sid, rst in self.restricted_sets.items():
            r.restricted_sets[sid] = {"sParts": rst["sParts"],
                                      "prs": [rst["prs"][0], rst["prs"][1]],
                                      "rids": set(rst["rids"])}

        return r

    def recomputeQuery(self, side, data=None, restrict=None):
        return self.queries[side].recompute(side, data, restrict)

    def invLiteralsSide(self, side, ex_anon=False):
        return self.queries[side].invLiterals(ex_anon)

    def invLiterals(self, ex_anon=False):
        return [self.invLiteralsSide(0, ex_anon), self.invLiteralsSide(1, ex_anon)]

    def invColsSide(self, side, ex_anon=False):
        return self.queries[side].invCols(ex_anon)

    def invCols(self, ex_anon=False):
        return [self.invColsSide(0, ex_anon), self.invColsSide(1, ex_anon)]

    def getNormalized(self, data=None, side=None):
        if side is not None:
            sides = [side]
        else:
            sides = [0, 1]
        queries = [self.queries[side] for side in [0, 1]]
        c = [False, False]
        for side in sides:
            queries[side], c[side] = self.queries[side].algNormalized()
        if c[0] or c[1]:
            red = Redescription.fromQueriesPair(queries, data)
            # check that support is same
            # if self.supports() != red.supports():
            #     print("ERROR ! SUPPORT CHANGED WHEN NORMALIZING...")
            #     pdb.set_trace()
            return red, True
        else:
            return self, False

    def minusOneLiteralSupps(self, data, restrict=None):
        sm_lits = [{}, {}]
        supp_sets_org = [None, None, None, None]
        for side in [0, 1]:
            sm = self.query(side).recompute(side, data, restrict, sm_lits[side])
            supp_sets_org[side] = sm[0]
            supp_sets_org[2+side] = sm[1]

        minuses = []
        for side in [0, 1]:
            supp_sets = [set(s) for s in supp_sets_org]
            for minus in self.query(side).minusOneLiteralSupps(side, data, restrict, sm_lits[side]):
                supp_sets[side] = minus["sm_q"][0]
                supp_sets[2+side] = minus["sm_q"][1]
                minus["sparts"] = SParts(data.getSSetts(), data.nbRows(), supp_sets)
                minus["acc"] = minus["sparts"].acc()
                minus["side"] = side
                minuses.append(minus)
        return minuses

    def getPruned(self, data):
        r, modr = self, False
        minuses = self.minusOneLiteralSupps(data)
        best_pos = numpy.argmax([m["acc"] for m in minuses])
        if minuses[best_pos]["acc"] > self.score():
            best = minuses[best_pos]
            queries = [None, None]
            queries[best["side"]] = best["q"]
            queries[1-best["side"]] = self.query(1-best["side"]).copy()
            cand = Redescription(queries[0], queries[1], best["sparts"])
            r, modb = cand.getPruned(data)
            modr = True
        return r, modr

    def normalize(self, data=None):
        return self.getNormalized(data)

    def recompute(self, data):
        (nsuppL, missL) = self.recomputeQuery(0, data)
        (nsuppR, missR) = self.recomputeQuery(1, data)
#        print(self.disp())
#        print(' '.join(map(str, nsuppL)) + ' \t' + ' '.join(map(str, nsuppR)))
        self.sParts = SParts(data.getSSetts(), data.nbRows(), [nsuppL, nsuppR, missL, missR])
        self.prs = [self.queries[0].proba(0, data), self.queries[1].proba(1, data)]
        if self.hasCondition():
            qC = self.getQueryC()
            supp_cond, miss_cond = qC.recompute(-1, data)
            self.setCondition(qC, supp_cond)
        self.dict_supp_info = None
        WithEVals.recompute(self, data)

    def check(self, data):
        result = 0
        details = None
        if self.supports() is not None:  # TODO: sparts
            (nsuppL, missL) = self.recomputeQuery(0, data)
            (nsuppR, missR) = self.recomputeQuery(1, data)

            details = (len(nsuppL.symmetric_difference(self.supports().supp(0))) == 0,
                       len(nsuppR.symmetric_difference(self.supports().supp(1))) == 0,
                       len(missL.symmetric_difference(self.supports().miss(0))) == 0,
                       len(missR.symmetric_difference(self.supports().miss(1))) == 0)
            result = 1
            for detail in details:
                result *= detail
        return (result, details)

    def normalizeSides(self):
        if self.supports().getSSetts().isDataSymmetric():
            if sorted(self.query(1).invCols()) < sorted(self.query(0).invCols()):
                self.swapSides()

    def swapSides(self):
        self.queries[0], self.queries[1] = self.queries[1], self.queries[0]
        self.lAvailableCols[0], self.lAvailableCols[1] = self.lAvailableCols[1], self.lAvailableCols[0]
        self.sParts.swapSides()
        self.dict_supp_info = {}
        self.extras["track"].append((-1, -1))

    def hasMissing(self):
        return self.supports().hasMissing()

    def getStatus(self):
        return self.extras["status"]

    def setStatus(self, status):
        self.extras["status"] = status

    def getEnabled(self, details={}):
        return 1*(self.getStatus() > 0)

    def isEnabled(self, details={}):
        return self.getStatus() > 0

    def flipEnabled(self):
        self.setStatus(-self.getStatus())

    def setEnabled(self):
        self.setStatus(1)

    def setDisabled(self):
        self.setStatus(-1)

    def setDiscarded(self):
        self.setStatus(-2)

    # GET FIELDS INFO INVOLVING ADDITIONAL DETAILS (PRIMARILY FOR SIREN)
    def getQueryU(self, side, details={}):
        q = self.query(side)
        if q is None:
            return ""
        names = None
        fmts = None
        if details is not None:
            if "names" in details:
                names = details["names"][side]
            if "fmts" in details:
                fmts = details["fmts"][side]
        return q.styledDisp(names=names, fmts=fmts, style="utf")

    def getQueryLU(self, details={}):
        return self.getQueryU(0, details)

    def getQueryRU(self, details={}):
        return self.getQueryU(1, details)

    def getQueriesU(self, details={}):
        return self.getQueryLU(details) + "---" + self.getQueryRU(details)

    track_fields = ["code", "cid", "hash", "uid"]
    def makeTrackTuple(self, side, op, cid, literal, init=False):
        ### tracking history of building the queries, as a list of (x, c) pairs
        ### x is even if extending with an OR, odd otherwise
        ### x is negative if extending on the left (0), positive on the right (1)
        ### x >= if column part of initialization
        ### c is the id of the column
        return ((2*side-1) * ((1+2*init) + op.isOr()), cid, hash(literal), self.getUid())

    def initTrack(self):
        track = []
        for side in [0, 1]:
            for lit in self.queries[side].invLiterals():
                track.append(self.makeTrackTuple(side, self.queries[side].getOp(), lit.colId(), lit, init=True))        
        self.setTrack(track)

    
    @classmethod
    def filterTrackTuples(tcl, track=[], details={}):
        if details.get("filter-init") is not None:
            track = [t for t in track if (details["filter-init"] and abs(t[0]) >= 3) or (not details["filter-init"] and abs(t[0]) < 3)]
        if details.get("filter-or") is not None:
            track = [t for t in track if (details["filter-or"] and t[0]%2 == 0) or (not details["filter-or"] and t[0]%2 == 1)]
        if details.get("filter-side") is not None:
            track = [t for t in track if (details["filter-side"] and t[0] < 0) or (not details["filter-side"] and t[0] > 0)]
        return track

    def setTrack(self, track=[]):
        self.extras["track"] = track

    def appendTrack(self, t, p=None):
        if p is None:
            self.extras["track"].append(t)
        else:
            self.extras["track"].insert(p, t)

    def extendTrack(self, track=[]):
        self.extras["track"].extend(track)

    def getTrack(self, details={}):
        if details is not None and any([re.match("filter-", k) for k in details.keys()]):            
            track = self.filterTrackTuples([t for t in self.extras["track"]], details)
            tpos = [1]
            if details.get("tfield") in self.track_fields:
                tpos = [self.track_fields.index(details.get("tfield"))]
            elif type(details.get("tfield")) is int and 0 <= details.get("tfield") < len(self.track_fields):
                tpos = [details.get("tfield")]
            elif type(details.get("tfield")) is list and all([t in self.track_fields for t in details.get("tfield")]):
                tpos = [self.track_fields.index(t) for t in details.get("tfield")]
            elif type(details.get("tfield")) is list and all([0 <= t < len(self.track_fields) for t in details.get("tfield")]):
                tpos = details.get("tfield")
            if details.get("aim", None) == "list" or details.get("format", None) == "str":
                return ";".join([":".join(["%s" % t[ti] for ti in tpos]) for t in track])
            return [tuple(t[ti] for ti in tpos) for t in track]
                
        elif details is not None and (details.get("aim", None) == "list" or details.get("format", None) == "str"):
            return ";".join(["%s:%s" % (t[0], t[1]) for t in self.getTrack()])
        else:
            return self.extras["track"]

    def getTrackStr(self):
        return self.getTrack({"format": "str"})

    def getTrackInitStr(self):
        return self.getTrack({"format": "str", "filter-init": True, "tfield": ["cid", "hash"]})
    def getTrackExtStr(self):
        return self.getTrack({"format": "str", "filter-init": False, "filter-or": True})

    
    def getSortAble(self, details={}):
        if details is not None and details.get("aim") == "sort":
            return (self.getStatus(), self.getUid())
        return ""

    # def getCohesion(self, details={}):
    #     return self.getDetail("cohesion")
    # def getCohesionNat(self, details={}):
    #     return self.getDetail("cohesion_nat")

    def getTypeParts(self, details={}):
        return self.supports().getTypeParts()

    def getMethodPVal(self, details={}):
        return self.supports().getMethodPVal()

    # RESTRICTED SETS

    def setRestrictedSuppSets(self, data, supp_sets=None):
        resets_ids = []
        if supp_sets is None:
            if data.hasLT():
                supp_sets = data.getLT()
            else:
                supp_sets = {ACTIVE_RSET_ID: data.nonselectedRows()}
        for sid, sset in supp_sets.items():
            old_rids = None
            if sid in self.restricted_sets:
                old_rids = self.restricted_sets[sid]["rids"]
            if len(sset) == 0:
                if sid not in self.restricted_sets or len(self.restricted_sets[sid]["rids"]) > 0:
                    resets_ids.append(sid)
                self.restricted_sets[sid] = {"sParts": None,
                                             "prs": None,
                                             "rids": set()}

            elif sid not in self.restricted_sets or self.restricted_sets[sid]["rids"] != sset:
                (nsuppL, missL) = self.recomputeQuery(0, data, sset)
                (nsuppR, missR) = self.recomputeQuery(1, data, sset)
                if len(missL) + len(missR) > 0:
                    rsParts = SParts(data.getSSetts(), sset, [nsuppL, nsuppR, missL, missR])
                else:
                    rsParts = SParts(data.getSSetts(), sset, [nsuppL, nsuppR])

                self.restricted_sets[sid] = {"sParts": rsParts,
                                             "prs": [self.queries[0].proba(0, data, sset),
                                                     self.queries[1].proba(1, data, sset)],
                                             "rids": set(sset)}
                resets_ids.append(sid)

        if len(resets_ids) > 0:
            self.dict_supp_info = None
            return True
        return False

    def getRSet(self, details={}):
        if type(details) is dict:
            rset_id = details.get("rset_id")
        else:
            rset_id = details
        if rset_id is not None:
            if rset_id == "cond" and self.hasCondition():
                return {"sParts": self.getSupportsC()}
            elif rset_id in self.restricted_sets:
                return self.restricted_sets[rset_id]
            # elif rset_id == "all" or rset_id == ACTIVE_RSET_ID:
            else:
                return {"sParts": self.supports()}
            return None
        # elif ACTIVE_RSET_ID in self.restricted_sets:
        #     return self.restricted_sets[ACTIVE_RSET_ID]
        else:
            return {"sParts": self.supports()}

    def getRSetParts(self, details={}):
        rset = self.getRSet(details)
        if rset is not None:
            return rset.get("sParts")

    def getRSetABCD(self, details={}):
        ssp = self.getRSetParts(details)
        if ssp is not None:
            return ssp.get("sParts").getVectorABCD(force_list=True, rest_ids=ssp.get("rids"))

    def getAccRatio(self, details={}):
        if details is not None and (details.get("rset_id_num") in self.restricted_sets
                                    or details.get("rset_id_den") in self.restricted_sets):
            acc_num = self.getRSetParts(details.get("rset_id_num")).acc()
            acc_den = self.getRSetParts(details.get("rset_id_den")).acc()
            return tool_ratio(acc_num, acc_den)
        return 1.

    def getLenRatio(self, details={}):
        if details is not None and (details.get("rset_id_num") in self.restricted_sets
                                    or details.get("rset_id_den") in self.restricted_sets):
            len_num = self.getRSetParts(details.get("rset_id_num")).lenI()
            len_den = self.getRSetParts(details.get("rset_id_den")).lenI()
            return tool_ratio(len_num, len_den)
        return 1.

    def getAcc(self, details={}):
        return self.getRSetParts(details).acc()

    def getPVal(self, details={}):
        return self.getRSetParts(details).pVal()

    def getLenP(self, details={}):
        if "part_id" in details:
            return self.getRSetParts(details).lenP(details["part_id"])
        return -1

    def getLenI(self, details={}):
        return self.getRSetParts(details).lenI()

    def getLenU(self, details={}):
        return self.getRSetParts(details).lenU()

    def getLenL(self, details={}):
        return self.getRSetParts(details).lenL()

    def getLenR(self, details={}):
        return self.getRSetParts(details).lenR()

    def getLenO(self, details={}):
        return self.getRSetParts(details).lenO()

    def getLenN(self, details={}):
        return self.getRSetParts(details).lenN()

    def getLenA(self, details={}):
        return self.getRSetParts(details).lenA()

    def getLenB(self, details={}):
        return self.getRSetParts(details).lenB()

    def getSuppI(self, details={}):
        return self.getRSetParts(details).suppI()

    def getSuppU(self, details={}):
        return self.getRSetParts(details).suppU()

    def getSuppL(self, details={}):
        return self.getRSetParts(details).suppL()

    def getSuppR(self, details={}):
        return self.getRSetParts(details).suppR()

    def getSuppO(self, details={}):
        return self.getRSetParts(details).suppO()

    def getSuppN(self, details={}):
        return self.getRSetParts(details).suppN()

    def getSuppA(self, details={}):
        return self.getRSetParts(details).suppA()

    def getSuppB(self, details={}):
        return self.getRSetParts(details).suppB()

    def getArea(self, which="I", rset_id=None, details={}):
        if which == "L" or which == "R" or which == "A" or which == "B":
            if which == "A":
                side, which_supp = "LHS", "I"
            elif which == "B":
                side, which_supp = "RHS", "I"
            else:
                side, which_supp = which+"HS", which
            lq = self.getQueryProp("len", which="q", rset_id=side, details=details)
            lsupp = self.getSPartsProp("area", which=which_supp, rset_id=rset_id, details=details)
            return lq*lsupp
        if which == "I" or which == "U":
            lqLHS = self.getQueryProp("len", which="q", rset_id="LHS", details=details)
            lqRHS = self.getQueryProp("len", which="q", rset_id="RHS", details=details)
            lsupp = self.getSPartsProp("area", which=which, rset_id=rset_id, details=details)
            return (lqLHS+lqRHS)*lsupp
        if which == "P":
            lqLHS = self.getQueryProp("len", which="q", rset_id="LHS", details=details)
            lqRHS = self.getQueryProp("len", which="q", rset_id="RHS", details=details)
            lsuppLHS = self.getSPartsProp("area", which="L", rset_id=rset_id, details=details)
            lsuppRHS = self.getSPartsProp("area", which="R", rset_id=rset_id, details=details)
            return lqLHS*lsuppLHS+lqRHS*lsuppRHS
        return -1

    def getQueryProp(self, what, which=None, rset_id=None, details={}):
        if Query.hasPropWhat(what) and rset_id in HAND_SIDE:
            if type(HAND_SIDE[rset_id]) is int:
                q = self.query(HAND_SIDE[rset_id])
                if q is not None:
                    dts = {"side": HAND_SIDE[rset_id]}
                    if details is not None:
                        dts.update(details)
                    return q.getProp(what, which, dts)
            else:
                sides = HAND_SIDE[rset_id] or [0, 1]
                xps = self.getQueryProp(what, which, rset_id=sides[0], details=details)
                for side in sides[1:]:
                    if type(xps) is set:
                        xps.update(self.getQueryProp(what, which, rset_id=side, details=details))
                    else:
                        xps += self.getQueryProp(what, which, rset_id=side, details=details)
                return xps

    def getSPartsProp(self, what, which=None, rset_id=None, details={}):
        if SParts.hasPropWhat(what):  # info from supp parts
            rset_parts = self.getRSetParts(rset_id)
            if rset_parts is None:
                return None
            prp = rset_parts.getProp(what, which)
            if what == "supp" or what == "set":
                return mapSuppNames(prp, details)
            return prp

    def getRidsProp(self, what, which=None, rset_id=None, details={}):
        if rset_id is not None and which == self.which_rids:  # ids details for folds subsets
            rset_ids = self.getRestrictedRids(rset_id)
            if rset_ids is None:
                return None
            if what == "len" or what == "card":
                return len(rset_ids)
            elif what == "supp" or what == "set":
                return mapSuppNames(rset_ids, details)
            elif what == "perc":
                return tool_ratio(100.*len(rset_ids), self.nbRows())
            elif what == "ratio":
                return tool_ratio(len(rset_ids), self.nbRows())

    def getProp(self, what, which=None, rset_id=None, details={}):
        if what == "extra":
            return self.getExtra(which, details)
        if what == "area":
            return self.getArea(which, rset_id, details)
        if Query.hasPropWhat(what) and rset_id in HAND_SIDE:
            return self.getQueryProp(what, which, rset_id, details)
        if rset_id is not None and which == self.which_rids:  # ids details for folds subsets
            return self.getRidsProp(what, which, rset_id, details)
        if SParts.hasPropWhat(what):  # info from supp parts
            return self.getSPartsProp(what, which, rset_id, details)
        else:  # other redescription info
            return self.getPropD(what, details=details)

    def getExpProp(self, exp, details={}):
        ws = self.getRP().getPrimitiveWs(exp)
        if ws[0] is not None:
            return self.getProp(ws[0], ws[1], ws[2], details)

    ###############################################
    # METHODS FOR COMPARING TWO REDESCRIPTIONS

    def compare(self, other, details={}):
        if self.score(details) > other.score(details):
            return Redescription.diff_score
        elif self.score(details) == other.score(details):
            cpair = Query.comparePair(self.queries[0], self.queries[1], other.queries[0], other.queries[1])
            if cpair == 0 and (self.hasCondition() or other.hasCondition()):
                if self.hasCondition() and other.hasCondition():
                    return self.getQueryC().compare(other.getQueryC())/10.
                elif other.hasCondition():
                    return -.5
                return .5
            return cpair
        else:
            return -Redescription.diff_score

    def sameRectangles(self, other, details={}):
        return self.supp(0, details) == other.supp(0, details) and self.supp(1, details) == other.supp(1, details) and \
            self.invColsSide(0) == other.invColsSide(0) and self.invColsSide(1) == other.invColsSide(1)
    # def sameRectangles(self, other, details={}):
    #     return self.supp(0, details) == other.supp(0, details) and self.supp(1, details) != other.supp(1, details) and \
    #       self.invColsSide(0) != other.invColsSide(0) and self.invColsSide(1) != other.invColsSide(1)

    def interArea(self, other, side, details={}):
        if other is not None:
            return len(other.supp(side, details) & self.supp(side, details)) * len(other.invColsSide(side) & self.invColsSide(side))
        return 0

    def unionArea(self, other, side, details={}):
        if other is not None:
            return len(other.supp(side, details) | self.supp(side, details)) * len(other.invColsSide(side) | self.invColsSide(side))
        return 0

    def overlapAreaSide(self, other, side, details={}):
        if len(other.invColsSide(side) & self.invColsSide(side)) == 0:
            return 0
        areaU = self.unionArea(other, side, details)
        return tool_ratio(self.interArea(other, side, details), areaU)

    def overlapAreaTotal(self, other, details={}):
        areaUL = self.unionArea(other, 0, details)
        areaUR = self.unionArea(other, 1, details)
        return tool_ratio(self.interArea(other, 0) + self.interArea(other, 1), areaUL+areaUR)

    def overlapAreaMax(self, other, details={}):
        return max(self.overlapAreaSide(other, 0, details), self.overlapAreaSide(other, 1, details))

    def overlapRows(self, other, details={}):
        if other is not None:
            return tool_ratio(len(other.getSuppI(details) & self.getSuppI(details)), min(other.getLenI(details), self.getLenI(details)))
        return 0

    def overlapCols(self, other, side, details={}):
        if other is not None:
            return tool_ratio(len(other.invColsSide(side) & self.invColsSide(side)), min(len(other.invColsSide(side)), len(self.invColsSide(side))))
        return 0
    
    def oneSideIdentical(self, other, details={}):
        return self.queries[0] == other.queries[0] or self.queries[1] == other.queries[1]

    def bothSidesIdentical(self, other, details={}):
        return self.queries[0] == other.queries[0] and self.queries[1] == other.queries[1]

    def equivalent(self, other, details={}):
        return abs(self.compare(other, details)) < Query.diff_balance

    def superceding(self, other, details={}):
        return (self.oneSideIdentical(other) and not self.equivalent(other)) or self.bothSidesIdentical(other)

    def getOverlapPair(self, other, how="ratio", which="I", rset_id=None, details={}):
        rows = False
        cols = False
        if type(how) is tuple:
            rows = how[0] == "rows"
            cols = how[0] == "cols"
            how = how[1]
        if which == "L" or which == "R" or which == "A" or which == "B":
            if which == "A":
                side, which_supp = "LHS", "I"
            elif which == "B":
                side, which_supp = "RHS", "I"
            else:
                side, which_supp = which+"HS", which
            if rows:
                own_sq, other_sq = set([1]), set([1])
            else:
                own_sq = self.getQueryProp("Cset", which="q", rset_id=side, details=details)
                other_sq = other.getQueryProp("Cset", which="q", rset_id=side, details=details)
            if cols:
                own_ssupp, other_ssupp = set([1]), set([1])
            else:
                own_ssupp = self.getSPartsProp("set", which=which_supp, rset_id=rset_id, details=details)
                other_ssupp = other.getSPartsProp("set", which=which_supp, rset_id=rset_id, details=details)
            if how == "inter":
                return len(own_sq & other_sq) * len(own_ssupp & other_ssupp)
            elif how == "union":
                return len(own_sq | other_sq) * len(own_ssupp | other_ssupp)
            elif how == "ratio":
                areaI = len(own_sq & other_sq) * len(own_ssupp & other_ssupp)
                areaU = len(own_sq | other_sq) * len(own_ssupp | other_ssupp)
                return tool_ratio(areaI, areaU)
            elif how == "fractmin":
                areaI = len(own_sq & other_sq) * len(own_ssupp & other_ssupp)
                areaU = min(len(own_sq)*len(own_ssupp), len(other_sq)*len(other_ssupp))
                return tool_ratio(areaI, areaU)

        if which == "I" or which == "U":
            if rows:
                own_sqLHS, own_sqRHS, other_sqLHS, other_sqRHS, = set([1]), set([1]), set([1]), set([1])
            else:
                own_sqLHS = self.getQueryProp("Cset", which="q", rset_id="LHS", details=details)
                own_sqRHS = self.getQueryProp("Cset", which="q", rset_id="RHS", details=details)
                other_sqLHS = other.getQueryProp("Cset", which="q", rset_id="LHS", details=details)
                other_sqRHS = other.getQueryProp("Cset", which="q", rset_id="RHS", details=details)
            if cols:
                own_ssupp, other_ssupp = set([1]), set([1])
            else:
                own_ssupp = self.getSPartsProp("set", which=which, rset_id=rset_id, details=details)
                other_ssupp = other.getSPartsProp("set", which=which, rset_id=rset_id, details=details)
            if how == "inter":
                return (len(own_sqLHS & other_sqLHS)+len(own_sqRHS & other_sqRHS))*len(own_ssupp & other_ssupp)
            elif how == "union":
                return (len(own_sqLHS | other_sqLHS)+len(own_sqRHS | other_sqRHS))*len(own_ssupp | other_ssupp)
            elif how == "ratio":
                areaI = (len(own_sqLHS & other_sqLHS)+len(own_sqRHS & other_sqRHS))*len(own_ssupp & other_ssupp)
                areaU = (len(own_sqLHS | other_sqLHS)+len(own_sqRHS | other_sqRHS))*len(own_ssupp | other_ssupp)
                return tool_ratio(areaI, areaU)
            elif how == "fractmin":
                areaI = (len(own_sqLHS & other_sqLHS)+len(own_sqRHS & other_sqRHS))*len(own_ssupp & other_ssupp)
                areaU = min((len(own_sqLHS)+len(own_sqRHS))*len(own_ssupp),
                            (len(other_sqLHS)+len(other_sqRHS))*len(other_ssupp))
                return tool_ratio(areaI, areaU)

        if which == "P":
            if rows:
                own_sqLHS, own_sqRHS, other_sqLHS, other_sqRHS = set([1]), set([1]), set([1]), set([1])
            else:
                own_sqLHS = self.getQueryProp("Cset", which="q", rset_id="LHS", details=details)
                own_sqRHS = self.getQueryProp("Cset", which="q", rset_id="RHS", details=details)
                other_sqLHS = other.getQueryProp("Cset", which="q", rset_id="LHS", details=details)
                other_sqRHS = other.getQueryProp("Cset", which="q", rset_id="RHS", details=details)
            if cols:
                own_ssuppLHS, own_ssuppRHS, other_ssuppLHS, other_ssuppRHS = set([1]), set([1]), set([1]), set([1])
            else:
                own_ssuppLHS = self.getSPartsProp("set", which="L", rset_id=rset_id, details=details)
                own_ssuppRHS = self.getSPartsProp("set", which="R", rset_id=rset_id, details=details)
                other_ssuppLHS = other.getSPartsProp("set", which="L", rset_id=rset_id, details=details)
                other_ssuppRHS = other.getSPartsProp("set", which="R", rset_id=rset_id, details=details)
            if how == "inter":
                return len(own_sqLHS & other_sqLHS)*len(own_ssuppLHS & other_ssuppLHS) + \
                    len(own_sqRHS & other_sqRHS)*len(own_ssuppRHS & other_ssuppRHS)
            elif how == "union":
                return len(own_sqLHS | other_sqLHS)*len(own_ssuppLHS | other_ssuppLHS) + \
                    len(own_sqRHS | other_sqRHS)*len(own_ssuppRHS | other_ssuppRHS)
            elif how == "ratio":
                areaI = len(own_sqLHS & other_sqLHS)*len(own_ssuppLHS & other_ssuppLHS) + \
                    len(own_sqRHS & other_sqRHS)*len(own_ssuppRHS & other_ssuppRHS)
                areaU = len(own_sqLHS | other_sqLHS)*len(own_ssuppLHS | other_ssuppLHS) + \
                    len(own_sqRHS | other_sqRHS)*len(own_ssuppRHS | other_ssuppRHS)
                return tool_ratio(areaI, areaU)
            elif how == "fractmin":
                areaI = len(own_sqLHS & other_sqLHS)*len(own_ssuppLHS & other_ssuppLHS) + \
                    len(own_sqRHS & other_sqRHS)*len(own_ssuppRHS & other_ssuppRHS)
                areaU = min(len(own_sqLHS)*len(own_ssuppLHS)+len(own_sqRHS)*len(own_ssuppRHS),
                            len(other_sqLHS)*len(other_ssuppLHS)+len(other_sqRHS)*len(other_ssuppRHS))
                return tool_ratio(areaI, areaU)
        return -1

    def getPairProp(self, other, what, which=None, rset_id=None, details={}):
        if what in Redescription.pair_map_how:
            # xps = self.getOverlapPair(other, Redescription.pair_map_how[what], which, rset_id, details)
            # print(">>>", self.disp())
            # print("<<<", other.disp())
            # x = self.overlapRows(other)
            # print(what, which, rset_id, xps, x)
            # return xps
            return self.getOverlapPair(other, Redescription.pair_map_how[what], which, rset_id, details)

        if rset_id is not None and rset_id != "all":
            if details is None:
                details = {"rset_id": rset_id}
            else:
                details.update({"rset_id": rset_id})

        if what in Redescription.info_pair_what_dets:  # other redescription info
            methode = eval(Redescription.info_pair_what_dets[what])
            if callable(methode):
                return methode(other, details)
        elif what in Redescription.info_pair_what_side and which in Redescription.pair_map_side:  # other redescription info
            methode = eval(Redescription.info_pair_what_side[what])
            if callable(methode):
                return methode(other, Redescription.pair_map_side[which], details)

    def getExpPairProp(self, other, exp, details={}):
        ws = self.getRP().getPairPrimitiveWs(exp)
        if ws[0] is not None:
            return self.getPairProp(other, ws[0], ws[1], ws[2], details)
   ###############################################

# PRINTING AND PARSING METHODS
    # FROM HERE ALL PRINTING AND READING
    def red2strLegacy(self):
        return '%s terms:\t (%i): %s\t%s\t%s\t%s' % (self.getNbAvC({"format": "str"}), len(self), self.dispQueries(sep=" "), self.dispStats(sep=" "), self.dispSuppL(sep=" "), self.getTrack({"format": "str"}))

    def __str__(self):
        return "%s%s\t%s\t%s\t%s\tlengthQs:%d nbAvC:%s track:%s" % (self.class_letter, self.getUid(), self.dispQueries(sep="\t"), self.dispStats(sep=" "), self.dispSuppL(sep=" "), len(self), self.getNbAvC({"format": "str"}), self.getTrack({"format": "str"}))

    # def __str__(self):
    #     return "%s\t%s\t%s" % (self.query(0), self.query(1), self.score())

    def dispQueries(self, names=[None, None, None], sep='\t', fmts=[None, None, None]):
        sides = [0, 1]
        if self.hasCondition():
            sides.append(-1)
        return sep.join(["q%s:%s" % (side, self.query(side).disp(names=names[side], fmts=fmts[side])) for side in sides])

    def dispStats(self, sep="\t"):
        supp_str = ""
        if self.supports() is not None:
            supp_str += self.supports().dispStats(sep)
        if self.hasCondition():
            supp_str += sep+"COND:"
            if self.getSupportsC() is not None:
                supp_str += self.getSupportsC().dispStats(sep)
        return supp_str

    def dispSupp(self, sep="\t"):
        supp_str = ""
        if self.supports() is not None:
            supp_str += self.supports().dispSupp(sep)
        if self.hasCondition():
            supp_str += sep+"COND:"
            if self.getSupportsC() is not None:
                supp_str += self.getSupportsC().dispSupp(sep)
        return supp_str

    def dispSuppL(self, sep="\t"):
        supp_str = ""
        if self.supports() is not None:
            supp_str += self.supports().dispSuppL(sep)
        if self.hasCondition():
            supp_str += sep+"COND:"
            if self.getSupportsC() is not None:
                supp_str += self.getSupportsC().dispSuppL(sep)
        return supp_str

    def prepareQuery(self, side, details={}, named=False):
        q = self.query(side)
        if q is None:
            return ""
        style = details.get("style")
        names = None
        fmts = None
        if details is not None:
            if "names" in details:
                names = details["names"][side]
            if "fmts" in details:
                fmts = details["fmts"][side]
        if style is None:
            return q.disp(names=names, fmts=fmts)
        else:
            return q.styledDisp(names=names, fmts=fmts, style=style)

    def prepareQueryRHS(self, details):
        return self.prepareQuery(1, details)

    def prepareQueryLHS(self, details):
        return self.prepareQuery(0, details)

    def prepareQueryCOND(self, details):
        return self.prepareQuery(-1, details)

    def disp(self, names=[None, None, None], row_names=None, with_fname=False, rid="", nblines=1, delim="", last_one=False, list_fields="basic", modifiers={}, style=None, sep=None, rp=None, fmts=[None, None, None]):
        return self.getRP(rp).disp(self, names, row_names, with_fname, rid, nblines, delim, last_one, list_fields, modifiers, style, sep, fmts)

    @classmethod
    def parse(tcl, stringQueries, stringSupp=None, data=None, rp=None, sid=None):
        if data is not None and data.hasNames():
            names = data.getNames()
        else:
            names = [None, None]
        (queryL, queryR, lpartsList) = tcl.getRP(rp).parseQueries(stringQueries, names=names)
        supportsS = None
        if data is not None and stringSupp is not None and type(stringSupp) == str and re.search('\t', stringSupp):
            supportsS = SParts.parseSupport(stringSupp, data.nbRows(), data.getSSetts())
        if sid is not None:
            lpartsList["sid"] = sid
        return tcl.initParsed(queryL, queryR, lpartsList, data, supportsS)

    @classmethod
    def prepareRid(tcl, rid=None, sid=None):
        if type(rid) is dict:
            sid = rid.get("sid")
            if "all:rid:" in rid:
                rid = rid.get("all:rid:")
            elif "rid" in rid:
                rid = rid.get("rid")
            else:
                rid = None
        if rid is None:
            return tcl.generateNextUid()

        try:
            rid = rid.lstrip(tcl.class_letter)
            rid = int(rid)
        except (AttributeError, TypeError, ValueError):
            pass
        if sid is not None:
            rid = "%s.%s" % (rid, sid)
            parts = rid.split(".")
            if len(parts) > 2 and parts[-1] == parts[-2]:
                rid = ".".join(parts[:-1])
            try:
                rid = float(rid)
            except (TypeError, ValueError):
                pass
        return rid

    @classmethod
    def initParsed(tcl, queryL, queryR, lpartsList={}, data=None, supportsS=None):
        status_enabled = None
        iid = tcl.prepareRid(lpartsList)
        if "all:extra:status" in lpartsList:
            status_enabled = int(lpartsList.pop("all:extra:status"))
        if "status_enabled" in lpartsList:
            status_enabled = int(lpartsList.pop("status_enabled"))
        r = None
        if supportsS is not None:
            r = tcl(queryL, queryR, supportsS.supparts(), data.nbRows(), [set(), set()], [queryL.proba(0, data), queryR.proba(1, data)], data.getSSetts(), iid=iid)

            for key, v in lpartsList.items():
                tv = RedProps.getEVal(r, key)
                if tv != v:
                    raise Warning("Something wrong in the supports ! (%s: %s ~ %s)\n" % (key, v, tv))

        if r is None:
            r = tcl(queryL, queryR, iid=iid)
            if data is not None:
                r.recompute(data)
            else:
                r.cache_evals = lpartsList

        if "queryCOND" in lpartsList:
            qC = lpartsList["queryCOND"]
            supp_cond = None
            if data is not None:
                supp_cond, miss_cond = qC.recompute(-1, data)
            r.setCondition(qC, supp_cond)
        if r is not None and status_enabled is not None:
            r.setStatus(status_enabled)
        return r


if __name__ == '__main__':
    # print(Redescription.exp_details.keys())
    from classData import Data
    import sys

    rep = "/home/egalbrun/short/vaalikone_FILES/"
    rep = "/home/egalbrun/short/raja_time/"
    data = Data([rep+"data_LHS_lngtid.csv", rep+"data_RHS_lngtid.csv", {}, "NA"], "csv")

    rp = Redescription.getRP()
    # filename = rep+"redescriptions.csv"
    filename = rep+"tid_test.queries"
    filep = open(filename, mode='r')
    reds, _ = rp.parseRedList(filep, data)
    with open("/home/egalbrun/short/tmp_queries.txt", mode='w') as fp:
        # pdb.set_trace()
        fp.write(rp.printRedList(reds, missing=True))
        ## fp.write(rp.printTexRedList(reds, names = [data.getNames(0), data.getNames(1)], nblines=3, standalone=True))
        ## fields=[-1, "CUST:XX=q0:containsC:0", "Lnb_queryLHS", "Lset_queryLHS", "Lnb_queryRHS", "Lset_queryRHS", "containsAND_queryRHS", "containsOR_queryRHS"]
    exit()
    # rep = "/home/galbrun/TKTL/redescriptors/data/vaalikone/"
    # data = Data([rep+"vaalikone_profiles_test.csv", rep+"vaalikone_questions_test.csv", {}, "NA"], "csv")

    # reds = []
    # with open("../../bazar/queries.txt") as f:
    #     for line in f:
    #         if len(line.strip().split("\t")) >= 2:
    #             try:
    #                 tmpLHS = Query.parse(line.strip().split("\t")[0], data.getNames(0))
    #                 tmpRHS = Query.parse(line.strip().split("\t")[1], data.getNames(1))
    #             except:
    #                 continue
    #             r = Redescription.fromQueriesPair([tmpLHS, tmpRHS], data)
    #             reds.append(r)

    # with open("../../bazar/queries_list2.txt", mode='w') as f:
    #     f.write(rp.printRedList(reds))

    # with open("../../bazar/queries_list2.txt", mode='r') as f:
    #     reds, _ = rp.parseRedList(f, data)

    # for red in reds:
    #     print(red.disp())
