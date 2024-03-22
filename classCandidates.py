import re
import os.path

try:
    from classProps import WithProps
    from classSParts import SSizes
    from classQuery import Op, Term, Literal, Query
    from classRedescription import Redescription
except ModuleNotFoundError:
    from .classProps import WithProps
    from .classSParts import SSizes
    from .classQuery import Op, Term, Literal, Query
    from .classRedescription import Redescription

import pdb

LIT_SEP = " AND "
# grep "class .*:\|def .*:" classCandidates.py > functions_classCandidates.py


class ExtensionWarning(Warning):
    def __init__(self, value):
        self.value = value

    def __str__(self):
        return repr(self.value)


def clpToStr(clp, preff="", suff=""):  # return empty string to skip
    if clp is not None:
        return preff+str([list(clp[i]) for i in [0, -1]])+suff
        # return preff+str([list(c) for c in clp])+suff
    # return ""


def clpToRepr(clp):
    if clp is not None:
        return str([list(c) for c in clp])


#######################################
# SORTING METHODS FOR INITIAL PAIRS
DEF_SCORE = float("-inf")
# drop_set is assumed to be a subset of all_keys,
# so "drop_set.symmetric_difference(all_keys)" is in effect the same as set(all_keys).difference(drop_set)


def fifo_sort(all_keys, store, drop_set=set()):
    return sorted(drop_set.symmetric_difference(all_keys), reverse=True)


def filo_sort(all_keys, store, drop_set=set()):
    return sorted(drop_set.symmetric_difference(all_keys))


def overall_sort(all_keys, store, drop_set=set()):
    return sorted(drop_set.symmetric_difference(all_keys), key=lambda x: store.cmpKeyCand(store[x]))


def alternate_sort(all_keys, store, drop_set=set()):
    def_cv = store.cmpKeyCand(None)
    def_rank = -1  # len(all_keys)+1
    best_sides = [{}, {}]
    loc_data = {}
    for k in all_keys:
        loc_data[k] = {"tail_brk": store.cmpKeyCand(store[k])}
        for side in [0, 1]:
            col = store[k].getCid(side)
            if col is not None and col >= 0:
                if col in best_sides[side]:
                    best_sides[side][col].append(k)
                else:
                    best_sides[side][col] = [k]

    for side in [0, 1]:
        for col in best_sides[side]:
            best_sides[side][col].sort(key=lambda x: loc_data[x]["tail_brk"], reverse=True)
        for c, vs in best_sides[side].items():
            for pp, v in enumerate(vs):
                loc_data[v]["rank_%d" % side] = pp
    sort_dt = sorted([(loc_data[k].get("tail_brk", def_cv) - max(loc_data[k].get("rank_0", def_rank), loc_data[k].get("rank_1", def_rank)), k)
                      for k in drop_set.symmetric_difference(all_keys)])
    return [k[-1] for k in sort_dt]


SORT_METHODS = {"overall": overall_sort,
                "alternate": alternate_sort,
                "fifo": fifo_sort,
                "filo": filo_sort}
DEFAULT_METHOD = overall_sort


###############################################
# Generic Candidate class and Stores where to keep them

class Candidate(WithProps):
    what = "Candidate"
    info_what_extra = set()

    def __init__(self):
        WithProps.__init__(self)

    def getPropD(self, what, details={}, default=None):
        if what in self.info_what_extra and what in self.extras:
            return self.extras[what]
        v = WithProps.getPropD(self, what, details, default)
        if what in self.info_what_extra:
            self.extras[what] = v
        return v

    def getLits(self):
        if self.getSide() == 1:
            return (None, self.lit)
        return (self.lit, None)

    def getLit(self, side=None):
        if side is None or side == self.getSide():
            return self.lit

    def dispLit(self, side=None):
        lit = self.getLit(side)
        if lit is None:
            return "-"
        elif type(lit) is list:
            return LIT_SEP.join(["%s" % ll for ll in lit])
        return "%s" % lit

    def getCid(self, side=None):
        lit = self.getLit(side)
        if isinstance(lit, Literal):
            return lit.colId()

    def getAcc(self):
        return -1

    def isNeg(self, side=None):
        lit = self.getLit(side)
        if isinstance(lit, Literal):
            return self.getLit().isNeg()
        return False

    # condition

    def getCondition(self):
        return None

    def hasCondition(self):
        return self.getCondition() is not None

    def getCidC(self):
        if self.hasCondition():
            return self.getCondition().getCid()

    def getLitC(self):
        if self.hasCondition():
            return self.getCondition().getLit()

    def dispLitC(self):
        if self.hasCondition():
            return self.getCondition().dispLit()

    def getAccC(self):
        if self.hasCondition():
            return self.getCondition().getAcc()
        return self.getAcc()

    def mkRed(self, data, red=None):
        return Redescription.fromInitialPair(self.getLits(), data, self.getLitC())

    def checkRed(self, red, thres=0.0001):
        sp_KO = False
        msg = ""
        if self.getSSizes() is not None:
            org_sp = str([l for l in self.getSSizes().lparts()])
            new_sp = str([l for l in red.supports().lparts()])
            msg += "spart %s %s %s" % (org_sp, new_sp, (org_sp != new_sp)*"!!!")
            sp_KO = org_sp != new_sp
            if sp_KO:
                print(self)
                print(msg)
                pdb.set_trace()

        acc_KO = (red.getAcc() - self.getAcc())**2 > thres
        accC_KO = red.hasCondition() and (red.getAcc("cond") - self.getAccC())**2 > thres
        if sp_KO or acc_KO or accC_KO:
            msg += "\nacc     %s %s %s" % (red.getAcc(), self.getAcc(), acc_KO*"!!!")
            if red.hasCondition():
                msg += "\nacc cond %s %s %s" % (red.getAcc("cond"), self.getAccC(), accC_KO*"!!!")
            return msg

    # comparisons
    def toKey(self):
        return None

    def __eq__(self, other):
        return isinstance(other, self.__class__) and self.toKey() == other.toKey()

    def __ne__(self, other):
        return not isinstance(other, self.__class__) or self.toKey() != other.toKey()

    def __lt__(self, other):
        return isinstance(other, self.__class__) and self.toKey() < other.toKey()

    def __le__(self, other):
        return isinstance(other, self.__class__) and self.toKey() <= other.toKey()

    def __gt__(self, other):
        return not isinstance(other, self.__class__) or self.toKey() > other.toKey()

    def __ge__(self, other):
        return not isinstance(other, self.__class__) or self.toKey() >= other.toKey()


###############################################
# Extensions candidates
class ExtCand(Candidate):

    what = "Extension"
    info_what = {"exp": "self.dispExp()",
                 "lit": "self.getLit()", "cid": "self.getCid()",
                 "acc": "self.getAcc()", "supp_ratio": "self.getSuppRatio()",
                 "toBlue": "self.getVarBlue()", "toRed": "self.getVarRed()"}
    info_what_args = {"score": "self.getScore", "pValQ": "self.pValQuery", "pValR": "self.pValRed"}
    info_what_extra = set(["score", "pValQ", "pValR"])

    def __init__(self, ssetts, adv, clp, lit_parts, condition=None):
        Candidate.__init__(self)
        # adv is a tuple: acc, var_num, var_den, contri, fix_num, fix_den
        # clp is a triple [lparts,[ lmiss,] lout, lin], where each is a list of counts in support sets
        self.condition = condition
        self.adv = adv
        # side, op, tmp_neg, lit = lit_parts
        self.setLit(lit_parts[-1], lit_parts[0])
        if len(lit_parts) == 4:
            self.op = lit_parts[1]
            self.through_neg = lit_parts[2]
        else:
            self.op = False
            self.through_neg = False
        self.clp = clp
        self.ssetts = ssetts
        # self.setClpSSizes(clp, ssetts)

    def setClpSSizes(self, clp=None, ssetts=None, prs=[-1, -1]):
        if not hasattr(self, "clp") or clp is not None:
            self.clp = clp
        self.ssizes = None
        if ssetts is not None and self.clp is not None:
            N, lsupports = ssetts.clpToLSupports(self.clp, side=self.getSide(), op=self.getOp())
            self.ssizes = SSizes(ssetts, N, lsupports, prs)

    def getSSizes(self):
        return None  # self.ssizes

    def getSSetts(self):
        return self.ssetts
        # if self.ssizes is not None:
        #     return self.ssizes.getSSetts()

    def toKey(self):
        return tuple(self.adv)+(self.side, self.op)+self.lit.toKey()

    def setLit(self, lit, side=0):
        self.lit = lit
        self.side = side

    def getOp(self):
        return self.op

    def getSide(self):
        return self.side

    def getAcc(self):
        return self.adv[0]

    def getVarBlue(self):
        return self.adv[1]

    def getVarRed(self):
        return self.adv[2]

    def getClp(self):
        return self.clp

    def setCondition(self, condition):
        self.condition = condition

    def getCondition(self):
        return self.condition

    def mkRed(self, data, red=None):
        if red is not None:
            supp = data.supp(self.getSide(), self.getLit())
            miss = data.miss(self.getSide(), self.getLit())
            red_ext = red.mkRed(data, self.getSide(), self.getOp(), self.getLit(), supp, miss)
            if self.hasCondition():
                litC = self.getLitC()
                if type(litC) is list:
                    # if len(litC) > 1: pdb.set_trace()
                    qC = Query(OR=False, buk=litC)
                else:
                    qC = Query(buk=[litC])
                supp_cond, miss_cond = qC.recompute(-1, data)
                red_ext.setCondition(qC, supp_cond)
            return red_ext

    def __repr__(self):
        lit = repr(self.getLit())
        clp = clpToRepr(self.clp)
        condition = ""
        if self.hasCondition():
            condition = ", condition=" + repr(self.getCondition())
        return f"{self.__class__.__name__}(None, {self.adv}, {clp}, ({self.side}, {self.op}, {self.through_neg}, {lit}){condition})"

    def __str__(self):
        if self.getSide() == -1:
            lit = self.dispLit()
        else:
            lit = "(%d, %s, %s)" % (self.getSide(), Op(self.getOp()), self.dispLit())
        if any([c != -1 for c in self.adv[1:]]):
            adv = str(self.adv)
        else:
            adv = "%f" % self.adv[0]
        clp = clpToStr(self.clp, preff="\t")
        tmp = ("%s: %s -> %s%s" %
               (self.what, lit, adv, clp))
        if self.hasCondition():
            tmp += "\n\t"+str(self.getCondition())
        return tmp

    def dispExp(self):
        strPieces = ["", "", ""]
        strPieces[self.getSide()] = "%s %s" % (Op(self.getOp()), self.dispLit())
        if self.hasCondition():
            # score_of = self.getCondition()
            strPieces[2] = self.getCondition().dispLit()
        return "* %20s <==> * %20s %20s" % tuple(strPieces)

    def getScore(self, N, prs=None, coeffs=None, base_acc=0, **kargs):
        return self.getWScore(N, prs, coeffs, base_acc, **kargs)

    def getWScore(self, N, prs=None, coeffs=None, base_acc=0, **kargs):
        # HERE: HOW TO SCORE WITH CONDITION?
        sc = coeffs["impacc"]*self.impacc(base_acc) \
            + coeffs["rel_impacc"]*self.relImpacc(base_acc) \
            + self.pValRedScore(N, prs, coeffs) \
            + self.pValQueryScore(N, prs, coeffs)
        if False:  # self.hasCondition():
            sc += self.getCondition().getScore(N, prs, coeffs, base_acc, **kargs)
        return sc

    def relImpacc(self, base_acc=0):
        if base_acc != 0:
            return (self.adv[0] - base_acc)/base_acc
        else:
            return self.adv[0]

    def impacc(self, base_acc=0):
        return (self.adv[0] - base_acc)

    def pValQueryScore(self, N, prs=None, coeffs=None, **kargs):
        if coeffs is None:
            return self.pValQuery(N, prs)
        elif coeffs["pval_query"] < 0:
            return coeffs["pval_query"] * self.pValQuery(N, prs)
        elif coeffs["pval_query"] > 0:
            return -coeffs["pval_fact"]*(coeffs["pval_query"] < self.pValQuery(N, prs))
        else:
            return 0

    def pValRedScore(self, N, prs, coeffs=None, **kargs):
        if coeffs is None:
            return self.pValRed(N, prs)
        elif coeffs["pval_red"] < 0:
            return coeffs["pval_red"] * self.pValRed(N, prs)
        elif coeffs["pval_red"] > 0:
            return -coeffs["pval_fact"]*(coeffs["pval_red"] < self.pValRed(N, prs))
        else:
            return 0

    def pValQuery(self, N=0, prs=None, **kargs):
        if self.getSSetts() is not None:
            return self.getSSetts().pValQueryCand(self.side, self.op, self.isNeg(), self.clp, N, prs)
        return 0

    def pValRed(self, N=0, prs=None, **kargs):
        if self.getSSetts() is not None:
            return self.getSSetts().pValRedCand(self.side, self.op, self.isNeg(), self.clp, N, prs)
        return 0


class CondCand(ExtCand):
    what = "Condition"

    def setLit(self, lit, side=-1):
        self.lit = lit
        self.side = -1


###############################################
# Initial candidates
class TermCand(Candidate):

    what = "Term"
    info_what = {"exp": "self.dispExp()",
                 "side": "self.getSide()", "lit": "self.getLit()", "cid": "self.getCid()",
                 "supp_len": "self.getSuppLen()", "supp_ratio": "self.getSuppRatio()"}

    def __init__(self, ssetts, side, lit, supp_len=-1, N=0):
        Candidate.__init__(self)
        self.setLit(lit, side)
        self.supp_len = supp_len
        self.N = N

    def setLit(self, lit, side=0):
        self.lit = lit
        self.side = side
        # if isinstance(self.lit, Lit) and self.cid != self.lit.colId():
        #     self.cid = self.lit.colId()

    def getSide(self):
        return self.side

    def getAcc(self):
        return 0

    def getSuppLen(self):
        return self.supp_len

    def getSuppRatio(self):
        if self.N == 0:
            return self.supp_len
        else:
            return self.supp_len/self.N

    def makeTrg(self, data, tcols=None):
        supp = numpy.zeros(data.nbRows(), dtype=bool)
        supp[list(data.col(self.getSide(), self.getCid()).suppLiteral(self.getLit()))] = 1
        if tcols is not None:
            involved = [tcols[self.getSide()].get(data.getMatLitK(self.getSide(), self.getLit(), bincats=True))]
        else:
            involved = []
        return {"side": self.getSide(), "target": supp, "involved": involved, "src": self.getLit()}

    def dispExp(self):
        strPieces = ["", "", ""]
        strPieces[self.getSide()] = self.getLit()
        return "* %20s <==> * %20s %20s" % tuple(strPieces)

    def __repr__(self):
        lit = repr(self.getLit())
        return f"{self.__class__.__name__}(None, {self.side}, {lit}, {self.supp_len}, {self.N})"

    def __str__(self):
        if self.getSide() == 1:
            return "%s: (---, %s) -> %s" % (self.what, self.getLit(), self.getSuppLen())
        else:
            return "%s: (%s, ---) -> %s" % (self.what, self.getLit(), self.getSuppLen())


class PairCand(Candidate):

    what = "Pair"
    info_what = {"exp": "self.dispExp()",
                 "litL": "self.getLit(0)", "litR": "self.getLit(1)",
                 "cidL": "self.getCid(0)", "cidR": "self.getCid(1)",
                 "acc": "self.getAcc()",
                 "toBlue": "self.getVarBlue()", "toRed": "self.getVarRed()"}
    info_what_args = {"score": "self.getScore", "pValR": "self.pValRed"}
    info_what_extra = set(["score", "pValR"])

    def __init__(self, ssetts, adv, clp, litL, litR, negL=False, negR=False, swap=False, condition=None):
        Candidate.__init__(self)
        self.swapped = swap
        if swap:
            self.lits = [litR, litL]
            self.thrg_neg = [negR, negL]
        else:
            self.lits = [litL, litR]
            self.thrg_neg = [negL, negR]
        if type(adv) is float:
            self.acc = adv
        else:
            self.acc = adv[0]
        self.condition = condition
        self.setClpSSizes(clp, ssetts)

    def setClpSSizes(self, clp=None, ssetts=None, prs=None):
        if not hasattr(self, "clp") or clp is not None:
            self.clp = clp
        self.ssizes = None
        if ssetts is not None and self.clp is not None:
            N, lsupports = ssetts.clpToLSupports(self.clp, pair=True, swap=self.swapped)
            self.ssizes = SSizes(ssetts, N, lsupports, prs)

    def getSSizes(self):
        return self.ssizes

    def getSSetts(self):
        if self.ssizes is not None:
            return self.ssizes.getSSetts()

    def getClp(self):
        return self.clp

    def getAcc(self):
        return self.acc

    def setLit(self, lit, side=0):
        self.lits[side] = lit

    def getLits(self):
        return self.lits

    def getLit(self, side=0):
        return self.lits[side]

    def setCondition(self, condition):
        self.condition = condition

    def getCondition(self):
        return self.condition

    # def swapSides(self):
    #     self.lits = [self.lits[1], self.lits[0]]
    #     self.thrg_neg = [self.thrg_neg[1], self.thrg_neg[0]]

    def fillCurrentElement(self, element=None):
        if element is not None:
            if len(element) > 1 and isinstance(element[0], Term):
                if type(self.lits[element[1]]) is bool:
                    self.setLit(Literal(self.lits[element[1]], element[0]), element[1])
                else:
                    self.setLit(Literal(False, element[0]), element[1])

    def doublyNegatedLit(self, side):
        return isinstance(self.lits[side], Literal) and self.lits[side].isNeg() != self.thrg_neg[side]

    def dispExp(self):
        strPieces = [self.dispLit(0), self.dispLit(1), ""]
        if self.hasCondition():
            # score_of = self.getCondition()
            strPieces[2] = self.getCondition().dispLit()
        return "* %20s <==> * %20s %20s" % tuple(strPieces)

    def getScore(self, N, prs=None, coeffs=None, base_acc=0, **kargs):
        return self.getWScore(N, prs, coeffs, base_acc, **kargs)

    def getWScore(self, N, prs=None, coeffs=None, base_acc=0, **kargs):
        coeffs = kargs.get("coeffs")
        if coeffs is None:
            return self.getAcc()
        # HERE: HOW TO SCORE WITH CONDITION?
        sc = max(coeffs["impacc"], coeffs["rel_impacc"])*self.getAcc() + self.pValRedScore(**kargs)
        # + coeffs["rel_impacc"]*self.relImpacc(base_acc)
        # + self.pValQueryScore(N, prs, coeffs)
        if False:  # self.hasCondition():
            sc += self.getCondition().getScore(**kargs)
        return sc

    def pValRedScore(self, **kargs):
        coeffs = kargs.get("coeffs")
        if coeffs is None:
            return self.pValRed(**kargs)
        elif coeffs["pval_red"] < 0:
            return coeffs["pval_red"] * self.pValRed(**kargs)
        elif coeffs["pval_red"] > 0:
            return -coeffs["pval_fact"]*(coeffs["pval_red"] < self.pValRed(**kargs))
        else:
            return 0

    def pValRed(self, **kargs):
        if self.getSSizes() is not None:
            return self.getSSizes().pVal()
        return -1

    def __repr__(self):
        litL = repr(self.getLit(0))
        litR = repr(self.getLit(1))
        clp = clpToRepr(self.clp)
        condition = ""
        if self.hasCondition():
            condition = ", condition="+repr(self.getCondition())
        return f"{self.__class__.__name__}(None, {self.acc}, {clp}, {litL}, {litR}, {self.thrg_neg[0]}, {self.thrg_neg[1]}{condition})"

    def __str__(self):
        condition = ""
        if self.hasCondition():
            condition = "\t+ " + str(self.getCondition())
        clp = clpToStr(self.clp, preff="\t")
        # return "%s: (%s, %s) -> %s%s%s" % (self.what, self.getLit(0), self.getLit(1), self.getAcc(), clp, condition)
        return "%s: (%s, %s) -> %s\t%s%s%s" % (self.what, self.getLit(0), self.getLit(1), self.getAcc(), self.pValRed(), clp, condition)

    @classmethod
    def pattCand(tcl):
        ip_patt = tcl.what+":\s+\((?P<litL>.*), (?P<litR>.*)\)\s+->\s+(?P<acc>\d\.\d+)"
        cond_patt = "\t\+ Condition:\s+(?P<litC>.*)\s+->\s+(?P<accC>[\(\)0-9\., ]+)"
        return ip_patt+"(\t(?P<clp>[\[\]0-9, ]*))?"+"(?P<cond>"+cond_patt+"(\t(?P<clpC>[\[\]0-9, ]*))?)?"


###############################################
# Storage of candidates
class DictStore(dict):

    def __init__(self, d={}, track_changes=True):
        dict.__init__(self, d)
        self.change_keys = set(self.keys())
        self.track_changes = track_changes
        self.current_elem = None

    def setCurrentElement(self, elem):
        self.current_elem = elem

    def getCurrentElement(self):
        return self.current_elem

    def clear(self):
        dict.clear(self)
        self.change_keys = set()
        self.current_elem = None

    def setTracking(self, track_changes=False):
        self.track_changes = track_changes

    def ready_items(self):
        if self.track_changes:
            return [(k, v) for (k, v) in self.items() if k not in self.change_keys]
        return self.items()

    def ready_keys(self):
        if self.track_changes:
            return [k for k in self.keys() if k not in self.change_keys]
        return self.keys()

    def __setitem__(self, key, value):
        try:
            dict.__setitem__(self, key, value)
        except:
            raise
        else:
            self.change_keys.add(key)

    def __delitem__(self, key):
        try:
            dict.__delitem__(self, key)
        except:
            raise
        else:
            self.change_keys.discard(key)

    def toNextRound(self):
        self.change_keys = set()
        self.current_elem = None

    def currentRoundKeys(self):
        return self.change_keys

    def currentRoundIter(self):
        for k in self.change_keys:
            yield self[k]

    def nextPos(self, cand):
        return None

    def add(self, cand):
        if self.keepCand(cand):
            pos = self.nextPos(cand)
            if pos is not None and ((pos not in self) or (self.cmpKeyCand(cand) > self.cmpKeyCand(self[pos]))):
                self[pos] = cand
                return pos


class ListStore(DictStore):
    # this is actually a dictionary with incremental keys

    def __init__(self, d=[]):
        DictStore.__init__(self, enumerate(d))
        self.next_pos = len(self)

    def nextPos(self, cand):
        return self.next_pos

    def add(self, cand):
        if self.keepCand(cand):
            pos = self.nextPos(cand)
            if pos is not None:
                self[pos] = cand
                self.next_pos += 1
                return pos


class XStore(object):

    @classmethod
    def getStoredClass(tcl):
        return Candidate

    @classmethod
    def getStoredClassName(tcl):
        return tcl.getStoredClass().__name__

    @classmethod
    def getStoredClassWhat(tcl):
        return tcl.getStoredClass().what

    # operates on a candidate store that provides items method
    # only store the necessary constraints
    def __init__(self, N=0, ssetts=None, constraints=None, current=None):
        self.N = N
        self.ssetts = ssetts
        self.initConstraints(constraints)
        self.setCurrent(current)

    def initConstraints(self, constraints):
        self.constraints = {}

    def setCurrent(self, current=None):
        self.c_vals = {}

    def copyCurrent(self):
        return self.c_vals

    def emptyCopy(self):
        return self.__class__(self.N, self.ssetts, self.constraints, self.copyCurrent())

    def getSSetts(self):
        return self.ssetts

    # Match addPair/addExt to corresponding init

    def mkCondCand(self, *cand_pieces):
        return CondCand(self.getSSetts(), *cand_pieces[1:])  # Drop extra nego args, needed for pairs, not extensions

    def mkExtCand(self, *cand_pieces):
        return ExtCand(self.getSSetts(), *cand_pieces[1:])  # Drop extra nego args, needed for pairs, not extensions

    def mkRedFromCand(self, cand, data):
        return cand.mkRed(data)

    def getFoundReds(self, data):
        # checking  >= self.constraints["init_minscore"] is done when adding the candidate
        # red.restrictAvailable(side, self.c_vals["avC"][side]) or setFull is done by mkRedFromCand
        return [self.mkRedFromCand(cand, data) for (pos, cand) in self.items()]

    def getFoundCands(self):
        return [cand for (pos, cand) in self.items()]

    # For evaluation and comparison of candidates

    def prepDetails(self):
        return {"N": self.N, "coeffs": self.constraints.get("score_coeffs")}

    def getCmpCriterion(self):
        return None

    def getCmpDefault(self):
        return DEF_SCORE

    def getKeepCriterion(self):
        return None

    def getKeepMinThres(self):
        return 0

    def keepCand(self, cand):
        if isinstance(cand, self.getStoredClass()):
            if self.getKeepCriterion() is None:  # unconditional
                return True
            kv = cand.getPropD(self.getKeepCriterion(), details=self.prepDetails(), default=None)
            return kv is not None and kv >= self.getKeepMinThres()
        return False

    def cmpKeyCand(self, cand):
        if isinstance(cand, self.getStoredClass()):
            return cand.getPropD(self.getCmpCriterion(), details=self.prepDetails(), default=self.getCmpDefault())
        return self.getCmpDefault()

    def betterCand(self, cand, ccmp):
        if self.getCmpCriterion() is None:  # unconditional
            return True
        return self.cmpKeyCand(cand) > self.cmpKeyCand(ccmp)

    # For printing content details
    @classmethod
    def getDispProps(tcl):
        return []

    def dispSpec(self):
        return ""

    def dispHeader(self, sep="\t"):
        return sep.join([dsp_def % lbl for (lbl, prop, fmt, dsp_def) in self.getDispProps()])

    def dispCand(self, cand, sep="\t"):
        # if isinstance(cand, Candidate):
        details = self.prepDetails()
        pieces = []
        for (lbl, prop, fmt, dsp_def) in self.getDispProps():
            v = cand.getPropD(prop, details=details, default=None)
            if v is None:
                pieces.append(dsp_def % "")
            else:
                pieces.append(fmt % v)
        return sep.join(pieces)

    def __str__(self):
        dsp = "%s:%s\n" % (self.__class__.__name__, self.dispSpec())
        dsp += self.dispHeader()
        for k, cand in self.items():
            dsp += "\n\t%s" % self.dispCand(cand)
        return dsp


###############################################
# Storage of extensions candidates
class ExtsStore(XStore):

    @classmethod
    def getStoredClass(tcl):
        return ExtCand

    def initConstraints(self, constraints):
        self.constraints = {}
        if type(constraints) is dict:  # probably copying
            self.constraints["score_coeffs"] = constraints.get("score_coeffs")
            self.constraints["min_impr"] = constraints.get("min_impr")
            self.constraints["max_var"] = tuple(constraints.get("max_var"))
        elif constraints is not None:  # an actual constraint instance
            self.constraints["score_coeffs"] = constraints.getCstr("score_coeffs")
            self.constraints["min_impr"] = constraints.getCstr("min_impr")
            self.constraints["max_var"] = (constraints.getCstr("max_var", side=0), constraints.getCstr("max_var", side=1))
        else:
            self.constraints["score_coeffs"] = None
            self.constraints["min_impr"] = 0
            self.constraints["max_var"] = (-1, -1)

    def setCurrent(self, current=None):
        if type(current) is dict:
            self.current = None
            self.c_vals = dict(current)
        elif current is not None:  # an actual redescription
            self.current = current
            self.c_vals = {"base_acc": self.current.getAcc(),
                           "prs": self.current.probas()}
        else:
            self.current = None
            self.c_vals = {"base_acc": 0, "prs": None}

    def copyCurrent(self):
        return self.current

    def mkRedFromCand(self, cand, data):
        red = cand.mkRed(data, self.current)
        red.setFull(self.constraints["max_var"])
        return red

    def prepDetails(self):
        return {"N": self.N, "coeffs": self.constraints["score_coeffs"],
                "prs": self.c_vals["prs"], "base_acc": self.c_vals["base_acc"]}

    def getCmpCriterion(self):
        return "score"

    def getCmpDefault(self):
        return DEF_SCORE

    def getKeepCriterion(self):
        return "score"

    def getKeepMinThres(self):
        return self.constraints["min_impr"]

    @classmethod
    def getDispProps(tcl):
        return [("\t  %20s        %20s        %20s" % ("LHS extension", "RHS extension", "Condition"), "exp", "%s", "%70s"),
                ("score", "score", "%+1.7f", "%8s"), ("Jacc", "acc", "%1.7f", "%8s"),
                ("pVal Q", "pValQ", "%1.7f", "%8s"), ("pVal R", "pValR", "%1.7f", "%8s"),
                ("to Blue", "toBlue", "% 5i", "%5s"), ("to Red", "toRed", "% 5i", "%5s")]

    def dispSpec(self):
        return "\nRedescription: %s" % self.current


class ExtsBatch(DictStore, ExtsStore):

    def __init__(self, N=0, ssetts=None, constraints=None, current=None):
        ExtsStore.__init__(self, N, ssetts, constraints, current)
        DictStore.__init__(self)

    def addPair(self, score, litL, litR, negL=False, negR=False):
        # print("--- ExtsBatch addPair", score, litL, litR)
        # does not store pairs
        pass

    def addExt(self, *cand_pieces):  # nego, adv, clp, lit_parts
        return self.add(self.mkExtCand(*cand_pieces))

    def clear(self, current=None):
        DictStore.clear(self)
        ExtsStore.setCurrent(self, current)

    def nextPos(self, cand):
        if isinstance(cand, ExtCand):
            return (cand.getSide(), cand.getOp())


###############################################
# Storage of initial candidates
class InitsStore(XStore):

    @classmethod
    def getStoredClass(tcl):
        return Candidate

    def initConstraints(self, constraints):
        self.constraints = {}
        if type(constraints) is dict:
            self.constraints["score_coeffs"] = constraints.get("score_coeffs")
            self.constraints["init_minscore"] = constraints.get("init_minscore")
        elif constraints is not None:
            self.constraints["score_coeffs"] = constraints.getCstr("score_coeffs")
            self.constraints["init_minscore"] = constraints.getCstr("init_minscore")
        else:
            self.constraints["score_coeffs"] = None
            self.constraints["init_minscore"] = 0

    def setCurrent(self, current=None):
        if type(current) is dict:  # probably copying
            if current.get("avC") is not None:
                self.c_vals = {"avC": [list(current["avC"][side]) if current["avC"][side] is not None else None for side in [0, 1]]}
            else:
                self.c_vals = None
        elif current is not None:  # an actual redescription
            self.c_vals = {"avC": [current.lAvailableCols[side] for side in [0, 1]]}
        else:
            self.c_vals = None

    def mkRedFromCand(self, cand, data):
        red = cand.mkRed(data)
        if self.c_vals is not None:
            for side in [0, 1]:
                red.restrictAvailable(side, self.c_vals["avC"][side])
        return red

    def dispSpec(self):
        return " (%s >= %s)" % (self.getKeepCriterion(), self.getKeepMinThres())


class PairsBatch(ListStore, InitsStore):

    @classmethod
    def getStoredClass(tcl):
        return PairCand

    def __init__(self, N=0, ssetts=None, constraints=None, current=None):
        InitsStore.__init__(self, N, ssetts, constraints, current)
        ListStore.__init__(self)

    def getCmpCriterion(self):
        return "score"

    def getCmpDefault(self):
        return DEF_SCORE

    def getKeepCriterion(self):
        return "score"

    def getKeepMinThres(self):
        return self.constraints["init_minscore"]

    @classmethod
    def getDispProps(tcl):
        return [("\t  %20s        %20s        %20s" % ("LHS literal", "RHS literal", "Condition"), "exp", "%s", "%70s"),
                ("score", "score", "%+1.7f", "%8s"), ("Jacc", "acc", "%1.7f", "%8s")]

    def addPairCand(self, cand):
        # print("--- InitialPairs addPair", adv, litL, litR)
        cand.fillCurrentElement(self.getCurrentElement())
        return self.add(cand)
    
    def addPair(self, adv, clp, litL, litR, negL=False, negR=False, side=1):
        # print("--- InitialPairs addPair", adv, litL, litR)
        cand = PairCand(self.getSSetts(), adv, clp, litL, litR, negL, negR, swap=(side == 0))
        self.addPairCand(cand)
        
    def addExt(self, nego, adv, clp, lit_parts):
        # only one literal provided, and indication whether the other is negated (nego), will be filled in later, and flipped around if needed
        # print("--- InitialPairs addExt", nego, adv, clp, lit_parts)
        cand = PairCand(self.getSSetts(), adv, clp, nego, lit_parts[-1], nego, lit_parts[2], swap=(lit_parts[0] == 0))
        cand.fillCurrentElement(self.getCurrentElement())
        self.add(cand)

    def clear(self, current=None):
        ListStore.clear(self)
        InitsStore.setCurrent(self, current)

    # Post-processing
    ########################################################
    def removeNegDuplicates(self):
        dneg, oth, drop = ([], [], [])
        for ck in self.currentRoundKeys():
            if self[ck].doublyNegatedLit(0) or self[ck].doublyNegatedLit(1):
                dneg.append(ck)
            else:
                oth.append(ck)
        for ck in dneg:
            i = 0
            while i < len(oth):
                if self[oth[i]].getLit(0) == self[ck].getLit(0) and \
                   self[oth[i]].getLit(1) == self[ck].getLit(1):
                    drop.append(self[ck])
                    del self[ck]
                    i = len(oth)
                i += 1
        # if len(drop) > 0:
        #     print("DROP", drop)
        return drop


class InitialCands(object):
    # to use, must implement getStoredClass and getCmpCriterion class methods

    def __init__(self, sort_meth="overall", max_out=-1):
        self.max_out = max_out
        if sort_meth in SORT_METHODS:
            self.sort_meth_name = sort_meth
            self.sort_meth = SORT_METHODS[sort_meth]
        else:
            self.sort_meth_name = "default"
            self.sort_meth = DEFAULT_METHOD
        self.list_out = []
        self.drop_set = set()
        self.sorted_ids = None

    def msgLoaded(self):
        return "Loaded %d %s, will try at most %d (ordered by %s, %s)" % \
            (len(self), self.getStoredClassName(), self.getMaxOut(), self.getCmpCriterion(), self.sort_meth_name)

    def msgFound(self):
        already = ""
        if self.getNbOut() > 0:
            already = " tried %d before testing all," % self.getNbOut()
        return "Found %d %s,%s will try at most %d (ordered by %s, %s)" % \
            (len(self), self.getStoredClassName(), already, self.getMaxOut(), self.getCmpCriterion(), self.sort_meth_name)

    @classmethod
    def pattFound(tcl):
        return "Found \d+ %s," % tcl.getStoredClassName()

    def resetSort(self):
        self.sorted_ids = None

    def clear(self):
        self.list_out = []
        self.drop_set = set()
        self.resetSort()

    def setMaxOut(self, n):
        self.max_out = n

    def getNbOut(self):
        return len(self.list_out)

    def getMaxOut(self):
        return self.max_out

    def getRemainOut(self):
        if self.max_out == -1:
            return self.max_out
        else:
            return self.max_out - self.getNbOut()

    def exhausted(self):
        return (self.max_out > -1) and (self.getNbOut() >= self.max_out)

    # Adding and retrieving candidates
    def _sort(self):
        if self.sorted_ids is None:
            self.sorted_ids = self.sort_meth(self.ready_keys(), self, self.drop_set)

    # cond allows to provide a method that must evaluate to True on the candidate, else it is skipped
    def getNext(self, cond=None):
        if len(self) > 0 and not self.exhausted():
            self._sort()
            if len(self.sorted_ids) > 0:
                nid = self.sorted_ids.pop()
                self.drop_set.add(nid)
                if cond is not None:
                    while not cond(self[nid]) and len(self.sorted_ids) > 0:
                        nid = self.sorted_ids.pop()
                        self.drop_set.add(nid)
                    # exhausted candidates, no good one
                    if not cond(self[nid]) and len(self.sorted_ids) == 0:
                        return None
                self.list_out.append(nid)
                return self[nid]
        return None

    # cond allows to provide a method that must evaluate to True on the candidate, else it is skipped
    def getNextRed(self, data, cond=None):
        cand = self.getNext(cond)
        if cand is not None:
            return self.mkRedFromCand(cand, data)


class InitialTerms(ListStore, InitsStore, InitialCands):
    @classmethod
    def getStoredClass(tcl):
        return TermCand

    def __init__(self, N=0, ssetts=None, constraints=None, current=None, sort_meth="overall", max_out=-1):
        InitsStore.__init__(self, N, ssetts, constraints, current)
        ListStore.__init__(self)
        InitialCands.__init__(self, sort_meth, max_out)

    def getCmpCriterion(self):
        return "supp_ratio"

    def getCmpDefault(self):
        return -1

    def getKeepCriterion(self):
        return "supp_ratio"

    def getKeepMinThres(self):
        return self.constraints["init_minscore"]

    @classmethod
    def getDispProps(tcl):
        return [("\t  %20s        %20s        %20s" % ("LHS literal", "RHS literal", ""), "exp", "%s", "%70s"),
                ("Support ratio", "supp_ratio", "%+1.7f", "%8s")]

    def clear(self, current=None):
        ListStore.clear(self)
        InitialCands.clear(self)
        InitsStore.setCurrent(self, current)

    def add(self, cand):
        nid = ListStore.add(self, cand)
        if nid is not None:
            self.resetSort()
        return nid

    def setExploredDone(self):
        self.setTracking(False)

#######################################
# List of initial candidates storing pairs with tracking of what goes in and out, saving and explore list for parallel computation


class InitialPairs(PairsBatch, InitialCands):

    def __init__(self, N=0, ssetts=None, constraints=None, current=None, sort_meth="overall", max_out=-1, save_filename=None):
        PairsBatch.__init__(self, N, ssetts, constraints, current)
        InitialCands.__init__(self, sort_meth, max_out)
        self.setExploreList()
        self.setSaveFilename(save_filename)

    # Adding and retrieving candidates
    ########################################################
    def add(self, cand):
        nid = PairsBatch.add(self, cand)
        if nid is not None:
            self.resetSort()
            self.resetSaved()
        return nid

    def clear(self, current=None):
        PairsBatch.clear(self)
        InitialCands.clear(self)
        self.setExploreList()
        self.resetSaved()

    # For saving to/loading from file
    ########################################################

    def setSaveFilename(self, save_filename):
        self.save_filename = None
        self.saved = True
        if type(save_filename) is str and len(save_filename) > 0:
            self.save_filename = save_filename
            self.saved = False

    def resetSaved(self):
        self.save = False

    def getSaveFilename(self):
        return self.save_filename

    def saveToFile(self):  # TODO: also save condition
        if self.save_filename is not None and not self.saved:
            try:
                with open(self.save_filename, "w") as f:
                    f.write("\n".join([repr(k) + "\t" + repr(v) for k, v in self.items()])+"\n")
                    done = self.getExploredDone()
                    if done is None:
                        self.saved = True
                    else:
                        f.write("DONE: "+" ".join(["%d-%d" % d for d in done])+"\n")
                return True
            except IOError:
                pass
        return False

    def loadFromFile(self, data=None):
        loaded = False
        done = set()
        if self.save_filename is not None and os.path.isfile(self.save_filename):
            with open(self.save_filename) as f:
                loaded, done = self.loadFromLogFile(f, data)
            if not loaded:
                with open(self.save_filename) as f:
                    loaded, done = self.loadFromStoreFile(f, data)
            if done is None:
                self.setTracking(False)
        return loaded, done

    def loadFromStoreFile(self, f, data=None):
        done = None
        for line in f:
            if re.match("DONE:", line):
                done = set([tuple(map(int, d.split("-"))) for d in line.strip().split(" ")[1:]])
            else:
                parts = line.strip().split("\t")
                if len(parts) == 2:
                    cand = eval(parts[1])
                    cand.setClpSSizes(ssetts=self.getSSetts())
                    if data is None or (data.col(0, cand.getCid(0)).isEnabled() and data.col(1, cand.getCid(1)).isEnabled()):
                        self.add(cand)
        return True, done

    def loadFromLogFile(self, f, data=None):
        log_patt = "\[\[log@[0-9]+\s+\]\]\s+"
        cand_patt = log_patt+self.getStoredClass().pattCand()
        end_patt = log_patt+self.pattFound()

        done = set()
        for line in f:
            if re.match(end_patt, line):
                return len(done) > 0, None
            tmp = re.match(cand_patt, line)
            if tmp is not None:
                condition = None
                if tmp.group("cond") is not None:
                    if tmp.group("clpC") is not None:
                        clpC = eval(tmp.group("clpC"))
                    else:
                        clpC = None
                    condition = CondCand(self.getSSetts(), adv=eval(tmp.group("accC")), clp=clpC, lit_parts=[Literal(l) for l in tmp.group("litC").split(LIT_SEP)])
                if tmp.group("clp") is not None:
                    clp = eval(tmp.group("clp"))
                else:
                    clp = None
                cand = PairCand(self.getSSetts(), adv=eval(tmp.group("acc")), clp=clp, litL=Literal(tmp.group("litL")), litR=Literal(tmp.group("litR")), condition=condition)
                done.add((cand.getCid(0), cand.getCid(1)))
                if data is None or (data.col(0, cand.getCid(0)).isEnabled() and data.col(1, cand.getCid(1)).isEnabled()):
                    self.add(cand)
        return len(done) > 0, done

    # For parallel computation of initial candidates
    ########################################################

    def setExploreList(self, explore_list=[], pointer=-1, batch_size=0, done=set()):
        if done is None:
            done = set()
        self.explore_initc = {"list": explore_list, "pointer": pointer, "batch_size": batch_size, "done": done}

    def addExploredPair(self, pair):
        if self.explore_initc["done"] is not None:
            self.explore_initc["done"].add(pair)

    def getExploreList(self):
        return self.explore_initc["list"]

    def getExploredDone(self):
        return self.explore_initc["done"]

    def setExploredDone(self):
        self.setTracking(False)
        self.explore_initc["done"] = None
        self.explore_initc["pointer"] = -1

    def getExplorePointer(self):
        return self.explore_initc["pointer"]

    def setExplorePointer(self, pointer=-1):
        self.explore_initc["pointer"] = pointer

    def incrementExplorePointer(self):
        self.explore_initc["pointer"] += 1

    def getExploreNextBatch(self, pointer=None, bsize=None):
        if pointer is None:
            pointer = self.explore_initc["pointer"]
        if bsize is None:
            bsize = self.explore_initc["batch_size"]
        return self.explore_initc["list"][pointer*bsize:(pointer+1)*bsize]


def initCands(pt, data, constraints, save_filename=None):
    if pt == "T":
        return InitialTerms(data.nbRows(), constraints.getSSetts(), constraints, None,
                            constraints.getCstr("pair_sel"), constraints.getCstr("max_inits"))
    else:
        return InitialPairs(data.nbRows(), constraints.getSSetts(), constraints, None,
                            constraints.getCstr("pair_sel"), constraints.getCstr("max_inits"), save_filename=save_filename)
