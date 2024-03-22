from functools import reduce
from scipy.special import gammaln
from scipy.stats import binom
import numpy
import random
import re
import pdb

# A and B can be either True, False, or None (missing)


def land(A, B):
    if A is None and B == False:
        return B
    return A and B


def lor(A, B):
    if A is None and B == False:
        return None
    return A or B

def cmp_eq(a, b):
    return a == b

def cmp_neq(a, b):
    return a != b

def cmp_lower(a, b):
    if b is not None:
        if a is None:
            return True
        return a < b
    return False


def cmp_greater(a, b):
    if a is not None:
        if b is None:
            return True
        return a > b
    return False


def cmp_leq(a, b):
    if a is not None:
        if b is None:
            return False
        return a <= b
    return True


def cmp_geq(a, b):
    if b is not None:
        if a is None:
            return False
        return a >= b
    return True

cmp_fun_map = {"<": cmp_lower, ">": cmp_greater,
               "=": cmp_eq, "<>": cmp_neq,
               "<=": cmp_leq, ">=": cmp_geq}


def cmp_vals(a, b):
    if a == b:
        return 0
    if a < b:
        return -1
    return 1


def cmp_lists(A, B):
    # comparing two sets
    # 1. compare sizes of sets
    if len(A) < len(B):
        return -2
    elif len(A) > len(B):
        return 2
    # 2. compare elements of sets
    for a, b in zip(A, B):
        if a < b:
            return -1
        elif a > b:
            return 1
    return 0


def cmp_sets(A, B):
    # comparing two sets
    # 1. compare sizes of sets
    if len(A) < len(B):
        return -2
    elif len(A) > len(B):
        return 2
    # 2. compare elements of sets
    for a, b in zip(sorted(A), sorted(B)):
        if a < b:
            return -1
        elif a > b:
            return 1
    return 0


def cmp_listsets(AA, BB):
    # comparing two lists of sets
    # 1. compare length of lists
    if len(AA) < len(BB):
        return -3
    elif len(AA) > len(BB):
        return 3
    # 2. compare sizes of sets
    for A, B in zip(AA, BB):
        if len(A) < len(B):
            return -2
        elif len(A) > len(B):
            return 2
    # 3. compare elements of sets
    for a, b in zip(AA, BB):
        for a, b in zip(sorted(A), sorted(B)):
            if a < b:
                return -1
            elif a > b:
                return 1
    return 0


def cmp_reclists(A, B):
    # comparing two sets
    # 1. compare sizes of sets
    if len(A) < len(B):
        return -3
    elif len(A) > len(B):
        return 3
    # 2. compare elements of sets
    for a, b in zip(A, B):
        if type(a) is list:
            if type(b) is list:
                tmp = cmp_reclists(a, b)
                if tmp != 0:
                    return tmp
            else:
                return -2
        else:
            if type(b) is list:
                return 2
            else:
                tmp = cmp_vals(a, b)
                if tmp != 0:
                    return tmp
    return 0


def tool_ratioXX(num, den):
    if num is None or den is None:
        return None
    if den == 0:
        if num > 0:
            return float("Inf")
        else:
            return 0.
    else:
        return num/den


def tool_ratio(num, den):
    if den == 0:
        if num > 0:
            return float("Inf")
        else:
            return 0.
    try:
        return num/den
    except TypeError:
        return None


def tool_hypergeomPMF(k, M, n, N):
    tot, good = M, n
    bad = tot - good
    return numpy.exp(gammaln(good+1) - gammaln(good-k+1) - gammaln(k+1) + gammaln(bad+1)
                     - gammaln(bad-N+k+1) - gammaln(N-k+1) - gammaln(tot+1)
                     + gammaln(tot-N+1) + gammaln(N+1))
# same as the following but numerically more precise
# return comb(good,k) * comb(bad,N-k) / comb(tot,N)


def tool_pValOver(kInter, nbRows, suppL, suppR, lU=None):
    # probability that two sets of these size have intersection equal or larger than kInter
    if lU is None:
        return sum([tool_hypergeomPMF(k, nbRows, suppL, suppR) for k in range(kInter, min(suppL, suppR)+1)])
    else:
        return sum([tool_hypergeomPMF(k, nbRows, suppL, suppR) for k in range(0, lU - kInter+1)])


def tool_pValSupp(nbRows, supp, pr, lU=None):
    # probability that an itemset with marginal probability pr has support equal or larger than supp
    if lU is None:
        return 1-binom.cdf(supp-1, nbRows, pr)
    else:
        return binom.cdf(lU-supp, nbRows, pr)


status = [(True, False), (False, True), (True, True), (False, False),
          (True, None), (None, True), (False, None), (None, False), (None, None)]
labels_status = {True: "1", False: "0", None: "?"}
labelsm_status = {True: "x", False: "o", None: "m"}
labelsu_status = {True: "x", False: "o", None: "?"}
# labelsu_status = {True:u"\u2081", False:u"\u2080", None:u"\u2098"}
# WITH UNICODE
sym_status = labelsu_status
# ## WITHOUT UNICODE
# sym_status = labels_status

labels_iom = {True: "into", False: "out", None: "imiss"}


class SSetts(object):

    labels = ["E%s%s" % (labelsm_status[slhs], labelsm_status[srhs]) for (slhs, srhs) in status]
    labels_sparts = ["E%s%s" % (labels_status[slhs], labels_status[srhs])
                     for (slhs, srhs) in status]
    sym_sparts = [sym_status[slhs]+sym_status[srhs] for (slhs, srhs) in status]

    # define the numerical values for the parts
    for i, l in enumerate(labels):
        exec("%s = %d" % (l, i))
        exec("global %s; %s = %d" % (l, l, i))

    io_labels = [(0, "tot"), (1, "imiss"), (-2, "out"), (-1, "into")]
    io_labels_dict = dict(io_labels)
    io_ord = dict([(v, k) for (k, v) in enumerate([0, -1, -2, 1])])
    for i, l in io_labels:
        exec("%s = %d" % (l, i))
        exec("global %s; %s = %d" % (l, l, i))

    # Setters to filter the parts on assignment

    def tSing(self, s):
        return self.bottom <= s <= self.top

    def tPair(self, p):
        return self.tSing(p[1]) and (p[0] != self.index_drop)

    # def tTrip(self, t):
    #     return self.tPair(p[0]) and self.tSing(t[1])

    def fSings(self, l):
        return list(filter(self.tSing, l))

    def fPairs(self, l):
        return list(filter(self.tPair, l))

    # def fTrip(self, l):
    #     return list(filter(self.tTrip, l))

    @property
    def IDS_inter(self):
        return self._IDS_inter

    @IDS_inter.setter
    def IDS_inter(self, x):
        self._IDS_inter = self.fSings(x)

    @property
    def IDS_uncovered(self):
        return self._IDS_uncovered

    @IDS_uncovered.setter
    def IDS_uncovered(self, x):
        self._IDS_uncovered = self.fSings(x)

    @property
    def IDS_diff(self):
        return self._IDS_diff

    @IDS_diff.setter
    def IDS_diff(self, x):
        self._IDS_diff = self.fSings(x)

    @property
    def IDS_suppL(self):
        return self._IDS_suppL

    @IDS_suppL.setter
    def IDS_suppL(self, x):
        self._IDS_suppL = self.fSings(x)

    @property
    def IDS_fixnum(self):
        return self._IDS_fixnum

    @IDS_fixnum.setter
    def IDS_fixnum(self, x):
        self._IDS_fixnum = [self.fPairs(x[0]), self.fPairs(x[1])]

    @property
    def IDS_varnum(self):
        return self._IDS_varnum

    @IDS_varnum.setter
    def IDS_varnum(self, x):
        self._IDS_varnum = [self.fPairs(x[0]), self.fPairs(x[1])]

    @property
    def IDS_fixden(self):
        return self._IDS_fixden

    @IDS_fixden.setter
    def IDS_fixden(self, x):
        self._IDS_fixden = [self.fPairs(x[0]), self.fPairs(x[1])]

    @property
    def IDS_varden(self):
        return self._IDS_varden

    @IDS_varden.setter
    def IDS_varden(self, x):
        self._IDS_varden = [self.fPairs(x[0]), self.fPairs(x[1])]

    @property
    def IDS_out(self):
        return self._IDS_out

    @IDS_out.setter
    def IDS_out(self, x):
        self._IDS_out = [self.fPairs(x[0]), self.fPairs(x[1])]

    @property
    def IDS_cont(self):
        return self._IDS_cont

    @IDS_cont.setter
    def IDS_cont(self, x):
        self._IDS_cont = [self.fPairs(x[0]), self.fPairs(x[1])]

    @property
    def IDS_nsupp(self):
        return self._IDS_nsupp

    @IDS_nsupp.setter
    def IDS_nsupp(self, x):
        self._IDS_nsupp = [self.fPairs(x[0]), self.fPairs(x[1])]

    @property
    def IDS_supp(self):
        return self._IDS_supp

    @IDS_supp.setter
    def IDS_supp(self, x):
        self._IDS_supp = self.fSings(x)

    @property
    def IDS_miss(self):
        return self._IDS_miss

    @IDS_miss.setter
    def IDS_miss(self, x):
        self._IDS_miss = self.fSings(x)

    @property
    def IDS_negated(self):
        return self._IDS_negated

    @IDS_negated.setter
    def IDS_negated(self, x):
        self._IDS_negated = [self.fSings(x[0]), self.fSings(x[1])]

    @property
    def IDS_init(self):
        return self._IDS_init

    @IDS_init.setter
    def IDS_init(self, x):
        self._IDS_init = self.fSings(x)

    def rSings(self, l):
        return ", ".join([self.labels[s] for s in sorted(l)])

    def rPairs(self, l):
        return ", ".join(["(%s, %s)" % (self.io_labels_dict[p[0]], self.labels[p[1]]) for p in sorted(l, key=lambda x: (self.io_ord[x[0]], x[1]))])

    def dSings(self, l):
        return ", ".join([self.labels[s] for s in l])

    def dPairs(self, l):
        return ", ".join(["%5s_%s" % (self.io_labels_dict[p[0]], self.labels[p[1]]) for p in l])

    def dispPartsDef(self):
        dsp = "\n".join(["%sself.IDS_%s = [%s]" % (8*" ", lbl, self.rSings(parts)) for (lbl, parts)
                         in [("inter", self.IDS_inter), ("uncovered", self.IDS_uncovered),
                             ("diff", self.IDS_diff), ("suppL", self.IDS_suppL)]])
        dsp += "\n\n"
        dsp += "\n".join(["%sself.IDS_%s = [[%s],\n%s[%s]]" % (8*" ", lbl, self.rPairs(parts[0]), 27*" ", self.rPairs(parts[1])) for lbl, parts
                          in [("fixnum", self.IDS_fixnum), ("varnum", self.IDS_varnum),
                              ("fixden", self.IDS_fixden), ("varden", self.IDS_varden),
                              ("out", self.IDS_out), ("cont", self.IDS_cont), ("nsupp", self.IDS_nsupp)]])
        return dsp

    def __str__(self):
        dsp = "SSetts: %s, %s" % (self.getTypeParts(), self.getMethodPVal())
        dsp += "\n".join(["%10s:\t%s" % (lbl, self.dSings(parts)) for (lbl, parts)
                          in [("inter", self.IDS_inter), ("uncovered", self.IDS_uncovered),
                              ("diff", self.IDS_diff), ("suppL", self.IDS_suppL)]])
        dsp += "\n--- AND\n"
        dsp += "\n".join(["%10s:\t%s" % (lbl, self.dPairs(parts[0])) for lbl, parts
                          in [("fixnum", self.IDS_fixnum), ("varnum", self.IDS_varnum),
                              ("fixden", self.IDS_fixden), ("varden", self.IDS_varden),
                              ("out", self.IDS_out), ("cont", self.IDS_cont), ("nsupp", self.IDS_nsupp)]])
        dsp += "\n--- OR\n"
        dsp += "\n".join(["%10s:\t%s" % (lbl, self.dPairs(parts[1])) for lbl, parts
                          in [("fixnum", self.IDS_fixnum), ("varnum", self.IDS_varnum),
                              ("fixden", self.IDS_fixden), ("varden", self.IDS_varden),
                              ("out", self.IDS_out), ("cont", self.IDS_cont), ("nsupp", self.IDS_nsupp)]])
        return dsp

    # (Exo, Eox, Exx, Eoo[, Exm, Emx, Eom, Emo, Emm])
    # (tot,[ imiss,] out, into)

    # Truth table
    # A  B   A and B   A or B
    # -----------------------
    # o  x      o        x
    # m  o      o        m
    # x  m      m        x
    # o  o      o        o
    # m  m      m        m
    # x  x      x        x
    # o  m      o        m
    # x  o      o        x
    # m  x      m        x

    # extending A op X
    #  op   :   AND (False)  :   OR (True)
    #  X    :   x    o    m  :   x    o    m
    # ------:----------------:-----------------
    # Exo   :  xo   oo   mo  :  ---- xo ----
    # Eox   :  ---- ox ----  :  xx   ox   mx
    # Exx   :  xx   ox   mx  :  ---- xx ----
    # Eoo   :  ---- oo ----  :  xo   oo   mo
    # Exm   :  xm   om   mm  :  ---- xm ----
    # Emx   :  mx   ox   mx  :  xx   mx   mx
    # Eom   :  ---- om ----  :  xm   om   mm
    # Emo   :  mo   oo   mo  :  xo   mo   mo
    # Emm   :  mm   om   mm  :  xm   mm   mm

    # methods for initializing parts without/with missing values
    meths_init_part_ids = {True: {}, False: {}}

    ##############################################################
    # No missing values        J = |Exx| / |Exx|+|Exo|+|Eox|
    ##############################################################
    def init_part_ids_none(self):
        # TO COMPUTE ACCURACY after building
        # Part of intersection, appear in the numerator
        self.IDS_inter = [Exx]
        # Complement of union, do not appear in the denominator
        self.IDS_uncovered = [Eoo]
        # Symmetric difference, in the denominator, not numerator
        self.IDS_diff = [Exo, Eox]
        # Side specific
        self.IDS_suppL = [Exo, Exx]

        # TO COMPUTE ADVANCE while building, INDEXED BY OPERATOR (0: AND, 1: OR)
        # Parts in numerator, independent of X
        self.IDS_fixnum = [[], [(tot, Exx)]]
        # Parts in numerator, dependent of X
        self.IDS_varnum = [[(into, Exx)], [(into, Eox)]]
        # Parts in denominator, independent of X
        self.IDS_fixden = [[(tot, Eox), (tot, Exx)], [(tot, Exo), (tot, Eox), (tot, Exx)]]
        # Parts in denominator, dependent of X
        self.IDS_varden = [[(into, Exo)], [(into, Eoo)]]
        # Parts left uncovered (OUT), (always dependent of X)
        self.IDS_out = [[(tot, Eoo), (out, Exo)], [(out, Eoo)]]
        # Parts in contribution (CONT), (always dependent of X)
        # Contribution: AND entities moved from diff to out, OR: entities moved from diff to inter
        self.IDS_cont = [[(out, Exo)], [(into, Eox)]]
        # Parts in the new support of the extended query
        self.IDS_nsupp = [[(into, Exo), (into, Exx)], [(tot, Exo), (tot, Exx), (into, Eox), (into, Eoo)]]
    meths_init_part_ids[False]["none"] = init_part_ids_none

    ##############################################################
    # REJECTIVE       J = |Exx| / |Exx|+|Exo|+|Eox|
    ##############################################################
    def init_part_ids_rejective(self):
        self.IDS_inter = [Exx]
        self.IDS_uncovered = [Eoo]
        self.IDS_diff = [Exo, Eox]
        self.IDS_suppL = [Exo, Exx]

        self.IDS_fixnum = [[],
                           [(tot, Exx)]]
        self.IDS_varnum = [[(into, Exx)],
                           [(into, Eox), (into, Emx)]]
        self.IDS_fixden = [[(tot, Eox)],
                           [(tot, Exo), (tot, Exx)]]
        self.IDS_varden = [[(into, Exo), (into, Exx), (out, Exx), (out, Emx)],
                           [(into, Eox), (into, Eoo), (into, Emx), (into, Emo), (out, Eox)]]
        self.IDS_out = [[(tot, Eoo), (out, Exo), (out, Emo)],
                        [(out, Eoo)]]
        self.IDS_cont = [[(out, Exo)],
                         [(into, Eox)]]
        self.IDS_nsupp = [[(into, Exo), (into, Exx)],
                          [(tot, Exo), (tot, Exx), (into, Eox), (into, Eoo), (into, Emx), (into, Emo)]]
    meths_init_part_ids[True]["rejective"] = init_part_ids_rejective
    meths_init_part_ids[True]["grounded"] = init_part_ids_rejective

    ##############################################################
    # OPTIMISTIC      J = |Exx|+|Exm|+|Emx|+|Emm| /
    #                     |Exo|+|Eox|+|Exx|+|Exm|+|Emx|+|Emm|
    ##############################################################
    def init_part_ids_optimistic(self):
        self.IDS_inter = [Exx, Exm, Emx, Emm]
        self.IDS_uncovered = [Eoo, Eom, Emo]
        self.IDS_diff = [Exo, Eox]
        self.IDS_suppL = [Exo, Exx, Exm, Emx, Emm]

        self.IDS_fixnum = [[],
                           [(tot, Exx), (tot, Exm), (tot, Emx), (tot, Emm)]]
        self.IDS_varnum = [[(into, Exx), (into, Exm), (into, Emx), (into, Emm), (imiss, Exx), (imiss, Exm), (imiss, Emx), (imiss, Emm)],
                           [(into, Eox), (into, Eom), (imiss, Eox), (imiss, Eom)]]
        self.IDS_fixden = [[(tot, Eox), (tot, Exx), (tot, Emx)],
                           [(tot, Exo), (tot, Eox), (tot, Exx), (tot, Exm), (tot, Emx), (tot, Emm)]]
        self.IDS_varden = [[(into, Exo), (into, Exm), (into, Emm), (imiss, Exm), (imiss, Emm)],
                           [(into, Eoo), (into, Eom), (into, Emo), (imiss, Eom)]]
        self.IDS_out = [[(tot, Eoo), (tot, Eom), (tot, Emo), (out, Exo), (out, Exm), (out, Emm), (imiss, Exo)],
                        [(out, Eoo), (out, Eom), (out, Emo), (imiss, Eoo), (imiss, Emo)]]
        self.IDS_cont = [[(out, Exo), (imiss, Exo)],
                         [(into, Eox), (imiss, Eox)]]
        self.IDS_nsupp = [[(into, Exo), (into, Exx), (into, Exm), (into, Emx), (into, Emm), (imiss, Exx), (imiss, Exm), (imiss, Emx), (imiss, Emm)],
                          [(tot, Exo), (tot, Exx), (tot, Exm), (tot, Emx), (tot, Emm), (into, Eox), (into, Eoo), (into, Eom), (into, Emo), (imiss, Eox), (imiss, Eom)]]
    meths_init_part_ids[True]["optimistic"] = init_part_ids_optimistic

    ##############################################################
    # PESSIMISTIC     J = |Exx| /
    #                     |Exo|+|Eox|+|Exx|+|Exm|+|Emx|+|Eom|+|Emo|+|Emm|
    ##############################################################
    # --- corrected (oct.17) from    J= |Exx|   /
    # ---            |Exo|+|Eox|+|Exx|+|Eom|+|Emo|+|Emm|
    def init_part_ids_pessimistic(self):
        self.IDS_inter = [Exx]
        self.IDS_uncovered = [Eoo]
        self.IDS_diff = [Exo, Eox, Exm, Emx, Eom, Emo, Emm]
        self.IDS_suppL = [Exo, Exx, Exm, Emo]

        self.IDS_fixnum = [[],
                           [(tot, Exx)]]
        self.IDS_varnum = [[(into, Exx)],
                           [(into, Eox), (into, Emx)]]
        self.IDS_fixden = [[(tot, Eox), (tot, Exx), (tot, Exm), (tot, Emx), (tot, Eom), (tot, Emm)],
                           [(tot, Exo), (tot, Eox), (tot, Exx), (tot, Exm), (tot, Emx), (tot, Eom), (tot, Emo), (tot, Emm)]]
        self.IDS_varden = [[(into, Exo), (into, Emo), (imiss, Exo), (imiss, Emo)],
                           [(into, Eoo), (imiss, Eoo)]]
        self.IDS_out = [[(tot, Eoo), (out, Exo), (out, Emo)],
                        [(out, Eoo)]]
        self.IDS_cont = [[(out, Exo), (out, Emo)],
                         [(into, Eox), (into, Emx)]]
        self.IDS_nsupp = [[(into, Exo), (into, Exx), (into, Exm), (into, Emo), (imiss, Exo), (imiss, Emo)],
                          [(tot, Exo), (tot, Exx), (tot, Exm), (tot, Emo), (into, Eox), (into, Eoo), (into, Emx), (into, Eom), (into, Emm), (imiss, Eoo)]]
    meths_init_part_ids[True]["pessimistic"] = init_part_ids_pessimistic

    ##############################################################
    # POSITIVE        J = |Exx|+|Exm|+|Emx|+|Emm| /
    #                     |Exo|+|Eox|+|Exx|+|Exm|+|Emx|+|Eom|+|Emo|+|Emm|
    ##############################################################
    def init_part_ids_positive(self):
        self.IDS_inter = [Exx, Exm, Emx, Emm]
        self.IDS_uncovered = [Eoo]
        self.IDS_diff = [Exo, Eox, Eom, Emo]
        self.IDS_suppL = [Exo, Exx, Exm, Emx, Emo, Emm]

        self.IDS_fixnum = [[],
                           [(tot, Exx), (tot, Exm), (tot, Emx), (tot, Emm)]]
        self.IDS_varnum = [[(into, Exx), (into, Exm), (into, Emx), (into, Emm), (imiss, Exx), (imiss, Exm), (imiss, Emx), (imiss, Emm)],
                           [(into, Eox), (into, Eom), (imiss, Eox), (imiss, Eom)]]
        self.IDS_fixden = [[(tot, Eox), (tot, Exx), (tot, Exm), (tot, Emx), (tot, Eom), (tot, Emm)],
                           [(tot, Exo), (tot, Eox), (tot, Exx), (tot, Exm), (tot, Emx), (tot, Eom), (tot, Emo), (tot, Emm)]]
        self.IDS_varden = [[(into, Exo), (into, Emo), (imiss, Exo), (imiss, Emo)],
                           [(into, Eoo), (imiss, Eoo)]]
        self.IDS_out = [[(tot, Eoo), (out, Exo), (out, Emo)],
                        [(out, Eoo)]]
        self.IDS_cont = [[(out, Exo), (out, Emo)],
                         [(into, Eox), (into, Eom), (imiss, Eox), (imiss, Eom)]]
        self.IDS_nsupp = [[(into, Exo), (into, Exx), (into, Exm), (into, Emx), (into, Emo), (into, Emm), (imiss, Exo), (imiss, Exx), (imiss, Exm), (imiss, Emx), (imiss, Emo), (imiss, Emm)],
                          [(tot, Exo), (tot, Exx), (tot, Exm), (tot, Emx), (tot, Emo), (tot, Emm), (into, Eox), (into, Eoo), (into, Eom), (imiss, Eox), (imiss, Eoo), (imiss, Eom)]]
    meths_init_part_ids[True]["positive"] = init_part_ids_positive

    ##############################################################
    # NEGATIVE           J = |Exx| / |Exo|+|Eox|+|Exx|+|Exm|+|Emx|
    ##############################################################
    def init_part_ids_negative(self):
        self.IDS_inter = [Exx]
        self.IDS_uncovered = [Eoo, Eom, Emo, Emm]
        self.IDS_diff = [Exo, Eox, Exm, Emx]
        self.IDS_suppL = [Exo, Exx, Exm]

        self.IDS_fixnum = [[],
                           [(tot, Exx)]]
        self.IDS_varnum = [[(into, Exx)],
                           [(into, Eox), (into, Emx)]]
        self.IDS_fixden = [[(tot, Eox), (tot, Exx), (tot, Emx)],
                           [(tot, Exo), (tot, Eox), (tot, Exx), (tot, Exm), (tot, Emx)]]
        self.IDS_varden = [[(into, Exo), (into, Exm)],
                           [(into, Eoo), (into, Eom), (into, Emo), (into, Emm)]]
        self.IDS_out = [[(tot, Eoo), (tot, Eom), (tot, Emo), (tot, Emm), (out, Exo), (out, Exm), (imiss, Exo), (imiss, Exm)],
                        [(out, Eoo), (out, Eom), (out, Emo), (out, Emm), (imiss, Eoo), (imiss, Eom), (imiss, Emo), (imiss, Emm)]]
        self.IDS_cont = [[(out, Exo), (out, Exm), (imiss, Exo), (imiss, Exm)],
                         [(into, Eox), (into, Emx)]]
        self.IDS_nsupp = [[(into, Exo), (into, Exx), (into, Exm)],
                          [(tot, Exo), (tot, Exx), (tot, Exm), (into, Eox), (into, Eoo), (into, Emx), (into, Eom), (into, Emo), (into, Emm)]]
    meths_init_part_ids[True]["negative"] = init_part_ids_negative
    meths_init_part_ids[True]["basic"] = init_part_ids_negative

    ##############################################################
    # EXCLU (special) J = |Exo|+|Eox| / |Exo|+|Eox|+|Exx|
    ##############################################################
    def init_part_ids_exclu(self):
        self.IDS_inter = [Exo, Eox]
        self.IDS_uncovered = [Eoo]
        self.IDS_diff = [Exx]
        self.IDS_suppL = [Exo]

        self.IDS_fixnum = [[(tot, Eox)],
                           [(tot, Exo)]]
        self.IDS_varnum = [[(into, Exo), (out, Exx), (out, Emx)],
                           [(into, Eoo), (into, Emo), (out, Eox)]]
        self.IDS_fixden = [[(tot, Eox)],
                           [(tot, Exo), (tot, Exx)]]
        self.IDS_varden = [[(into, Exo), (into, Exx), (out, Exx), (out, Emx)],
                           [(into, Eox), (into, Eoo), (into, Emx), (into, Emo), (out, Eox)]]
        self.IDS_out = [[(tot, Eoo), (out, Exo), (out, Emo)],
                        [(out, Eoo)]]
        self.IDS_cont = [[(out, Exx)],
                         [(into, Emo), (out, Eox)]]
        self.IDS_nsupp = [[(into, Exo), (into, Exx)],
                          [(tot, Exo), (tot, Exx), (tot, Exm), (into, Eox), (into, Eoo), (into, Emx), (into, Eom), (into, Emo), (into, Emm)]]
    meths_init_part_ids["exclu"] = init_part_ids_exclu

    # without missing values, none, negative, optimistic and pessimistic are equivalent
    defaults_init_part_ids = {False: "none", True: "rejective"}
    ##############################################################

    map_label_part = dict([(s, p) for (p, s) in enumerate(labels)])
    map_status_part = dict([(s, p) for (p, s) in enumerate(status)])
    @classmethod
    def mapStatusToSPart(tcl, status):
        return tcl.map_status_part.get(status, -1)

    @classmethod
    def mapStatusToSPart(tcl, status):
        return tcl.map_status_part.get(status, -1)

    def __init__(self, has_missing, type_parts=None, methodpVal="Marg", data_symmetric=False):
        self.type_parts = None
        self.resetPartsIds(has_missing, type_parts)
        self.setMethodPVal(methodpVal)
        self.data_symmetric = data_symmetric

    def isDataSymmetric(self):
        return self.data_symmetric

    def getTypeParts(self):
        return self.type_parts

    def getMethodPVal(self):
        return self.methodpVal

    def reset(self, has_missing, type_parts=None, methodpVal=None, data_symmetric=False):
        if type_parts is not None:
            self.resetPartsIds(has_missing, type_parts)
        if methodpVal is not None:
            self.setMethodPVal(methodpVal)
        self.data_symmetric = data_symmetric

    def getAssigns(self, op, side):
        return self._assigns[op, side]

    def setAssigns(self, has_missing):
        if has_missing:
            states = [True, False, None]
        else:
            states = [True, False]
        assigns = {(False, 0): [], (True, 0): [], (False, 1): [], (True, 1): []}
        for A in states:
            for B in states:
                for X in states:
                    lfrom = "(%s, E%s%s)" % (labels_iom[X], labelsm_status[A], labelsm_status[B])
                    assigns[(False, 0)].append("(%s, E%s%s)" % (lfrom, labelsm_status[land(A, X)], labelsm_status[B]))
                    assigns[(True, 0)].append("(%s, E%s%s)" % (lfrom, labelsm_status[lor(A, X)], labelsm_status[B]))
                    assigns[(False, 1)].append("(%s, E%s%s)" % (lfrom, labelsm_status[A], labelsm_status[land(B, X)]))
                    assigns[(True, 1)].append("(%s, E%s%s)" % (lfrom, labelsm_status[A], labelsm_status[lor(B, X)]))

        self._assigns = {}
        for k, vs in assigns.items():
            self._assigns[k] = eval("["+", ".join(vs)+"]")

    def resetPartsIds(self, has_missing, type_parts=None):
        if type_parts not in self.meths_init_part_ids[has_missing]:
            type_parts = self.defaults_init_part_ids[has_missing]

        if self.type_parts == type_parts:
            return

        # indexes from the parts when looking from the right (A=L, B=R) or the left (A=R,B=L)
        self.side_index = [eval(",".join([self.labels[i] for i in [0, 1, 2, 3, 4, 5, 6, 7, 8]])),
                           eval(",".join([self.labels[i] for i in [1, 0, 2, 3, 5, 4, 7, 6, 8]]))]

        # indexes for the intersections with parts
        # (into: part inter X_True, out: part inter X_False, miss: part inter X_Missing, tot: total part = into + out + miss)
        # indexed for the intersections with parts when considering positive or negative X
        self.neg_index = [[tot, imiss, out, into], [tot, imiss, into, out]]

        # (Exo, Eox, Exx, Eoo, Exm, Emx, Eom, Emo, Emm) = range(9)
        self.last_nonmiss = Eoo

        if has_missing:
            self.bottom = Exo
            self.top = Emm
            self.index_drop = None
        else:
            self.bottom = Exo
            self.top = Eoo
            self.index_drop = imiss

        self.setAssigns(has_missing)
        self.meths_init_part_ids[has_missing][type_parts](self)

        # TO COMPUTE SUPPORTS, no index, common to all part types with missing values
        self.IDS_supp = [Exx, Exo, Exm]
        self.IDS_miss = [Emx, Emo, Emm]
        # indexes swaping when negating one side (0: negating A, 1: negating B)
        self.IDS_negated = [[Eoo, Exx, Eox, Exo, Eom, Emx, Exm, Emo, Emm],
                            [Exx, Eoo, Exo, Eox, Exm, Emo, Eom, Emx, Emm]]

        # initializing part counts
        self.IDS_init = [Eox, Eoo, Eom]

        self.IDS_num = []
        self.IDS_den = []
        for op in [0, 1]:
            self.IDS_num.append([x[1] for x in self.IDS_varnum[op]+self.IDS_fixnum[op]])
            self.IDS_den.append([x[1] for x in self.IDS_varden[op]+self.IDS_fixden[op]])
        self.type_parts = type_parts
    # return part label

    def getLabels(self):
        return self.labels

    def getLabel(self, id):
        return self.labels[id]

    # return the index corresponding to part_id when looking from given side
    def partId(self, part_id, side=0):
        return self.side_index[side][part_id]

    # return the index corresponding to part_id when negating given side
    def negatedPartId(self, part_id, side=0):
        return self.IDS_negated[side][part_id]

    # return the index corresponding to inout and possible negation

    def inOutId(self, inout_id, neg=0):
        return self.neg_index[neg][inout_id]

    # sums the values in parts that correspond to part_id indexes given in parts_id
    # parts_id can be
    # * a list of pairs (inout, part_id), inout are then ignored
    # * a list of values part_id

    def sumPartsId(self, side, parts_id, parts):
        if type(parts) == int:
            return 1*(parts in [self.partId(part_id[1], side) for part_id in parts_id])
        elif len(parts_id) > 0:
            if type(parts_id[0]) == int:
                ids = parts_id
            elif len(parts_id[0]) == 2:
                (inout, ids) = zip(*parts_id)
            else:
                ids = []
            return sum([parts[self.partId(part_id, side)] for part_id in set(ids)])
        return 0

    def suppPartRange(self):
        return range(self.bottom, self.top+1)

    def suppPartRangeNoMiss(self):
        return range(self.bottom, self.last_nonmiss+1)

    # sums the values in parts that correspond to inout and part_id indexes given in parts_id
    # parts_id must be
    # * a list of pairs (inout, part_id)

    def sumPartsIdInOut(self, side, neg, parts_id, parts):
        return sum([parts[self.inOutId(part_id[0], neg)][self.partId(part_id[1], side)] for part_id in parts_id])

    # return parts reordered to match the new indexes of parts corresponding to negation of given side

    def negateParts(self, side, parts):
        return [parts[self.negatedPartId(p, side)] for p in range(len(parts))]

    def negatePartsV(self, side, parts):
        return parts[self.IDS_negated[side][:len(parts)]]

    # sets the method to compute p-values
    def setMethodPVal(self, methodpVal="Marg"):
        try:
            self.methodpVal = methodpVal.capitalize()
            eval("self.pVal%sQueryCand" % (self.methodpVal))
            eval("self.pVal%sRedCand" % (self.methodpVal))
            # self.pValQueryCand = eval("self.pVal%sQueryCand" % (self.methodpVal))
            # self.pValRedCand = eval("self.pVal%sRedCand" % (self.methodpVal))
        except AttributeError:
            raise Exception("Oups method to compute the p-value does not exist !")

    def pValRedCand(self, side, op, neg, lParts, N, prs=None, method=""):
        meth = eval("self.pVal%sRedCand" % (self.methodpVal))
        return meth(side, op, neg, lParts, N, prs)

    def pValQueryCand(self, side, op, neg, lParts, N, prs=None):
        meth = eval("self.pVal%sQueryCand" % (self.methodpVal))
        return meth(side, op, neg, lParts, N, prs)
        # return 0 # self.pValSuppQueryCand(side, op, neg, lParts, N, prs)

    # query p-value using support probabilities (binomial), for candidates
    def pValSuppQueryCand(self, side, op, neg, lParts, N, prs=None):
        if prs is None:
            return 0
        else:
            lInter = self.sumPartsId(side, self.IDS_supp, lParts[self.inOutId(self.into, neg)])
            lX = sum(lParts[self.inOutId(self.into, neg)])
            if op:
                return 1-tool_pValSupp(N, lInter, prs[side] + lX/N - prs[side]*lX/N)
            else:
                return tool_pValSupp(N, lInter, prs[side]*lX/N)

    # query p-value using marginals (binomial), for candidates
    def pValMargQueryCand(self, side, op, neg, lParts, N, prs=None):
        if prs is None:
            return 0
        else:
            lInter = self.sumPartsId(side, self.IDS_supp, lParts[self.inOutId(self.into, neg)])
            lsupp = self.sumPartsId(side, self.IDS_supp, lParts[self.inOutId(self.tot, neg)])
            lX = sum(lParts[self.inOutId(self.into, neg)])
            if op:
                # return 1-tool_pValSupp(N, lInter, lsupp*lX/(N*N))
                vv = 1-tool_pValSupp(N, lInter, lsupp*lX/(N*N))
            else:
                # return tool_pValSupp(N, lInter, lsupp*lX/(N*N))
                vv = tool_pValSupp(N, lInter, lsupp*lX/(N*N))
            # print("---- pVal Marg Query", lInter, lsupp, lX, N, lsupp*lX/(N*N), "-->", vv)
            return vv

    # query p-value using support sizes (hypergeom), for candidates
    def pValOverQueryCand(self, side, op, neg, lParts, N, prs=None):
        if prs is None:
            return 0
        else:
            lInter = self.sumPartsId(side, self.IDS_supp, lParts[self.inOutId(self.into, neg)])
            lsupp = self.sumPartsId(side, self.IDS_supp, lParts[self.inOutId(self.tot, neg)])
            lX = sum(lParts[self.inOutId(self.into, neg)])
            if op:
                return 1-tool_pValOver(lInter, N, lsupp, lX)
            else:
                return tool_pValOver(lInter, N, lsupp, lX)

    # redescription p-value using support probabilities (binomial), for candidates

    def pValSuppRedCand(self, side, op, neg, lParts, N, prs=None):
        lInter = self.sumPartsIdInOut(side, neg, self.IDS_fixnum[op] + self.IDS_varnum[op], lParts)
        lUnion = None
        if self.type_parts == "exclu":
            lUnion = self.sumPartsIdInOut(side, neg, self.IDS_fixden[op] + self.IDS_varden[op], lParts)

        lX = sum(lParts[self.inOutId(self.into, neg)])
        # if self.pValOut: pdb.set_trace()
        if prs is None:
            lO = self.sumPartsId(1-side, self.IDS_supp, lParts[self.inOutId(self.tot, neg)])
            return tool_pValSupp(N, lInter, (lO*lX)/(N*N), lU=lUnion)
        elif op:
            return tool_pValSupp(N, lInter, prs[1-side]*(prs[side] + lX/N - prs[side]*lX/N), lU=lUnion)
        else:
            return tool_pValSupp(N, lInter, prs[1-side]*(prs[side] * lX/N), lU=lUnion)

    # redescription p-value using marginals (binomial), for candidates

    def pValMargRedCand(self, side, op, neg, lParts, N, prs=None):
        lInter = self.sumPartsIdInOut(side, neg, self.IDS_fixnum[op] + self.IDS_varnum[op], lParts)
        lUnion = None
        if self.type_parts == "exclu":
            lUnion = self.sumPartsIdInOut(
                side, neg, self.IDS_fixden[op] + self.IDS_varden[op], lParts)

        lO = self.sumPartsId(1-side, self.IDS_supp, lParts[self.inOutId(self.tot, neg)])
        lS = self.sumPartsIdInOut(side, neg, self.IDS_nsupp[op], lParts)
        return tool_pValSupp(N, lInter, (lO*lS)/(N*N), lU=lUnion)

    # redescription p-value using support sizes (hypergeom), for candidates
    def pValOverRedCand(self, side, op, neg, lParts, N, prs=None):
        lInter = self.sumPartsIdInOut(side, neg, self.IDS_fixnum[op] + self.IDS_varnum[op], lParts)
        lUnion = None
        if self.type_parts == "exclu":
            lUnion = self.sumPartsIdInOut(
                side, neg, self.IDS_fixden[op] + self.IDS_varden[op], lParts)

        lO = self.sumPartsId(1-side, self.IDS_supp, lParts[self.inOutId(self.tot, neg)])
        lS = self.sumPartsIdInOut(side, neg, self.IDS_nsupp[op], lParts)
        return tool_pValOver(lInter, N, lO, lS, lU=lUnion)

    # initialize parts counts
    # default count for every part is zero
    # pairs contains a list of (part_id, value)
    # if value is non negative, the count of part_id is set to that value
    # if value is negative, the count of part_id is set to - value - sum of the other parts set so far
    def makeLPartsXX(self, pairs=[], side=0):
        lp = [0 for i in range(self.top+1)]
        for (part_id, val) in pairs:
            if self.partId(part_id, side) < len(lp):
                if val < 0:
                    tmp = sum(lp)
                    lp[self.partId(part_id, side)] = -val - tmp
                else:
                    lp[self.partId(part_id, side)] = val
            else:
                if val > 0:
                    raise Exception("Some missing data where there should not be any!")
        return lp

    def makeLParts(self, counts=None, pids=None, side=0):
        lp = numpy.zeros(self.top+1, dtype=int)
        if counts is not None:
            if pids is None:
                pids = self.getInitPartIds(side)
                lp[pids] = counts[:len(pids)]  # cut off the counts if no missing parts
            else:
                lp[pids] = counts
        return lp

    def clpToLSupports(self, clp, side=0, op=True, pair=False, swap=False):
        # retrieves the sizes of parts for SSizes from list of sizes clp (counts of parts from candidate, built by greedy mining algos)
        lsupports = numpy.zeros(self.top+1, dtype=int)
        if pair:
            for cfrom, lto in self.getAssigns(1, True):
                lsupports[self.partId(lto, side=1*swap)] += clp[cfrom[0]][self.partId(cfrom[1], side=1*swap)]
        else:
            for cfrom, lto in self.getAssigns(side, op):
                lsupports[lto] += clp[cfrom[0]][cfrom[1]]
        if self.top == self.last_nonmiss:
            return lsupports.sum(), list(lsupports[:-1])  # no missing values, drop Eoo, implicit
        else:
            return lsupports.sum(), list(lsupports)

    def getInitPartIds(self, side=0):
        # side == 0 -> Exo, Eoo, Emo
        # side == 1 -> Exo, Eoo, Emo
        return [self.partId(pid, side) for pid in self.IDS_init]


class SSizes(object):

    class_letter = "z"
    # PROPS WHAT
    info_what = {"acc": "self.acc()", "pval": "self.pVal()"}
    props_what = ["len", "card", "perc", "ratio", "area"]
    Pwhat_match = "(" + "|".join(list(info_what.keys()) + props_what) + ")"
    @classmethod
    def hasPropWhat(tcl, what):
        return re.match(tcl.Pwhat_match, what) is not None

    # PROPS WHICH
    sets_letters = "PIULROABN"
    Pwhich_match = "(" + "|".join(["["+sets_letters+"]"] + list(SSetts.map_label_part.keys())) + ")"
    @classmethod
    def hasPropWhich(tcl, which):
        return re.match(tcl.Pwhich_match, which) is not None

    props_stats = [("acc", None), ("len", "I"), ("pval", None)]

    @classmethod
    def prepare_lsupports_sizes(tcl, lsupports, N):  # retrieves the sizes of parts from list of sizes, which must be of length 3 or 9
        if type(lsupports) == list and len(lsupports) == 9 and lsupports[8] + lsupports[7] + lsupports[6] + lsupports[5] + lsupports[4] == 0:
            lsupports = lsupports[0:3]  # all missing empty -> cut off and proceed

        # two lsupports: interpreted as (suppL, suppR) -> cannot deduce size of Exx
        # three lsupports: interpreted as (Exo, Eox, Exx)
        if type(lsupports) == list and len(lsupports) == 3:
            missing = False
            sSizes = [lsupports[0], lsupports[1], lsupports[2]]
        # four lsupports: interpreted as (suppL, suppR, missL, missR) -> cannot deduce size of parts
        # nine lsupports: interpreted as (Exo, Eox, Exx, Eoo, Exm, Emx, Eom, Emo, Emm)
        elif type(lsupports) == list and len(lsupports) == 9:
            missing = True
            sSizes = [support for support in lsupports]
        # else:  not valid
        else:
            missing = False
            sSizes = None
        return missing, sSizes

    @classmethod
    def prepare_supports_sizes(tcl, supports, N):  # retrieves the sizes of parts from list of sets, which must be of length 2, 3, 4 or 9
        if type(supports) == list and len(supports) == 4 and len(supports[2]) + len(supports[3]) == 0:
            supports = supports[0:2]
        elif type(supports) == list and len(supports) == 9 and len(supports[8]) + len(supports[7]) + len(supports[6]) + len(supports[5]) + len(supports[4]) == 0:
            supports = supports[0:3]

        # two supports: interpreted as (suppL, suppR)
        if type(supports) == list and len(supports) == 2:
            (suppL, suppR) = supports
            missing = False
            sSizes = [len(suppL - suppR),
                      len(suppR - suppL),
                      len(suppL & suppR)]
        # three supports: interpreted as (Exo, Eox, Exx)
        elif type(supports) == list and len(supports) == 3:
            missing = False
            sSizes = [len(supports[0]), len(supports[1]), len(supports[2])]
        # four supports: interpreted as (suppL, suppR, missL, missR)
        elif type(supports) == list and len(supports) == 4:
            missing = True
            (suppL, suppR, missL, missR) = supports
            sSizes = [len(suppL - suppR - missR),
                      len(suppR - suppL - missL),
                      len(suppL & suppR),
                      len(set(range(N)) - suppL - suppR - missL - missR),
                      len(suppL & missR),
                      len(suppR & missL),
                      len(missR - suppL - missL),
                      len(missL - suppR - missR),
                      len(missL & missR)]
        # nine supports: interpreted as (Exo, Eox, Exx, Eoo, Exm, Emx, Eom, Emo, Emm)
        elif type(supports) == list and len(supports) == 9:
            missing = True
            sSizes = [len(support) for support in supports]
        # else: not valid
        else:
            missing = False
            sSizes = None
        return missing, sSizes

    def __init__(self, ssetts, N, lsupports, prs=[1, 1]):
        # init from dict_info
        self.ssetts = ssetts
        if type(N) == dict:
            sdict = N
            self.missing = False
            self.sSizes = [set() for i in range(len(self.ssetts.getLabels()))]
            self.prs = [-1, -1]
            self.N = 0
            for i, supp_key in enumerate(self.ssetts.getLabels()):
                if supp_key in sdict and type(sdict[supp_key]) is int:
                    if i > 3 and sdict[supp_key] > 0:
                        self.missing = True
                    self.sSizes[i] = sdict.pop(supp_key)

            if 'pr_0' in sdict:
                self.prs[0] = sdict.pop('pr_0')
            if 'pr_1' in sdict:
                self.prs[1] = sdict.pop('pr_1')
            if 'N' in sdict:
                self.N = sdict.pop('N')
            if not self.missing:
                del self.sSizes[4:]

        else:
            if type(N) is set:
                self.N = len(N)
                bk = N
            else:
                self.N = N
                bk = None

            self.vect = None
            self.missing, self.sSizes = self.prepare_lsupports_sizes(lsupports, self.N)
            if self.sSizes is None:
                self.sSizes = [0 for i in range(self.ssetts.top+1)]
                bk = None

            if bk is not None:
                if len(self.sSizes) == 3:
                    self.sSizes.append(bk)
                else:
                    self.sSizes[self.ssetts.Eoo] = bk
                for si, sp in enumerate(self.sSizes):
                    if si != self.ssetts.Eoo:
                        self.sSizes[self.ssetts.Eoo] -= sp
            if prs is None:
                self.prs = [self.lenSupp(0)/self.N, self.lenSupp(1)/self.N]
            else:
                self.prs = prs

    def swapSides(self):
        _swapped_ids = [self.ssetts.partId(i, side=1) for i in range(len(self.sSizes))]
        self.sSizes = [self.sSizes[x] for x in _swapped_ids]
        self.prs[0], self.prs[1] = self.prs[1], self.prs[0]

    @classmethod
    def privatize_support(tcl, ssetts, N, supports, prs=[1, 1], budget=0):
        has_missing, lsupports = tcl.prepare_supports_sizes(supports, N)
        lsupports = [numpy.clip(round(x+numpy.random.laplace(scale=1/budget)), 0, N) for x in lsupports]
        return SSizes(ssetts, N, lsupports, prs)

    def copy(self):
        return SSizes(self.ssetts, self.N, self.sSizes, prs=self.prs if self.prs is None else list(self.prs))

    def nbStored(self):
        return len(self.sSizes)

    def __eq__(self, other):
        return isinstance(other, SSizes) and self.N == other.N and cmp_lists(self.sSizes, other.sSizes) == 0

    def __ne__(self, other):
        return not isinstance(other, SSizes) or self.N != other.N or cmp_lists(self.sSizes, other.sSizes) != 0
    # !! if not the same length or not the same total, set lists are not comparable

    def __lt__(self, other):
        if isinstance(other, SSizes) and self.N == other.N:
            c = cmp_lists(self.sSizes, other.sSizes)
            return c > -3 and c < 0
        return False

    def __le__(self, other):
        if isinstance(other, SSizes) and self.N == other.N:
            c = cmp_lists(self.sSizes, other.sSizes)
            return c > -3 and c <= 0
        return False

    def __gt__(self, other):
        if isinstance(other, SSizes) and self.N == other.N:
            c = cmp_lists(self.sSizes, other.sSizes)
            return c < 3 and c > 0
        return False

    def __ge__(self, other):
        if isinstance(other, SSizes) and self.N == other.N:
            c = cmp_lists(self.sSizes, other.sSizes)
            return c < 3 and c >= 0
        return False

    def getTypeParts(self):
        return self.ssetts.getTypeParts()

    def getMethodPVal(self):
        return self.ssetts.getMethodPVal()

    def proba(self, side):
        if self.prs is None:
            return -1
        return self.prs[side]

    def pVal(self):
        try:
            return eval("self.pVal%s()" % self.getMethodPVal())
        except AttributeError:
            raise Exception("Oups method to compute the p-value does not exist !")

    def getSSetts(self):
        return self.ssetts

    def nbRows(self):
        return self.N

    def toDict(self, with_Eoo=False):
        sdict = {}
        for i in range(self.nbStored()):
            sdict["card_" + self.ssetts.getLabel(i)] = self.lpart(i)
            sdict["perc_" + self.ssetts.getLabel(i)] = self.lpart(i) * 100. / self.N
        if with_Eoo:
            sdict["card_" + self.ssetts.getLabel(SSetts.Eoo)] = self.lpart(SSetts.Eoo)
            sdict["perc_" + self.ssetts.getLabel(SSetts.Eoo)] = self.lpart(SSetts.Eoo) * 100. / self.N
        for side in [0, 1]:
            if self.prs is not None and self.prs[side] != -1:
                sdict["pr_" + str(side)] = self.prs[side]
        sdict["N"] = self.N
        for info_key, info_meth in self.info_what.items():
            sdict[info_key] = eval(info_meth)
        return sdict

    # contains missing values
    def hasMissing(self):
        return self.missing

    # return copy of the probas
    def probas(self):
        if self.prs is not None:
            return list(self.prs)

    def lpart(self, part_id, side=0):
        pid = self.ssetts.partId(part_id, side)
        if pid < self.nbStored():
            return self.sSizes[pid]
        elif part_id == self.ssetts.Eoo:
            return self.N - self.sSizes[0] - self.sSizes[1] - self.sSizes[2]
        else:
            return 0

    def lparts(self, side=0):
        return [self.lpart(i, side) for i in range(self.ssetts.top+1)]

    def nbParts(self):
        return self.ssetts.top+1

    def lparts_union(self, ids, side=0):
        return sum([self.lpart(i, side) for i in ids])

    def lenSupp(self, side=0):
        return self.lparts_union(self.ssetts.IDS_supp, side)

    def lenNonSupp(self, side=0):
        return self.N - self.lenSupp(side) - self.lenMiss(side)

    def lenMiss(self, side=0):
        if not self.missing:
            return 0
        else:
            return self.lparts_union(self.ssetts.miss_ids, side)

    # LENGTHS
    # corresponding lengths
    def lenSide(self, side):
        return self.lparts_union(self.ssetts.IDS_suppL, side)

    def lenP(self, i, side=0):
        return self.lpart(i, side)

    def lenD(self, side=0):
        return self.lparts_union(self.ssetts.IDS_diff, side)

    def lenI(self, side=0):
        return self.lparts_union(self.ssetts.IDS_inter, side)

    def lenU(self, side=0):
        return self.lenD(side)+self.lenI(side)

    def lenL(self, side=0):
        return self.lenSide(0)

    def lenR(self, side=0):
        return self.lenSide(1)

    def lenO(self, side=0):
        return self.lparts_union(self.ssetts.IDS_uncovered, side)

    def lenA(self, side=0):
        return self.lparts_union(self.ssetts.IDS_suppL, side)

    def lenB(self, side=0):
        return self.lparts_union(self.ssetts.IDS_suppL, 1-side)

    def lenN(self, side=0):
        if self.nbStored() == 4:
            return self.lparts_union(range(4), side)
        else:
            return self.N

    # def lenD(self, side=0):
    #     return self.lparts_union(self.ssetts.IDS_diff, side)
    # def lenI(self, side=0):
    #     return self.lparts_union(self.ssetts.IDS_inter, side)
    # def lenU(self, side=0):
    #     return self.lparts_union(self.ssetts.IDS_inter+self.ssetts.IDS_diff, side)
    #     return self.suppI(side) | self.suppD(side)
    # def lenL(self, side=0):
    #     return self.lenSide(0)
    # def lenR(self, side=0):
    #     return self.lenSide(1)
    # def lenO(self, side=0):
    #     return self.lparts_union(self.ssetts.IDS_uncovered, side)

    def getProp(self, what, which=None):
        if what in self.info_what:
            return eval(self.info_what[what])
        wt = what
        if what == "card" or what == "area":
            wt = "len"
        methode = eval("self.%s" % wt)
        if callable(methode):
            return methode(which)

    def len(self, which="I"):
        if which in SSetts.map_label_part:
            return self.lenP(SSetts.map_label_part[which])
        elif which in self.sets_letters:
            return eval("self.len%s()" % which)

    def ratio(self, which="I"):
        return tool_ratio(self.len(which), self.nbRows())

    def perc(self, which="I"):
        return tool_ratio(self.len(which), self.nbRows()/100.)

    # accuracy
    def acc(self, side=0):
        lenI = self.lenI(side)
        return tool_ratio(lenI, lenI+self.lenD(side))

    # redescription p-value using support probabilities (binomial), for redescriptions
    def pValSupp(self):
        if self.prs == [-1, -1] or self.N == -1:
            return -1
        elif self.lenSupp(0)*self.lenSupp(1) == 0:
            return 1.
        else:
            lUnion = self.lenU() if self.getTypeParts() == "exclu" else None
            return tool_pValSupp(self.N, self.lenI(), self.prs[0]*self.prs[1], lU=lUnion)

    # redescription p-value using marginals (binomial), for redescriptions
    def pValMarg(self):
        if self.N == -1:
            return -1
        elif self.lenSupp(0)*self.lenSupp(1) == 0:
            return 1.
        else:
            lUnion = self.lenU() if self.getTypeParts() == "exclu" else None
            return tool_pValSupp(self.N, self.lenI(), (self.lenSupp(0)*self.lenSupp(1))/(numpy.int64(self.N)*numpy.int64(self.N)), lU=lUnion)

    # redescription p-value using support sizes (hypergeom), for redescriptions
    def pValOver(self):
        if self.N == -1:
            return -1
        elif self.lenSupp(0)*self.lenSupp(1) == 0:
            return 1.
        else:
            lUnion = self.lenI() if self.getTypeParts() == "exclu" else None
            return tool_pValOver(self.lenI(), self.N, self.lenSupp(0), self.lenSupp(1), lU=lUnion)

    # update support probabilities
    @classmethod
    def updateProba(tcl, prA, prB, OR):
        if type(prA) == int and prA == -1:
            return prB
        elif OR:
            return prA + prB - prA*prB
        else:
            return prA*prB

    # update support probabilities
    @classmethod
    def updateProbaMass(tcl, prs, OR):
        if len(prs) == 1:
            return prs[0]
        elif OR:
            return reduce(lambda x, y: x+y-x*y, prs)
        else:
            return numpy.prod(prs)

    # PRINTING
    ##############
    # def __str__(self):
    #         s = "|"
    #         r = "||\n|"
    #         if self.missing: up_to = self.ssetts.Emm
    #         else: up_to = self.ssetts.Eoo
    #         for i in range(up_to+1):
    #             s += "|%s" % (3*self.ssetts.getLabel(i))
    #             r += "| % 4i " % self.lpart(i,0)
    #         return s+r+"||"

    def __str__(self):
        return "SUPPORT:" + self.dispSuppL(sep=" ")

    def dispSuppL(self, sep="\t"):
        return sep.join(["card_" + self.ssetts.getLabel(i)+":" + str(self.lpart(i)) for i in range(self.ssetts.top+1)])

    def dispStats(self, sep="\t"):
        return sep.join(["%s%s:%s" % (what, which or "", self.getProp(what, which)) for (what, which) in self.props_stats])


class SParts(SSizes):

    class_letter = "s"
    # PROPS WHAT
    props_what = ["set", "supp"]+SSizes.props_what
    Pwhat_match = "(" + "|".join(list(SSizes.info_what.keys()) + props_what) + ")"
    @classmethod
    def prepare_supports_parts(tcl, supports, N):  # retrieves the support parts from list of sets, which must be of length 2, 3, 4 or 9
        if type(supports) == list and len(supports) == 4 and len(supports[2]) + len(supports[3]) == 0:
            supports = supports[0:2]
        elif type(supports) == list and len(supports) == 9 and len(supports[8]) + len(supports[7]) + len(supports[6]) + len(supports[5]) + len(supports[4]) == 0:
            supports = supports[0:3]

        # two supports: interpreted as (suppL, suppR)
        if type(supports) == list and len(supports) == 2:
            (suppL, suppR) = supports
            missing = False
            sParts = [set(suppL - suppR),
                      set(suppR - suppL),
                      set(suppL & suppR)]
        # three supports: interpreted as (Exo, Eox, Exx)
        elif type(supports) == list and len(supports) == 3:
            missing = False
            sParts = [set(supports[0]), set(supports[1]), set(supports[2])]
        # four supports: interpreted as (suppL, suppR, missL, missR)
        elif type(supports) == list and len(supports) == 4:
            missing = True
            (suppL, suppR, missL, missR) = supports
            sParts = [set(suppL - suppR - missR),
                      set(suppR - suppL - missL),
                      set(suppL & suppR),
                      set(range(N)) - suppL - suppR - missL - missR,
                      set(suppL & missR),
                      set(suppR & missL),
                      set(missR - suppL - missL),
                      set(missL - suppR - missR),
                      set(missL & missR)]
        # nine supports: interpreted as (Exo, Eox, Exx, Eoo, Exm, Emx, Eom, Emo, Emm)
        elif type(supports) == list and len(supports) == 9:
            missing = True
            sParts = [set(support) for support in supports]
        # else: not valid
        else:
            missing = False
            sParts = None
        return missing, sParts

    def __init__(self, ssetts, N, supports, prs=[1, 1]):
        # sParts is a partition of the rows (Eoo is not explicitely stored when there are no missing values)
        # init from dict_info
        self.ssetts = ssetts
        if type(N) == dict:
            sdict = N
            self.missing = False
            self.sParts = [set() for i in range(len(self.ssetts.getLabels()))]
            self.prs = [-1, -1]
            self.N = 0
            for i, supp_key in enumerate(self.ssetts.getLabels()):
                if supp_key in sdict:
                    if i > 3 and len(sdict[supp_key]) > 0:
                        self.missing = True
                    self.sParts[i] = set(sdict.pop(supp_key))

            if "pr_0" in sdict:
                self.prs[0] = sdict.pop("pr_0")
            if "pr_1" in sdict:
                self.prs[1] = sdict.pop("pr_1")
            if "N" in sdict:
                self.N = sdict.pop("N")
            if not self.missing:
                del self.sParts[4:]
        else:
            if type(N) is set:
                self.N = len(N)
                bk = N
            else:
                self.N = N
                bk = None

            self.vect = None
            self.missing, self.sParts = self.prepare_supports_parts(supports, self.N)
            if self.sParts is None:
                self.sParts = [set() for i in range(self.ssetts.top+1)]
                bk = None

            if bk is not None:
                if len(self.sParts) == 3:
                    self.sParts.append(set(bk))
                else:
                    self.sParts[self.ssetts.Eoo] = set(bk)
                for si, sp in enumerate(self.sParts):
                    if si != self.ssetts.Eoo:
                        self.sParts[self.ssetts.Eoo] -= sp
            if prs is None:
                self.prs = [self.lenSupp(0)/self.N, self.lenSupp(1)/self.N]
            else:
                self.prs = prs

    def swapSides(self):
        # _swapped_ids = [self.ssetts.partId(i, side=1) for i in range(len(self.sParts))]
        # pdb.set_trace()
        self.sParts = [self.sParts[self.ssetts.partId(i, side=1)] for i in range(len(self.sParts))]
        self.prs[0], self.prs[1] = self.prs[1], self.prs[0]

    def privatize(self, budget=0):
        return self.privatize_support(self.ssetts, self.N, self.sParts, prs=list(self.prs), budget=budget)

    def copy(self):
        return SParts(self.ssetts, self.N, self.sParts, prs=list(self.prs))

    def nbStored(self):
        return len(self.sParts)

    def __eq__(self, other):
        return isinstance(other, SParts) and self.N == other.N and cmp_listsets(self.sParts, other.sParts) == 0

    def __ne__(self, other):
        return not isinstance(other, SParts) or self.N != other.N or cmp_listsets(self.sParts, other.sParts) != 0
    # !! if not the same length or not the same total, set lists are not comparable

    def __lt__(self, other):
        if isinstance(other, SParts) and self.N == other.N:
            c = cmp_listsets(self.sParts, other.sParts)
            return c > -3 and c < 0
        return False

    def __le__(self, other):
        if isinstance(other, SParts) and self.N == other.N:
            c = cmp_listsets(self.sParts, other.sParts)
            return c > -3 and c <= 0
        return False

    def __gt__(self, other):
        if isinstance(other, SParts) and self.N == other.N:
            c = cmp_listsets(self.sParts, other.sParts)
            return c < 3 and c > 0
        return False

    def __ge__(self, other):
        if isinstance(other, SParts) and self.N == other.N:
            c = cmp_listsets(self.sParts, other.sParts)
            return c < 3 and c >= 0
        return False

    def toDict(self, with_Eoo=False):
        sdict = {}
        for i in range(self.nbStored()):
            sdict[self.ssetts.getLabel(i)] = self.part(i)
            sdict["card_" + self.ssetts.getLabel(i)] = self.lpart(i)
            sdict["perc_" + self.ssetts.getLabel(i)] = self.lpart(i) * 100. / self.N
        if with_Eoo:
            sdict[self.ssetts.getLabel(SSetts.Eoo)] = self.part(SSetts.Eoo)
            sdict["card_" + self.ssetts.getLabel(SSetts.Eoo)] = self.lpart(SSetts.Eoo)
            sdict["perc_" + self.ssetts.getLabel(SSetts.Eoo)] = self.lpart(SSetts.Eoo) * 100. / self.N
        for side in [0, 1]:
            if self.prs[side] != -1:
                sdict["pr_" + str(side)] = self.prs[side]
        sdict["N"] = self.N
        for info_key, info_meth in SParts.info_what.items():
            sdict[info_key] = eval(info_meth)
        return sdict

    # return support (used to create new instance of SParts)
    def supparts(self):
        return self.sParts

    # return new instance of SParts corresponding to negating given side
    def negate(self, side=0):
        if self.missing:
            return SParts(self.ssetts, self.N, self.ssetts.negateParts(side, self.sParts))
        else:
            self.sParts.append(self.part(self.ssetts.Eoo))
            n = self.ssetts.negateParts(side, self.sParts)
            return SParts(self.ssetts, self.N, n[0:-1])

    def part(self, part_id, side=0):
        pid = self.ssetts.partId(part_id, side)
        if pid < self.nbStored():
            return self.sParts[pid]
        elif part_id == self.ssetts.Eoo:
            return set(range(self.N)) - self.sParts[0] - self.sParts[1] - self.sParts[2]
        else:
            return set()

    def lpart(self, part_id, side=0):
        pid = self.ssetts.partId(part_id, side)
        if pid < self.nbStored():
            return len(self.sParts[pid])
        elif part_id == self.ssetts.Eoo:
            return self.N - len(self.sParts[0]) - len(self.sParts[1]) - len(self.sParts[2])
        else:
            return 0

    def parts(self, side=0):
        return [self.part(i, side) for i in range(self.ssetts.top+1)]

    def parts4M(self, side=0):
        if self.missing:
            return [self.part(i, side) for i in range(self.ssetts.Eoo+1)]+[set().union(*[self.part(i, side) for i in range(self.ssetts.Eoo+1, self.ssetts.top+1)])]
        else:
            return self.parts(side)

    def partInterX(self, suppX, part_id, side=0):
        pid = self.ssetts.partId(part_id, side)
        if pid < self.nbStored():
            return set(suppX & self.sParts[pid])
        elif part_id == self.ssetts.Eoo:
            return set(suppX - self.sParts[0] - self.sParts[1] - self.sParts[2])
        else:
            return set()

    def lpartInterX(self, suppX, part_id, side=0):
        pid = self.ssetts.partId(part_id, side)
        if pid < self.nbStored():
            return len(suppX & self.sParts[pid])
        elif part_id == self.ssetts.Eoo:
            return len(suppX - self.sParts[0] - self.sParts[1] - self.sParts[2])
        else:
            return 0

    def partsInterX(self, suppX, side=0):
        return [self.partInterX(suppX, i, side) for i in range(self.ssetts.top+1)]
        
    def lpartsInterX(self, suppX, side=0):
        if self.missing:
            return [self.lpartInterX(suppX, i, side) for i in range(self.ssetts.top+1)]
        else:
            la = self.lpartInterX(suppX, self.ssetts.Exo, side)
            lb = self.lpartInterX(suppX, self.ssetts.Eox, side)
            lc = self.lpartInterX(suppX, self.ssetts.Exx, side)
            tmp = [la, lb, lc, len(suppX) - la - lb - lc]
            for i in range(len(tmp), self.ssetts.top+1):
                tmp.append(0)
            return tmp

    def part_union(self, ids, side=0):
        union = set()
        for i in ids:
            union |= self.part(i, side)
        return union

    def supp(self, side=0):
        return self.part_union(self.ssetts.IDS_supp, side)

    def nonSupp(self, side=0):
        if not self.missing:
            return set(range(self.N)) - self.supp(side)
        else:
            return self.part_union(set(range(self.ssetts.top+1)) - set(self.ssetts.IDS_supp + self.ssetts.IDS_miss), side)

    def miss(self, side=0):
        if not self.missing:
            return set()
        else:
            return self.part_union(self.ssetts.IDS_miss, side)

    # SUPPORTS
    def suppSide(self, side):
        return self.part_union(self.ssetts.IDS_suppL, side)

    def suppP(self, i, side=0):
        return self.part(i, side)

    def suppD(self, side=0):
        return self.part_union(self.ssetts.IDS_diff, side)

    def suppI(self, side=0):
        return self.part_union(self.ssetts.IDS_inter, side)

    def suppU(self, side=0):
        return self.part_union(self.ssetts.IDS_inter+self.ssetts.IDS_diff, side)

    def suppL(self, side=0):
        return self.suppSide(0)

    def suppR(self, side=0):
        return self.suppSide(1)

    def suppO(self, side=0):
        return self.part_union(self.ssetts.IDS_uncovered, side)

    def suppA(self, side=0):
        return self.part_union(self.ssetts.IDS_suppL, side)

    def suppB(self, side=0):
        return self.part_union(self.ssetts.IDS_suppL, 1-side)

    def suppN(self, side=0):
        if self.nbStored() == 4:
            return self.part_union(range(4), side)
        else:
            return set(range(self.N))

    def getProp(self, what, which=None):
        if what in SParts.info_what:
            return eval(SParts.info_what[what])
        wt = what
        if what == "card" or what == "area":
            wt = "len"
        elif what == "supp":
            wt = "set"
        methode = eval("self.%s" % wt)
        if callable(methode):
            return methode(which)

    def set(self, which="I"):
        if which in SSetts.map_label_part:
            return self.suppP(SSetts.map_label_part[which])
        elif which in self.sets_letters:
            return eval("self.supp%s()" % which)

    # moves the instersection of supp with part with index id_from to part with index id_to
    def moveInter(self, side, id_from, id_to, supp):
        self.sParts[self.ssetts.partId(id_to, side)] |= (self.sParts[self.ssetts.partId(id_from, side)] & supp)
        self.sParts[self.ssetts.partId(id_from, side)] -= supp

    def moveInterAllOut(self, side, supp):
        out_id = self.getSSetts().Eoo
        lparts = []
        for i in range(self.nbStored()):
            lparts.append(len(self.sParts[i] & supp))
            if i != out_id:
                self.sParts[i].difference_update(supp)
            else:
                self.sParts[i].update(supp)
        if len(lparts) == out_id:
            lparts.append(len(supp)-sum(lparts))
        return lparts

    # update supports and probabilities resulting from appending X to given side with given operator

    def update(self, side, OR, suppX, missX=None):
        self.vect = None
        union = None
        self.prs[side] = SParts.updateProba(self.prs[side], len(suppX)/self.N, OR)

        if not self.missing and (type(missX) == set and len(missX) > 0):
            self.missing = True
            if self.nbStored() == 3:
                self.sParts.append(set(range(self.N)) - self.sParts[0] - self.sParts[1] - self.sParts[2])
            else:
                union = set(self.sParts[0] | self.sParts[1] | self.sParts[2] | self.sParts[3])
            self.sParts.extend([set(), set(), set(), set(), set()])

        if self.missing and self.ssetts.top > self.ssetts.Eoo:
            if OR:  # OR
                ids_from_to_supp = [(self.ssetts.Eox, self.ssetts.Exx), (self.ssetts.Eoo, self.ssetts.Exo),
                                    (self.ssetts.Emx, self.ssetts.Exx), (self.ssetts.Emo, self.ssetts.Exo),
                                    (self.ssetts.Eom, self.ssetts.Exm), (self.ssetts.Emm, self.ssetts.Exm)]
                for (id_from, id_to) in ids_from_to_supp:
                    self.moveInter(side, id_from, id_to, suppX)

                if (type(missX) == set and len(missX) > 0):
                    ids_from_to_miss = [(self.ssetts.Eox, self.ssetts.Emx), (self.ssetts.Eoo, self.ssetts.Emo),
                                        (self.ssetts.Eom, self.ssetts.Emm)]
                    for (id_from, id_to) in ids_from_to_miss:
                        self.moveInter(side, id_from, id_to, missX)

            else:  # AND
                if (type(missX) == set and len(missX) > 0):
                    suppXB = set(range(self.N)) - suppX - missX
                else:
                    suppXB = set(range(self.N)) - suppX
                ids_from_to_suppB = [(self.ssetts.Exo, self.ssetts.Eoo), (self.ssetts.Exx, self.ssetts.Eox),
                                     (self.ssetts.Exm, self.ssetts.Eom), (self.ssetts.Emx, self.ssetts.Eox),
                                     (self.ssetts.Emo, self.ssetts.Eoo), (self.ssetts.Emm, self.ssetts.Eom)]
                for (id_from, id_to) in ids_from_to_suppB:
                    self.moveInter(side, id_from, id_to, suppXB)

                if (type(missX) == set and len(missX) > 0):
                    ids_from_to_miss = [(self.ssetts.Exo, self.ssetts.Emo), (self.ssetts.Exx, self.ssetts.Emx),
                                        (self.ssetts.Exm, self.ssetts.Emm)]
                    for (id_from, id_to) in ids_from_to_miss:
                        self.moveInter(side, id_from, id_to, missX)

        else:
            if OR:  # OR
                self.sParts[self.ssetts.partId(self.ssetts.Exo, side)] |= (suppX
                                                                           - self.sParts[self.ssetts.partId(self.ssetts.Eox, side)]
                                                                           - self.sParts[self.ssetts.partId(self.ssetts.Exx, side)])
                self.sParts[self.ssetts.partId(self.ssetts.Exx, side)] |= (suppX
                                                                           & self.sParts[self.ssetts.partId(self.ssetts.Eox, side)])
                self.sParts[self.ssetts.partId(self.ssetts.Eox, side)] -= suppX

            else:  # AND
                self.sParts[self.ssetts.partId(self.ssetts.Eox, side)] |= (self.sParts[self.ssetts.partId(self.ssetts.Exx, side)]
                                                                           - suppX)
                self.sParts[self.ssetts.partId(self.ssetts.Exx, side)] &= suppX
                self.sParts[self.ssetts.partId(self.ssetts.Exo, side)] &= suppX
        if union is not None:
            self.sParts[self.ssetts.Eoo] = union - self.sParts[self.ssetts.Exx] - self.sParts[self.ssetts.Eox] - self.sParts[self.ssetts.Exo]

    # computes vector ABCD (vector containg for each row the index of the part it belongs to)
    def makeVectorABCD(self, force_list=False, rest_ids=None):
        if self.vect is None or (force_list and type(self.vect) is not list):
            if self.nbStored() == 4 and not force_list:
                # svect = {}
                self.vect = {}
                for partId in range(self.nbStored()):
                    for i in self.sParts[partId]:
                        self.vect[i] = partId
            else:
                self.vect = [self.ssetts.Eoo for i in range(self.N)]
                map_rest = {}
                if rest_ids is not None:
                    map_rest = dict([(vvv, vvi) for (vvi, vvv) in enumerate(sorted(rest_ids))])
                for partId in range(self.nbStored()):
                    for i in self.sParts[partId]:
                        self.vect[map_rest.get(i, i)] = partId

    def getVectorABCD(self, force_list=False, rest_ids=None):
        self.makeVectorABCD(force_list, rest_ids)
        if type(self.vect) is dict:
            return None
        return list(self.vect)

    # returns the index of the part the given row belongs to, vectorABCD need to have been computed
    def partRow(self, row):
        return self.vect[row]

    # return the index of the part the given row belongs to
    # or the intersection of the mode of X with the different parts if row == -1, vectorABCD need to have been computed
    def lpartsRow(self, row, X=None):
        lp = None
        if row == -1 and X is not None:
            if self.missing:
                lp = [len(X.interMode(self.sParts[i])) for i in range(self.ssetts.top+1)]
            else:
                lp = [0 for i in range(self.nbParts())]
                lp[0] = len(X.interMode(self.sParts[0]))
                lp[1] = len(X.interMode(self.sParts[1]))
                lp[2] = len(X.interMode(self.sParts[2]))
                lp[3] = X.lenMode() - lp[0] - lp[1] - lp[2]
        elif row is not None:
            lp = self.vect[row]
        return lp

    # PRINTING
    ##############
    def dispSupp(self, sep="\t"):
        supportStr = ""
        for i in sorted(self.supp(0)):
            supportStr += "%i " % i
        supportStr += sep
        for i in sorted(self.supp(1)):
            supportStr += "%i " % i
        if self.missing:
            supportStr += sep
            for i in sorted(self.miss(0)):
                supportStr += "%i " % i
            supportStr += sep
            for i in sorted(self.miss(1)):
                supportStr += "%i " % i
        return supportStr

    # compute the resulting support and missing when combining X and Y with given operator
    @classmethod
    def partsSuppMiss(tcl, OR, XSuppMiss, YSuppMiss):
        if XSuppMiss is None:
            return YSuppMiss
        elif YSuppMiss is None:
            return XSuppMiss
        elif OR:
            supp = set(XSuppMiss[0] | YSuppMiss[0])
            miss = set(XSuppMiss[1] | YSuppMiss[1]) - supp
        else:
            miss = set((XSuppMiss[1] & YSuppMiss[1]) | (XSuppMiss[1] & YSuppMiss[0]) | (YSuppMiss[1] & XSuppMiss[0]))
            supp = set(XSuppMiss[0] & YSuppMiss[0])
        return (supp, miss)

    @classmethod
    def partsSuppMissMass(tcl, OR, SuppMisses):
        if len(SuppMisses) == 1:
            return SuppMisses[0]
        elif len(SuppMisses) > 1:
            if OR:
                supp = reduce(set.union, [X[0] for X in SuppMisses])
                miss = reduce(set.union, [X[1] for X in SuppMisses]) - supp
            else:
                supp = reduce(set.intersection, [X[0] for X in SuppMisses])
                miss = reduce(set.intersection, [X[0].union(X[1]) for X in SuppMisses]) - supp
            return (supp, miss)

    # Make binary out of supp set
    @classmethod
    def suppVect(tcl, N, supp, val=1):
        vect = None
        if 2*len(supp) < N:
            st = supp
            v = val
            if val == 1:
                vect = numpy.zeros(N)
            else:
                vect = numpy.ones(N)
        else:
            st = set(range(N)) - supp
            v = 1-val
            if val == 0:
                vect = numpy.zeros(N)
            else:
                vect = numpy.ones(N)
        for i in st:
            vect[i] = v
        return vect

    @classmethod
    def parseSupport(tcl, stringSupp, N, ssetts):
        partsSupp = stringSupp.rsplit("\t")
        if len(partsSupp) == 2:
            return tcl(ssetts, N, [tcl.parseSupportPart(partsSupp[0]), tcl.parseSupportPart(partsSupp[1])])
        elif len(partsSupp) == 4:
            return tcl(ssetts, N, [tcl.parseSupportPart(partsSupp[0]), tcl.parseSupportPart(partsSupp[1]),
                                   tcl.parseSupportPart(partsSupp[2]), tcl.parseSupportPart(partsSupp[3])])
        return None

    @classmethod
    def parseSupportPart(tcl, string):
        nsupp = set()
        for i in string.strip().rsplit():
            try:
                nsupp.add(int(i))
            except TypeError as detail:
                raise Exception("Unexpected element in the support: %s\n" % i)
        return nsupp


if __name__ == "__main__":
    pass
    # sss = {}
    # ks = ["none", "rejective", "optimistic", "pessimistic", "positive", "negative"]  # , "exclu"]
    # for k in ks:
    #     sss[(k, False)] = str(SSetts(False, k))
    #     sss[(k, True)] = str(SSetts(True, k))
    #     # print("\n### %s (no miss)" % k, SSetts(False, k))
    # for ki in range(len(ks)):
    #     if sss[(ks[ki], False)] == sss[(ks[ki], True)]:
    #         print(" --- ", ks[ki])
    #     for kj in range(ki):
    #         if sss[(ks[ki], False)] == sss[(ks[kj], False)]:
    #             print("no miss", ks[ki], ks[kj])
    #         if sss[(ks[ki], True)] == sss[(ks[kj], True)]:
    #             print("missing", ks[ki], ks[kj])

    # for k in ks[:1]:
    #     # print("\n### %s (missing)" % k, SSetts(True, k))
    #     print("\n\n    def init_part_ids_%s(self):\n%s" % (k, SSetts(False, k).dispPartsDef()))

    # print("\n")
    # for k in ks[1:]:
    #     # print("\n### %s (missing)" % k, SSetts(True, k))
    #     print("\n\n    def init_part_ids_%s(self):\n%s" % (k, SSetts(True, k).dispPartsDef()))

    #     ssetts = SSetts()
    #     spartsA = SParts(ssetts, 200, [range(50), range(50, 100), range(100, 150)])
    #     spartsB = SParts(ssetts, 300, [range(50), range(50, 100), range(100, 150)])
    #     spartsC = SParts(ssetts, 200, [range(10), range(10, 100), range(100, 150)])
    #     assert(spartsA==spartsA)
    #     assert(not spartsA>=spartsB); assert(not spartsA<spartsB) # not comparable
    #     assert(not spartsA<spartsC); assert(spartsA>spartsC) # comparable
