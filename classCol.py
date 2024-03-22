import os.path
import time
import copy
import re
import numpy

try:
    from toolXML import BOOL_MAP
    from classQuery import Op, Term, AnonTerm, ConstantTerm, BoolTerm, CatTerm, NumTerm, Literal, TimeTools, NA_str_c
    from classSParts import tool_ratio, cmp_vals, cmp_lists
    from classProps import WithEVals, VarProps, mapSuppNames, ACTIVE_RSET_ID
    import bellman
    # import misc_codes
except ModuleNotFoundError:
    from .toolXML import BOOL_MAP
    from .classQuery import Op, Term, AnonTerm, ConstantTerm, BoolTerm, CatTerm, NumTerm, Literal, TimeTools, NA_str_c
    from .classSParts import tool_ratio, cmp_vals, cmp_lists    
    from .classProps import WithEVals, VarProps, mapSuppNames, ACTIVE_RSET_ID
    from . import bellman
    # from . import misc_codes

import pdb

NA_num = numpy.nan
NA_bool = -1
NA_cat = -1

MODE_VALUE = 0
MODE_RATIO = 0.1


def associate_term_class(col_class, term_class=None):
    col_class.assoc_term = term_class
    if term_class is not None:
        col_class.type_id = term_class.type_id
        col_class.type_letter = term_class.type_letter
        col_class.type_name = term_class.type_name


def getValSpec(param, i):
    out = param
    if type(param) is list:
        out = param[i]
    elif type(param) is dict:
        if i in param:
            out = param[i]
        elif None in param:
            out = param[None]
    return out


def mapSuppNames(supp, details={}):
    if details.get("named", False) and "row_names" in details:
        return [details["row_names"][t] for t in supp]
    return supp


class DataError(Exception):
    def __init__(self, value):
        self.value = value

    def __str__(self):
        return repr(self.value)


class ColM(WithEVals):

    class_letter = "v"
    ##############################################################
    ##############################################################
    # PROPS WHAT
    info_what_dets = {}
    info_what_mask = {"set": "set"}
    info_what = {"type": "self.getType()", "missing": "self.getMissInfo()",
                 "density": "self.getDensity()",
                 "categories": "self.getCategories()",
                 "min": "self.getMin()", "max": "self.getMax()"}
    Pwhat_match = "(" + "|".join(["[a-zA-Z]+", "extra"]+list(info_what.keys()) +
                                 list(info_what_dets.keys())+list(info_what_mask.keys())) + ")"

    info_which_mask = {"vect": "self.getVect", "vals": "self.getVals"}
    # PROPS WHICH
    Pwhich_match = "(" + "|".join(list(info_which_mask.keys())+[WithEVals.which_rids]) + ")"

    RP = None
    @classmethod
    def setupRP(tcl, fields_fns=None):
        elems_typs = [("v", ColM)]
        VarProps.setupProps(ColM, elems_typs)
        tcl.RP = VarProps(fields_fns)

    # filtering
    def getProp(self, what, which=None, rset_id=None, details={}):
        if what == "extra":
            return self.getExtra(which, details)
        if which == self.which_rids:  # ids details for folds subsets
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
        elif which in self.info_which_mask:
            rset_ids = self.getRestrictedRids(rset_id, as_list=True)
            methode_which = eval(self.info_which_mask[which])
            if callable(methode_which):
                x = methode_which(rset_ids, details)
                if what in self.info_what_mask:
                    methode = eval(self.info_what_mask[what])
                else:
                    methode = eval(what)
                if callable(methode):
                    return methode(x)
        else:
            return self.getPropD(what, details=details)

    def getExpProp(self, exp, details={}):
        ws = self.getRP().getPrimitiveWs(exp)
        if ws[0] is not None:
            return self.getProp(ws[0], ws[1], ws[2], details)

    def setRestrictedSuppSets(self, data, supp_sets=None):
        resets_ids = []
        if supp_sets is None:
            if data.hasLT():
                supp_sets = data.getLT()
            else:
                supp_sets = {ACTIVE_RSET_ID: sorted(data.nonselectedRows())}
        for sid, sset in supp_sets.items():
            if sid not in self.restricted_sets or self.restricted_sets[sid]["rids"] != sset:
                self.restricted_sets[sid] = {"rids": sset}
                resets_ids.append(sid)

        if len(resets_ids) > 0:
            return True
        return False
    ##############################################################
    ##############################################################

    width = 0
    typespec_placeholder = "<!-- TYPE_SPECIFIC -->"
    NA = NA_bool
    NA_specimen_str = ["na", "nan", "-", "-1"]
    anon_term = AnonTerm

    @classmethod
    def getAnonTermClass(tcl):
        return tcl.anon_term

    @classmethod
    def getAssocTermClass(tcl):
        return tcl.assoc_term

    @classmethod
    def initSums(tcl, N):
        return [0 for i in range(N)]

    @classmethod
    def parseList(tcl, listV, indices=None, force=False, native_missing_check=True):
        return None

    @classmethod
    def fromVect(tcl, vect_data, prec=None, enabled=True):
        return None

    def __init__(self, N=-1, nmiss=set()):
        WithEVals.__init__(self)
        if nmiss is None:
            nmiss = set()
        self.N = N
        self.missing = nmiss
        self.infofull = {"in": (-1, True), "out": (-1, True)}
        self.vect = None
        self.extras.update({"status": 1, "side": -1, "id": -1, "gid": -1})

    def cmpCol(self, other):
        if other is None or not isinstance(other, ColM):
            return 1
        else:
            return cmp_vals(self.getId(), other.getId())

    def cmpType(self, other):
        if other is None or not isinstance(other, ColM):
            return 1
        else:
            return cmp_vals(self.typeId(), other.typeId())

    def cmpVals(self, other):
        if other is None or not isinstance(other, ColM):
            return 3
        if self.N == other.N:
            return cmp_lists(self.getVect(), other.getVect())
        elif self.N < other.N:
            return -2
        return 2

    def getUid(self):
        return (self.getSide(), self.getId())

    def setUid(self, iid=None):
        pass
        # raise Warning("Col UID is ready-only!")

    def setSideIdName(self, side, cid, name=None):
        self.setId(cid)
        self.extras["side"] = side
        if name is None:
            self.extras["name"] = Term.pattVName % cid
        else:
            self.extras["name"] = name

    def simpleBool(self):
        return False

    def nbRows(self):
        return self.N

    def rows(self):
        return set(range(self.nbRows()))

    def setGroupId(self, gid=-1):
        self.extras["gid"] = gid

    def getGroupId(self):
        return self.extras.get("gid", -1)

    def hasGroup(self):
        return self.getGroupId() > -1

    def setId(self, nid):
        self.extras["id"] = nid

    def hasMissing(self):
        return self.missing is not None and len(self.missing) > 0

    def nbMissing(self):
        if self.missing is not None:
            return len(self.missing)
        return 0

    def valToStr(self, val):
        if val == self.NA:
            return NA_str_c
        return val

    def areDataEquiv(self, vA, vB):
        return vA == vB

    def getPrec(self, details={}):
        return 0

    def getTimePrec(self, details={}):
        return -1

    def getFmt(self, details={}):
        return {"prec": self.getPrec(details), "time_prec": self.getTimePrec(details)}

    def density(self):
        return 1.0

    def minGap(self):
        return 0

    def isDense(self, thres=None):
        if thres is None:
            thres = 0.5
        return self.density() > thres

    def getName(self, details={}):
        if self.hasName():
            return self.extras["name"]
        else:
            return Term.pattVName % self.getId()

    def hasName(self):
        return self.extras.get("name") is not None

    def getSide(self, details={}):
        return self.extras.get("side")

    def getId(self, details={}):
        return self.extras.get("id")

    def getAnonTerm(self):
        return self.getAnonTermClass()(self.getId(), self.typeId())

    def upSumsRows(self, sums_rows):
        pass

    def sumCol(self):
        return 0

    def numEquiv(self, v):
        try:
            return int(v)
        except:
            pass
        return self.NA

    def mkVector(self):
        self.vect = numpy.ones(self.N, dtype=numpy.int)*self.NA

    def getVector(self, bincats=False, nans=None):
        if self.vect is None:
            self.mkVector()
        if self.hasMissing() and nans is not None and \
                not ((numpy.isnan(nans) and numpy.isnan(self.NA)) or nans == self.NA):  # Not the same nan...
            tmp = numpy.array(self.vect, dtype=numpy.float, copy=True)
            tmp[tmp == self.NA] = nans
            return tmp
        return self.vect

    def getSortAble(self, details={}):
        if details.get("aim") == "sort":
            return (self.getEnabled(), self.getId())
        return ""

    def getType(self, details={}):
        return self.type_name

    def getDensity(self, details={}):
        return None

    def getCategories(self, details={}):
        return None

    def getMin(self, details={}):
        return None

    def getMax(self, details={}):
        return None

    def getMissInfo(self, details={}):
        return "%d%%: %d" % (100*self.nbMissing()/float(self.N), self.nbMissing())

    def getOrd(self):
        return []

    def getRange(self):
        return dict([(k, v) for (v, k) in enumerate(self.getOrd())])

    # def getCohesion(self, details={}):
    #     return "%1.4f" % self.cohesion
    def getVect(self, mask=None, details={}):
        vect = self.getVector()
        if mask is not None:
            return vect[mask]
        return vect

    def getVals(self, mask=None, details={}):
        return [self.getValFromNum(v) for v in self.getVect(mask, details)]

    def typeId(self):
        return self.type_id

    def miss(self):
        return self.missing

    def negSuppTerm(self, term):
        return self.rows() - self.suppTerm(term) - self.miss()

    def suppLiteral(self, literal):  # TODO:literal has a time constraint added!
        if isinstance(literal, Term):  # It's a term, not a literal
            return self.suppTerm(literal)
        elif isinstance(literal, Literal):
            if literal.isNeg():
                return self.negSuppTerm(literal.getTerm())
            else:
                return self.suppTerm(literal.getTerm())

    def lMiss(self):
        return len(self.miss())

    def lSuppTerm(self, term):  # column specific implementation
        return 0

    def lNegSuppTerm(self, term):
        return self.nbRows() - self.lsuppTerm(term) - self.lMiss()

    def lSuppLiteral(self, literal):
        if isinstance(literal, Term):  # It's a literal, not a term
            return self.lsuppTerm(literal)
        elif isinstance(literal, Literal):
            if literal.isNeg():
                return self.lNegSuppTerm(literal.getTerm())
            else:
                return self.lsuppTerm(literal.getTerm())

    def getEnabled(self, details={}):
        return self.extras["status"]

    def isEnabled(self, details={}):
        return self.getEnabled() > 0

    def flipEnabled(self):
        self.extras["status"] = 1-self.extras["status"]

    def setStatus(self, status=True):
        self.extras["status"] = int(status)

    def setEnabled(self):
        self.extras["status"] = 1

    def setDisabled(self):
        self.extras["status"] = 0

    def __str__(self):
        act = ""
        if not self.isEnabled():
            act = " (OFF)"
        if self.hasGroup():
            act += " [gid=%d]" % self.getGroupId()
        # return "%s variable %i %s%s, %d missing values" % (self.getType(), self.getId(), self.getName().encode('ascii', 'replace'), act, self.lMiss())
        return "%s variable %i %s%s, %d missing values" % (self.getType(), self.getId(), self.getName(), act, self.lMiss())

    def suppInBounds(self, min_in=-1, min_out=-1):
        return (self.infofull["in"][1] and self.infofull["out"][1])

    def usable(self, min_in=-1, min_out=-1, checkable=True):
        return self.suppInBounds(min_in, min_out) and (not checkable or self.isEnabled())


associate_term_class(ColM, Term)

# For evaluating and getting support of constant term (colId = -1)


class ConstantColM(ColM):

    def __init__(self, N=-1):
        ColM.__init__(self, N)

    def suppTerm(self, term):
        if isinstance(term, ConstantTerm) and term.truthEval():
            return set(range(self.N))
        return set()

    def lsuppTerm(self, term):
        if isinstance(term, ConstantTerm) and term.truthEval():
            return self.N
        return set()


class BoolColM(ColM):
    width = -1
    values_eq = {True: 1, False: 0}
    NA = NA_bool

    values = BOOL_MAP

    @classmethod
    def parseList(tcl, listV, indices=None, force=False, native_missing_check=True):
        if type(listV) is list and len(set(listV).difference([True, False, None])) == 0:
            miss = set([i for (i, v) in enumerate(listV) if v is None])
            trues = set([i for (i, v) in enumerate(listV) if v == True])
            return BoolColM(trues, len(listV), miss)
        miss = set()
        if force:
            if type(listV) is list:
                miss = set([i for (i, v) in enumerate(listV) if v is None])
                listV = set([i for (i, v) in enumerate(listV)
                             if v is not None and BoolColM.values.get(v.lower(), True)])
            elif type(listV) is not set:
                tt = set()
                ok = True
                for idx, v in listV.items():
                    try:
                        if float(v) != 0:
                            tt.add(idx)
                    except ValueError:
                        ok = False
                if ok:
                    listV = tt
        if type(listV) is set:
            if type(indices) is int:
                trues = set(indices)
                N = indices
            elif type(indices) is dict:
                trues = set([indices.get(i, None) for i in listV])
                trues.discard(None)
                miss = set([indices.get(i, None) for i in miss])
                miss.discard(None)
                N = max(indices.values())+1
            else:
                raise ValueError('Sparse requires indices')
            return BoolColM(trues, N, miss)
        if indices is None:
            indices = dict([(v, v) for v in range(len(listV))])
        trues = set()
        miss = set()
        if type(listV) is dict:
            ttt = set(listV.keys()).intersection(indices.keys())
        else:
            ttt = [i for i in indices.keys() if i < len(listV)]

        val_nonbool = set([listV[i].lower() for i in ttt if listV[i] is not None]
                          ).difference(BoolColM.values.keys())
        val_na = val_nonbool.intersection(BoolColM.NA_specimen_str)
        na_v = BoolColM.NA
        if len(val_na) == 1:
            na_v = val_na.pop()
        elif len(val_nonbool) > 0:
            return None

        for i in ttt:
            j = indices[i]
            if listV[i] is None or listV[i].lower() == na_v:
                miss.add(j)
            else:
                v = listV[i].lower()
                if v not in BoolColM.values:
                    return None
                elif BoolColM.values[v]:
                    trues.add(j)

        return BoolColM(trues, max(indices.values())+1, miss)

    @classmethod
    def fromVect(tcl, vect_data, prec=None, enabled=True):
        tmp = vect_data
        if prec is not None:
            tmp = numpy.around(tmp, prec)

        col = tcl(set([i for (i, v) in enumerate(tmp) if v != 0]), len(tmp))
        if not enabled:
            col.flipEnabled()
        return col

    def toList(self, sparse=False, fill=False, cats_int_fmt=False):
        if sparse:
            t = int(True)
            tmp = [(n, t) for n in self.hold]+[(n, self.NA) for n in self.missing]
            if fill and self.N-1 not in self.hold and self.N-1 not in self.missing:
                tmp.append((self.N-1, int(False)))
            return tmp
        else:
            # return map(self.valToStr, self.getVector())
            return self.getVector()

    def density(self):
        if self.N == self.nbMissing():
            return 0.0
        else:
            return self.sumCol()/float(self.N-self.nbMissing())

    def minGap(self):
        return 1.

    def getInitTerms(self, minIn=0, minOut=0, productivity="medium"):
        if self.sumCol() >= minIn and self.N-(self.sumCol()+self.nbMissing()) >= minOut:
            return [(self.getAssocTermClass()(self.getId()), self.sumCol())]
        else:
            return []

    def simpleBool(self):
        return not self.hasMissing() and self.density() > 0

    def __str__(self):
        return ColM.__str__(self) + (", %i Trues" % (self.lTrue()))

    def upSumsRows(self, sums_rows):
        for i in self.hold:
            sums_rows[i] += 1

    def sumCol(self):
        return len(self.hold)

    def getOrd(self):
        return [False, True]

    def getNbValues(self):
        return 2

    def mkVector(self):
        self.vect = numpy.ones(self.N, dtype=numpy.int)*self.numEquiv(False)
        self.vect[list(self.missing)] = self.NA
        self.vect[list(self.hold)] = self.numEquiv(True)

    def getDensity(self, details={}):
        if self.N > 0:
            return self.density()
            # return "%1.4f" % self.density()
        return 0

    def __init__(self, ncolSupp=set(), N=-1, nmiss=set()):
        ColM.__init__(self, N, nmiss)
        self.hold = ncolSupp
        self.missing -= self.hold
        if self.N == len(self.missing):
            print("Boolean variable contains only missing values, this is suspect!..")
        elif len(self.hold) == 0:
            print("Boolean variable contains only False entries, this is suspect!..")
        elif len(self.hold) == (self.N - len(self.missing)):
            print("Boolean variable contains only True entries, this is suspect!..")

    def subsetCol(self, row_ids=None):
        if row_ids is None:
            hold = set(self.hold)
            miss = set(self.missing)
            N = self.nbRows()
        else:
            miss = set()
            hold = set()
            N = sum([len(news) for news in row_ids.values()])
            for old in self.missing.intersection(row_ids.keys()):
                miss.update(row_ids[old])
            for old in self.hold.intersection(row_ids.keys()):
                hold.update(row_ids[old])
        tmp = BoolColM(hold, N, miss)
        tmp.extras = self.copyExtras()
        tmp.infofull = {"in": tuple(self.infofull["in"]), "out": tuple(self.infofull["out"])}
        return tmp

    def getValue(self, rid, pref=None):
        if self.vect is None:
            if rid in self.missing:
                return self.NA

            if pref == "bnum":
                return int(rid in self.hold)
            return rid in self.hold
        else:
            if pref == "bnum":
                return self.vect[rid]
            return BoolColM.values.get(self.vect[rid], self.NA)

    def getNumValue(self, rid):
        return int(self.getValue(rid))

    def getValFromNum(self, n):
        return n == 1

    def supp(self):
        return self.hold

    def suppTerm(self, term):
        if term is not None and term.isAnon():
            return set()
        return set(self.hold)

    def lsuppTerm(self, term):
        if term is not None and term.isAnon():
            return 0
        return len(self.hold)

    def lTrue(self):
        return self.sumCol()

    def lFalse(self):
        return self.nbRows() - self.lTrue() - self.lMiss()

    def suppInBounds(self, min_in=-1, min_out=-1):
        if self.infofull["in"][0] != min_in:
            self.infofull["in"] = (min_in, self.lTrue() >= min_in)
        if self.infofull["out"][0] != min_out:
            self.infofull["out"] = (min_out, self.lFalse() >= min_out)
        return (self.infofull["in"][1] and self.infofull["out"][1])


associate_term_class(BoolColM, BoolTerm)


class CatColM(ColM):
    width = 1
    NA = NA_cat
    CATS_INT_FMT = "c%d"

    def getLiteralCat(self, neg, best_cat, allw_neg=True):
        if CatTerm.is_multi_cat(best_cat):
            in_cats = set()
            out_cats = set(self.cats())
            for c in best_cat:
                if self.isCatId(c):
                    in_cats.add(self.getCatForId(c))
                    out_cats.discard(self.getCatForId(c))
                elif c in out_cats:
                    in_cats.add(c)
                    out_cats.discard(c)                    
            if len(in_cats) == 0 or len(in_cats) == self.nbCats():  # unbounded
                return None
            if allw_neg and (len(out_cats) < len(in_cats)): ## actually fewer categories out, best negate and list out
                return Literal(not neg, self.getAssocTermClass()(self.getId(), out_cats))
            else:
                return Literal(neg, self.getAssocTermClass()(self.getId(), in_cats))
        if self.isCatId(best_cat):
            return Literal(neg, self.getAssocTermClass()(self.getId(), self.getCatForId(best_cat)))
        elif best_cat in self.cats():
            return Literal(neg, self.getAssocTermClass()(self.getId(), best_cat))
        return None # best_cat is not among the categories of this column

    def __init__(self, ncolSupp={}, N=-1, nmiss=set()):
        ColM.__init__(self, N, nmiss)
        self.sCats = ncolSupp
        self.ord_cats = sorted(self.sCats.keys())
        self.cards = sorted([(cat, self.lsuppCat(cat))
                                 for cat in self.cats()], key=lambda x: x[1])
        if self.N == len(self.missing):
            print("Categorical variable contains only missing values, this is suspect!..")
        elif len(self.ord_cats) == 1:
            print("Categorical variable has only one category %s, this is suspect!.." % (self.ord_cats))
        elif len(self.ord_cats) == 0:
            print("Categorical variable has no categories, this is suspect!..")

    @classmethod
    def initSums(tcl, N):
        return [{} for i in range(N)]

    @classmethod
    def parseList(tcl, listV, indices=None, force=False, native_missing_check=True):
        if indices is None:
            indices = dict([(v, v) for v in range(len(listV))])
        cats = {}
        miss = set()
        if type(listV) is dict:
            ttt = set(listV.keys()).intersection(indices.keys())
        else:
            ttt = [i for i in indices.keys() if i < len(listV)]
        for i in ttt:
            j = indices[i]
            v = listV[i]
            if v is None or (native_missing_check and v == str(tcl.NA)):
                miss.add(j)
            else:
                if v in cats:
                    cats[v].add(j)
                else:
                    cats[v] = set([j])
        if len(cats) > 0 or force:
            return tcl(cats, max(indices.values())+1, miss)
        else:
            return None

    @classmethod
    def fromVect(tcl, vect_data, prec=0, enabled=True):
        tmp = vect_data
        if prec is not None:
            tmp = numpy.around(tmp, prec)
        cats = {}
        for (j, v) in enumerate(tmp):
            if v in cats:
                cats[v].add(j)
            else:
                cats[v] = set([j])

        col = tcl(cats, N=len(tmp))
        if not enabled:
            col.flipEnabled()
        return col

    def toList(self, sparse=False, fill=False, cats_int_fmt=None):
        ### cats_int_fmt: formatting categories indices?
        ##         None or string not containing "%" -> using names,
        ##         string containing % -> using formatting string,
        ##         True -> using default cats_int_fmt,
        ##         else -> using integers
        if sparse:
            dt = [(i, NA_str_c) for i in self.missing]
            for ci, cat in enumerate(self.ord_cats):
                if cats_int_fmt is not None:
                    if cats_int_fmt == True:
                        dt.extend([(i, self.CATS_INT_FMT % ci) for i in self.sCats[cat]])
                    elif type(cats_int_fmt) is str:
                        if "%" in cats_int_fmt:
                            dt.extend([(i, cats_int_fmt % ci) for i in self.sCats[cat]])
                        else:
                            dt.extend([(i, cat) for i in self.sCats[cat]])
                    else:
                        dt.extend([(i, ci) for i in self.sCats[cat]])
                else:
                    dt.extend([(i, cat) for i in self.sCats[cat]])
            return dt
        else:
            if cats_int_fmt is not None:
                if cats_int_fmt == True:
                    return [self.CATS_INT_FMT % v if v >= 0 else NA_str_c for v in self.getVector()]
                elif type(cats_int_fmt) is str:
                    if "%" in cats_int_fmt:
                        return [cats_int_fmt % v if v >= 0 else NA_str_c for v in self.getVector()]
                    else:
                        return [self.ord_cats[v] if v >= 0 else NA_str_c for v in self.getVector()]
                else:
                    return [v if v >= 0 else NA_str_c for v in self.getVector()]
            else:
                return [self.ord_cats[v] if v >= 0 else NA_str_c for v in self.getVector()]

    def mkVector(self):
        self.vect = numpy.ones(self.N, dtype=numpy.int)*self.numEquiv(self.NA)
        for v, cat in enumerate(self.ord_cats):
            self.vect[list(self.sCats[cat])] = v

    def getVector(self, bincats=False, nans=None):
        if bincats:  # binarize the categories, i.e return a matrix rather than vector
            vect = numpy.zeros((self.N, self.nbCats()), dtype=numpy.int)
            for v, cat in enumerate(self.ord_cats):
                vect[list(self.sCats[cat]), v] = 1
            return vect

        if self.vect is None:
            self.mkVector()
        if self.hasMissing() and nans is not None and \
                not ((numpy.isnan(nans) and numpy.isnan(self.NA)) or nans == self.NA):  # Not the same nan...
            tmp = numpy.array(self.vect, dtype=numpy.float, copy=True)
            tmp[tmp == self.NA] = nans
            return tmp
        return self.vect

    def minGap(self):
        return 1

    def getInitTerms(self, minIn=0, minOut=0, productivity="medium"):
        terms = []
        for cat in self.cats():
            if len(self.sCats[cat]) >= minIn and self.N-(len(self.sCats[cat])+self.nbMissing()) >= minOut:
                terms.append((self.getAssocTermClass()(self.getId(), cat), len(self.sCats[cat])))
        return terms

    def __str__(self):
        return ColM.__str__(self) + (", %i categories" % self.nbCats()) #+ (" "+ ", ".join(["%s: %d" % (k, len(self.sCats[k])) for k in self.ord_cats]))

    def getOrd(self):
        return list(self.ord_cats)

    def getNbValues(self):
        return self.nbCats()

    def getCategories(self, details={}):
        if self.nbCats() < 5:
            return ("%d [" % self.nbCats()) + ', '.join(["%s:%d" % (catL, len(self.sCats[catL])) for catL in self.cats()]) + "]"
        else:
            return ("%d [" % self.nbCats()) + ', '.join(["%s:%d" % (catL, len(self.sCats[catL])) for catL in self.cats()[:3]]) + "...]"

    def upSumsRows(self, sums_rows):
        for cat, rows in self.sCats.items():
            for i in rows:
                sums_rows[i][cat] = sums_rows[i].get(cat, 0)+1

    def sumCol(self):
        return dict(self.cards)

    def getCatForVal(self, v, missing_str=None):
        if v != self.NA:
            try:
                vint = int(v)
                return self.ord_cats[vint]
            except:
                pass
        if missing_str is not None:
            return missing_str
        return self.NA

    def getValue(self, rid, pref=None):
        if pref == "cnum":
            return self.getNumValue(rid)
        elif pref is not None and self.getNumValue(rid) >= 0:
            if pref == True:
                return self.CATS_INT_FMT % self.getNumValue(rid)
            elif type(pref) is str:
                if "%" in pref:
                    return pref % self.getNumValue(rid)
                else:
                    return self.getCatForVal(self.getNumValue(rid))
            else:
                return self.getNumValue(rid)
        return self.getCatForVal(self.getNumValue(rid))

    def getNumValue(self, rid):
        if self.vect is None:
            self.getVector()
        if rid < len(self.vect):
            return self.vect[rid]
        else:
            return self.NA

    def numEquiv(self, v):
        if v == "#LOW#":
            return -0.5  # 0
        if v == "#HIGH#":
            return -0.5  # len(self.ord_cats)-1

        try:
            return self.ord_cats.index(v)
        except:
            return self.NA

    def subsetCol(self, row_ids=None):
        if row_ids is None:
            scats = dict(self.sCats)
            miss = set(self.missing)
            N = self.nbRows()
        else:
            miss = set()
            scats = {}
            N = sum([len(news) for news in row_ids.values()])
            for old in self.missing.intersection(row_ids.keys()):
                miss.update(row_ids[old])
            for cat, rs in self.sCats.items():
                scats[cat] = set()
                for old in rs.intersection(row_ids.keys()):
                    scats[cat].update(row_ids[old])
        tmp = CatColM(scats, N, miss)
        tmp.extras = self.copyExtras()
        tmp.infofull = {"in": tuple(self.infofull["in"]), "out": tuple(self.infofull["out"])}
        return tmp

    def modeCat(self):
        return self.cards[-1][0]
    
    def isCatId(self, cat):
        return CatTerm.is_int_cat(cat) and -1 <= cat < self.nbCats()
    def getCatForId(self, cat):
        return self.ord_cats[cat]
    
    def getValFromNum(self, n):
        if 0 <= n < self.nbCats():
            return self.ord_cats[int(n)]
        return self.NA

    def iter_cats(self):
        return list(self.sCats.items())

    def cats(self):
        return self.ord_cats

    def nbCats(self):
        return len(self.ord_cats)
    
    def getSCat(self, cat): ## no safety checks
        if CatTerm.is_int_cat(cat):
            return self.sCats[self.ord_cats[cat]]
        else:
            return self.sCats[cat]
                
    def suppCat(self, cat):
        try:
            return set(self.getSCat(cat))
        except (TypeError, KeyError) as e:
            supp = set()
            if CatTerm.is_multi_cat(cat):
                cc = cat
            else:
                cc = [cat]
            for c in cc:
                if self.isCatId(c):
                    supp.update(self.sCats[self.getCatForId(c)])
                else:
                    supp.update(self.sCats.get(c, set()))
            return supp
                
    def lsuppCat(self, cat):
        try:
            return len(self.getSCat(cat))
        except (TypeError, KeyError) as e:
            lsupp = 0
            if CatTerm.is_multi_cat(cat):
                cc = cat
            else:
                cc = [cat]
            for c in cc:
                if self.isCatId(cat):
                    lsupp += len(self.sCats[self.getCatForId(c)])
                else:
                    lsupp += len(self.sCats.get(c, []))
            return lsupp

    def suppTerm(self, term):
        if term.isAnon():
            return set()
        return self.suppCat(term.cat)

    def lsuppTerm(self, term):
        if term.isAnon():
            return set()
        return self.lsuppCat(term.cat)

    def suppInBounds(self, min_in=-1, min_out=-1):
        if self.infofull["in"][0] != min_in:
            self.infofull["in"] = (
                min_in, (self.cards[-1][1] >= min_in or (self.nbRows() - self.cards[0][1]) >= min_in))
        if self.infofull["out"][0] != min_out:
            self.infofull["out"] = (min_out, (self.cards[-1][1] >=
                                              min_out or (self.nbRows() - self.cards[0][1]) >= min_out))
        return (self.infofull["in"][1] and self.infofull["out"][1])


associate_term_class(CatColM, CatTerm)


class NumColM(ColM):
    width = 0
    NA = NA_num
    
    p_patt = "^-?\d+(?P<dec>(\.\d+)?)$"
    # alt_patt = "^[+-]?\d+.?\d*(?:[Ee][-+]\d+)?$"
    alt_patt = "^[+-]?\d+\.?\d*(?:[Ee][-+]?\d+)?$"
    prec_patt = "^[+-]?\d+(?P<dec>\.\d+)?([Ee](?P<esgn>[-+])?(?P<epw>\d+))?$"
    @classmethod
    def parseVal(tcl, v, j=0, vals=None, miss=set(), prec=None, exclude=False, matchMiss=False):
        if (matchMiss is not False and v == matchMiss) or v == str(tcl.NA):
            miss.add(j)
            return v, prec
        if type(v) != str:
            v = "%s" % v
        tmatch = re.match(tcl.p_patt, v)
        if not tmatch:
            atmatch = re.match(tcl.alt_patt, v)
            if not atmatch:
                if matchMiss is False:
                    miss.add(j)
                return v, prec
            sfv = v
        else:
            sfv = str(float(v))

        pprec = 0
        mtch = re.match(tcl.prec_patt, sfv)
        if mtch is not None:
            if mtch.group("dec") is not None:
                pprec = len(mtch.group("dec"))-1
            if mtch.group("epw") is not None:
                if mtch.group("esgn") is not None and mtch.group("esgn") == "-":
                    pprec += int(mtch.group("epw"))
                else:
                    pprec -= int(mtch.group("epw"))
        if prec is None or pprec > prec:
            prec = pprec

        val = float(v)
        if vals is not None and (exclude is False or val != exclude):
            vals.append((val, j))
        return val, prec

    @classmethod
    def parseList(tcl, listV, indices=None, force=False, native_missing_check=True):
        prec = None
        if indices is None:
            indices = dict([(v, v) for v in range(len(listV))])
        miss = set()
        vals = []
        N = max(indices.values())+1
        if type(listV) is dict:
            ttt = set(listV.keys()).intersection(indices.keys())
        else:
            ttt = [i for i in indices.keys() if i < len(listV)]
        for i in ttt:
            j = indices[i]
            val, prec = tcl.parseVal(listV[i], j, vals, miss, prec, matchMiss=None)
        if len(vals) > 0 and (len(vals) + len(miss) == N or type(listV) is dict):
            return tcl(vals, N, miss, prec)
        elif force:
            # pdb.set_trace()
            return tcl(vals, N, miss, prec, force=True)
        else:
            return None

    @classmethod
    def fromVect(tcl, vect_data, prec=None, enabled=True):
        tmp = vect_data
        if prec is not None:
            tmp = numpy.around(tmp, prec)

        col = tcl([(v, i) for (i, v) in enumerate(tmp)], N=len(tmp), prec=prec)
        if not enabled:
            col.flipEnabled()
        return col

    def toList(self, sparse=False, fill=False, cats_int_fmt=False):
        if self.isDense():  # and not self.hasMissing():
            if sparse:
                return list(enumerate(self.getVector()))
            else:
                return self.getVector()
        else:
            tmp = dict([(i, v) for (v, i) in self.sVals])
            if self.nbRows()-1 not in tmp and fill:
                tmp[self.nbRows()-1] = tmp[-1]
            if sparse:
                if -1 in tmp:
                    tmp.pop(-1)
                return tmp.items()
            else:
                return [tmp.get(i, tmp[-1]) for i in range(self.nbRows())]

    def getInitTerms(self, minIn=0, minOut=0, productivity="medium"):
        terms = []
        low_idx, hi_idx = (len(self.sVals), 0)
        if self.lenNonMode() < minIn:
            low_idx, hi_idx = (-1, -1)
        # if self.lenMode() >= minIn and self.lenNonMode() >= minOut:
        elif self.lenMode() > 0:  # and self.lenNonMode() >= minIn: #minOut:
            # MAKE TERMS OUT OF LOWER than mode and GREATER than mode
            idx = self.modeIdx()  # index of the mode in sVals
            low_idx, hi_idx = (idx-1, idx+1)
            while low_idx > 0 and self.sVals[low_idx][0] == self.sVals[idx][0]:
                low_idx -= 1
            while hi_idx < len(self.sVals) and self.sVals[hi_idx][0] == self.sVals[idx][0]:
                hi_idx += 1
            if low_idx >= 0 and low_idx+1 >= minIn and (self.nbRows() - (low_idx+1)) >= minOut:
                terms.append((self.getAssocTermClass()(
                    self.getId(), float("-Inf"), self.sVals[low_idx][0]), low_idx+1))

            if hi_idx < len(self.sVals) and (len(self.sVals)-hi_idx) >= minIn \
                    and (self.nbRows() - (len(self.sVals)-hi_idx)) >= minOut:
                terms.append((self.getAssocTermClass()(
                    self.getId(), self.sVals[hi_idx][0], float("Inf")), len(self.sVals)-hi_idx))

        if self.lenMode() == 0 or low_idx == 0:  # every non mode is above mode
            # MAKE TERMS top half
            split_idx = (hi_idx+len(self.sVals))//2
            while split_idx < len(self.sVals) and self.sVals[split_idx][0] == self.sVals[split_idx-1][0]:
                split_idx += 1
            if split_idx < len(self.sVals) and (len(self.sVals)-split_idx) >= minIn \
                    and (self.nbRows() - (len(self.sVals)-split_idx)) >= minOut:
                terms.append((self.getAssocTermClass()(
                    self.getId(), self.sVals[split_idx][0], float("Inf")), len(self.sVals)-split_idx))

        if self.lenMode() == 0 or hi_idx == len(self.sVals):  # every non mode is below mode
            # MAKE TERMS bottom half
            split_idx = low_idx//2
            while split_idx > 0 and self.sVals[split_idx][0] == self.sVals[split_idx+1][0]:
                split_idx -= 1
            if split_idx >= 0 and split_idx+1 >= minIn and (self.nbRows() - (split_idx+1)) >= minOut and \
                    (len(terms) == 0 or (len(self.sVals)-terms[-1][1]) != split_idx+1):
                terms.append((self.getAssocTermClass()(self.getId(), float(
                    "-Inf"), self.sVals[split_idx][0]), split_idx+1))

        if not self.hasMoreInMode():
            # print("--- %s" % self)
            max_agg = [max(4*minIn, (self.nbRows() - minOut)/4),
                       min(4*minIn, (self.nbRows() - minOut)/2),
                       2*minIn]
            # print("Generated %d terms before levels\t%s" % (len(terms), max_agg))
            # tmp_t = len(terms)

            upTo = 2
            seen = set()
            if productivity == "high":
                upTo = len(max_agg)
            elif productivity == "low":
                upTo = 1
                seen = None

            for lvl in range(upTo):
                if lvl == 0 or max_agg[lvl] < max_agg[lvl-1]:
                    tt = self.collapseBuckets(max_agg[lvl])
                    for i in range(len(tt[0])):
                        count, lowb, upb = (len(tt[0][i]), tt[1][i], tt[-2][i])
                        if lvl == 0 or (lowb, upb) not in seen:
                            if seen is not None:
                                seen.add((lowb, upb))
                            if count >= minIn and (self.nbRows() - count) >= minOut:
                                if i == 0:
                                    lowb = float("-Inf")
                                elif i == len(tt[0])-1:
                                    upb = float("Inf")
                                terms.append((self.getAssocTermClass()(
                                    self.getId(), lowb, upb), count))

                                # print("Found term\t", terms[-1][0], terms[-1][1], count)
                    # print("Generated %d terms at level %d (%d agg buckets)" % (len(terms)-tmp_t, lvl, len(tt[0])))
                    # tmp_t = len(terms)
        # if len(terms) == 0:
        #     print("Nothing found %s" % self)
        return terms

    def getRoundThres(self, thres, which):
        # return thres ### NO ROUNDING, many digits...
        i = 0
        while i < len(self.sVals) - 1 and self.sVals[i][0] < thres:
            i += 1
        if which == "high":
            # print(thres, which, self.sVals[i-1][0])
            return self.sVals[i-1][0]
        else:
            # print(thres, which, self.sVals[i][0])
            return self.sVals[i][0]

    def __str__(self):
        return ColM.__str__(self) + (", %i values not in mode" % self.lenNonMode())

    def upSumsRows(self, sums_rows):
        for (v, i) in self.sVals:
            sums_rows[i] += v
        if self.mode[0] == 1:
            for i in set(range(self.N)) - self.mode[1]:
                sums_rows[i] += self.sVals[-1]
        if self.mode[0] == -1:
            for i in self.mode[1]:
                sums_rows[i] += self.sVals[-1]

    def sumCol(self):
        tt = 0
        if len(self.sVals) > 0:
            tt = sum(zip(*self.sVals)[0])
        # Add mode values, one has already been counted
        if self.mode[0] == 1:
            tt += (self.N - len(self.mode[1]) - 1)*self.sVals[-1]
        if self.mode[0] == -1:
            tt += (len(self.mode[1]) - 1)*self.sVals[-1]
        return tt

    def getValue(self, rid, pref=None):
        if self.vect is None:
            self.getVector()
        if type(self.vect) is dict:
            return self.vect.get(rid, self.vect[-1])
        else:
            self.getVector()
            return self.vect[rid]

    def valToStr(self, val):
        #     if (numpy.isnan(val) and numpy.isnan(self.NA)) or \
        #            val == self.NA:
        #         return NA_str_c
        #     return val
        # def valToTimeStr(self, val):
        if (numpy.isnan(val) and numpy.isnan(self.NA)) or \
                val == self.NA:
            return NA_str_c
        if self.getTimePrec() is not None:
            return TimeTools.format_time(val, time_prec=self.getTimePrec())
        return val

    def getNumValue(self, rid):
        return self.getValue(rid)

    def areDataEquiv(self, vA, vB):
        if "uvals" not in self.cache:
            self.makeCountTools()
        # print("DataEquiv? %s vs. %s [%s]" % (vA, vB, list(enumerate(self.cache["uvals"][:10]))))
        if vA == vB:
            return True
        icurrent, iAe, iBe, iAg, iBg = (0, None, None, None, None)
        while icurrent < len(self.cache["uvals"]) and (iAg is None or iBg is None):
            if iAg is None and self.cache["uvals"][icurrent] > vA:
                iAg = icurrent
            if iBg is None and self.cache["uvals"][icurrent] > vB:
                iBg = icurrent
            if iAe is None and self.cache["uvals"][icurrent] >= vA:
                iAe = icurrent
            if iBe is None and self.cache["uvals"][icurrent] >= vB:
                iBe = icurrent
            icurrent += 1
        # print("iAe=%s iBe=%s iAg=%s iBg=%s" % (iAe, iBe, iAg, iBg))
        ans = (False, False)
        if iAe != iAg and iBe != iBg:  # if both values are in the data
            if iAe == iBe and iAg == iBg:  # if they are equal
                # This should not happen, as we tested for equality of values at the start
                ans = (True, True)
        elif iAe == iAg and iBe == iBg:  # if neither value is in the data
            if iAe == iBe:  # and iAg == iBg: # if they are in the same bin
                ans = (True, True)
        elif iAe == iAg:  # and iBe!=iBg: # if B is in the data but not A
            ans = ((iAe == iBg), (iAe == iBe))
        elif iBe == iBg:  # and iAe!=iAg: # if A is in the data but not B
            ans = ((iBe == iAg), (iBe == iAe))
        # print("A in data %s\tB in data %s\tA ~ B %s" % (iAe!=iAg, iBe!=iBg, ans))
        return ans[0], ans[1], iAe != iAg, iBe != iBg

    def numEquiv(self, v):
        try:
            tmp = float(v)
            if tmp < self.getMin():
                tmp = self.getMin()
            elif tmp > self.getMax():
                tmp = self.getMax()
            return tmp
        except:
            pass
        return self.NA

    def getValFromNum(self, n):
        return n

    def minGap(self):
        if self.vect is None:
            self.mkVector()
        return numpy.min(numpy.diff(numpy.unique(self.vect[numpy.isfinite(self.vect)])))

    def mkVector(self):
        if self.isDense():
            self.vect = numpy.ones(self.N)*self.NA
        else:
            # mode rows -> default value is zero
            self.vect = numpy.zeros(self.N)
            self.vect[list(self.missing)] = self.NA

        if len(self.sVals) > 0:
            vals, ids = [], []
            if self.mode[1] is not None:
                tmp = [(vv, ii) for (vv, ii) in self.sVals if ii != -1]
                if len(tmp) > 0:
                    vals, ids = zip(*tmp)
            else:
                vals, ids = zip(*self.sVals)
            self.vect[list(ids)] = vals

    def getRange(self, details={}):
        return (self.getMin(details), self.getMax(details))

    def getMin(self, details={}):
        if len(self.sVals) > 0:
            return self.sVals[0][0]
        return MODE_VALUE  # DEBUG

    def getMax(self, details={}):
        if len(self.sVals) > 0:
            return self.sVals[-1][0]
        return MODE_VALUE

    def getNbValues(self):
        if "uvals" not in self.cache:
            self.makeCountTools()
        return len(self.cache["uvals"])

    def compPrec(self, details={}):
        for (v, i) in self.sVals:
            if self.prec is None or len(str(v % 1))-2 > self.prec:
                self.prec = len(str(v % 1))-2

    def getPrec(self, details={}):  # Nb of decimal digits
        if self.prec is None:
            self.compPrec()
        return max(0, self.prec)

    def getOMagn(self, details={}):  # Order of magnitude
        if self.prec is None:
            self.compPrec()
        return -self.prec

    def compTimePrec(self, details={}):
        ii = [si for si, ss in enumerate(zip(
            *[re.split("[:,-]", TimeTools.format_time(v, time_prec=1)) for (v, i) in self.sVals])) if len(set(ss)) > 1]
        if len(ii) > 1:
            self.tprec = TimeTools.range_to_time_prec(numpy.min(ii), numpy.max(ii))
        else:
            self.tprec = TimeTools.TIME_FMT  # + ", %H:%M"

    def getTimePrec(self, details={}):
        if TimeTools.isTimeVarName(self.getName()) and self.tprec is None:
            self.compTimePrec()
        return self.tprec

    def __init__(self, ncolSupp=[], N=-1, nmiss=set(), prec=None, force=False, mode=None):
        ColM.__init__(self, N, nmiss)
        self.cache = {}
        self.prec = prec
        self.tprec = None
        self.sVals = ncolSupp
        self.buk_data = {}
        if mode is None:
            self.sVals.sort()
            self.setMode(force)
            # self.setMode(force=True) # DEBUG
        else:  # specially for subsetCol, values already sorted
            self.mode = mode
        if self.N == len(self.missing):
            print("Numerical variable contains only missing values, this is suspect!..")

    def subsetCol(self, row_ids=None):
        mode_rids = None
        hasMode = False
        if row_ids is None:
            svals = [(v, i) for (v, i) in self.sVals]
            miss = set(self.missing)
            N = self.nbRows()
            if self.mode[1] is not None:
                hasMode = True
                mode_rids = set(self.mode[1])
                mode_idx = self.mode[2]
        else:
            miss = set()
            svals = []
            N = sum([len(news) for news in row_ids.values()])
            for old in self.missing.intersection(row_ids.keys()):
                miss.update(row_ids[old])

            if self.mode[1] is not None:
                mode_rids = set()
                hasMode = True
                for old in self.mode[1].intersection(row_ids.keys()):
                    mode_rids.update(row_ids[old])

            for v, old in self.sVals:
                if hasMode and old == -1:
                    mode_idx = len(svals)
                    svals.append((MODE_VALUE, -1))
                else:
                    svals.extend([(v, new) for new in row_ids.get(old, [])])

        mode = (MODE_VALUE, None, None)
        if hasMode:
            mode = (self.mode[0], mode_rids, mode_idx)

        tmp = NumColM(svals, N, miss, self.prec, mode=mode)
        tmp.extras = self.copyExtras()
        tmp.infofull = {"in": tuple(self.infofull["in"]), "out": tuple(self.infofull["out"])}
        return tmp

    def setMode(self, force=False):
        # The mode is indicated by a special entry in sVals with row id -1,
        # all rows which are not listed in either sVals or missing are assumed to take that value, equal to MODE_VALUE
        tmpV = [(v, i) for (v, i) in self.sVals if v != MODE_VALUE]
        if force or (len(self.sVals)+self.nbMissing() > 0 and len(tmpV)+self.nbMissing() != self.N
                     and (len(self.sVals)+self.nbMissing() < self.N or len(self.sVals) - len(tmpV) > MODE_RATIO*self.N)):
            # LAST CONDITIONS: either is already sparse or turning to sparse would bring gain (more than MODE_RATIO*N of values equal to mode)
            self.sVals = tmpV
            if len(self.sVals) > 0:
                vs, rids = zip(*self.sVals)
                rids = set(rids)
            else:
                rids = set()
            if -1 in rids:  # should have been removed from tmpV
                raise DataError("Error reading real values, non default mode value!")
            if len(rids) != len(self.sVals):
                raise DataError("Error reading real values, multiple values for a row!")
            # find where to insert mode in sVals
            mode_idx = 0
            while mode_idx < len(self.sVals) and self.sVals[mode_idx][0] < MODE_VALUE:
                mode_idx += 1
            self.sVals.insert(mode_idx, (MODE_VALUE, -1))

            if 2*len(rids) > self.N:  # less than half of rows in mode -> record ids of rows outside mode
                self.mode = (-1, set(range(self.N)) - rids - self.missing, mode_idx)
            else:
                self.mode = (1, rids, mode_idx)
        else:  # MODE unused
            self.mode = (0, None, None)

    def density(self):
        if self.mode[0] != 0:
            if self.mode[0] == 1:
                return len(self.mode[1])/float(self.N)
            else:
                return 1-len(self.mode[1])/float(self.N)
        return 1.0

    def isDense(self, thres=None):
        if self.mode[0] != 0:
            if thres is None:
                return False
            else:
                return self.density() > thres
        return True

    def interNonMode(self, suppX):
        if self.mode[0] == -1:
            return suppX - self.mode[1] - self.miss()
        elif self.mode[0] == 1:
            return suppX & self.mode[1]
        else:
            return suppX - self.miss()

    def interMode(self, suppX):
        if self.mode[0] == 1:
            return suppX - self.mode[1] - self.miss()
        elif self.mode[0] == -1:
            return suppX & self.mode[1]
        else:
            return set()

    def hasMoreInMode(self):
        return self.mode[0] == 1

    def lenNonMode(self):
        if self.mode[0] == -1:
            return self.nbRows() - len(self.mode[1]) - len(self.miss())
        elif self.mode[0] == 1:
            return len(self.mode[1])
        else:
            return self.nbRows() - len(self.miss())

    def lenMode(self):
        if self.mode[0] == 1:
            return self.nbRows() - len(self.mode[1]) - len(self.miss())
        elif self.mode[0] == -1:
            return len(self.mode[1])
        else:
            return 0

    def nonModeSupp(self):
        if self.mode[0] == -1:
            return set(range(self.nbRows())) - self.mode[1] - self.miss()
        elif self.mode[0] == 1:
            return self.mode[1]
        else:
            return set(range(self.nbRows()))-self.miss()

    def modeSupp(self):
        if self.mode[0] == 1:
            return set(range(self.nbRows())) - self.mode[1] - self.miss()
        elif self.mode[0] == -1:
            return self.mode[1]
        else:
            return set()

    def modeIdx(self):  # index of the mode in sVals, None if no mode
        return self.mode[2]

    def suppInBounds(self, min_in=-1, min_out=-1):
        if self.infofull["in"][0] != min_in:
            self.infofull["in"] = (min_in, self.lenNonMode() >= min_in)
        if self.infofull["out"][0] != min_out:
            self.infofull["out"] = (min_out, self.lenNonMode() >= min_out)
        return (self.infofull["in"][1] or self.infofull["out"][1])

    @classmethod
    def buk_ind_maxes(tcl, buckets):  # in case of collapsed bucket the threshold is different
        if len(buckets) > 3 and buckets[3] is not None:
            return 3
        return 1

    @classmethod
    def buk_excl_bi(tcl, buckets):  # in case of collapsed or tail buckets, index of the bucket around mode value, should be left alone 
        if len(buckets) > 4:
            return buckets[4]
        return None

    def hasBuckets(self, which=None):
        return which in self.buk_data
    
    def buckets(self, which=None, params=None):
        params = self.prepareBucketsParams(which, params)
        if not self.hasBuckets(which) or (params is not None and self.buk_data[which].get("params") != params):
            self.buk_data[which] = {"buks": self.makeBucketsMore(which, params), "params": params}
        return self.buk_data[which]["buks"]

    def prepareBucketsParams(self, which, params):
        if params is None:
            return None
        if params.get("nbb") is not None:
            params["max_agg"] = self.nbRows()/float(params.pop("nbb"))
        if which == "collapsed":
            params["checknext"] = True
            params["base_buckets"] = "tails"
        return params

    def makeBuckets(self):
        if self.sVals[0][1] != -1:
            bucketsSupp = [set([self.sVals[0][1]])]
        else:
            bucketsSupp = [set()]
        bucketsVal = [self.sVals[0][0]]
        bukMode = None
        for (val, row) in self.sVals:
            if row == -1:
                if val != bucketsVal[-1]:  # should be ...
                    bucketsVal.append(val)
                    bucketsSupp.append(set())
                bukMode = len(bucketsVal)-1
            else:
                if val == bucketsVal[-1]:
                    bucketsSupp[-1].add(row)
                else:
                    bucketsVal.append(val)
                    bucketsSupp.append(set([row]))
        return (bucketsSupp, bucketsVal, bukMode)
    
    def bucketWithEdges(self, bin_edges, separate_mode=True, standalone_mode=False):
        mode_buk, standalone_buk = (None, None)
        colB_supp = []
        colB_min = []
        colB_max = []
        current_bi = 0

        for (val, row) in self.sVals:
            if row == -1 and separate_mode: # put mode values in a bucket of their own
                if len(colB_max) == 0 or colB_max[-1] != val:  # should be ...
                    colB_min.append(val)
                    colB_max.append(val)
                    colB_supp.append(set())
                mode_buk = len(colB_min)-1
                standalone_buk = len(colB_min)-1 if standalone_mode else None
            else:
                if len(colB_min) == 0 or (current_bi < len(bin_edges) and val >= bin_edges[current_bi]) or mode_buk == len(colB_min)-1:
                    ## NEW BUCKET NEEDED
                    while current_bi < len(bin_edges) and val >= bin_edges[current_bi]:
                        current_bi += 1
                    colB_min.append(val)
                    colB_max.append(val)
                    if row == -1: # consider mode values like other values
                        colB_supp.append(set(self.modeSupp()))
                    else:
                        colB_supp.append(set([row]))
                else:
                    colB_max[-1] = val                    
                    if row == -1: # consider mode values like other values
                        colB_supp[-1].update(self.modeSupp())
                    else:
                        colB_supp[-1].add(row)
        return (colB_supp, colB_min, mode_buk, colB_max, standalone_buk)
    
    def similarHeightBuckets(self, nb_buckets=25, separate_mode=False):
        nb_agg = self.nbRows()/nb_buckets
        if self.lenMode() > nb_agg and nb_buckets > 1 and not separate_mode:
            return self.collapseBuckets((self.nbRows()-self.lenMode())/(nb_buckets-1), checknext=True)
        else:
            return self.collapseBuckets(nb_agg, checknext=True)
            
    def equalWidthBuckets(self, nb_buckets=25, separate_mode=False):
        return self.bucketWithEdges(numpy.linspace(self.sVals[0][0], self.sVals[-1][0], num=nb_buckets, endpoint=False), separate_mode=separate_mode)
        
    def bellmanBuckets(self, nb_buckets=25, separate_mode=False, bellman_criterion=None):
        if separate_mode:
            assign, bounds, cost, opt_k = bellman.bin_segments([x[0] for x in self.sVals], nb_buckets, bellman_criterion)
        else:
            assign, bounds, cost, opt_k = bellman.bin_segments(self.getVector(), nb_buckets, bellman_criterion)
        return self.bucketWithEdges(bounds, separate_mode=separate_mode)
    
    def tailsBuckets(self, lower_tail_agg=0, upper_tail_agg=1, base_buckets=None):
        if self.hasBuckets(base_buckets):
            tmp = self.buckets(base_buckets)
        else:
            tmp = self.buckets()
        if lower_tail_agg == 0 or upper_tail_agg == 0:
            return tmp
        nbs = [0, 0]
        counts = [len(x) for x in tmp[0]]
        if tmp[2] is not None:
            counts[tmp[2]] += self.lenMode()
        for i, v in enumerate([lower_tail_agg, upper_tail_agg]):
            # v==0: don't cut, v==-1: cut off entirely
            if v == 0 or v == -1:
                nbs[i] = v
            elif v <= 0:
                if v < -1:  # nb of distinct values
                    nbs[i] = -v
                else:  # fraction of distinct values
                    nbs[i] = int(-v*len(tmp[0]))
            else:
                if v >= 1:  # number of distinct rows
                    nbv = v
                else:  # fraction of distinct rows, quantile
                    nbv = v*(self.nbRows()-self.nbMissing())
                # iterate in reverse order for upper tail
                ccounts = numpy.cumsum(counts[::(-i*2+1)])
                while nbs[i] < len(ccounts) and ccounts[nbs[i]] < nbv:
                    nbs[i] += 1

        bottom_mid = max(0, nbs[0])
        top_mid = min(len(tmp[0]), len(tmp[0])-nbs[1]-1)
        if bottom_mid >= top_mid:
            return tmp

        mode_buk = tmp[2]
        if tmp[2] is not None:
            if tmp[2] > top_mid:  # mode is above the merge middle
                mode_buk = tmp[2]-(top_mid-bottom_mid)+1
            elif tmp[2] >= bottom_mid:  # mode is in the merge middle
                mode_buk = bottom_mid

        colB_supp = []
        colB_min = []
        colB_max = []
        # the lower tail
        colB_supp.extend(tmp[0][:bottom_mid])
        colB_min.extend(tmp[1][:bottom_mid])
        colB_max.extend(tmp[1][:bottom_mid])
        # the aggregated part
        standalone_buk = len(colB_supp)
        colB_supp.append(set([]).union(*tmp[0][bottom_mid:top_mid]))
        colB_min.append(tmp[1][bottom_mid])
        colB_max.append(tmp[1][top_mid-1])
        # the lower tail
        colB_supp.extend(tmp[0][top_mid:])
        colB_min.extend(tmp[1][top_mid:])
        colB_max.extend(tmp[1][top_mid:])
        # print(self, colB_min[standalone_buk], colB_max[standalone_buk], len(colB_supp[standalone_buk]))
        # if mode_buk is not None:
        #     print("Mode", colB_min[mode_buk], colB_max[mode_buk], colB_supp[mode_buk])
        return (colB_supp, colB_min, mode_buk, colB_max, standalone_buk)

    def collapseBuckets(self, max_agg, base_buckets=None, checknext=False):
        # collapsing from low to up, could do reverse...
        if self.hasBuckets(base_buckets):
            tmp = self.buckets(base_buckets)
        else:
            tmp = self.buckets()
        tmp_supp = set([])
        bucket_min = tmp[1][0]
        colB_supp = []
        colB_min = []
        colB_max = []
        bukMode = None
        max_id = self.buk_ind_maxes(tmp)
        standalone_buk = tmp[2]
        if len(tmp) > 4 and tmp[4] is not None:
            standalone_buk = tmp[4]
        new_standalone_buk = None
        # colB_max= [None]
        for i in range(len(tmp[1])):
            # would exceed aggregated size
            if len(tmp_supp) > max_agg or (checknext and i > 0 and len(tmp[0][i]) > max_agg) or \
                    (standalone_buk is not None and ((standalone_buk == i and i > 0)  or standalone_buk == i-1)):  # if there is mode bucket, leave it alone (but there is nothing to store if it's the first bucket, i.e. standalone_buk == i == 0)
                colB_supp.append(tmp_supp)
                colB_min.append(bucket_min)
                colB_max.append(tmp[max_id][i-1])
                bucket_min = tmp[1][i]
                tmp_supp = set([])
            tmp_supp.update(tmp[0][i])
            if tmp[2] is not None and tmp[2] == i:
                bukMode = len(colB_supp)
            if standalone_buk is not None and standalone_buk == i:
                new_standalone_buk = len(colB_supp)

        colB_supp.append(tmp_supp)
        colB_min.append(bucket_min)
        colB_max.append(tmp[max_id][-1])
        # colB_max[0] = colB_max[1]
        return (colB_supp, colB_min, bukMode, colB_max, new_standalone_buk)

    BUCKET_METHODS = {"collapsed": collapseBuckets,
                      "tails": tailsBuckets,
                      "bellman": bellmanBuckets,
                      "equal-width": equalWidthBuckets,
                      "similar-height": similarHeightBuckets}
    
    def makeBucketsMore(self, which, params):
        if params is None:
            params = {}
        if which in self.BUCKET_METHODS:
            return self.BUCKET_METHODS[which](self, **params)
            # bucks = self.BUCKET_METHODS[which](self, **params)
            # if bucks[3] is not None:
            #     print("Col", self.getId(), " buckets:", list(zip(bucks[1], bucks[3])))
            # else:
            #     print("Col", self.getId(), " buckets:", bucks[1])
            # return bucks
        return self.makeBuckets()

    def dispBuckets(self, buckets, sep="\n\t"):
        bstr = ["## %d buckets" % len(buckets[1])]
        max_id = self.buk_ind_maxes(buckets)
        for bi, min_v in enumerate(buckets[1]):
            bstr.append("[%s, %s] %s%s" % (min_v, buckets[max_id][bi],
                                           self.lenMode() if buckets[2] == bi else len(buckets[0][bi]),
                                           "M"*(buckets[2] == bi) + "A"*(len(buckets) > 3 and buckets[-1] == bi)))
        return sep.join(bstr)
            
    
    def suppTerm(self, term):
        if term.isAnon():
            return set()
        suppIt = set()
        for (val, row) in self.sVals:
            if val > term.upb:
                return suppIt
            elif val >= term.lowb:
                if row == -1:
                    suppIt.update(self.modeSupp())
                else:
                    suppIt.add(row)
        return suppIt

    def lsuppTerm(self, term):
        if term.isAnon():
            return 0
        lsuppIt = 0
        for (val, row) in self.sVals:
            if val > term.upb:
                return lsuppIt
            elif val >= term.lowb:
                if row == -1:
                    lsuppIt += self.lenMode()
                else:
                    suppIt += 1
        return suppIt

    def makeCountTools(self):
        vs, rids = zip(*self.sVals)
        map_rids = -numpy.ones(self.N+1, dtype=int)  # extra entry for mode, just in case
        map_rids[list(rids)] = numpy.arange(len(rids))
        self.cache["boundaries"] = numpy.concatenate(
            [[0], numpy.where(numpy.diff(vs) != 0)[0]+1, [len(vs)]])
        # vxs = numpy.array(vs)
        # print([vxs[self.cache["boundaries"][bi]:self.cache["boundaries"][bi+1]]
        #        for bi in range(len(self.cache["boundaries"])-1)])
        self.cache["uvals"] = [vs[bi] for bi in self.cache["boundaries"][:-1]]
        self.cache["map_rids"] = map_rids

    def makeCountsMatrix(self, supports):
        if "uvals" not in self.cache:
            self.makeCountTools()
        map_rids = self.cache["map_rids"]
        C = numpy.zeros((len(self.sVals)+1, supports.nbParts()), dtype=int)  # extra entry for mode
        for spi in range(C.shape[1]):
            spset = supports.part(spi)
            C[map_rids[list(spset)], spi] = 1
            C[map_rids[-1], spi] = len(self.modeSupp() & spset)
        # if spi+2 == C.shape[1]:  # last part, Eoo, is implicit
        #     C[:, -1] = 1*(numpy.sum(C[:, :-1], axis=1) == 0)
        #     C[map_rids[-1], -1] = self.lenMode() - numpy.sum(C[map_rids[-1], :-1])
        counts = numpy.vstack([numpy.sum(C[self.cache["boundaries"][bi]:self.cache["boundaries"][bi+1], :], axis=0)
                               for bi in range(len(self.cache["boundaries"])-1)])
        return counts

    def findCutpoints(self, ssetts, side, supports, op, counts=None):
        if counts is None:
            counts = self.makeCountsMatrix(supports)
        cutpoints = []  # list of [low boundary id, up boundary id, count_var_num, count_var_den]
        ids_var_num = [ssetts.partId(part_id, side) for (io, part_id) in ssetts.IDS_varnum[op]]
        ids_var_den = [ssetts.partId(part_id, side) for (io, part_id) in ssetts.IDS_varden[op]]

        counts_var_num = numpy.sum(counts[:, ids_var_num], axis=1)
        counts_var_den = numpy.sum(counts[:, ids_var_den], axis=1)

        for i in range(counts.shape[0]):
            if counts_var_num[i]+counts_var_den[i] > 0:
                if len(cutpoints) > 0 and \
                    ((cutpoints[-1][2] + counts_var_num[i]) == 0 or
                     (cutpoints[-1][3] + counts_var_den[i]) == 0):  # block remains "pure"
                    cutpoints[-1][1] = i
                    cutpoints[-1][2] += counts_var_num[i]
                    cutpoints[-1][3] += counts_var_den[i]
                else:
                    cutpoints.append([i, i, counts_var_num[i], counts_var_den[i]])  # initialize
        return cutpoints, counts

    def getLiteralNum(self, neg, val_range):
        if val_range[-1] in ["uvals", "cutpoints"]:
            return self.getLiteralCutpoint(neg, val_range[:-1])
        elif self.hasBuckets(val_range[-1]):
            buckets = self.buckets(val_range[-1])
            bUp = self.buk_ind_maxes(buckets)
            return self.getLiteralBuk(neg, buckets[1], val_range[:-1], buckets[bUp])
        else:
            return Literal(neg, self.getAssocTermClass()(self.getId(), val_range[0], val_range[-1]))
    
    def getLiteralBuk(self, neg, buk_op, bound_ids, buk_op_top=None):
        same_top = False
        if buk_op_top is None:
            buk_op_top = buk_op
            same_top = True
        if (bound_ids[0] is None and bound_ids[1] is None) or \
           (bound_ids[0] == 0 and bound_ids[1] == len(buk_op)-1) or \
                (same_top and bound_ids[0] is not None and bound_ids[1] is not None and bound_ids[0] > bound_ids[1]):
            return None  # unbounded
        elif bound_ids[0] is None or bound_ids[0] == 0:
            if neg:
                lowb = buk_op[bound_ids[1]+1]
                upb = float('Inf')
                n = False
            else:
                lowb = float('-Inf')
                upb = buk_op_top[bound_ids[1]]
                n = False
        elif bound_ids[1] is None or bound_ids[1] == len(buk_op)-1:
            if neg:
                lowb = float('-Inf')
                upb = buk_op_top[bound_ids[0]-1]
                n = False
            else:
                lowb = buk_op[bound_ids[0]]
                upb = float('Inf')
                n = False
        else:
            lowb = buk_op[bound_ids[0]]
            upb = buk_op_top[bound_ids[1]]
            n = neg
        return Literal(n, self.getAssocTermClass()(self.getId(), lowb, upb))

    def getLiteralCutpoint(self, neg, bound_ids):
        if "uvals" not in self.cache:
            raise DataError("Cutpoints boundaries of unidentified origin!")
        return self.getLiteralBuk(neg, self.cache["uvals"], bound_ids)


associate_term_class(NumColM, NumTerm)

# if __name__ == '__main__':
#     vals = ["+4.9091", "+4.909100", "+4.9e4", "4.9e-4", "-4e-4", "+4.9e34", "4.9e-34", "-4e-34"]
#     for v in vals:
#         print(v, NumColM.parseVal(v))
