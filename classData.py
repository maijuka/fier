import os.path
import time
import re
from itertools import chain
import numpy

try:
    from classContent import Container, ContentCollection
    from classCol import DataError, ColM, ConstantColM, BoolColM, CatColM, NumColM
    from classDataExtension import DataExtension
    from classQuery import Op, Term, Literal, Query, NA_str_c
    from classSParts import SSetts
    from toolICList import ICList
    import csv_reader
except ModuleNotFoundError:
    from .classContent import Container, ContentCollection
    from .classCol import DataError, ColM, ConstantColM, BoolColM, CatColM, NumColM
    from .classDataExtension import DataExtension
    from .classQuery import Op, Term, Literal, Query, NA_str_c
    from .classSParts import SSetts
    from .toolICList import ICList
    from . import csv_reader

import pdb

FORCE_WRITE_DENSE = False


def all_subclasses(cls):
    return cls.__subclasses__() + [g for s in cls.__subclasses__()
                                   for g in all_subclasses(s)]


class ExtensionPool(object):

    extensions_map = {}
    extras_map = {}
    filenames_map = {}
    for tcl in all_subclasses(DataExtension):
        key = tcl.getKey()
        extensions_map[key] = {"class": tcl}
        extensions_map[key].update(tcl.getExtensionDetails())
        for ek in tcl.getExtrasKeys():
            extras_map[ek] = key
        for fk in tcl.getFilenamesKeys():
            filenames_map[fk] = key

    def __init__(self, data):
        self.data = data
        self.extensions = {}

    def getExtensionClass(self, key):
        if key in self.extensions_map:
            return self.extensions_map[key]["class"]

    def initExtension(self, key, filenames=None, params=None, details={}):
        tcl = self.getExtensionClass(key)
        if tcl is not None:
            ext = tcl(self.data, filenames, params, details)
            self.extensions[key] = ext
            return ext

    def getExtension(self, key):
        return self.extensions.get(key)

    def delExtension(self, key):
        if key in self.extensions:
            del self.extensions[key]

    def hasActiveExtension(self, key):
        return key in self.extensions

    def getActiveExtensionKeys(self):
        return self.extensions.keys()

    def getAvailableExtensionKeys(self):
        return self.extensions_map.keys()

    def computeExtras(self, item, extras=None, details={}):
        if extras is None:
            extras = [k for (k, kk) in self.extras_map.items() if kk in self.extensions]
        extra_values = dict([(e, None) for e in extras])
        extensions_to_extras = {}
        for extra_key in extras:
            if extra_key in self.extras_map:
                extension_key = self.extras_map[extra_key]
                if extension_key in self.extensions:
                    if extension_key not in extensions_to_extras:
                        extensions_to_extras[extension_key] = []
                    extensions_to_extras[extension_key].append(extra_key)
        for (extension_key, extra_keys) in extensions_to_extras.items():
            extra_values.update(self.extensions[extension_key].computeExtras(item, extra_keys, details))
        return extra_values

    def loadExtensions(self, ext_keys=[], filenames=None, params=None, details={}):
        if len(ext_keys) == 0:  # prevent loading
            return

        if ext_keys is None or "(auto)" in ext_keys:
            if filenames is None:
                filenames = self.getExtensionsFilesDict()
                eks = self.extensions_map.keys()
            else:
                eks = set()
                for fk, filename in filenames.items():
                    if fk in self.filenames_map:
                        eks.add(self.filenames_map[fk])
        else:
            eks = ext_keys
            if filenames is None:
                filenames = self.getExtensionsFilesDict()

        params_collect = {}
        for ek in eks:
            ext = self.initExtension(ek, filenames, params, details)
            params_collect.update(ext.getParams())
        return params_collect

    def saveExtensions(self, filenames=None, details={}):
        if filenames is None:
            filenames = self.getExtensionsActiveFilesDict()
        eks = set()
        for fk, filename in filenames.items():
            if self.filenames_map.get(fk) in self.extensions_map:
                eks.add(self.filenames_map[fk])
        for ek in eks:
            ext = self.extensions[ek]
            ext.doWithFiles("save", filenames, details)

    def getExtensionsFilesDict(self, exts=None):
        fdict = {}
        if exts is None:
            exts = self.extensions_map.keys()
        for ek in exts:
            if ek in self.extensions_map:
                fdict.update(self.extensions_map[ek]["class"].getFilesDict())
        return fdict

    def getExtensionsActiveFilesDict(self):
        fdict = {}
        for ek, ext in self.extensions.items():
            fdict.update(ext.getActiveFilesDict())
        return fdict


class RowE(object):

    def __init__(self, rid, data):
        self.rid = rid
        self.data = data

    def getValue(self, side, col=None):
        if col is None:
            if side.get("aim", None) == "sort":
                t = self.data.getNumValue(side["side"], side["col"], self.rid)
                return {BoolColM.NA: None, CatColM.NA: None, NumColM.NA: None}.get(t, t)
            elif side.get("aim", None) == "row":
                return self.data.getNumValue(side["side"], side["col"], self.rid)
            else:
                return self.data.getValue(side["side"], side["col"], self.rid)
        else:
            return self.data.getValue(side, col, self.rid)

    def getEnabled(self, details={}):
        if self.rid not in self.data.selectedRows():
            return 1
        else:
            return 0

    def isEnabled(self, details={}):
        return self.getEnabled(details) > 0

    def flipEnabled(self):
        if self.rid in self.data.selectedRows():
            self.data.removeSelectedRow(self.rid)
        else:
            self.data.addSelectedRow(self.rid)

    def setEnabled(self):
        self.data.removeSelectedRow(self.rid)

    def setDisabled(self):
        self.data.addSelectedRow(self.rid)

    def getId(self, details={}):
        return self.rid

    def getUid(self, details={}):
        return self.getId()

    def getRName(self, details={}):
        return self.data.getRName(self.rid)


TYPES_SMAP = {}
for c in ColM.__subclasses__():
    TYPES_SMAP[c.type_id] = c
    TYPES_SMAP[str(c.type_id)] = c
    TYPES_SMAP[c.type_letter] = c


class ContainerVars(Container):
    def getShortStr(self):
        return "Vars #%s (%d)" % (self.getUid(), len(self))


class Data(ContentCollection):

    enabled_codes = {(0, 0): "F", (1, 1): "T", 0: "F", 1: "T", (0, 1): "L", (1, 0): "R"}
    enabled_codes_rev_simple = {"F": 0, "T": 1}
    enabled_codes_rev_double = {"F": (0, 0), "T": (1, 1), "L": (0, 1), "R": (1, 0)}
    var_types = [None, BoolColM, CatColM, NumColM]
    all_types_map = dict([(None, None)]+[(v.type_letter, v) for v in var_types[1:]])
    real_types = var_types[1:]
    real_types_name_to_id = dict([(c.type_name, c.type_id) for c in real_types])
    NA_str_def = NA_str_c
    NA_str = NA_str_c
    NA_str = NA_str_c
    compat_diff = True

    @classmethod
    def setCompat(tcl, var_compat=True):
        if type(var_compat) is bool:
            tcl.compat_diff = var_compat
        elif type(var_compat) is str:
            tcl.compat_diff = var_compat.startswith("d")

    @classmethod
    def isCompatD(tcl):
        return tcl.compat_diff

    @classmethod
    def getNamesTids(tcl):
        return [(c.type_name, c.type_id) for c in tcl.real_types]

    @classmethod
    def getTidForName(tcl, name):
        if type(name) is list:
            return [tcl.real_types_name_to_id.get(n) for n in name]
        return tcl.real_types_name_to_id.get(name)

    @classmethod
    def isTypeId(tcl, tid, name, default_accept=False):
        if tid is None:
            return default_accept
        if type(name) is list:
            return tid in tcl.getTidForName(name)
        return tcl.getTidForName(name) == tid

    @classmethod
    def getColClassForName(tcl, name):
        if type(name) is list:
            return [TYPES_SMAP.get(tcl.real_types_name_to_id.get(n)) for n in name]
        return TYPES_SMAP.get(tcl.real_types_name_to_id.get(name))

    container_class = ContainerVars

    def __init__(self, cols=[[], []], N=0, coords=None, rnames=None, single_dataset=False):
        ContentCollection.__init__(self)
        self.single_dataset = single_dataset
        self.fold = None
        self.as_array = [None, None, None]
        self.selected_rows = set()
        self.condition_dt = None
        self.extensions = ExtensionPool(self)
        if type(N) == int:
            data_cols = cols
            self.N = N
            self.rnames = rnames
        elif type(N) == str:
            try:
                data_cols, self.N, coords, self.rnames, self.selected_rows, self.condition_dt, self.single_dataset, Data.NA_str = readDNCFromCSVFiles(cols)
            except DataError:
                data_cols, self.N, coords, self.rnames = [[], []], 0, None, None
                raise
        else:
            print("Input non recognized!")
            data_cols, self.N, coords, self.rnames = [[], []], 0, None, None
            raise
        self.setCoords(coords)

        dt_cs = list(enumerate(data_cols))
        if self.condition_dt is not None and len(self.condition_dt["cols"]) > 0:
            dt_cs.append((-1, self.condition_dt["cols"]))
        self.initLists(dt_cs)
        # self.initLists(data_cols)
        self.ssetts = SSetts(self.hasMissing(), data_symmetric=self.isSymmetric())

    def sameAs(self, other, side=None, cid=None):
        if not isinstance(other, Data):
            return False
        if cid is None:
            if self.getSides(side) == other.getSides(side) and \
                    all([len(self.colsSide(sside)) == len(other.colsSide(sside)) for sside in self.getSides(side)]):
                for sside in self.getSides(side):
                    for colA, colB in zip(self.colsSide(sside), other.colsSide(sside)):
                        if colA.cmpType(colB) != 0 or colA.cmpVals(colB) != 0:
                            return False
                return True
        else:
            colA = self.col(side, cid)
            colB = other.col(side, cid)
            if colA is None and colB is None:
                return True
            elif colA is not None and colB is not None:
                return colA.cmpType(colB) == 0 and colA.cmpVals(colB) == 0
        return False

    def recomputeCols(self, side=None, cid=None):
        if cid is None:
            for sside in self.getSides(side):
                for col in self.colsSide(sside):
                    try:
                        col.recompute(self)
                    except KeyError:
                        raise Exception("Error when recomputing column! (%s,%s)" % (side, cid))
                        # pdb.set_trace()
                        # col.recompute(self)
        else:
            col = self.col(side, cid)
            if col is not None:
                col.recompute(self)

    def computeExtras(self, item, extras=None, details={}):
        return self.extensions.computeExtras(item, extras, details)

    def initExtension(self, key, filenames=None, params=None, details={}):
        self.extensions.initExtension(key, filenames, params, details)

    def getExtension(self, key):
        return self.extensions.getExtension(key)

    def delExtension(self, key):
        self.extensions.delExtension()

    def hasActiveExtension(self, key):
        return self.extensions.hasActiveExtension(key)

    def getActiveExtensionKeys(self):
        return self.extensions.getActiveExtensionKeys()

    def getAvailableExtensionKeys(self):
        return self.extensions.getAvailableExtensionKeys()

    def saveExtensions(self, filenames=None, details={}):
        self.extensions.saveExtensions(filenames, details)

    def loadExtensions(self, ext_keys=[], filenames=None, params=None, details={}):
        return self.extensions.loadExtensions(ext_keys, filenames, params, details)

    def getExtensionsFilesDict(self, exts=None):
        return self.extensions.getExtensionsFilesDict(exts)

    def getExtensionsActiveFilesDict(self):
        return self.extensions.getExtensionsActiveFilesDict()

    def isSingleD(self):
        return self.single_dataset

    def isSymmetric(self):
        return self.isSingleD() and all([self.col(0, col.getId()).getEnabled() == self.col(1, col.getId()).getEnabled() for col in self.colsSide(0)])

    def setColsSide(self, side, cols):
        self.clearList(side)
        for c in cols:
            self.appendCol(c, side)

    def appendCol(self, col, side=0):
        self.addItem(col, side)

    def addItem(self, item, trg_lid=None, trg_pos=-1):
        ContentCollection.addItem(self, item, trg_lid, trg_pos)
        item.recompute(self)

    def colsSide(self, side):
        # if side == -1 and self.condition_dt is not None:
        #     return self.condition_dt.get("cols", [])
        if side in self.containers:
            return [self.items[iid] for iid in self.containers[side].naturalOrder()]  # keep the original order when enumerating columns
        return []

    def col(self, side, literal):
        cid = None
        ccs = self.items
        if False:  # side == -1 and self.condition_dt is not None:
            ccs = self.condition_dt["cols"]
            if (type(literal) is int) or (isinstance(literal, numpy.integer)):
                cid = literal
            elif (isinstance(literal, Term) or isinstance(literal, Literal)):
                cid = literal.colId()
            if cid is not None and cid > len(ccs):
                err_cid = cid
                cid = None
                raise DataError("This columns does not exist! CID=%s" % str(err_cid))

        elif side in self.containers:
            if (type(literal) is int) or (isinstance(literal, numpy.integer)):
                cid = (side, literal)
            elif (isinstance(literal, Term) or isinstance(literal, Literal)):
                cid = (side, literal.colId())

        if cid is not None:

            if cid[1] == -1:  # constant column
                return ConstantColM(self.N)

            elif cid[1] > self.getLen(side) or cid not in ccs:
                err_cid = cid
                cid = None
                raise DataError("This columns does not exist! CID=%s" % str(err_cid))

            else:
                if (isinstance(literal, Term) or isinstance(literal, Literal)) and not literal.isAnon() and literal.typeId() != ccs[cid].typeId():
                    err_cid = cid
                    cid = None
                    raise DataError("The type of literal does not match the type of the corresponding variable (on side %s col %s type %s ~ lit %s type %s [%s])!" % (side, err_cid, ccs[err_cid].typeId(), literal.colId(), literal.typeId(), literal))
                else:
                    return ccs[cid]

    def getElement(self, iid):
        return self.col(iid[0], iid[1])

    def getSides(self, side=None):
        if side is not None:
            return [side]
        return self.getList()

    def keys(self):
        kys = []
        for side in self.getSides():
            for col in self.colsSide(side):
                kys.append((side, col.getId()))
        return kys

    def getValue(self, side, col, rid):
        return self.col(side, col).getValue(rid)

    def getNumValue(self, side, col, rid):
        return self.col(side, col).getNumValue(rid)

    def replaceSideFromMatrix(self, mat_data, prec=None, vnames=None, side=0, enabled=True, vtypes=None):
        if self.nbRows() != mat_data.shape[1]:
            raise DataError("Side replacement with different number of rows is prohibited! (%d vs. %d)" % (self.nbRows(), mat_data.shape[1]))
        if vnames is None and self.nbCols(side) == mat_data.shape[0]:
            vnames = ["RND_%s" % n for n in self.getNames(side)]
        cols = prepareSideFromMatrix(mat_data, prec=prec, vnames=vnames, side=side, enabled=enabled, vtypes=vtypes)
        back = self.colsSide(side)
        self.setColsSide(side, cols[side])
        return back

    # TODO: STOP
    def hasMissing(self, side=None):
        for sside in self.getSides(side):
            for c in self.colsSide(sside):
                if c.hasMissing():
                    return True
        return False

    def getAllTypes(self, side=None):
        typs = []
        for sside in self.getSides(side):
            typs.extend([col.typeId() for col in self.colsSide(sside)])
        return set(typs)

    def getCommonType(self, side):
        s = set([col.type_letter for col in self.colsSide(side)])
        if len(s) == 1:
            return s.pop()
        return None

    def arePairTypesIn(self, cidA, cidB, sideA=0, sideB=1, tset=[]):
        if tset is not None and len(tset) > 0:
            if self.col(sideA, cidA).type_id > self.col(sideB, cidB).type_id:
                k = "%s%s" % (self.col(sideB, cidB).type_letter, self.col(sideA, cidA).type_letter)
            else:
                k = "%s%s" % (self.col(sideA, cidA).type_letter, self.col(sideB, cidB).type_letter)
            return k in tset
        return False

    def hasGroups(self, side=None):
        if side is None:
            return any([self.hasGroups(sside) for sside in self.getSides(side)])
        return any([col.hasGroup() for col in self.colsSide(side)])

    def diffGid(self, gidA, gidB):
        return gidA == -1 or gidB == -1 or gidA != gidB

    def sameGid(self, gidA, gidB):
        return gidA == -1 or gidB == -1 or gidA == gidB

    def areGroupCompat(self, cidA, cidB, sideA=0, sideB=1):
        gidA = self.col(sideA, cidA).getGroupId()
        gidB = self.col(sideB, cidB).getGroupId()
        if self.isCompatD():
            return self.diffGid(gidA, gidB)
        return self.sameGid(gidA, gidB)

    def upColsCompat(self, cAv, side, cids, other_side=False):
        if cAv is None or len(cAv) == 0:
            return
        seen_groups = set([-1])
        if other_side:
            ss = 1-side
        else:
            ss = side
            cAv.difference_update(cids)

        for cid in cids:
            gid = self.col(side, cid).getGroupId()
            if gid not in seen_groups:
                seen_groups.add(gid)
                if self.isCompatD():
                    cAv.difference_update([c for c in cAv if self.col(ss, c).getGroupId() == gid])
                else:
                    cAv.intersection_update([c for c in cAv if self.col(ss, c).getGroupId() == gid or self.col(ss, c).getGroupId() == -1])
        return cAv

    def getSSetts(self):
        return self.ssetts

    def resetSSetts(self, parts_type=None, method_pval=None):
        return self.ssetts.reset(self.hasMissing(), parts_type, method_pval, self.isSymmetric())

    def getRName(self, rid):
        if self.rnames is not None and rid < len(self.rnames):
            return self.rnames[rid]
        return "#%d" % rid

    def getRNames(self):
        if self.rnames is not None:
            return list(self.rnames)
        return ["#%d" % rid for rid in self.rows()]

    def getStats(self, group=None):
        if group is None:
            # Group all columns from both side together
            group = []
            for side in [0, 1]:
                group.extend([(side, i) for i in range(data.nbCols(side))])
        elif type(group) == int and group in [0, 1]:
            # Group all columns from that side together
            side = group
            group = [(side, i) for i in range(data.nbCols(side))]

        sums_rows = [None for t in Data.var_types]
        sums_cols = []
        details = []
        for side, col in group:
            tid = self.col(side, col).typeId()
            if sums_rows[tid] is None:
                sums_rows[tid] = self.col(side, col).initSums(self.N)
            self.col(side, col).upSumsRows(sums_rows[tid])
            sums_cols.append(self.col(side, col).sumCol())
            details.append((side, col, tid))
        return sums_rows, sums_cols, details

    def getLitsToBinaryMatrix(self, lits_dets):
        mat = numpy.zeros((self.nbRows(), len(lits_dets)), dtype=bool)
        for li, (side, literal) in enumerate(lits_dets):
            supp = self.supp(side, literal)
            mat[list(supp), li] = 1
        return mat

    @classmethod
    def getMatrixCols(tcl, cols, nb_rows=None, bincats=False, nans=None):
        if len(cols) > 0:
            if nb_rows is None:
                nb_rows = cols[0].nbRows()
            return numpy.hstack([col.getVector(bincats, nans).reshape((nb_rows, -1)) for col in cols]).T
        return numpy.array([])

    def getMatLitK(self, side, lit, bincats=False):
        term, cid,  off = lit, lit, 0
        if isinstance(lit, Literal):
            term = lit.getTerm()
        if isinstance(term, Term):
            cid = term.colId()
            col = self.col(side, cid)
            if bincats and col is not None and col.typeId() == 2:
                for cc in term.getCat():
                    off = col.numEquiv(cc)
        return (side, cid, off)

    def getMatrix(self, side_cols=None, store=True, types=None, only_able=False, bincats=False, nans=None):
        compare_cols = None
        if store and self.as_array[0] == (side_cols, types, only_able, bincats):
            compare_cols = sorted(self.as_array[1][-1].keys())

        if store:
            self.as_array[0] = (side_cols, types, only_able, bincats)

        if types is None:
            types = [BoolColM.type_id, CatColM.type_id, NumColM.type_id]

        if side_cols is None:
            side_cols = [(side, None) for side in self.getSides()]

        mcols = {}
        details = []
        off = 0
        mat = None
        tcols = []
        for side, col in side_cols:
            if col is None:
                tcols = [c.getId() for c in self.colsSide(side)]
            else:
                tcols = [col]
            tcols = [c for c in tcols if self.col(side, c).typeId() in types and (not only_able or self.col(side, c).isEnabled())]
            if len(tcols) > 0:
                for col in tcols:
                    bids = [0]
                    if bincats and self.col(side, col).typeId() == 2:
                        bids = range(self.col(side, col).nbCats())
                    mcols[(side, col)] = len(details)
                    for bid in bids:
                        mcols[(side, col, bid)] = off
                        off += 1
                    details.append({"side": side, "col": col, "type": self.col(side, col).typeId(), "name": self.col(side, col).getName(), "enabled": self.col(side, col).getEnabled(), "bincats": bids})

        if compare_cols is not None and compare_cols == sorted(mcols.keys()):
            return self.as_array[1]

        cols = [self.col(d["side"], d["col"]) for d in details]
        mat = self.getMatrixCols(cols, nb_rows=self.nbRows(), bincats=bincats, nans=nans)
        if store:
            self.as_array[1] = (mat, details, mcols)
        return mat, details, mcols

    # FOLDS
    #######################

    def getFold(self, nbsubs=10, coo_dim=None, grain=10., force=False):
        if coo_dim is not None and \
            not ((self.isGeospatial() and coo_dim < 0 and abs(coo_dim)-1 < len(self.getCoords())) or
                 (coo_dim > 0 and coo_dim < self.nbCols(0)+self.nbCols(1)+1)):
            coo_dim = None

        if (self.fold is None) or (self.fold["source"] != "auto") \
                or self.fold["parameters"].get("nbsubs", None) != nbsubs \
                or self.fold["parameters"].get("coo_dim", None) != coo_dim \
                or self.fold["parameters"].get("grain", None) != grain:
            if coo_dim is None:
                vals = None
                grain = None
            elif self.isGeospatial() and coo_dim < 0 and abs(coo_dim)-1 < len(self.getCoords()):
                vals = self.getCoordPoints()[:, abs(coo_dim)-1]
            else:  # is in the variables
                if coo_dim-1 >= self.nbCols(0):
                    col = self.col(1, coo_dim-self.nbCols(0)-1)
                else:
                    col = self.col(0, coo_dim-1)
                vals = col.getVector()

            self.fold = {"source": "auto",
                         "parameters": {"coo_dim": coo_dim, "grain": grain, "nbsubs": nbsubs},
                         "folds": self.rsubsets_fold(nbsubs, vals, grain)}
            skeys = ["%d" % i for (i, v) in enumerate(self.fold['folds'])]
            self.fold["fold_ids"] = dict([(v, k) for (k, v) in enumerate(skeys)])
        return self.fold['folds']

    def dropLT(self):
        if self.fold is not None:
            if "lt_ids" in self.fold:
                del self.fold["lt_ids"]
                del self.fold["lt_sids"]

    def assignLT(self, learn_sids, test_sids):
        if self.fold is None:
            return
        rids = {"learn": set(), "test": set()}
        for (which, sids) in [("learn", learn_sids), ("test", test_sids)]:
            for sid in sids:
                if sid in self.fold["fold_ids"]:
                    rids[which].update(self.fold["folds"][self.fold["fold_ids"][sid]])
        self.fold["lt_ids"] = rids
        self.fold["lt_sids"] = {"learn": learn_sids, "test": test_sids}

    def hasFolds(self):
        return self.fold is not None

    def hasAutoFolds(self):
        return self.fold is not None and self.fold["source"] == "auto"

    def hasLT(self):
        return self.fold is not None and "lt_ids" in self.fold

    def getLT(self):
        if self.hasLT():
            return self.fold["lt_ids"]
        else:
            return {}

    def getLTsids(self):
        if self.hasLT():
            return self.fold["lt_sids"]
        else:
            return {}

    def getFoldsInfo(self):
        return self.fold

    def extractFolds(self, side, colid):
        folds = None
        if self.isTypeId(self.col(side, colid).typeId(), "Categorical"):
            col = self.col(side, colid)
            col.setDisabled()
            folds = dict([(re.sub("^F:", "", f), set(fsupp)) for (f, fsupp) in col.iter_cats()])
            skeys = sorted(folds.keys())
            self.fold = {"source": "data",
                         "parameters": {"side": side, "colid": colid, "colname": col.getName()},
                         "fold_ids": dict([(v, k) for (k, v) in enumerate(skeys)]),
                         "folds": [folds[k] for k in skeys]}
        elif self.isTypeId(self.col(side, colid).typeId(), "Boolean"):
            col = self.col(side, colid)
            col.setDisabled()
            folds = {"True": set(col.supp()), "False": set(col.negSuppTerm(None))}
            skeys = sorted(folds.keys())
            self.fold = {"source": "data",
                         "parameters": {"side": side, "colid": colid, "colname": col.getName()},
                         "fold_ids": dict([(v, k) for (k, v) in enumerate(skeys)]),
                         "folds": [folds[k] for k in skeys]}

        return folds

    def getFoldsStats(self, side, colid):
        folds = numpy.array(self.col(side, colid).getVector())
        counts_folds = 1.*numpy.bincount(folds)
        nb_folds = len(counts_folds)
        return {"folds": folds, "counts_folds": counts_folds, "nb_folds": nb_folds}

    def findCandsFolds(self, strict=False):
        folds = self.getColsByName(csv_reader.FOLDS_PATT)
        if strict:
            return folds
        more = self.getColsMoreFolds(folds)
        return folds + more

    def getColsByName(self, pattern):
        results = []
        for side in self.getSides():
            for col in self.colsSide(side):
                if re.search(pattern, col.getName()):
                    results.append((side, col.getId()))
        return results

    def getColsMoreFolds(self, folds=[], nbcats=[1, 6], inclB=True):
        results = []
        for side in self.getSides():
            for col in self.colsSide(side):
                (sito, ci) = (side, col.getId())
                if ((sito, ci) not in folds) and (
                        (self.isTypeId(col.typeId(), "Categorical") and col.nbCats() > nbcats[0] and col.nbCats() < nbcats[1])
                        or (self.isTypeId(col.typeId(), "Boolean") and inclB)):
                    results.append((sito, ci))
        return results

    # creating folds subsets
    def rsubsets_fold(self, nbsubs=10, fold_vals=None, grain=10.):
        # uv, uids = numpy.unique(numpy.mod(numpy.floor(self.getCoords()[0]*grain),nbsubs), return_inverse=True)
        # return [set(numpy.where(uids==uv[i])[0]) for i in range(len(uv))]
        if fold_vals is not None:
            uv, uids = numpy.unique(fold_vals, return_inverse=True)
            if len(uv) > nbsubs:
                nv = numpy.floor(fold_vals*grain)
                uv, uids = numpy.unique(nv, return_inverse=True)
                sizes = [(len(uv)//nbsubs, nbsubs - len(uv) % nbsubs), (len(uv)//nbsubs+1, len(uv) % nbsubs)]
                maps_to = numpy.hstack([[i]*sizes[0][0] for i in range(sizes[0][1])]+[[i+sizes[0][1]]*sizes[1][0] for i in range(sizes[1][1])])
                numpy.random.shuffle(maps_to)
                subsets_ids = [set() for i in range(nbsubs)]
                for i in range(len(uv)):
                    subsets_ids[maps_to[i]].update(numpy.where(uids == i)[0])
            else:
                subsets_ids = [set(numpy.where(uids == i)[0]) for i in range(len(uv))]
        else:
            sizes = [(self.nbRows()//nbsubs, nbsubs - self.nbRows() % nbsubs), (self.nbRows()//nbsubs+1, self.nbRows() % nbsubs)]
            maps_to = numpy.hstack([[i]*sizes[0][0] for i in range(sizes[0][1])]+[[i+sizes[0][1]]*sizes[1][0] for i in range(sizes[1][1])])
            numpy.random.shuffle(maps_to)
            subsets_ids = [set() for i in range(nbsubs)]
            for i in range(nbsubs):
                subsets_ids[i].update(numpy.where(maps_to == i)[0])
        return subsets_ids

    # old version for reremi run subfolds
    def get_LTfold(self, row_idsT):
        row_idsL = set(range(self.nbRows()))
        row_idsT = row_idsL.intersection(row_idsT)
        row_idsL.difference_update(row_idsT)
        return self.subset(row_idsL), self.subset(row_idsT)

    def addFoldsCol(self, subsets=None, sito=1):
        suff = "cust"
        if subsets is None and self.fold is not None:
            if self.fold["source"] == "auto":
                subsets = dict([(k, self.fold["folds"][kk]) for (k, kk) in self.fold["fold_ids"].items()])
                suff = "%dc%sg%s" % (self.nbCols(sito), (self.fold["parameters"]["coo_dim"] or "N"), (self.fold["parameters"]["grain"] or "N"))
        if type(subsets) is list:
            subsets = dict(enumerate(subsets))
        if subsets is not None and type(subsets) is dict:
            self.addVecCol(dict([("F:%s" % i, s) for (i, s) in subsets.items()]), csv_reader.FOLDS_PREF+suff, sito)
            return suff
    #######################

    def addVecCol(self, vec, name=None, sito=1, force_type=None):
        vname = "v%d" % self.nbCols(sito)
        if name is not None:
            # MAKE SURE THIS IS COHERENT WITH FOLDS_PATT
            vname = vname+":"+name
        col = None
        if force_type == "clusters" or force_type == "support":
            ks = numpy.unique(vec)
            lbls = []
            if force_type == "support":
                lbls = SSetts.labels_sparts
            if numpy.max(ks) > len(lbls):
                lbls = ["S%d" % i for i in range(numpy.max(ks)+1)]
            parts = {}
            for k in ks:
                parts[lbls[k]] = set(numpy.where(vec == k)[0])
            vec, force_type = (parts, None)

        if force_type == "values" and isinstance(vec, ColM):
            if name is not None:
                vname = name
            col = vec
        elif force_type is not None:
            cclass = self.getColClassForName(force_type)
            col = cclass.parseList(vec, force=True)
        elif type(vec) is set:
            col = BoolColM(vec, self.nbRows())
        elif type(vec) is dict:
            col = CatColM(vec, self.nbRows())
        else:
            col = CatColM(vec, self.nbRows())

        if col is not None:
            col.setSideIdName(sito, self.nbCols(sito), vname)
            self.appendCol(col, sito)

    def subset(self, row_ids=None, shuff_side=None):
        coords = None
        rnames = None
        if row_ids is None:
            N = self.nbRows()
        else:
            if type(row_ids) is set:
                row_ids = sorted(row_ids)
            if type(row_ids) is list:
                row_ids = dict([(r, [ri]) for (ri, r) in enumerate(row_ids)])
            N = sum([len(news) for news in row_ids.values()])
        if self.rnames is not None:
            if row_ids is None:
                rnames = list(self.rnames)
            else:
                rnames = ["-" for i in range(N)]
                for old, news in row_ids.items():
                    for new in news:
                        rnames[new] = self.rnames[old]
        if self.coords is not None:
            if row_ids is None:
                coords = [list(self.coords[i]) for i in range(len(self.coords))]
            else:
                maps_to = numpy.array([0 for i in range(N)])
                for old, news in row_ids.items():
                    maps_to[news] = old
                coords = [[self.coords[i][j] for j in maps_to] for i in range(len(self.coords))]

        cols = [[], []]
        for side in [0, 1]:
            sss = None
            if shuff_side is None or (side in shuff_side):
                if type(shuff_side) is dict and shuff_side[side] is not None:
                    sss = shuff_side[side]
                    if type(sss) is list:
                        sss = dict([(v, [i]) for i, v in enumerate(sss)])
                else:
                    sss = row_ids

            for col in self.colsSide(side):
                if sss is not None:
                    tmp = col.subsetCol(sss)
                else:  # just copy
                    tmp = col.subsetCol()
                tmp.setSideIdName(side, len(cols[side]))
                cols[side].append(tmp)
        return Data(cols, N, coords, rnames, self.isSingleD())

    def shuffle_sides(self, shuff_sides):
        shuffled_rids = range(self.nbRows())
        numpy.random.shuffle(shuffled_rids)
        shuffled_d = dict([(v, [i]) for i, v in enumerate(shuffled_rids)])
        return self.subset(shuffled_d, shuff_sides), shuffled_d

    def hasSelectedRows(self):
        return len(self.selected_rows) > 0

    def selectedRows(self):
        return self.selected_rows

    def getVizRows(self, details={}):
        if details is not None and details.get("rset_id", None) in self.getLT():
            if len(self.selected_rows) == 0:
                return set(self.getLT()[details["rset_id"]])
            return set(self.getLT()[details["rset_id"]]) - self.selected_rows
        return self.nonselectedRows()

    def getUnvizRows(self, details={}):
        if details is not None and details.get("rset_id", None) in self.getLT():
            return self.selected_rows.union(*[s for (k, s) in self.getLT().items() if k != details["rset_id"]])
        return self.selected_rows

    def nonselectedRows(self):
        return self.rows() - self.selected_rows

    def addSelectedRow(self, rid):
        self.selected_rows.add(rid)

    def removeSelectedRow(self, rid):
        self.selected_rows.discard(rid)

    def getRows(self):
        return [RowE(i, self) for i in range(self.nbRows())]

    def getRow(self, i):
        if i >= 0 and i < self.nbRows():
            return RowE(i, self)

    def getInfo(self):
        miss_str = ""
        if self.hasMissing():
            miss_str = " (some entries missing)"
        return "Data: %s x %s%s, using %s Jaccard index" % (self.rowsInfo(), self.colsInfo(), miss_str, self.getSSetts().getTypeParts())

    def __str__(self):
        dts_str = ""
        if self.hasMissing():
            dts_str += "M"
        if self.isSymmetric():
            dts_str += "S"
        if len(dts_str) > 0:
            dts_str = " ("+dts_str+")"
        return "Data: %s x %s%s" % (self.rowsInfo(), self.colsInfo(), dts_str)
        # if self.nbRowsEnabled() == self.nbRows() and \
        #   self.nbColsEnabled(0) == self.nbCols(0) and self.nbColsEnabled(1) == self.nbCols(1):
        #     return "%i x %i+%i data" % ( self.nbRows(), self.nbCols(0), self.nbCols(1))
        # return "%i(+%i) x %i(+%i)+%i(+%i) data" \
        #   % ( self.nbRowsEnabled(), self.nbRowsDisabled(),
        #       self.nbColsEnabled(0), self.nbColsDisabled(0),
        #       self.nbColsEnabled(1), self.nbColsDisabled(1))

    def colsInfo(self):
        if self.nbColsEnabled(0) == self.nbCols(0) and self.nbColsEnabled(1) == self.nbCols(1):
            return "%i+%i" % (self.nbCols(0), self.nbCols(1))
        return "%i(+%i)+%i(+%i)" \
            % (self.nbColsEnabled(0), self.nbColsDisabled(0),
               self.nbColsEnabled(1), self.nbColsDisabled(1))

    def rowsInfo(self):
        if self.nbRowsEnabled() == self.nbRows():
            return "%i" % self.nbRows()
        return "%i(+%i)" % (self.nbRowsEnabled(), self.nbRowsDisabled())

    def writeCSV(self, outputs, thres=0.1, full_details=False, inline=False, cats_int_fmt=None):
        if FORCE_WRITE_DENSE:
            thres = 0.
        # check whether some row name is worth storing
        rids = {}
        if self.rnames is not None:
            rids = dict(enumerate([prepareRowName(rname, i, self) for i, rname in enumerate(self.rnames)]))
        elif len(self.selectedRows()) > 0:
            rids = dict(enumerate([prepareRowName(i+1, i, self) for i in range(self.N)]))
        mean_denses = [numpy.mean([col.density() for col in self.colsSide(0)]),
                       numpy.mean([col.density() for col in self.colsSide(1)])]
        # FORCE SPARSE mean_denses = [0,0]
        # FORCE DENSE  mean_denses = [1,1]
        argmaxd = 0
        if mean_denses[0] < mean_denses[1]:
            argmaxd = 1
        if mean_denses[1-argmaxd] > thres:  # BOTH SIDES ARE DENSE
            styles = {argmaxd: {"meth": "dense", "details": True, "cats_int_fmt": cats_int_fmt},
                      1-argmaxd: {"meth": "dense", "details": full_details, "cats_int_fmt": cats_int_fmt}}
        elif mean_denses[argmaxd] > thres:  # ONE SIDE IS DENSE
            methot = "triples"
            if not self.hasDisabledCols(1-argmaxd) and not self.hasGroups(1-argmaxd) and sum([not col.simpleBool() for col in self.colsSide(1-argmaxd)]) == 0:
                methot = "pairs"
            styles = {argmaxd: {"meth": "dense", "details": True, "cats_int_fmt": cats_int_fmt},
                      1-argmaxd: {"meth": methot, "details": full_details, "inline": inline, "cats_int_fmt": cats_int_fmt}}
        else:  # BOTH SIDES ARE SPARSE
            simpleBool = [sum([not col.simpleBool() for col in self.colsSide(0)]) == 0,
                          sum([not col.simpleBool() for col in self.colsSide(1)]) == 0]
            if self.isGeospatial() or len(rids) > 0:
                if not simpleBool[1-argmaxd]:  # is not only boolean so can have names and coords
                    methot = "pairs"
                    if self.hasDisabledCols(argmaxd) or self.hasGroups(argmaxd) or not simpleBool[argmaxd]:
                        methot = "triples"
                    styles = {argmaxd: {"meth": methot, "details": full_details, "cats_int_fmt": cats_int_fmt},
                              1-argmaxd: {"meth": "triples", "details": True, "inline": inline, "cats_int_fmt": cats_int_fmt}}
                else:  # otherwise argmax has it
                    methot = "pairs"
                    if self.hasDisabledCols(1-argmaxd):
                        methot = "triples"
                    styles = {argmaxd: {"meth": "triples", "details": True, "inline": inline, "cats_int_fmt": cats_int_fmt},
                              1-argmaxd: {"meth": methot, "details": full_details, "cats_int_fmt": cats_int_fmt}}
            else:
                styles = {argmaxd: {"meth": "pairs", "details": full_details, "cats_int_fmt": cats_int_fmt},
                          1-argmaxd: {"meth": "pairs", "details": full_details, "cats_int_fmt": cats_int_fmt}}
                for side in [0, 1]:
                    if not simpleBool[side] or len(cids[side]) > 0:
                        styles[side]["meth"] = "triples"
                        styles[side]["inline"] = inline

        meths = {"pairs": self.writeCSVSparsePairs, "triples": self.writeCSVSparseTriples, "dense": self.writeCSVDense}
        # meths = {"pairs": self.writeCSVDense, "triples": self.writeCSVDense, "dense": self.writeCSVDense}
        sides = [0, 1]
        if self.isSingleD() and (len(outputs) == 1 or outputs[0] == outputs[1] or outputs[1] is None):
            sides = [0]
            styles[0]["details"] = True
            styles[0]["single_dataset"] = True
        for side in sides:
            # check whether some column name is worth storing
            cids = {}
            if sum([col.getName() != col.getId() or not col.isEnabled() or col.hasGroupId() for col in self.colsSide(side)]) > 0:
                type_smap = None
                if full_details and styles[side]["meth"] == "dense":
                    type_smap = {}
                cids = dict(enumerate([prepareColumnName(col, type_smap) for col in self.colsSide(side)]))
                meth = meths[styles[side].pop("meth")]
            with open(outputs[side], "w") as fp:
                csvf = csv_reader.start_out(fp)
                meth(side, csvf, rids=rids, cids=cids, **styles[side])

    def writeCSVDense(self, side, csvf, rids={}, cids={}, details=True, inline=False, single_dataset=False, cats_int_fmt=None):
        discol = []
        groupscol = []
        if self.hasDisabledCols(side) or (single_dataset and self.hasDisabledCols()):
            discol.append(csv_reader.ENABLED_COLS[0])
            # ADD VALUES IN THE ENABLED COL TO FILL FOR SPECIAL COLS
            if details and self.hasSelectedRows():
                discol.append(0)
            if details and self.isGeospatial():
                discol.append(0)
                discol.append(0)
            if details and self.isConditional() and not self.isGeoConditional():
                discol.extend([0 for i in self.getColsC()])

            for col in self.colsSide(side):
                cid = col.getId()
                if single_dataset:
                    discol.append(Data.enabled_codes[(self.col(0, cid).getEnabled(), self.col(1, cid).getEnabled())])
                else:
                    discol.append(Data.enabled_codes[col.getEnabled()])

        if self.hasGroups(side) or (single_dataset and self.hasGroups()):
            groupscol.append(csv_reader.GROUPS_COLS[0])
            # ADD VALUES IN THE ENABLED COL TO FILL FOR SPECIAL COLS
            if details and self.hasSelectedRows():
                groupscol.append(-1)
            if details and self.isGeospatial():
                groupscol.append(-1)
                groupscol.append(-1)
            if details and self.isConditional() and not self.isGeoConditional():
                groupscol.extend([-1 for i in self.getColsC()])

            for col in self.colsSide(side):
                groupscol.append(col.getGroupId())

        header = []
        colsC = None
        inc_ids = (details and len(rids) > 0) or len(discol) > 0 or len(groupscol) > 0
        if inc_ids:
            header.append(csv_reader.IDENTIFIERS[0])
        if details and self.hasSelectedRows():
            header.append(csv_reader.ENABLED_ROWS[0])
        if details and self.isGeospatial():
            header.append(csv_reader.LONGITUDE[0])
            header.append(csv_reader.LATITUDE[0])
        if details and self.isConditional() and not self.isGeoConditional():
            header.extend([colC.getName() for colC in self.getColsC()])

        for col in self.colsSide(side):
            col.getVector()
            cid = col.getId()
            if len(header) > 0 or len(cids) > 0:
                header.append(cids.get(cid, cid))
        letter = self.getCommonType(side)
        if letter is not None:
            if len(header) == 0:
                header.append("")
            header[-1] += " # type=%s" % letter
            # header.append("type=%s" % letter)

        if len(header) > 0:
            csv_reader.write_row(csvf, header)
        if len(discol) > 0:
            csv_reader.write_row(csvf, discol)
        if len(groupscol) > 0:
            csv_reader.write_row(csvf, groupscol)

        for n in range(self.N):
            row = []
            if inc_ids:
                row.append(rids.get(n, n))
            if details and self.hasSelectedRows():
                row.append(Data.enabled_codes[n not in self.selectedRows()])
            if details and self.isGeospatial():
                row.append(":".join(map(str, self.coords[0][n])))
                row.append(":".join(map(str, self.coords[1][n])))

            if details and self.isConditional() and not self.isGeoConditional():
                for colC in self.getColsC():
                    if Data.isTypeId(colC.typeId(), "Categorical") and cats_int_fmt is not None:
                        row.append("%s" % colC.valToStr(colC.getValue(n, pref=cats_int_fmt)))
                    else:
                        row.append("%s" % colC.valToStr(col.getValue(n, pref="bnum")))                        
                    # if self.isTimeConditional():
                    #     row.append("TTA%s" % colC.valToTimeStr(v))
                    # else:
                    #     row.append("TTB%s" % colC.valToStr(v))

            for col in self.colsSide(side):
                if Data.isTypeId(col.typeId(), "Categorical") and cats_int_fmt is not None:
                    row.append("%s" % col.valToStr(col.getValue(n, pref=cats_int_fmt)))
                else:
                    row.append("%s" % col.valToStr(col.getValue(n, pref="bnum")))
                # row.append(col.valToStr(col.getValue(n, pref="bnum")))
            csv_reader.write_row(csvf, row)

    def writeCSVSparseTriples(self, side, csvf, rids={}, cids={}, details=True, inline=False, single_dataset=False, cats_int_fmt=None):
        header = [csv_reader.IDENTIFIERS[0], csv_reader.COLVAR[0], csv_reader.COLVAL[0]]
        letter = self.getCommonType(side)
        if letter is not None:
            header[-1] += " # type=%s" % letter
        csv_reader.write_row(csvf, header)
        if not inline:
            trids, tcids = {}, {}
        else:
            trids, tcids = rids, cids

        if details and self.isGeospatial():
            for n in range(self.N):
                csv_reader.write_row(csvf, [trids.get(n, n), csv_reader.LONGITUDE[0], ":".join(map(str,  self.coords[0][n]))])
                csv_reader.write_row(csvf, [trids.get(n, n), csv_reader.LATITUDE[0], ":".join(map(str,  self.coords[1][n]))])

        for n in self.selectedRows():
            csv_reader.write_row(csvf, [trids.get(n, n), csv_reader.ENABLED_ROWS[0], "F"])

        if self.hasDisabledCols(side) or (single_dataset and self.hasDisabledCols()):
            for col in self.colsSide(side):
                cid = col.getId()
                if single_dataset:
                    tmp = Data.enabled_codes[(self.col(0, cid).getEnabled(), self.col(1, cid).getEnabled())]
                else:
                    tmp = Data.enabled_codes[col.getEnabled()]
                if tmp != Data.enabled_codes[1]:
                    if inline:
                        csv_reader.write_row(csvf, [csv_reader.ENABLED_COLS[0], cids.get(cid, cid), tmp])
                    else:
                        csv_reader.write_row(csvf, [csv_reader.ENABLED_COLS[0], cid, tmp])

        if self.hasGroups(side) or (single_dataset and self.hasGroups()):
            for col in self.colsSide(side):
                cid = col.getId()
                if inline:
                    csv_reader.write_row(csvf, [csv_reader.GROUPS_COLS[0], cids.get(cid, cid), col.getGroupId()])
                else:
                    csv_reader.write_row(csvf, [csv_reader.GROUPS_COLS[0], cid, col.getGroupId()])

        fillR = False
        # if names are not written inline, add the entities names now, that serves to recovers the correct number of lines
        if details and len(rids) > 0 and not inline:
            for n in range(self.N):
                csv_reader.write_row(csvf, [n, -1, rids.get(n, n)])
        else:
            # otherwise it will need fill to recover the number of lines
            fillR = True

        for col in self.colsSide(side):
            ci = col.getId()
            fillC = False
            if not inline and len(cids) > 0:
                # if names are not written inline, add the variable's name now
                csv_reader.write_row(csvf, [-1, ci, cids.get(ci, ci)])
            else:
                fillC = True

            if ci == 0 and fillR:
                tmp = col.toList(sparse=True, fill=False, cats_int_fmt=cats_int_fmt)
                non_app = col.rows().difference(list(zip(*tmp))[0])
                for (n, v) in tmp:
                    csv_reader.write_row(csvf, [trids.get(n, n), tcids.get(ci, ci), v])
                for n in non_app:
                    csv_reader.write_row(csvf, [trids.get(n, n), tcids.get(ci, ci), 0])
            else:
                for (n, v) in col.toList(sparse=True, fill=fillR, cats_int_fmt=cats_int_fmt):
                    fillC = False
                    csv_reader.write_row(csvf, [trids.get(n, n), tcids.get(ci, ci), v])
                if fillC:
                    # Filling for column if it does not have any entry
                    csv_reader.write_row(csvf, [trids.get(0, 0), tcids.get(ci, ci), 0])

    # THIS FORMAT ONLY ALLOWS BOOLEAN WITHOUT COORS, IF NAMES THEY HAVE TO BE INLINE

    def writeCSVSparsePairs(self, side, csvf, rids={}, cids={}, details=True, inline=False, single_dataset=False, cats_int_fmt=None):
        # cats_int_fmt not relevant for this format, only listed for function signature compatibility
        # self.writeCSVSparseTriples(side, csvf, rids, cids, details, inline, single_dataset)
        # return
        header = [csv_reader.IDENTIFIERS[0], csv_reader.COLVAR[0]]
        letter = self.getCommonType(side)
        # if letter is not None: ### THIS CAN ONLY BE FULLY BOOLEAN...
        #     if len(header) == 0:
        #         header.append("")
        #     header[-1] += " # type=%s" % letter
        #     # header.append("type=%s" % letter)
        csv_reader.write_row(csvf, header)
        if not details:
            rids = {}
        for col in self.colsSide(side):
            ci = col.getId()
            for (n, v) in col.toList(sparse=True, fill=False, cats_int_fmt=cats_int_fmt):
                csv_reader.write_row(csvf, [rids.get(n, n), cids.get(ci, ci)])

    def disp(self):
        strd = str(self) + ":\n"
        strd += 'Left Hand side columns:\n'
        for col in self.colsSide(0):
            strd += "\t%s\n" % col
        strd += 'Right Hand side columns:\n'
        for col in self.colsSide(1):
            strd += "\t%s\n" % col
        return strd

    def rows(self):
        return set(range(self.N))

    def nbRows(self):
        return self.N

    def nbRowsEnabled(self):
        return self.N-len(self.selected_rows)

    def nbRowsDisabled(self):
        return len(self.selected_rows)

    def nbCols(self, side):
        return self.getLen(side)

    def nbColsEnabled(self, side):
        return len([c for c in self.colsSide(side) if c.isEnabled()])

    def nbColsDisabled(self, side):
        return len([c for c in self.colsSide(side) if not c.isEnabled()])

    def name(self, side, literal):
        return self.col(side, literal).getName()

    def supp(self, side, literal):
        return self.col(side, literal).suppLiteral(literal)

    def OutSupp(self, side, literal):
        return set(range(self.nbRows())) - self.supp(side, literal) - self.miss(side, literal)

    def miss(self, side, literal):
        return self.col(side, literal).miss()

    def literalSuppMiss(self, side, literal):
        return (self.supp(side, literal), self.miss(side, literal))

    def usableIds(self, min_in=-1, min_out=-1):
        return [[col.getId() for col in self.colsSide(0) if col.usable(min_in, min_out)],
                [col.getId() for col in self.colsSide(1) if col.usable(min_in, min_out)]]

    def getDisabledCols(self, side=None):
        dis = []
        for s in self.getSides(side):
            for col in self.colsSide(s):
                if not col.isEnabled():
                    dis.append((s, col.getId()))
        return dis

    def hasDisabledCols(self, side=None):
        for s in self.getSides(side):
            for col in self.colsSide(s):
                if not col.isEnabled():
                    return True
        return False

    def getIids(self, side=None):
        iids = []
        for s in self.getSides(side):
            iids.extend([(s, cc.getId()) for cc in self.colsSide(s)])
        return iids

    def isConditional(self):
        return self.condition_dt is not None

    def isTimeConditional(self):
        return self.condition_dt is not None and self.condition_dt.get("is_time", False)

    def isGeoConditional(self):
        return self.condition_dt is not None and self.condition_dt.get("is_geo", False)

    def getTimeCoord(self):
        if self.isTimeConditional() and len(self.condition_dt["cols"]) > 0:
            return self.condition_dt["cols"][0].getVect()

    def getTimeCoordCol(self):
        if self.isTimeConditional() and len(self.condition_dt["cols"]) > 0:
            return self.condition_dt["cols"][0]

    def getColsC(self):
        if self.isConditional():
            return self.condition_dt["cols"]

    def getTypeD(self):
        return {"geo": self.isGeospatial(), "time": self.isTimeConditional(), "cond": self.isConditional()}

    def isGeospatial(self):
        return self.coords is not None

    def hasGeoPoly(self):
        return self.coords is not None and (min([len(cs) for cs in self.coords[0]]) > 2)

    def getCoords(self):
        return self.coords

    def getCoordPoints(self):
        return self.coords_points

    def getCoordsExtrema(self):
        if self.isGeospatial():
            return [min(chain.from_iterable(self.coords[0])), max(chain.from_iterable(self.coords[0])), min(chain.from_iterable(self.coords[1])), max(chain.from_iterable(self.coords[1]))]
        return None

    def hasRNames(self):
        return self.rnames is not None

    def hasNames(self):
        for side in self.getSides():
            for col in self.colsSide(side):
                if col.hasName():
                    return True
        return False

    def getFmts(self, side=None):
        if side is None:
            return [[col.getFmt() for col in self.colsSide(side)] for side in self.getSides()]
        else:
            return [col.getFmt() for col in self.colsSide(side)]

    def getNames(self, side=None):
        if side is None:
            return [[col.getName() for col in self.colsSide(side)] for side in self.getSides()]
        else:
            return [col.getName() for col in self.colsSide(side)]

    def setNames(self, names):
        for side in self.getSides():
            try:
                if len(names[side]) == self.nbCols(side):
                    for i, col in enumerate(self.colsSide(side)):
                        col.name = names[side][i]
                else:
                    raise DataError('Number of names does not match number of variables!')
            except IndexError:
                raise DataError('Number of names does not match number of sides!')

    def setCoords(self, coords):
        # coords are NOT turned to a numpy array because polygons might have different numbers of points
        self.coords = None
        self.coords_points = None
        if coords is not None:
            if (len(coords) == 2 and len(coords[0]) == self.nbRows()):
                coords_tmp = coords
                coords_points_tmp = numpy.array([[coords[0][i][0], coords[1][i][0]] for i in range(len(coords[0]))])
                # check for duplicates and randomize
                ids_miss = numpy.where((coords_points_tmp[:, 1] == -361) & (coords_points_tmp[:, 0] == -361))[0]
                ids_pres = list(set(range(self.nbRows())).difference(ids_miss))
                # keys_cc = ["%s:%s" % (coords_points_tmp[v,0], coords_points_tmp[v,1]) for v in ids_pres]
                # pdb.set_trace()
                # if len(ids_pres) > len(set(keys_cc)):
                #     print(len(ids_pres), len(set(keys_cc)))
                miss_cc = (numpy.min(coords_points_tmp[ids_pres, 0]), numpy.min(coords_points_tmp[ids_pres, 1]))
                coords_points_tmp[ids_miss, 0] = miss_cc[0]
                coords_points_tmp[ids_miss, 1] = miss_cc[1]
                for cci in ids_miss:
                    coords_tmp[0][cci] = [miss_cc[0]]
                    coords_tmp[1][cci] = [miss_cc[1]]
                self.coords_points = coords_points_tmp
                self.coords = coords_tmp
            else:
                raise DataError('Number of coordinates does not match number of entities!')

    def prepareGeoCond(self):
        if self.isGeospatial():
            geo_cols = self.prepareGeoCols("cond_geo")
            self.condition_dt = {"cols": geo_cols, "is_geo": True}

    def prepareGeoCols(self, cname="coord_geo"):
        geo_cols = []
        if self.isGeospatial():
            coords_points = self.getCoordPoints()
            for ci in range(coords_points.shape[1]):
                ncolSupp = [(v, k) for (k, v) in enumerate(coords_points[:, ci])]
                col = NumColM(ncolSupp, N=coords_points.shape[0])
                col.setId(ci)
                col.side = -1
                col.name = cname + ("%d" % ci)
                geo_cols.append(col)
        return geo_cols

    def enableAll(self):
        for side in self.getSides():
            for i in range(self.nbCols(side)):
                self.col(side, i).setEnabled()  # reset all to active
        self.selected_rows = set()  # reset all to active

    def applyDisableMasks(self, mask_vars_LHS=None, mask_vars_RHS=None, mask_rows=None):
        m0 = self.applyDisableMaskCols(0, mask_vars_LHS)
        m1 = self.applyDisableMaskCols(1, mask_vars_RHS)
        mR = self.applyDisableMaskRows(mask_rows)
        return (m0, m1, mR)

    def applyDisableMaskCols(self, side, mask_vars=None):
        on, off = (set(), set())
        if mask_vars is not None and len(mask_vars) > 0:
            if ";" in mask_vars:
                mask_parts = mask_vars.strip().split(";")
            else:
                mask_parts = mask_vars.strip().split(",")
            if len(mask_parts) > 0:
                if all([re.match("^[0-9]*:[TF]$", p) for p in mask_parts]):
                    for p in mask_parts:
                        pps = p.split(":")
                        which, what = (int(pps[0]), self.enabled_codes_rev_simple.get(pps[1]))
                        if which < self.nbCols(side):
                            if what:
                                self.col(side, which).setEnabled()
                                on.add(which)
                            else:
                                self.col(side, which).setDisabled()
                                off.add(which)
                elif all([re.match("^[TF]$", p) for p in mask_parts]):
                    for which, p in enumerate(mask_parts):
                        what = self.enabled_codes_rev_simple.get(p)
                        if which < self.nbCols(side):
                            if what:
                                self.col(side, which).setEnabled()
                                on.add(which)
                            else:
                                self.col(side, which).setDisabled()
                                off.add(which)
                elif all([re.match("^[0-9]*$", p) for p in mask_parts[1:]]):
                    pfrom = 1
                    def_stat = True
                    def_set, other_set = (on, off)
                    if re.match("^[0-9]*$", mask_parts[0]):
                        pfrom = 0
                    elif mask_parts[0] == "T":
                        def_stat = False
                        def_set, other_set = (off, on)

                    for i in range(self.nbCols(side)):
                        self.col(side, i).setStatus(def_stat)  # reset all to active
                    for pp in mask_parts[pfrom:]:
                        which = int(pp)
                        if which < self.nbCols(side):
                            self.col(side, which).setStatus(not def_stat)
                            other_set.add(which)
                    def_set.update(set(range(self.nbCols(side))).difference(other_set))
        return on, off

    def applyDisableMaskRows(self, mask_rows=None):
        on, off = (set(), set())
        if mask_rows is not None and len(mask_rows) > 0:
            mask_parts = mask_rows.strip().split(",")
            if len(mask_parts) > 0:
                if all([re.match("^[0-9]*:[TF]$", p) for p in mask_parts]):
                    for p in mask_parts:
                        pps = p.split(":")
                        which, what = (int(pps[0]), self.enabled_codes_rev_simple.get(pps[1]))
                        if which < self.nbRows():
                            if what:
                                self.removeSelectedRow(which)
                                on.add(which)
                            else:
                                self.addSelectedRow(which)
                                off.add(which)
                elif all([re.match("^[TF]$", p) for p in mask_parts]):
                    for which, p in enumerate(mask_parts):
                        what = self.enabled_codes_rev_simple.get(p)
                        if which < self.nbRows():
                            if what:
                                self.removeSelectedRow(which)
                                on.add(which)
                            else:
                                self.addSelectedRow(which)
                                off.add(which)
                elif all([re.match("^[0-9]*$", p) for p in mask_parts]):
                    self.selected_rows = set()  # reset all to active
                    for pp in mask_parts:
                        which = int(pp)
                        if which < self.nbRows():
                            self.addSelectedRow(which)
                            off.add(which)
                    on = set(range(self.nbRows()).difference(off))
        return on, off


############################################################################
# READING METHODS
############################################################################
def readDNCFromCSVFiles(filenames):
    cols, N, coords, rnames, cond_col = [[], []], 0, None, None, None
    csv_params = {}
    other_params = {}
    single_dataset = False
    unknown_string = NA_str_c
    if len(filenames) >= 2:
        left_filename = filenames[0]
        right_filename = filenames[1]
        if len(filenames) >= 3:
            csv_params = filenames[2]
            if len(filenames) >= 4 and filenames[3] is not None:
                if type(filenames[3]) is dict:
                    other_params = filenames[3]
                else:
                    other_params["NA_str"] = filenames[3]
        try:
            tmp_data, single_dataset, unknown_str = csv_reader.importCSV(left_filename, right_filename, csv_params, other_params)
        except ValueError as arg:
            raise DataError('Data error reading csv: %s' % arg)
        # except csv_reader.CSVRError as arg:
        #     raise DataError(str(arg).strip("'"))
        cols, N, coords, rnames, disabled_rows, cond_col = parseDNCFromCSVData(tmp_data, single_dataset)

    return cols, N, coords, rnames, disabled_rows, cond_col, single_dataset, unknown_string


def prepareRowName(rname, rid=None, data=None):
    return "%s" % rname
    en = ""
    if rid is not None and data is not None and rid in data.selectedRows():
        en = "_"
    return "%s%s" % (en, rname)


def parseRowsNames(rnames):
    names = []
    for i, rname in enumerate(rnames):
        if rname is None:
            names.append("%d" % (i+1))
        else:
            names.append(rname)
    return names


def prepareColumnName(col, types_smap={}):
    return "%s" % col.getName()
    en = ""
    if not col.isEnabled():
        en = "_"
    if types_smap is None:
        return "%s%s" % (en, col.getName())
    else:
        return "%s[%s]%s" % (en, types_smap.get(col.typeId(), col.typeId()), col.getName())


def parseColumnName(name, types_smap={}):
    tmatch = re.match("^(\[(?P<type>[0-9])\])?(?P<name>.*)$", name)
    det = {"name": tmatch.group("name")}
    if tmatch.group("type") is not None and tmatch.group("type") in types_smap:
        det["type"] = types_smap[tmatch.group("type")]
    return name, det


def prepareColFromVector(vect_data, prec=None, enabled=True, vtype=None):
    if vtype in TYPES_SMAP:
        cClass = TYPES_SMAP[vtype]
    else:
        cClass = NumColM
    return cClass.fromVect(vect_data, prec, enabled)


def prepareSideFromMatrix(mat_data, prec=None, vnames=None, side=0, cols=None, enabled=True, vtypes=None):
    # like getMatrix, one row per variable
    if cols is None:
        cols = {side: []}
    for i in range(mat_data.shape[0]):
        col = prepareColFromVector(mat_data[i, :], prec=getValSpec(prec, i), enabled=getValSpec(enabled, i), vtype=getValSpec(vtypes, i))
        col.setSideIdName(side, len(cols[side]), getValSpec(vnames, i))
        cols[side].append(col)
    return cols


def parseDNCFromCSVData(csv_data, single_dataset=False):
    if single_dataset:
        sides = (0, 0)
    else:
        sides = (0, 1)
    type_ids_org = [CatColM, NumColM, BoolColM]
    types_smap = dict([(str(c.type_id), c) for c in type_ids_org])
    cols = [[], []]
    coords = None
    condition_dt = None
    single_dataset = False

    if csv_data.get("cond_col", None) is not None:
        cond_col = NumColM.parseList(csv_data["cond_col"], native_missing_check=False)
        cname = csv_reader.COND_COL[0]
        if csv_data.get("c_time", False):
            cname = csv_reader.COND_TIME
        cond_col.setSideIdName(-1, 0, cname)
        condition_dt = {"cols": [cond_col], "is_time": csv_data.get("c_time", False)}

    if csv_data.get("coord", None) is not None:
        try:
            tmp = list(zip(*csv_data["coord"]))
            coords = [tmp[1], tmp[0]]
            # coords = numpy.array([tmp[1], tmp[0]], dtype=object)
        except Exception:
            coords = None

    N = len(csv_data['data'][0]["order"])  # THE READER CHECKS THAT BOTH SIDES HAVE SAME SIZE
    if csv_data.get("ids", None) is not None and len(csv_data["ids"]) == N:
        rnames = parseRowsNames(csv_data["ids"])
    else:
        rnames = [Term.pattVName % n for n in range(N)]

    indices = [dict([(v, k) for (k, v) in enumerate(csv_data['data'][sides[0]]["order"])]),
               dict([(v, k) for (k, v) in enumerate(csv_data['data'][sides[1]]["order"])])]
    disabled_rows = set()

    for er in csv_reader.ENABLED_ROWS:
        for side in set(sides):
            if er in csv_data['data'][side]["headers"]:
                csv_data['data'][side]["headers"].remove(er)
                tmp = csv_data['data'][side]["data"].pop(er)
                if type(tmp) is dict:
                    tmp = tmp.items()
                else:
                    tmp = enumerate(tmp)
                for i, v in tmp:
                    if not Data.enabled_codes_rev_simple[v]:
                        disabled_rows.add(indices[side][i])

    for sito, side in enumerate(sides):
        for header in csv_data['data'][side]["headers"]:
            name, det = parseColumnName(header, types_smap)
            if len(name) == 0:
                continue
            values = csv_data['data'][side]["data"][name]
            col = None
            if Data.all_types_map.get(csv_data['data'][side]['type_all']) is not None:
                col = Data.all_types_map[csv_data['data'][side]['type_all']].parseList(values, indices[side], force=True, native_missing_check=False)
                if col is None:
                    print("DID NOT MANAGE GLOBAL TYPE FORCE PARSING...")
            if col is None:
                if "type" in det:
                    col = det["type"].parseList(values, indices[side], force=True, native_missing_check=False)
                    if col is None:
                        print("DID NOT MANAGE COL TYPE FORCE PARSING...")
            if col is None:
                type_ids = list(type_ids_org)
                while col is None and len(type_ids) >= 1:
                    col = type_ids.pop().parseList(values, indices[side], native_missing_check=False)

            if col is not None:
                if col.N == N:
                    if not det.get("enabled", True) or \
                            not Data.enabled_codes_rev_double.get(csv_data["data"][side][csv_reader.ENABLED_COLS[0]].get(name, None), (1, 1))[sito]:
                        col.flipEnabled()
                    if len(csv_data["data"][side][csv_reader.GROUPS_COLS[0]]) > 0:
                        gid = int(csv_data["data"][side][csv_reader.GROUPS_COLS[0]].get(name, -1))
                        col.setGroupId(gid)
                    # if not det.get("enabled", True) or \
                    col.setSideIdName(sito, len(cols[sito]), det.get("name", name))
                    cols[sito].append(col)
                else:
                    raise DataError('Unexpected number of rows in column %s on side %s! (%d instead of %d)' % (name, side, col.N, N))
            else:
                raise DataError('Unrecognized variable type for column %s on side %s!' % (name, side))
    return (cols, N, coords, rnames, disabled_rows, condition_dt)


def main():
    pass
    # rep = "/home/egalbrun/short/cats/"
    # nan="-nan"
    # data = Data([rep+"wagner_LHS"+nan+".csv", rep+"wagner_RHS"+nan+".csv", {}, "nan"], "csv")
    # # print("ORG", "\n", data.disp())
    # read_data = {}
    # for (n, p) in [("", {"thres": 0, "cats_int_fmt": None}),
    #              ("_cint", {"thres": 0, "cats_int_fmt": True}),
    #              ("_int", {"thres": 0, "cats_int_fmt": False}),
    #              ("_s", {"thres": 1, "cats_int_fmt": None}),
    #              ("_s-cint", {"thres": 1, "cats_int_fmt": True}),
    #              ("_s-int", {"thres": 1, "cats_int_fmt": False})]:
    #     data.writeCSV([rep+"out_LHS"+n+nan+".csv", rep+"out_RHS"+n+nan+".csv"], **p)
    #     read_data[n] = Data([rep+"out_LHS"+n+nan+".csv", rep+"out_RHS"+n+nan+".csv", {}, "nan"], "csv")
    #     print(n, "\n", read_data[n].disp())
    
    # pass
    # rep = "/home/egalbrun/short/cms/"
    # data = Data([rep+"toy500_LHS.csv", rep+"toy500_RHS.csv", {}, "nan"], "csv")

    # rep = "/home/egalbrun/TKTL/redescriptors/current/CSM_data/sample_2/"
    # data = Data([rep+"Patient_Claims_aggregated_Sample_2.csv", rep+"Patient_Claims_aggregated_Sample_2.csv", {"delimiter": ","}, "nan"], "csv")
    # print(data)

    # rep = "/home/egalbrun/TKTL/redescriptors/dp/test/dblp_pickedBB_empty/"
    # data = Data([rep+"data_LHS.csv", rep+"data_RHS.csv", {}, "nan"], "csv")

    # ADDING NEW VARIABLES
    # rep = "/home/egalbrun/TKTL/misc/ecometrics/china_paper/redescription_china/7_fossils/biotraits_focus+FOSS/"
    # data = Data([rep+"data_LHS.csv", rep+"data_RHS.csv", {}, "nan"], "csv")
    # vs = []
    # ### MAT (bio1) = 27 - 28.5 AL
    # vs.append(("modeled_bio1", 27 - 28.5*data.col(0, 3).getVector()))
    # # ### MINT (bio5) = 18 - 42.9 AL
    # vs.append(("modeled_bio5", 18 - 42.9*data.col(0, 3).getVector()))
    # # ### MAP (bio12) = 2491 - 289 HYP - 841 LOP
    # vs.append(("modeled_bio12", 2491 - 289.*data.col(0, 0).getVector()-841*data.col(0, 1).getVector()))
    # # ### NPP = 2601 - 144 HYP - 935 LOP
    # vs.append(("modeled_NPP", 2601 - 144.*data.col(0, 0).getVector()-935*data.col(0, 1).getVector()))
    # for (name, v) in vs:
    #     new_col = NumColM.parseList(v, force=True)
    #     new_col.setDisabled()
    #     new_col.setSideIdName(0, data.nbCols(0), name)
    #     data.appendCol(new_col, 0)
    # data.writeCSV([rep+"out_LHS.csv", rep+"out_RHS.csv"])
    # dataX = Data([rep+"out_LHS.csv", rep+"out_RHS.csv", {}, "nan"], "csv")
    # print(data)
    # print(dataX)
    # print(dataX.sameAs(data, side=0))
    # print(dataX.sameAs(data, side=1))

    # RANDOM SHUFFLING
    # rep = "/home/egalbrun/short/rpr/"
    # for (dL, dR) in [("LHS_50", "RHS_50"), ("RHS_50", "LHS_50")]: # ("THS", "RHS"),
    #     ## data = Data([rep+"EA_ethnoY.csv", rep+"EA_bioY.csv", {}, ""], "csv")
    #     # data = Data([rep+dt+"/data_LHSa.csv", rep+dt+"/data_RHSa.csv", {}, "NA"], "csv")
    #     data = Data([rep+"data_%s.csv" % dL, rep+"data_%s.csv" % dR, {"delimiter": ","}, "NA"], "csv")
    #     rdt, sd = data.shuffle_sides([0])
    #     pdb.set_trace()
    # print(data.cols[1][1].areDataEquiv(-21.88,-21.8))

    # ############# DATA EXTENSIONS
    # # rep_in = "/home/egalbrun/short/geo_small/"
    # # data_files = [rep_in+"ecoEu_teeth.csv", rep_in+"ecoEu_clim.csv"]
    # rep_in = "/home/egalbrun/x/"
    # data_files = [rep_in+"data_LHS.csv", rep_in+"data_RHS.csv"]
    # data = Data(data_files+[{}, ""], "csv")
    # data.initExtension("geoplus", params ={"wgrid_percentile": 90, "hgrid_percentile": 90})
    # print(data.computeExtras(data.col(0,0), extras=["cohesion"]))
    # data.getExtension("geoplus").prepAreasData(cells_colors=[0])
    # data.saveExtensions(details={"dir": rep_in})

    # # data = Data([rep_in+"ecoEu_teeth.csv", rep_in+"ecoEu_clim.csv", {}, ""], "csv")
    # # data.initExtension("geoplus")
    # # print(data.computeExtras(data.col(0,0), extras=["cohesion"]))
    # # data.saveExtensions(details={"dir": rep_in})

    # data2 = Data(data_files+[{}, ""], "csv")
    # data2.loadExtensions(ext_keys=["geoplus"], filenames={"extf_list_edges": rep_in+"list_edges.csv"})
    # print(data2.computeExtras(data2.col(0,0), extras=["cohesion"]))

    # OLD DATA FORMAT
    # data = Data(["/home/galbrun/dblp_data/filtered/conference_filtered.datnum",
    #              "/home/galbrun/dblp_data/filtered/coauthor_filtered.datnum",
    #              "/home/galbrun/dblp_data/filtered/conference_filtered.names",
    #              "/home/galbrun/dblp_data/filtered/coauthor_filtered.names",
    #              None,
    #              "/home/galbrun/dblp_data/filtered/coauthor_filtered.names"], "multiple")


if __name__ == '__main__':
    main()
