import os.path, time, re
import numpy

try:
    from classData import Data
    from classRedescription import Redescription
except ModuleNotFoundError:
    from .classData import Data
    from .classRedescription import Redescription

import pdb

def getPrec(counts):
    max_count = numpy.max(counts)
    dk = 1 / (max_count * (max_count-1.))
    prec_all = int(numpy.ceil(-numpy.log10(dk)))
    if prec_all < 0 or prec_all > 6:
        prec_all = None       
    return prec_all

def dotProduct(X, Y, Trids, sids=None):        
    Ysub = Y[:,Trids]
    Xsub  = X
    counts_sub = numpy.sum(X, axis=0)
    if sids is not None:
        Xsub, counts_sub = Xsub[:, sids], counts_sub[sids]
    counts_sub[counts_sub==0] = -1. ### will average to zero, just making sure we are not dividing by zero        
    return numpy.dot(Ysub, Xsub)/numpy.tile(counts_sub, (Ysub.shape[0], 1))

def swap_binary(mat, attempts=None):
    swapped = 0
    if attempts is None:
        attempts = numpy.sum(mat.shape)
    rows, cols = mat.nonzero()
    intra_maxrepeat = 10
    for i in range(attempts):
        ### update the indices
        posA = numpy.random.randint(rows.shape[0])
        posB = numpy.random.randint(rows.shape[0])
        nba = 0
        while (cols[posA] == cols[posB] or rows[posA] == rows[posB]) and nba < intra_maxrepeat:
            nba += 1
            posA = numpy.random.randint(rows.shape[0])
            posB = numpy.random.randint(rows.shape[0])
        if cols[posA] != cols[posB] and rows[posA] != rows[posB] and mat[rows[posB], cols[posA]] == 0 and mat[rows[posA], cols[posB]] == 0:
            swapped += 1
            mat[rows[posA], cols[posA]], mat[rows[posB], cols[posB]] = (0, 0)
            mat[rows[posB], cols[posA]], mat[rows[posA], cols[posB]] = (1, 1)
            # rows[posA], rows[posB], cols[posA], cols[posB] = rows[posA], rows[posB], cols[posA], cols[posB]
            cols[posA], cols[posB] = cols[posB], cols[posA]
    # print("SWAPPED", swapped, attempts)
    return mat


def prepareMatrices(data, dT):
    Tids = dT.getRNames()
    map_ids = dict([(v,k) for (k,v) in enumerate(Tids)])
    map_crows = []
    for i, ns in enumerate(data.getNames()):
        for ni, n in enumerate(ns):
            if n in map_ids:
                map_crows.append(((i, ni), map_ids[n]))
    dtd_sidecols, Trids = zip(*map_crows)
    X, xdets, xmcols = data.getMatrix(side_cols=dtd_sidecols)
    Y, ydets, ymcols = dT.getMatrix(side_cols=[(0,None)])
    return {"X": X, "xdets": xdets, "xmcols": xmcols,
            "Y": Y, "ydets": ydets, "ymcols": ymcols,
            "Trids": Trids}

def prepareCounts(X, data=None, count_vname="COUNTS", side=0):
    counts = numpy.sum(X, axis=0)
    if data is not None and len(data.getColsByName("^%s$" % count_vname)) == 0:
        col = prepareColFromVector(counts, prec=0)
        col.setSideIdName(side, data.nbCols(side), count_vname)        
        data.appendCol(col, side)
    return counts

def selectRids(data, select_red, select_union=False):
    sids = None
    if select_red is not None:
        r = Redescription.parse(select_red, data=data)
        if select_union:
            sids = sorted(r.getSuppU())
        else:
            sids = sorted(r.getSuppI())
        if len(sids) == 0:
            raise Warning("Selection empty!")
    return sids        

def prepare(data, dT, count_vname="COUNTS", vname_patt="MEAN_%s", prec_all=-1):
    store = prepareMatrices(data, dT)
    store["counts"] = prepareCounts(store["X"], data, count_vname, side=0)
    store["vnames"] = [vname_patt % s for s in dT.getNames(0)]
    if prec_all == -1:
        prec_all = getPrec(counts)
    store["prec"] = prec_all
    return store

def prepareAggData(data, dT, select_red=None, select_union=False, count_vname="COUNTS", vname_patt="MEAN_%s", prec_all=-1):
    store = prepare(data, dT, count_vname, vname_patt, prec_all)
    sids = selectRids(data, select_red, select_union)
        
    Z = dotProduct(store["X"], store["Y"], store["Trids"], sids)
    Dsub = data.subset(row_ids=sids)
    back = Dsub.replaceSideFromMatrix(Z, prec=store["prec"], vnames=store["vnames"], side=0)    
    counts_sub = store["counts"][sids] if sids is not None else store["counts"]
    col = Dsub.prepareColFromVector(counts_sub, prec=0, enabled=False)
    col.setSideIdName(0, Dsub.nbCols(0), count_vname)
    Dsub.appendCol(col, 0)
    return Dsub, sids, back, store

def prepareRndAggData(data, dT, rnd_meth="permute_LHS", select_red=None, select_union=False, count_vname="COUNTS", vname_patt="MEAN_%s", prec_all=-1):
    store = prepare(data, dT, count_vname, vname_patt, prec_all)
    sids = selectRids(data, select_red, select_union)
    Lsids = sids
    sdt = None
    X, Y, Trids = (store["X"], store["Y"], store["Trids"])

    occs_changed = False
    if "permute_traits" == rnd_meth:
        cids = numpy.random.permutation(Y.shape[1])
        Trids = cids[:len(Trids)]        
        store["random"] = {"meth": rnd_meth, "permute": cids, "nb": len(Trids)}
    elif re.match("permute_", rnd_meth):
        if sids is None:
            permuted_sids = range(data.nbRows())
        else:
            permuted_sids = list(sids)
        numpy.random.shuffle(permuted_sids)
        if re.search("permute_RHS", rnd_meth):
            sdt = {1: permuted_sids, 0: None}
        else:
            Lsids = permuted_sids
            sdt = {0: permuted_sids, 1: None}
        store["random"] = {"meth": rnd_meth, "permute": permuted_sids}
    elif "shuffle_traits" == rnd_meth:
        Y = Y.copy()
        rnds = []
        for j in range(Y.shape[0]):
            x = numpy.random.permutation(Y.shape[1])
            Y[j,:] = Y[j,x]
            rnds.append(x)
        store["random"] = {"meth": rnd_meth, "permute": rnds}
            
    elif "swaprnd_occs" == rnd_meth:
        occs_changed = True
        Xt = swap_binary(X.copy())
        store["random"] = {"meth": rnd_meth, "where": numpy.where(X!=Xt)}
        # print("Col-dist", numpy.sum((numpy.sum(X, axis=0) - numpy.sum(Xt, axis=0)) !=0))
        # print("Row-dist", numpy.sum((numpy.sum(X, axis=1) - numpy.sum(Xt, axis=1)) !=0))
        # print("Nb diffs", len(numpy.where(X!=Xt)[0]), len(numpy.where(X!=Xt)[0])/4.)
        X = Xt        

    Z = dotProduct(X, Y, Trids, Lsids)
    Dsub = data.subset(sids, sdt)
    counts_sub = store["counts"][Lsids] if sids is not None else store["counts"]    

    if occs_changed:
        Xp = numpy.append(X[:, Lsids], [counts_sub], axis=0)
        org = Dsub.replaceSideFromMatrix(Xp, side=0, vtypes={None: 1, X.shape[0]: 3})
        
    back = Dsub.replaceSideFromMatrix(Z, prec=store["prec"], vnames=store["vnames"], side=0)    
    col = Dsub.prepareColFromVector(counts_sub, prec=0, enabled=False)
    col.setSideIdName(0, Dsub.nbCols(0), count_vname)
    Dsub.appendCol(col, 0)
    return Dsub, sids, back, store


def prepareSubData(data, select_red=None, select_union=False):
    sids = selectRids(data, select_red, select_union)        
    Dsub = data.subset(row_ids=sids)
    return Dsub, sids

def prepareRndSubData(data, rnd_meth="-", select_red=None, select_union=False):
    sids = selectRids(data, select_red, select_union)
    shuffled = sids
    sdt = None
    if re.match("permute_", rnd_meth):
        if sids is None:
            sids = range(data.nbRows())
        shuffled = list(sids)
        numpy.random.shuffle(shuffled)
        sdt = {0: shuffled, 1: None}
        if re.search("permute_RHS", rnd_meth):
            sdt = {1: shuffled, 0: None}
    Dsub = data.subset(sids, sdt)
    return Dsub, shuffled


class RndFactory(object):


    def __init__(self, LHS_filename=None, RHS_filename=None, csv_params=None, unknown_string=None, traits_filename=None, org_data=None, traits_data=None):
        self.org_data = None
        self.traits_data = None
        self.store = {}
        if org_data is not None:
            self.setSides(org_data)
        elif LHS_filename is not None and RHS_filename is not None:
            self.readSides(LHS_filename, RHS_filename, csv_params, unknown_string)
        if traits_data is not None:
            self.setTraits(traits_data)
        elif traits_filename is not None:
            self.readTraits(traits_filename, csv_params, unknown_string)

    def setSides(self, data):
        self.org_data = data
            
    def setTraits(self, data):
        self.store = {}
        self.traits_data = data
            
    def readSides(self, LHS_filename, RHS_filename, csv_params=None, unknown_string=None):
        data = Data([LHS_filename, RHS_filename, csv_params, unknown_string], "csv")
        self.setSides(data)
            
    def readTraits(self, traits_filename, csv_params=None, unknown_string=None):
        traits_data = Data([traits_filename, None, csv_params, unknown_string], "csv")
        self.setTraits(traits_data)

    def setSeed(self, v):
        numpy.random.seed(v)
        
    def makeupData(self, with_traits=False, count_vname="COUNTS", vname_patt="MEAN_%s", prec_all=-1, select_red=None, select_union=False):
        if self.traits_data is not None and with_traits:
            Dsub, sids, back, store = prepareAggData(self.org_data, self.traits_data, select_red, select_union, count_vname, vname_patt, prec_all)
        else:
            Dsub, sids = prepareSubData(self.org_data, select_red, select_union)
            back = None
            store = {}
        return Dsub, sids, back, store

    def makeupRndData(self, rnd_meth="none", with_traits=False, count_vname="COUNTS", vname_patt="MEAN_%s", prec_all=-1, select_red=None, select_union=False):
        if rnd_meth == "none":
            return self.makeupData(with_traits, count_vname, vname_patt, prec_all, select_red, select_union)
        
        if self.traits_data is not None and with_traits:
            Dsub, sids, back, store = prepareRndAggData(self.org_data, self.traits_data, rnd_meth, select_red, select_union, count_vname, vname_patt, prec_all)
        else:
            Dsub, sids = prepareRndSubData(self.org_data, rnd_meth, select_red, select_union)
            back = None
            store = {}
        return Dsub, sids, back, store

            
#####################################################
#####################################################

def main():

    # select_red = "[NB_SPC > 2]\t[Z:AsiaFocus=F:in]"
    select_red = "[NB_SPC > 5]\t[bio1:TMeanY<25.0]"

    numpy.random.seed(1)
    csv_params = {}
    unknown_string = "nan"
    rep = "/home/egalbrun/TKTL/misc/ecometrics/redescription_teeth_continents/prepared_vX/"
    rep_out = "/home/egalbrun/Desktop/"

    rf = RndFactory(rep+"occurence_IUCN_all.csv", rep+"bio_IUCN_select.csv", csv_params, unknown_string)
    # Dsub, sids, _, _ = rf.makeupRndData(rnd_meth="permute_LHS") #, select_red="[]\t[Z:AsiaFocus=F:in]", select_union=True) #, prec_all=3)
    # Dsub.writeCSV([rep+"data_LHSv2yr1.csv", rep+"data_RHSv2yr1.csv"])
    
    rf.readTraits(rep+"traits_IUCN_all.csv", csv_params, unknown_string)
    # Dsub, sids, back, store = rf.makeupData(with_traits=True, count_vname="NB_SPC", select_red=select_red, prec_all=3)
    # print(Dsub)
    # suff = "orgX"
    # Dsub.writeCSV([rep_out+"data_LHSvT_"+suff+".csv", rep_out+"data_RHSvT_"+suff+".csv"])

    for rnd_meth in ["permute_RHS", "permute_LHS", "permute_traits", "shuffle_traits", "swaprnd_occs"]:
        for i in range(10):
            Dsub, sids, back, store = rf.makeupRndData(rnd_meth=rnd_meth, with_traits=True, count_vname="NB_SPC", select_red=select_red, prec_all=3)            
            
            print(Dsub)
            suff = "%s-%dY" % (rnd_meth, i)
            Dsub.writeCSV([rep_out+"data_LHSvT_"+suff+".csv", rep_out+"data_RHSvT_"+suff+".csv"])
            
if __name__ == '__main__':    
    main()

