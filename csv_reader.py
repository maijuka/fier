import sys
import re
import io
import numpy
import zipfile
import csv

try:
    from classQuery import ENTITY_MARK, Term, TimeTools, NA_str_c
    from toolICDict import first_value
except ModuleNotFoundError:
    from .classQuery import ENTITY_MARK, Term, TimeTools, NA_str_c
    from .toolICDict import first_value

import pdb

LATITUDE = ('lat', 'latitude', 'Lat', 'Latitude', 'lats', 'latitudes', 'Lats', 'Latitudes')
LONGITUDE = ('long', 'longitude', 'Long', 'Longitude', 'longs', 'longitudes', 'Longs', 'Longitudes')
IDENTIFIERS = ('id', 'identifier', 'Id', 'Identifier', 'ids', 'identifiers', 'Ids', 'Identifiers', 'ID', 'IDS')
COND_TIME = 'cond_time'
COND_COL = ('cond_var', 'cond_col', 'timeid', COND_TIME)

ENABLED_ROWS = ('enabled_row', 'enabled_rows')
ENABLED_COLS = ('enabled_col', 'enabled_cols')
GROUPS_COLS = ('groups_col', 'groups_cols')

COLVAR = ['cid', 'CID', 'cids', 'CIDS', 'variable', 'Variable', 'variables', 'Variables']
COLVAL = ['value', 'Value', 'values', 'Values']

FOLDS_PREF = "folds_subs_"
FOLDS_PATT = "^(v\d+:)?"+FOLDS_PREF


class CSVRError(Exception):
    def __init__(self, value):
        self.value = value

    def __str__(self):
        return repr(self.value)


def test_some_numbers(strgs):
    i = len(strgs)
    while i > 0:
        i -= 1
        try:
            float(strgs[i])
            return True
        except:
            pass
    return False


def test_all_numbers(strgs):
    i = len(strgs)
    while i > 0:
        i -= 1
        try:
            float(strgs[i])
        except:
            return False
    return True


class RowOrder(object):
    def __init__(self, LhasIds=False, RhasIds=False, Lids=None, Rids=None,
                 LhasTimes=False, RhasTimes=False, Ltimes=None, Rtimes=None,
                 LhasCoords=False, RhasCoords=False, Lcoords=None, Rcoords=None):        
        self.setKeys(LhasIds, RhasIds, Lids, Rids,
                     LhasTimes, RhasTimes, Ltimes, Rtimes,
                     LhasCoords, RhasCoords, Lcoords, Rcoords)
        self.orders = self.computeOrders()

    def resetKeys(self):
        self.orders = None
        self.Ns = [-1, -1]
        self.key_ids = None
        self.key_times = None
        self.key_coords = None

    def setKeyIds(self, Lids, Rids):
        self.key_ids = [Lids, Rids]

    def setKeyTimes(self, Ltimes, Rtimes):
        self.key_times = [Ltimes, Rtimes]

    def setKeyCoords(self, Lcoords, Rcoords):
        if len(Lcoords) == 2:
            Lkcoords = list(zip(*Lcoords))
        else:
            Lkcoords = Lcoords
        if len(Rcoords) == 2:
            Rkcoords = list(zip(*Rcoords))
        else:
            Rkcoords = Rcoords
        self.key_coords = [Lkcoords, Rkcoords]

    def hasKeyIds(self):
        return self.key_ids is not None

    def hasKeyTimes(self):
        return self.key_times is not None

    def hasKeyCoords(self):
        return self.key_coords is not None

    def getKeyId(self, side, row):
        if self.key_ids is not None and row < len(self.key_ids[side]):
            return [self.key_ids[side][row]]
        else:
            return []

    def getKeyTime(self, side, row):
        if self.key_times is not None and row < len(self.key_times[side]):
            return [self.key_times[side][row]]
        else:
            return []

    def getKeyCoord(self, side, row):
        if self.key_coords is not None and row < len(self.key_coords[side]):
            if len(self.key_coords[side][row]) == 2:
                return ["LAT"] + self.key_coords[side][row][0] + ["LONG"] +self.key_coords[side][row][1]
            return [self.key_coords[side][row]]
        else:
            return []

    def getKeyTuple(self, side, row):
        return tuple(self.getKeyId(side, row) + self.getKeyTime(side, row) + self.getKeyCoord(side, row))

    def getKeyStr(self, side, row):
        return "::".join(["%s" % v for v in self.getKeyId(side, row) + self.getKeyTime(side, row) + self.getKeyCoord(side, row)])

    def getAllKeysTuple(self, side):
        if self.Ns[side] > 0:
            return [self.getKeyTuple(side, row) for row in range(self.Ns[side])]
        return []

    def getAllKeysStr(self, side):
        if self.Ns[side] > 0:
            return [self.getKeyStr(side, row) for row in range(self.Ns[side])]
        return []

    def setKeys(self, LhasIds=False, RhasIds=False, Lids=None, Rids=None,
                LhasTimes=False, RhasTimes=False, Ltimes=None, Rtimes=None,
                LhasCoords=False, RhasCoords=False, Lcoords=None, Rcoords=None):
        self.resetKeys()
        if (LhasIds and RhasIds):
            self.setKeyIds(Lids, Rids)
        if (LhasTimes and RhasTimes):
            self.setKeyTimes(Ltimes, Rtimes)
        if (LhasCoords and RhasCoords):
            self.setKeyCoords(Lcoords, Rcoords)

    def validLen(self):
        Ns = [-1, -1]
        if self.hasKeyIds():
            for side in [0, 1]:
                Ns[side] = len(self.key_ids[side])
        if self.hasKeyTimes():
            for side in [0, 1]:
                if Ns[side] == -1:
                    Ns[side] = len(self.key_times[iside])
                elif Ns[side] != len(self.key_times[side]):
                    raise CSVRError('Incoherent length of keys! %d vs. %d (times on side %d)' % (Ns[side], len(self.key_times[side]), side))
        if self.hasKeyCoords():
            for side in [0, 1]:
                if Ns[side] == -1:
                    Ns[side] = len(self.key_coords[iside])
                elif Ns[side] != len(self.key_coords[side]):
                    raise CSVRError('Incoherent length of keys! %d vs. %d (coords on side %d)' % (Ns[side], len(self.key_coords[side]), side))
        if Ns[0] > 0 and Ns[1] > 0:
            self.Ns = Ns
            return True
        return False

    def hasOrders(self):
        return self.orders is not None

    def getOrders(self):
        return self.orders

    def computeOrders(self):
        if not self.validLen():
            return None

        # formatL = "::".join(["%s" for i in range(len(order_keys[0]))])
        # formatR = "::".join(["%s" for i in range(len(order_keys[1]))])
        # xLll = [formatL % p for p in zip(*order_keys[0])]
        # xRll = [formatR % p for p in zip(*order_keys[1])]
        Lll = self.getAllKeysStr(0)
        Rll = self.getAllKeysStr(1)
        # Rll = ["::".join(map(str, p)) for p in zip(*order_keys[1])]
        if len(set(Lll)) < len(Lll) or len(set(Rll)) < len(Rll):
            print('Those ids are no real ids, they are not unique!..')

        Lorder = sorted(range(len(Lll)), key=Lll.__getitem__)
        Rorder = sorted(range(len(Rll)), key=Rll.__getitem__)
        both = set(Lll).intersection(Rll)
        if len(both) == 0:
            raise CSVRError('Error while parsing the data, found no matching rows!')

        # Remove from Lorder and Rorder the parts that aren't in both
        i = 0
        while i < len(Lorder):
            if Lll[Lorder[i]] not in both:
                del Lorder[i]
            else:
                i += 1
        i = 0
        while i < len(Rorder):
            if Rll[Rorder[i]] not in both:
                del Rorder[i]
            else:
                i += 1
        return Lorder, Rorder


def start_out(fp):
    return csv.writer(fp, quoting=csv.QUOTE_MINIMAL)  # , delimiter=';', quotechar='"'
    # return csv.writer(fp, quoting=csv.QUOTE_NONNUMERIC) #, delimiter=';', quotechar='"'


def write_row(csvf, row_data):
    csvf.writerow(row_data)


def getFp(filename, write=False):
    fcl = False
    f = None
    if type(filename) is str:  # or type(filename) is unicode:
        if write:
            f = open(filename, 'w')
        else:
            f = open(filename, 'r')
        fcl = True
    elif isinstance(filename, io.TextIOBase) or isinstance(filename, zipfile.ZipExtFile):
        f = filename
    # else:
    #     pdb.set_trace()
    #     ### Because ZIPext files don't have a seek method...
    #     f = io.StringIO(filename.read())
    return f, fcl


def read_coords_csv(filename, csv_params={}, unknown_string=None):
    coords = []
    rnames = None
    f, fcl = getFp(filename)
    if f is not None:
        try:
            dialect = csv.Sniffer().sniff(f.read(2048))
        except Exception:
            dialect = "excel"
        f.seek(0)
        # header = csv.Sniffer().has_header(f.read(2048))
        # f.seek(0)
        csvreader = csv.reader(f, dialect=dialect, **csv_params)
        # Try to read headers
        head = [h for h in next(csvreader)]

        cpos = {}
        for i, h in enumerate(head):
            for clbls in [LATITUDE, LONGITUDE, IDENTIFIERS]:
                if h in clbls:
                    cpos[clbls[0]] = i

        if not (LATITUDE[0] in cpos and LONGITUDE[0] in cpos):
            return None, None
        cmax = max(cpos.values())
        if IDENTIFIERS[0] in cpos:
            rnames = []

        for row in csvreader:
            if re.match("\s*#", row[0]) or row[0] in ENABLED_ROWS+ENABLED_COLS+GROUPS_COLS:
                continue
            if len(row) < cmax+1:
                raise ValueError('number of columns does not match (is ' +
                                 str(len(row))+', should be at least' +
                                 str(cmax+1)+')')
            coords.append((float(row[cpos[LONGITUDE[0]]].strip()), float(row[cpos[LATITUDE[0]]].strip())))
            # coords.append((float(row[cpos[LATITUDE[0]]].strip()), float(row[cpos[LONGITUDE[0]]].strip())))
            if rnames is not None:
                tmp = row[cpos[IDENTIFIERS[0]]].strip()
                if tmp != type(tmp)(unknown_string):
                    rnames.append(tmp)
                else:
                    rnames.append(None)
    if fcl:
        f.close()
    # HERE DEBUG UTF-8
    return coords, rnames


def write_coords_csv(filename, coords, rnames=None, csv_params={}, unknown_string=None):
    f, fcl = getFp(filename, write=True)
    if f is not None:
        head = [LONGITUDE[0], LATITUDE[0]]
        if rnames is not None:
            head.append(IDENTIFIERS[0])
        f.write(",".join(head)+"\n")
        for i, coord in enumerate(coords):
            tow = ["%f" % c for c in coord]
            if rnames is not None:
                tow.append("%s" % rnames[i])
            f.write(",".join(tow)+"\n")
    if fcl:
        f.close()


def read_csv(filename, csv_params={}, unknown_string=None):
    type_all = None
    skipfirst = False
    f, fcl = getFp(filename)
    if f is not None:
        try:
            dialect = csv.Sniffer().sniff(f.read(2048))
        except Exception:
            dialect = "excel"
        f.seek(0)
        # header = csv.Sniffer().has_header(f.read(2048))
        # f.seek(0)
        csvreader = csv.reader(f, dialect=dialect, **csv_params)
        # Try to read headers
        head = [h for h in next(csvreader)]

        tmp = re.search('#*\s*"?type=(?P<type>\w)"?', head[-1])
        if tmp is not None:
            tt = head.pop().split("type=")[0].strip("# \"'")
            if len(tt) > 0:
                head.append(tt)
            type_all = tmp.group('type')
            if len(head) == 0:
                head = [h for h in next(csvreader)]
                skipfirst = True
        if test_some_numbers(head):
            # If we read a row with some numerical values, this was no header...
            head = [Term.pattVName % i for i in range(len(head))]
            data = dict(zip(head, [[] for i in range(len(head))]))
            f.seek(0)
        else:
            tmp = re.match('type=(?P<type>\w)', head[-1])
            if tmp is not None:
                head.pop()
                type_all = tmp.group('type')
            data = dict([(head[i].strip(), []) for i in range(len(head))])
            if len(data) != len(head):
                map_names = {}
                for i in range(len(head)):
                    if head[i] in map_names:
                        map_names[head[i]].append(i)
                        head[i] += "(duplicate#%d)" % (len(map_names[head[i]]))
                    else:
                        map_names[head[i]] = [i]
                data = dict([(head[i].strip(), []) for i in range(len(head))])
                # raise ValueError('some columns have the same name, this is a very bad idea...')
        no_of_columns = len(head)
        specials = {}
        for row in csvreader:
            if skipfirst:
                skipfirst = False
                continue
            if re.match("\s*#", row[0]):
                continue
            if len(row) == no_of_columns+1 and ((row[0] in GROUPS_COLS) or (row[0] in ENABLED_COLS)):
                # these would be for data not including ids
                specials[row.pop(0)] = len(first_value(data))
            if len(row) != no_of_columns:
                raise ValueError('number of columns does not match (is ' +
                                 str(len(row))+', should be ' +
                                 str(no_of_columns)+')')
            for i in range(len(row)):
                tmp = row[i].strip()
                if tmp != type(tmp)(unknown_string):
                    data[head[i]].append(tmp)
                else:
                    # print("Turned to None (in csv_reader)", tmp)
                    data[head[i]].append(None)
    if fcl:
        f.close()
    return head, data, type_all, specials


def parse_sparse(D, coord, ids, varcol, valcol, off=None):
    nids = None
    if varcol is None:
        # this is no sparse data...
        return D, coord, ids, False, False
    else:
        nll = set(ids)
        dictLL = {}
        special_rows = {}
        numerical_ids = True
        row_named = False
        sub = 1

        if valcol is not None:
            for i in nll:  # check for special rows, negative or token identifiers
                try:
                    if (i in IDENTIFIERS) or (i in ENABLED_COLS) or (i in GROUPS_COLS) or int(i) < 0:
                        special_rows[i] = {}
                except ValueError as ve:
                    pass
            nll.difference_update(special_rows.keys())

        if off is not None:
            if off+1 == len(nll):
                sub = 0
            else:
                sub = off

        for i in nll:
            # in general if numerical ids are provided that should be the row number
            # we expect rows to start at one ...
            try:
                adj_idx = int(i)-sub
            except ValueError as ve:
                # a row identifier cannot be parsed to int
                numerical_ids = False
                break
            if adj_idx < 0:
                # a row identifier is offset to negative
                numerical_ids = False
                break
            else:
                dictLL[i] = adj_idx

        if numerical_ids:
            if off is None and (-sub in dictLL.values()):
                # ... unless there was a zero
                dictLL = dict([(k, v+sub) for (k, v) in dictLL.items()])
            # if max(dictLL.values()) > 2*len(dictLL):
            if max(dictLL.values()) > 10*len(dictLL):
                print("Too large ids compared to number of rows (>10x)!...")
                numerical_ids = False
        if numerical_ids:
            maxn = max(dictLL.values())
            pattEName = ENTITY_MARK+"%0"+str(len("%d" % (maxn+1)))+"d"
            nll = [pattEName % (v+1) for v in range(maxn+1)]
            nids = nll
        if not numerical_ids:
            # if row identifiers are not numerical
            row_named = True
            special_rows = {}
            dictLL = dict([(v, k) for (k, v) in enumerate(sorted(nll))])
            if ids is not None:
                nids = [None for i in range(len(nll))]
                for ii, i in enumerate(ids):
                    nids[dictLL[ids[ii]]] = i

    nD = {'data': {}, 'headers': [], "sparse": True, "bool": False, ENABLED_COLS[0]: {}, GROUPS_COLS[0]: {}}
    if valcol is None:
        nD['bool'] = True
        # Turning the data from list of ids to sets, Boolean
        for rid, col in enumerate(D['data'][varcol]):
            if col in nD['headers']:
                nD['data'][col].add(dictLL[ids[rid]])
            else:
                nD['headers'].append(col)
                nD['data'][col] = set([dictLL[ids[rid]]])
    else:
        # Turning the data from list to dict row:value
        for rid, col in enumerate(D['data'][varcol]):
            if ids[rid] in special_rows:
                special_rows[ids[rid]][col] = D['data'][valcol][rid]
            elif col in nD['headers']:
                nD['data'][col][dictLL[ids[rid]]] = D['data'][valcol][rid]
            else:
                nD['headers'].append(col)
                nD['data'][col] = {dictLL[ids[rid]]: D['data'][valcol][rid]}

    # Retrieving names if any (column/row with index -1)
    # if contains_row_names:
    if "-1" in nD["headers"]:
        row_named = True
        nD["headers"].remove("-1")
        nids = [None for i in range(len(nll))]
        for ii, i in nD["data"].pop("-1").items():
            nids[ii] = i

    special_rids = list(special_rows.keys())
    for k in special_rids:
        if "-1" in special_rows[k]:
            special_rname = special_rows[k].pop("-1")
            if special_rname in special_rows:
                special_rows[k].update(special_rows[special_rname])
            special_rows[special_rname] = special_rows.pop(k)
            k = special_rname

    col_names = {}
    col_enabled = {}
    col_group = {}
    special_rids = list(special_rows.keys())
    for k in special_rids:
        if (k in IDENTIFIERS) or k == "-1":
            special_rows[k].update(col_names)
            col_names = special_rows.pop(k)
        elif k in ENABLED_COLS:
            special_rows[k].update(col_enabled)
            col_enabled = special_rows.pop(k)
        elif k in GROUPS_COLS:
            special_rows[k].update(col_group)
            col_group = special_rows.pop(k)

    if len(col_names) > 0:
        new_headers = []
        keysc = sorted(col_names.keys(), key=lambda x: int(x))

        col_enabled = dict([(col_names.get(str(c), c), v) for (c, v) in col_enabled.items()])
        col_group = dict([(col_names.get(str(c), c), v) for (c, v) in col_group.items()])

        for ko in keysc:
            cn = col_names[ko]
            if ko in nD["headers"]:
                nD["headers"].remove(ko)
                nD["data"][cn] = nD["data"].pop(ko)
            else:
                nD["data"][cn] = dict()
            new_headers.append(cn)
        nD["headers"] = new_headers + nD["headers"]

    nD[ENABLED_COLS[0]] = col_enabled
    nD[GROUPS_COLS[0]] = col_group
    # mapping the coordinates to the correct order
    ncoord = None
    hasc, coord_sp = has_coord_sparse(nD)
    if coord is not None and coord != (None, None):
        ncoord = [[None for i in range(len(dictLL))], [None for i in range(len(dictLL))]]
        for ii in range(len(coord[0])):
            ncoord[0][dictLL[ids[ii]]] = coord[0][ii]
            ncoord[1][dictLL[ids[ii]]] = coord[1][ii]
            if hasc and ((coord_sp[0].get(dictLL[ids[ii]]) != ncoord[0][dictLL[ids[ii]]] or
                          coord_sp[1].get(dictLL[ids[ii]]) != ncoord[1][dictLL[ids[ii]]])):
                raise CSVRError('Found incoherent coordinates! #%d sparse (%s, %s) vs. previous (%s, %s)'
                                % (dictLL[ids[ii]],
                                    coord_sp[0].get(dictLL[ids[ii]]), coord_sp[1].get(dictLL[ids[ii]]),
                                    ncoord[0][dictLL[ids[ii]]], ncoord[1][dictLL[ids[ii]]]))
    elif hasc:
        ncoord = [[coord_sp[0].get(i, None) for i in range(len(dictLL))], [coord_sp[1].get(i, None) for i in range(len(dictLL))]]

    return nD, ncoord, nids, hasc, row_named


def has_disabled_rows(D):
    # This is taken care of in the Data class
    hasDis = False
    dis = []
    for s in ENABLED_ROWS:
        if s in D['headers']:
            hasDis = True
            dis = [p for p in D['data'][s]]
            del D['data'][s]
            D['headers'].remove(s)
            break
    return (hasDis, dis)


def has_coord(D):
    hasCoord = False
    coord = (None, None)
    for s in LATITUDE:
        if s in D['headers']:
            for t in LONGITUDE:
                if t in D['headers']:
                    hasCoord = True
                    coord = ([list(map(float, (p or "-361").strip(" :").split(":"))) for p in D['data'][s]],
                             [list(map(float, (p or "-361").strip(" :").split(":"))) for p in D['data'][t]])
                    del D['data'][s]
                    del D['data'][t]
                    D['headers'].remove(s)
                    D['headers'].remove(t)
                    break
        if hasCoord:
            break

    return (hasCoord, coord)


def has_coord_sparse(D):
    hasCoord = False
    coord = (None, None)
    for s in LATITUDE:
        if s in D['headers']:
            for t in LONGITUDE:
                if t in D['headers']:
                    hasCoord = True
                    coord = (dict([(k, list(map(float, (v or "-361").strip(" :").split(":")))) for (k, v) in D['data'][s].items()]),
                             dict([(k, list(map(float, (v or "-361").strip(" :").split(":")))) for (k, v) in D['data'][t].items()]))
                    del D['data'][s]
                    del D['data'][t]
                    D['headers'].remove(s)
                    D['headers'].remove(t)
                    break
        if hasCoord:
            break
    return (hasCoord, coord)


def has_ids(D):
    hasIds = False
    ids = None
    for s in IDENTIFIERS:
        if s in D['headers']:
            if not hasIds:
                hasIds = True
                ids = D['data'][s]
            del D['data'][s]
            D['headers'].remove(s)
            break
    return (hasIds, ids)


def has_condition(D, other_params={}):
    date_fmt = {}
    for t, f in [("time_yearfirst", "yearfirst"), ("time_dayfirst", "dayfirst")]:
        if other_params.get(t) is not None:
            date_fmt[f] = other_params.get(t)
    col = None
    condIds = None
    isTime = False
    for s in COND_COL:
        if s in D['headers']:
            if col is None:
                col = D['data'][s]
                if len(set(col)):
                    condIds = numpy.copy(col)
                if re.search("time", s):
                    try:
                        col = ["%d" % TimeTools.parse_time(v, **date_fmt) for v in col]
                        isTime = True
                    except (TypeError, ValueError):
                        col = D['data'][s]
            del D['data'][s]
            D['headers'].remove(s)
            break
    return (col, condIds, isTime)


def get_discol(D, ids):
    discol = {}
    for s in ENABLED_COLS:
        if s in ids:
            if type(ids) is dict:
                p = ids.pop(s)
            else:
                p = ids.index(s)
                ids.remove(s)
            for c in D['headers']:
                discol[c] = D['data'][c].pop(p)
            break
    return discol


def get_groupcol(D, ids):
    groupcol = {}
    for s in GROUPS_COLS:
        if s in ids:
            if type(ids) is dict:
                p = ids.pop(s)
            else:
                p = ids.index(s)
                ids.remove(s)
            for c in D['headers']:
                groupcol[c] = D['data'][c].pop(p)
            break
    return groupcol


def is_sparse(D):
    colid = list(COLVAR)
    colv = list(COLVAL)
    varcol = None
    valcol = None
    while len(colid) > 0:
        s = colid.pop(0)
        if s in D['headers']:
            varcol = s
            colid = []
    while len(colv) > 0:
        s = colv.pop(0)
        if s in D['headers']:
            valcol = s
            colv = []
    return varcol, valcol


def get_size_hint(L, Lids, Lfull, R, Rids, Rfull):
    rtn = None
    if Lfull:
        nb_rowsL = len(first_value(L["data"]))
        if Rfull:
            nb_rowsR = len(first_value(R["data"]))
            if nb_rowsL == nb_rowsR:
                rtn = (nb_rowsL, 0, 0)
            else:
                rtn = (None, nb_rowsL, nb_rowsR)
        else:
            tt = list(map(int, set(Rids) - set(["-1"])))
            nrMinR, nrMaxR = min(tt), max(tt)
            if nb_rowsL == nrMaxR+1:
                rtn = (nb_rowsL, 0, 0)
            elif nb_rowsL == nrMaxR and nrMinR > 0:
                rtn = (nb_rowsL, 0, 1)
            else:
                rtn = (None, nb_rowsL, nrMaxR)

    else:
        tt = []
        for t in set(Lids):
            if t != "-1":
                try:
                    tt.append(int(t))
                except ValueError:
                    pass
        nrMinL, nrMaxL = min(tt), max(tt)

        if Rfull:
            nb_rowsR = len(first_value(R["data"]))
            if nb_rowsR == nrMaxL+1:
                rtn = (nb_rowsR, 0, 0)
            elif nb_rowsR == nrMaxL and nrMinL > 0:
                rtn = (nb_rowsR, 1, 0)
            else:
                rtn = (None, nrMaxL, nb_rowsR)

        else:
            tt = []
            for t in set(Rids):
                if t != "-1":
                    try:
                        tt.append(int(t))
                    except ValueError:
                        pass

            nrMinR, nrMaxR = min(tt), max(tt)
            if nrMaxR == nrMaxL:
                if nrMinR == 1 and nrMinL == 1:
                    rtn = (nrMaxR, 1, 1)
                else:
                    rtn = (nrMaxR+1, 0, 0)
            elif nrMaxR == nrMaxL+1 and nrMinR == 1:
                rtn = (nrMaxR, 0, 1)
            elif nrMaxR+1 == nrMaxL and nrMinL == 1:
                rtn = (nrMaxL, 1, 0)
            else:
                rtn = (None, nrMaxL, nrMaxR)
    return rtn


def row_order(L, R, other_params={}):
    # TODO catch the dense row containing info on enabled columns
    (LhasIds, Lids) = has_ids(L)
    (RhasIds, Rids) = has_ids(R)
    (Lvarcol, Lvalcol) = is_sparse(L)
    (Rvarcol, Rvalcol) = is_sparse(R)

    if Lvarcol is None:
        if LhasIds:
            L[ENABLED_COLS[0]] = get_discol(L, Lids)
            L[GROUPS_COLS[0]] = get_groupcol(L, Lids)
        elif len(L["specials"]) > 0:
            L[ENABLED_COLS[0]] = get_discol(L, L["specials"])
            L[GROUPS_COLS[0]] = get_groupcol(L, L["specials"])

    if Rvarcol is None:
        if RhasIds:
            R[ENABLED_COLS[0]] = get_discol(R, Rids)
            R[GROUPS_COLS[0]] = get_groupcol(R, Rids)
        elif len(R["specials"]) > 0:
            R[ENABLED_COLS[0]] = get_discol(R, R["specials"])
            R[GROUPS_COLS[0]] = get_groupcol(R, R["specials"])
    (Lcond_col, LcondIds, LisTime) = has_condition(L, other_params)
    (Rcond_col, RcondIds, RisTime) = has_condition(R, other_params)

    if not LhasIds and not RhasIds:
        if LcondIds is not None:
            LhasIds = True
            Lids = LcondIds
        if RcondIds is not None:
            RhasIds = True
            Rids = RcondIds
            
    (LhasCoord, Lcoord) = has_coord(L)
    (RhasCoord, Rcoord) = has_coord(R)

    nbr, offL, offR = get_size_hint(L, Lids, Lvarcol is None, R, Rids, Rvarcol is None)
    # if sizes don't match and either file had not header and could be sparse, try it...
    if nbr is None and not (LhasIds and RhasIds):
        sLids, sRids = (None, None)
        nbr_SF, nbr_FS, nbr_SS = (None, None, None)
        if not LhasIds and Lvarcol is None and len(L["headers"]) in [2, 3]:
            sLids = L["data"][L["headers"][0]]
            nbr_SF, offL_SF, offR_SF = get_size_hint(L, sLids, False, R, Rids, Rvarcol is None)

        if not RhasIds and Rvarcol is None and len(R["headers"]) in [2, 3]:
            sRids = R["data"][R["headers"][0]]
            nbr_FS, offL_FS, offR_FS = get_size_hint(L, Lids, Lvarcol is None, R, sRids, False)

            if sLids is not None:
                nbr_SS, offL_SS, offR_SS = get_size_hint(L, sLids, False, R, sRids, False)

        if (nbr_SF or nbr_SS) is not None:
            idsv = L["headers"].pop(0)
            LhasIds = True
            Lids = L["data"].pop(idsv)
            Lvarcol = L["headers"][0]
            if len(L["headers"]) == 2:
                Lvalcol = L["headers"][1]
        if (nbr_FS or nbr_SS) is not None:
            idsv = R["headers"].pop(0)
            RhasIds = True
            Rids = R["data"].pop(idsv)
            Rvarcol = R["headers"][0]
            if len(R["headers"]) == 2:
                Rvalcol = R["headers"][1]

        nbrTT, offLTT, offRTT = get_size_hint(L, Lids, Lvarcol is None, R, Rids, Rvarcol is None)
        if nbrTT is None and not (LhasIds and RhasIds):
            raise CSVRError('The two data sets are not of same size (%d ~ %d)' % (offL, offR))
        else:
            nbr, offL, offR = (nbrTT, offLTT, offRTT)

    # Sparse pairs storing with empty first row on sparse side
    if LhasIds and Lvarcol is not None:
        # if True:
        try:
            L, Lcoord, Lids, LhasCoord_sp, LhasIds = parse_sparse(L, Lcoord, Lids, Lvarcol, Lvalcol, offL)
            LhasCoord |= LhasCoord_sp
        except Exception as arg:
            raise CSVRError('Error while trying to parse sparse left hand side: %s' % arg)

    if RhasIds and Rvarcol is not None:
        try:
            R, Rcoord, Rids, RhasCoord_sp, RhasIds = parse_sparse(R, Rcoord, Rids, Rvarcol, Rvalcol, offR)
            RhasCoord |= RhasCoord_sp
        except IOError as arg:  # Exception as arg:
            raise CSVRError('Error while trying to parse sparse right hand side: %s' % arg)

    ro = RowOrder(LhasIds, RhasIds, Lids, Rids,
                  LisTime, RisTime, Lcond_col, Rcond_col,
                  LhasCoord, RhasCoord, Lcoord, Rcoord)

    if ro.hasOrders():
        Lorder, Rorder = ro.getOrders()

        coord = None
        try:
            # Order Lcoord according to Lorder
            if LhasCoord:
                coord = [(Lcoord[0][Lorder[i]], Lcoord[1][Lorder[i]]) for i in range(len(Lorder))]
            elif RhasCoord:
                coord = [(Rcoord[0][Rorder[i]], Rcoord[1][Rorder[i]]) for i in range(len(Rorder))]
        except Exception as arg:
            raise CSVRError('Error while trying to get the coordinates of data: %s' % arg)

        ids = None
        try:
            if LhasIds:
                ids = [Lids[Lorder[i]] for i in range(len(Lorder))]
            elif RhasIds:
                ids = [Rids[Rorder[i]] for i in range(len(Rorder))]
        except Exception as arg:
            raise CSVRError('Error while trying to get the ids of data: %s' % arg)

        # LDis = has_disabled_rows(L)
        # RDis = has_disabled_rows(R)

        cond_col = None
        CisTime = False
        try:
            # Order Lcoord according to Lorder
            if Lcond_col is not None:
                cond_col = [Lcond_col[Lorder[i]] for i in range(len(Lorder))]
                CisTime = LisTime
            elif Rcond_col is not None:
                cond_col = [Rcond_col[Rorder[i]] for i in range(len(Rorder))]
                CisTime = RisTime
        except Exception as arg:
            raise CSVRError('Error while trying to get the condition column of data: %s' % arg)

        return (L, R, Lorder, Rorder, coord, ids, cond_col, CisTime)

    else:
        # if not (LhasCoord or RhasCoord):
        #     # Neither has coordinates
        #     raise ValueError('At least one data file must have coordinates')
        # elif not (LhasCoord and RhasCoord):
        # Only one has coordinates (or none), do not re-order rows
        # TODO HERE PARSE SPARSE ALSO WITHOUT IDS ON BOTH SIDES
        if Lids is not None:  # e.g. from sparse
            nbrowsL = len(Lids)
        else:
            nbrowsL = len(first_value(L["data"]))

        if Rids is not None:  # e.g. from sparse
            nbrowsR = len(Rids)
        else:
            nbrowsR = len(first_value(R["data"]))

        # Sanity check
        if nbrowsR != nbrowsL:
            raise CSVRError('The two data sets without mapping ids are not of same size')

        data = L['data']
        head = L['headers']
        # extract the coordinates
        if LhasCoord:
            coord = list(zip(*Lcoord))
        elif RhasCoord:
            coord = list(zip(*Rcoord))
        else:
            coord = None

        if Lcond_col is not None:
            cond_col = Lcond_col
            CisTime = LisTime
        elif Rcond_col is not None:
            cond_col = Rcond_col
            CisTime = RisTime
        else:
            cond_col = None
            CisTime = False

        ids = None
        if LhasIds:  # and len(L["data"].values()[0]) == len(Lids):
            ids = Lids
        elif RhasIds:  # and if len(R["data"].values()[0]) == len(Rids):
            ids = Rids

        # LDis = has_disabled_rows(L)
        # RDis = has_disabled_rows(R)

        return (L, R, range(nbrowsL), range(nbrowsR), coord, ids, cond_col, CisTime)


def row_order_single(L, other_params={}):
    # TODO catch the dense row containing info on enabled columns
    (LhasIds, Lids) = has_ids(L)
    (Lvarcol, Lvalcol) = is_sparse(L)

    if Lvarcol is None:
        if LhasIds:
            L[ENABLED_COLS[0]] = get_discol(L, Lids)
            L[GROUPS_COLS[0]] = get_groupcol(L, Lids)
        elif len(L["specials"]) > 0:
            L[ENABLED_COLS[0]] = get_discol(L, L["specials"])
            L[GROUPS_COLS[0]] = get_groupcol(L, L["specials"])

    (Lcond_col, LcondIds, LisTime) = has_condition(L, other_params)
    if LcondIds is not None and not LhasIds:
        LhasIds = True
        Lids = LcondIds

    (LhasCoord, Lcoord) = has_coord(L)

    if LhasIds and Lvarcol is not None:
        if True:
            #        try:
            L, Lcoord, Lids, LhasCoord_sp, LhasIds = parse_sparse(L, Lcoord, Lids, Lvarcol, Lvalcol, off=0)
            LhasCoord |= LhasCoord_sp
        # except Exception as arg:
        #     raise CSVRError('Error while trying to parse sparse left hand side: %s' % arg)

        if LhasCoord:
            coord = list(zip(*Lcoord))
        else:
            coord = None
        ids = Lids

        if Lcond_col is not None:
            cond_col = Lcond_col
            CisTime = LisTime
        else:
            cond_col = None
            CisTime = False

        return (L, range(len(ids)), coord, ids, cond_col, CisTime)

    else:
        # if not (LhasCoord or RhasCoord):
        #     # Neither has coordinates
        #     raise ValueError('At least one data file must have coordinates')
        # elif not (LhasCoord and RhasCoord):
        # Only one has coordinates (or none), do not re-order rows
        # TODO HERE PARSE SPARSE ALSO WITHOUT IDS ON BOTH SIDES
        if Lids is not None:  # e.g. from sparse
            nbrowsL = len(Lids)
        else:
            nbrowsL = len(first_value(L["data"]))

        data = L['data']
        head = L['headers']
        # extract the coordinates
        if LhasCoord:
            coord = list(zip(*Lcoord))
        else:
            coord = None

        if LhasIds:  # and len(L["data"].values()[0]) == len(Lids):
            ids = Lids
        else:
            ids = None

        if Lcond_col is not None:
            cond_col = Lcond_col
            CisTime = LisTime
        else:
            cond_col = None
            CisTime = False

        # LDis = has_disabled_rows(L)

        return (L, range(nbrowsL), coord, ids, cond_col, CisTime)


def readTransCSV(trans_filename, csv_params={}, other_params={}):
    unknown_string = other_params.get("NA_str", NA_str_c)
    (Th, Td, Ttype, Tspecials) = read_csv(trans_filename, csv_params, unknown_string)
    T = {'data': Td, 'headers': Th, "sparse": False, "type_all": Ttype, ENABLED_COLS[0]: {}, GROUPS_COLS[0]: {}, "specials": Tspecials}

    (T, Torder, coord, ids, cond_col, CisTime) = row_order_single(T, other_params)
    T['order'] = Torder
    T["type_all"] = Ttype
    return T, ids


def importCSV(left_filename, right_filename, csv_params={}, other_params={}):
    single_dataset = (left_filename == right_filename) or (right_filename is None)
    unknown_string = other_params.get("NA_str", NA_str_c)
    try:
        (Lh, Ld, Ltype, Lspecials) = read_csv(left_filename, csv_params, unknown_string)
    except ValueError as arg:
        raise CSVRError("Error reading the left hand side data: %s" % arg)
    except csv.Error as arg:
        raise CSVRError("Error reading the left hand side data: %s" % arg)
    L = {'data': Ld, 'headers': Lh, "sparse": False, "type_all": Ltype, ENABLED_COLS[0]: {}, GROUPS_COLS[0]: {}, "specials": Lspecials}

    if single_dataset:
        (L, Lorder, coord, ids, cond_col, CisTime) = row_order_single(L, other_params)
        L['order'] = Lorder
        L["type_all"] = Ltype
        R = None
    else:
        try:
            (Rh, Rd, Rtype, Rspecials) = read_csv(right_filename, csv_params, unknown_string)
        except ValueError as arg:
            raise CSVRError("Error reading the right hand side data: %s" % arg)
        except csv.Error as arg:
            raise CSVRError("Error reading the right hand side data: %s" % arg)
        R = {'data': Rd, 'headers': Rh, "sparse": False, "type_all": Rtype, ENABLED_COLS[0]: {}, GROUPS_COLS[0]: {}, "specials": Rspecials}
        (L, R, Lorder, Rorder, coord, ids, cond_col, CisTime) = row_order(L, R, other_params)
        L['order'] = Lorder
        R['order'] = Rorder
        L["type_all"] = Ltype
        R["type_all"] = Rtype

    return {'data': (L, R), 'coord': coord, "ids": ids, 'cond_col': cond_col, "c_time": CisTime}, single_dataset, unknown_string


def print_out(data):
    keysL = sorted(data['data'][0]['headers'])
    keysR = sorted(data['data'][1]['headers'])

    line = "# ;"
    if data['ids'] is not None:
        line += "ID; "
    if data['coord'] is not None:
        line += "coord0; coord1; "
    if data['cond_col'] is not None:
        tstr = ""
        if data['c_time']:
            tstr = "(T)"
        line += "cond_col%s; " % tstr

    line += " | "
    line += "; ".join(["%s" % k for k in keysL])
    line += " || "
    line += "; ".join(["%s" % k for k in keysR])
    line += " |"
    print(line)

    for row in range(len(data['data'][0]['order'])):
        line = "%d; " % row
        if data['ids'] is not None:
            line += "%s; " % data['ids'][row]
        if data['coord'] is not None:
            line += "%s; %s; " % (data['coord'][row][0], data['coord'][row][1])
        if data['cond_col'] is not None:
            line += "%s; " % data['cond_col'][row]

        line += " | "
        if data["data"][0]["sparse"]:
            if data["data"][0]["bool"]:
                line += "; ".join(["%s" % int(data['data'][0]['order'][row] in data['data'][0]['data'][k]) for k in keysL])
            else:
                line += "; ".join(["%s" % data['data'][0]['data'][k].get(data['data'][0]['order'][row], "--") for k in keysL])
        else:
            line += "; ".join(["%s" % data['data'][0]['data'][k][data['data'][0]['order'][row]] for k in keysL])
        line += " || "
        if data["data"][1]["sparse"]:
            if data["data"][1]["bool"]:
                line += "; ".join(["%s" % int(data['data'][1]['order'][row] in data['data'][1]['data'][k]) for k in keysR])
            else:
                line += "; ".join(["%s" % data['data'][1]['data'][k].get(data['data'][1]['order'][row], "--") for k in keysR])
        else:
            line += "; ".join(["%s" % data['data'][1]['data'][k][data['data'][1]['order'][row]] for k in keysR])
        line += " |"
        print(line)


def main(argv=[]):

    rep = "/home/egalbrun/short/raja_small/"
    res, single, unknown_str = importCSV(rep+"data_RHSk.csv", rep+"data_RHSo.csv")
    print_out(res)


if __name__ == '__main__':
    main(sys.argv)
