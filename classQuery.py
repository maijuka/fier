import re
import random
import operator
import itertools
import numpy
import copy
import time
import datetime
from dateutil import parser
import time
import scipy.spatial.distance

try:
    from classSParts import SParts, SSetts, cmp_vals, cmp_reclists
    from redquery_parser import RedQueryParser
    from grako.exceptions import *  # @UnusedWildImport
except ModuleNotFoundError:
    from .classSParts import SParts, SSetts, cmp_vals, cmp_reclists
    from .redquery_parser import RedQueryParser
    from .grako.exceptions import *  # @UnusedWildImport
import pdb

NA_str_c = "nan"

ENTITY_MARK = 'n'
VARIABLE_MARK = 'v'
XPR_MARK = "="


def getNameCol(cid, names):
    try:
        return names[cid]
    except (IndexError, TypeError):
        raise Warning("Names does not contains this column cid=%d vs. len(names)=%d" % (cid, len(names)))
    return Term.pattVName % cid


def getFmtCol(cid, fmts):
    try:
        return fmts[cid]
    except (IndexError, TypeError):
        pass
    return {}
########################


TIME_FMT = "%Y-%b-%d,%H:%M:%S"


class TimeTools:

    MTCH_TIME = "_time$"
    TIME_FMT = TIME_FMT
    TIME_ELEMS = [TIME_FMT[i] for i in range(1, len(TIME_FMT), 3)]

    @classmethod
    def time_prec_to_range(tcl, time_prec):
        if time_prec is None:
            return (None, None)
        ii = [i for i, l in enumerate(tcl.TIME_ELEMS) if ("%"+l) in time_prec]
        return numpy.min(ii), numpy.max(ii)

    @classmethod
    def range_to_time_prec(tcl, tmin=None, tmax=None):
        if tmin is None:
            low = 0
        else:
            low = 3*tmin
        if tmax is None:
            up = len(tcl.TIME_FMT)
        else:
            up = 3*(tmax+1)-1
            if up >= len(tcl.TIME_FMT):
                up = len(tcl.TIME_FMT)
        return tcl.TIME_FMT[low:up]

    @classmethod
    def lower_time_prec(tcl, time_prec):
        rng = tcl.time_prec_to_range(time_prec)
        if rng[0] < rng[1]:
            return tcl.range_to_time_prec(rng[0], rng[1]-1)
        return tcl.range_to_time_prec(rng[0], rng[1])

    @classmethod
    def higher_time_prec(tcl, time_prec):
        rng = tcl.time_prec_to_range(time_prec)
        return tcl.range_to_time_prec(rng[0], rng[1]+1)

    @classmethod
    def isTimeVarName(tcl, var_name):
        return re.search(tcl.MTCH_TIME, var_name) is not None

    @classmethod
    def asTimeVar(tcl, cid, names, fmts):
        try:
            return type(names) == list and (tcl.isTimeVarName(getNameCol(cid, names)) or getFmtCol(cid, fmts).get("time_prec") is not None)
        except:
            pdb.set_trace()
            print(fmts)

    @classmethod
    def get_time_prec(tcl, time_struct):
        parts = [time_struct.tm_year, time_struct.tm_mon, time_struct.tm_mday, time_struct.tm_hour, time_struct.tm_min, time_struct.tm_sec]
        ii = [i for i, l in enumerate(parts) if l > 0]
        return tcl.range_to_time_prec(numpy.min(ii), numpy.max(ii))

    @classmethod
    def format_time(tcl, v, time_prec=None):
        time_struct = time.localtime(v)
        if time_prec is None:
            return "%s" % v
        if time_prec == 1:
            time_prec = TimeTools.TIME_FMT
        if time_prec == -1:
            time_prec = tcl.get_time_prec(time_struct)
            # if time_struct.tm_hour + time_struct.tm_min + time_struct.tm_sec == 0:
            #     time_prec = DATE_FMT
            # else:
            #     if time_struct.tm_year == 0:
            #         time_prec = 0
            #     else:
            #         time_prec = DATE_FMT + "," + TIME_FMT
        if time_prec == 0:
            return "%s" % datetime.timedelta(seconds=v)
        return time.strftime(time_prec, time_struct)

    @classmethod
    def parse_time(tcl, v, dayfirst=False, yearfirst=True):
        return time.mktime(parser.parse(v, dayfirst=dayfirst, yearfirst=yearfirst).timetuple())

# truth table (TT) manipulation tools


def foldRowsTT(tt):
    changed = False
    tts = numpy.argsort(numpy.abs(0.5*tt.shape[0] - tt.sum(axis=0)))
    for i in tts:
        cols = [j for j in range(tt.shape[1]) if j != i]
        idL = numpy.where(tt[:, i] == 1)[0]
        idR = numpy.where(tt[:, i] == 0)[0]
        ps = []
        if len(cols) == 0:  # No other column, always true
            if idL.size == idR.size == 1:
                ps = [(0, 0)]
        else:
            L = tt[idL, :]
            R = tt[idR, :]
            ps = zip(*numpy.where(scipy.spatial.distance.cdist(L[:, cols], R[:, cols], 'hamming') == 0))
        keep = numpy.ones(tt.shape[0], dtype=numpy.bool)
        for p in ps:
            changed = True
            tt[idL[p[0]], i] = -1
            keep[idR[p[1]]] = False
        tt = tt[keep, :]
    return tt, changed


def foldColsTT(tt):
    # org = tt.copy()
    idsO = []
    # idsM = []
    ps = zip(*numpy.where(scipy.spatial.distance.squareform((scipy.spatial.distance.pdist(tt, 'hamming')*tt.shape[1] == 2))))
    # while(len(ps)) > 0:
    for rows in ps:
        # rows = ps[0]
        cols = numpy.where(tt[rows[0], :] != tt[rows[1], :])[0]
        block = tt[rows, :][:, cols]
        pr, pc = numpy.where(block == -1)
        if len(pr) == 1 and pr[0] == 1:
            changed = True
            if pc[0] == 1:
                # tt[rows[0], cols[0]] = -1
                idsO.append((rows[0], cols[0]))
                # idsM.append((rows[1], cols[1]))
            else:
                # tt[rows[0], cols[1]] = -1
                idsO.append((rows[0], cols[1]))
                # idsM.append((rows[1], cols[0]))
    if len(idsO) > 0:
        rx, cx = zip(*idsO)
        tt[rx, cx] = -1
    return tt, len(idsO) > 0


def subsRowsTT(tt):
    tts = numpy.argsort(-numpy.sum(tt == -1, axis=1))
    keep = numpy.ones(tt.shape[0], dtype=numpy.bool)
    for row in tts:
        if keep[row]:
            cmask = tt[row, :] != -1
            ps = numpy.where(scipy.spatial.distance.cdist(tt[:, cmask], [tt[row, cmask]], 'hamming') == 0)[0]
            for p in ps:
                if p != row:
                    keep[p] = False
    return tt[keep, :], numpy.sum(~keep) > 0


def triplesRowsTT(tt):
    tts = numpy.argsort(-numpy.sum(tt == -1, axis=1))
    keep = numpy.ones(tt.shape[0], dtype=numpy.bool)
    pivots = []
    pairs = []
    for row in tts:
        # if row == 2:
        #     print(tt)
        #     pdb.set_trace()
        cmask = tt[row, :] == -1
        ps = zip(*numpy.where(scipy.spatial.distance.squareform((scipy.spatial.distance.pdist(tt[:, cmask] == -1, 'hamming') == 0))))
        for p in ps:
            if p[0] < p[1]:
                matches = tt[p[0], :] == tt[p[1], :]
                # MATCHING CONDITIONS:
                # a) where row == -1:
                # i) all -1 match, ii, there are non -1, ii) in all of them p[0] and p[1] are complements
                # b) where row != -1:
                # i) one value is -1 (if rows don't match), ii) the other (or both if match) is the value in row
                if (numpy.sum(cmask * matches) == 0 or numpy.max(tt[p, :][:, cmask * matches]) == -1) and \
                   numpy.sum(cmask * ~matches) > 0 and \
                   numpy.all(tt[p[0], cmask * ~matches] == (1-tt[p[1], cmask * ~matches])) and \
                   numpy.all(numpy.min(tt[p, :][:, ~cmask * ~matches], axis=0) == -1) and \
                   numpy.all(numpy.max(tt[p, :][:, ~cmask], axis=0) == tt[row, ~cmask]):
                    pivots.append(row)
                    pairs.append(p)

            # for p in ps:
            #     if p[0] < p[1]:
            #         matches = tt[p[0],:] == tt[p[1],:]
            #         if numpy.max(numpy.min(tt[p,:][:, ~(cmask + matches)],axis=0)) == -1 \
            #                and (numpy.sum(matches) == 0 or numpy.min(tt[p[0], matches]) > -1):
            #             pivots.append(row)
            #             pairs.append(p)

    if len(pivots) > 0:
        if len(pivots) > 1:
            print("More than one pivot", pivots, pairs)
            pdb.set_trace()
        keep[pivots] = False
        # print("Found triple")
    return tt[keep, :], numpy.sum(~keep) > 0


def simplerTT(tt):
    # nbrows = tt.shape[0]
    # nbnn = numpy.sum(tt>-1)
    # ttchanged = 0
    changed = [True, False, False, False]
    if tt.shape[0] > 0:
        # print("SIMPLER TT start\n", tt)
        while sum(changed) > 0:
            changed = [False, False, False, False]
            tt, changed[0] = foldRowsTT(tt)
            # if changed[0]:
            #     print("SIMPLER TT foldRows\n", tt)
            tt, changed[1] = foldColsTT(tt)
            # if changed[1]:
            #     print("SIMPLER TT foldCols\n", tt)
            tt, changed[2] = subsRowsTT(tt)
            # if changed[2]:
            #     print("SIMPLER TT subsRowsTT\n", tt)
            tt, changed[3] = triplesRowsTT(tt)
        #     if changed[3]:
        #         print("SIMPLER TT tripleRowsTT\n", tt)
        # print("SIMPLER TT done\n", tt)
    #         ttchanged += sum(changed)
    # if ttchanged > 0:
    #     print("SIMPLIFY FROM (%d, %d) TO (%d, %d)" % (nbrows, nbnn, tt.shape[0], numpy.sum(tt>-1)))
    return tt


def recurse_numeric(b, function, args={}):
    if type(b) is list:
        out = 0
        for bi, bb in enumerate(b):
            nargs = dict(args)
            if "trace" in args:
                nargs["trace"] = [bi]+args["trace"]
            out += recurse_numeric(bb, function, nargs)
        return out
    else:
        return function(b, **args)

# WARNING THIS DOES NOT RECURSE ON NEGATIONS


def recurse_list(b, function, args={}):
    if type(b) is list:
        out = []
        for bi, bb in enumerate(b):
            nargs = dict(args)
            if "trace" in args:
                nargs["trace"] = [bi]+args["trace"]
            tou = recurse_list(bb, function, nargs)
            if type(tou) is list:
                out.extend(tou)
            elif tou is not None:
                out.append(tou)
        return out
    elif isinstance(b, Literal):
        return function(b, **args)


def recurse_deep(b, function, args={}):
    if type(b) is list:
        out = []
        for bi, bb in enumerate(b):
            nargs = dict(args)
            if "trace" in args:
                nargs["trace"] = [bi]+args["trace"]
            tmp = recurse_deep(bb, function, nargs)
            if tmp is not None:
                out.append(tmp)
        return out
    else:
        return function(b, **args)


class SYM(object):

    SYMU_OR = '\u2228'
    SYMU_AND = '\u2227'
    SYMU_NOT = '\u00ac'
    SYMU_LEQ = '\u2264'
    SYMU_EIN = '\u2208'
    SYMU_NIN = '\u2209'
    SYMU_NEQ = '\u2260'

    SYMU_SETMIN = "\u2216"

    SYMU_ARRTOP = "\u2191"
    SYMU_ARRBOT = "\u2193"
    # SYMU_LEARN = '\u25e9'
    # SYMU_TEST = '\u25ea'
    # SYMU_LEARN = '\u25d6'
    # SYMU_TEST = '\u25d7'
    # SYMU_LEARN = '\u25d0'
    # SYMU_TEST = '\u25d1'
    SYMU_LEARN = '\u25d5'
    SYMU_TEST = '\u25d4'
    SYMU_RATIO = '\u2298'
    # SYMU_CROSS = '\u274c'
    SYMU_CROSS = '\u2715'
    SYMU_INOUT = '\u21d7'
    SYMU_OUTIN = '\u21d8'

    SYMO_OR = 'OR'
    SYMO_AND = 'AND'
    SYMO_NOT = 'NOT'
    SYMO_LEQ = '<'
    SYMO_EIN = 'IN'
    SYMO_NIN = 'NOT IN'
    SYMO_NEQ = '~'

    SYMO_SETMIN = "\\"

    SYMO_ARRTOP = "^"
    SYMO_ARRBOT = "v"
    SYMO_LEARN = '[l]'
    SYMO_TEST = '[t]'
    SYMO_RATIO = '[r]'
    SYMO_CROSS = 'X'
    SYMO_INOUT = '>>'
    SYMO_OUTIN = '<<'

    # WITH UNICODE
    SYM_OR = SYMU_OR
    SYM_AND = SYMU_AND
    SYM_NOT = SYMU_NOT
    SYM_LEQ = SYMU_LEQ
    SYM_EIN = SYMU_EIN
    SYM_NIN = SYMU_NIN
    SYM_NEQ = SYMU_NEQ

    SYM_SETMIN = SYMU_SETMIN

    SYM_ARRTOP = SYMU_ARRTOP
    SYM_ARRBOT = SYMU_ARRBOT
    SYM_LEARN = SYMU_LEARN
    SYM_TEST = SYMU_TEST
    SYM_RATIO = SYMU_RATIO
    SYM_CROSS = SYMU_CROSS
    SYM_INOUT = SYMU_INOUT
    SYM_OUTIN = SYMU_OUTIN

    # WITHOUT UNICODE
    # SYM_OR = SYMO_OR
    # SYM_AND = SYMO_AND
    # SYM_NOT = SYMO_NOT
    # SYM_LEQ = SYMO_LEQ
    # SYM_EIN = SYMO_EIN
    # SYM_NIN = SYMO_NIN
    # SYM_NEQ = SYMO_NEQ

    # SYM_SETMIN = SYMO_SETMIN

    # SYM_ARRTOP = SYMO_ARRTOP
    # SYM_ARRBOT = SYMO_ARRBOT
    # SYM_LEARN = SYMO_LEARN
    # SYM_TEST = SYMO_TEST
    # SYM_RATIO = SYMO_RATIO
    # SYM_CROSS = SYMO_CROSS
    # SYM_INOUT = SYMO_INOUT
    # SYM_OUTIN = SYMO_OUTIN


class DictStyleSymbols:

    def_style = "txt"
    dict_styles_syms = {def_style: {}}

    @classmethod
    def getTkKey(tcl, key, opening=0):
        if opening > 0:
            return "SYM_"+key+"_OPEN"
        elif opening < 0:
            return "SYM_"+key+"_CLOSE"
        return "SYM_"+key

    @classmethod
    def getKeyTk(tcl, tk):
        mtc = re.match("SYM_(.*)_OPEN", tk)
        if mtc:
            return mtc.group(1), 1
        mtc = re.match("SYM_(.*)_CLOSE", tk)
        if mtc:
            return mtc.group(1), -1
        mtc = re.match("SYM_(.*)", tk)
        if mtc:
            return mtc.group(1), 0
        return None, None
    
    @classmethod
    def getFmtKey(tcl, key, opening=0):
        return "%("+tcl.getTkKey(key, opening)+")s"

    @classmethod
    def hasKey(tcl, key, opening=0):
        return tcl.getTkKey(key, opening) in tcl.dict_styles_syms[tcl.def_style]

    @classmethod
    def pieceContainsKey(tcl, piece, key, opening=0):
        if opening == 0:
            return re.search("\("+tcl.getTkKey(key, opening)+"\)", piece)
        else:
            return re.search("\("+tcl.getTkKey(key, 1)+"\).*\("+tcl.getTkKey(key, -1)+"\)", piece)

    @classmethod
    def upgradeKey(tcl, piece, key, opening=0):
        if tcl.hasKey(key+"_BIG", opening) and tcl.pieceContainsKey(piece, key, opening):
            return key+"_BIG"
        return key
    
    @classmethod
    def enclose(tcl, piece, encl):
        if tcl.hasKey(encl, 1): # valid enclosing?
            if encl == "TEXT":
                mtc = re.match("\s*\$(.*)\$\s*$", piece)
                if mtc is not None:
                    return mtc.group(1)
            up_encl = tcl.upgradeKey(piece, encl, opening=1)            
            return tcl.getFmtKey(up_encl, 1) + piece + tcl.getFmtKey(up_encl, -1)
        return piece

    # @classmethod
    # def cleanEnclosed(tcl, piece, keys=None):
    #     if keys is None:
    #         all_keys = [tcl.getKeyTk(tk) for tk in tcl.dict_styles_syms[tcl.def_style]]
    #         keys = [c[0] for c in all_keys if c[1] == 1]
    #     old_piece = piece
    #     for key in keys:
    #         ## successives
    #         patt = "%\("+tcl.getTkKey(key, -1)+"\)s\s*%\("+tcl.getTkKey(key, 1)+"\)s"
    #         piece = re.sub(patt, " ", piece)
    #         ## empty
    #         patt = "%\("+tcl.getTkKey(key, 1)+"\)s\s*%\("+tcl.getTkKey(key, -1)+"\)s"
    #         piece = re.sub(patt, " ", piece)
    #     return piece

    
    @classmethod
    def prepare(tcl, tks, encl_all=[], encl_each=[], sep=",", protect=False):
        if type(tks) is str:
            tks = [tks]
        if type(encl_each) is str:
            encl_each = [encl_each]
        if type(encl_all) is str:
            encl_all = [encl_all]
        encl_all_pro = []
        if protect:
            encl_all_pro = ["PROTECT"]
        encl_tks = tks
        for encl in encl_each:
            encl_tks = [tcl.enclose(encl_tk, encl) for encl_tk in encl_tks]
        if tcl.hasKey(sep):
            sep = tcl.getFmtKey(sep)
        encl_tk = sep.join(encl_tks)
        for encl in encl_all+encl_all_pro:
            encl_tk = tcl.enclose(encl_tk, encl)
        return encl_tk
    
    @classmethod
    def getTermTk(tcl, tk, type_letter="", protect=True):
        if tcl.hasKey("TERM"+type_letter, 1):
            return tcl.prepare(tk, encl_all=["TERM"+type_letter], protect=protect)
        return tcl.prepare(tk, encl_all=["TERM"], protect=protect)
        
    
    @classmethod
    def addSymbol(tcl, tk, stym):
        if type(stym) is str:
            tcl.dict_styles_syms[tcl.def_style][tk] = stym
        elif type(stym) is dict:
            for (style, sym) in stym.items():
                if style not in tcl.dict_styles_syms:
                    tcl.dict_styles_syms[style] = {}
                tcl.dict_styles_syms[style][tk] = sym

        ## automatically add html symbol
        def_sym = tcl.dict_styles_syms[tcl.def_style][tk]                
        if tcl.dict_styles_syms.get("utf", {}).get(tk, def_sym) != def_sym:
            if "html" not in tcl.dict_styles_syms:
                tcl.dict_styles_syms["html"] = {}
            tcl.dict_styles_syms["html"][tk] = "".join(["&#%d;" % ord(c) if c!=" " else c for c in tcl.dict_styles_syms["utf"][tk]])
                            
    @classmethod
    def addKeySymbols(tcl, key, stym=None, stym_open=None, stym_close=None):
        if stym is not None:
            tcl.addSymbol(tcl.getTkKey(key), stym)
        if stym_open is not None:
            tcl.addSymbol(tcl.getTkKey(key, 1), stym_open)
        if stym_close is not None:
            tcl.addSymbol(tcl.getTkKey(key, -1), stym_close)

    @classmethod
    def addTermSymbols(tcl, type_letter, stym_open, stym_close):
        tcl.addKeySymbols("TERM"+type_letter, stym_open=stym_open, stym_close=stym_close)

                
    @classmethod
    def getStyleDict(tcl, style):
        d = dict(tcl.dict_styles_syms[tcl.def_style])
        d.update(tcl.dict_styles_syms.get(style, {}))
        return d

    @classmethod
    def getStyledTk(tcl, style, tk):
        if tk in tcl.dict_styles_syms.get(style, {}):
            return tcl.dict_styles_syms[style][tk]
        return tcl.dict_styles_syms[tcl.def_style].get(tk, "(?"+tk+")")

    @classmethod
    def replaceStyleTks(tcl, fstr, style):
        d = dict(tcl.dict_styles_syms[tcl.def_style])
        d.update(tcl.dict_styles_syms.get(style, {}))
        return fstr % d

    
    
DSS = DictStyleSymbols()
DSS.addKeySymbols("PROTECT", None, {"txt": "", "tex": "$", "utf": ""}, {"txt": "", "tex": "$", "utf": ""})
DSS.addKeySymbols("TEXT", None, {"txt": "", "tex": "\\text{", "utf": ""}, {"txt": "", "tex": "}", "utf": ""})
DSS.addKeySymbols("BRCK", None, {"txt": "(", "tex": "(", "utf": "("}, {"txt": ")", "tex": ")", "utf": ")"})
DSS.addKeySymbols("BRCK_BIG", None, {"txt": "(", "tex": "\\big(", "utf": "("}, {"txt": ")", "tex": "\\big)", "utf": ")"})
DSS.addKeySymbols("CURL", None, {"txt": "{", "tex": "\\{", "utf": "{"}, {"txt": "}", "tex": "\\}", "utf": "}"})
DSS.addKeySymbols("TERM", None, {"txt": "[", "tex": "[", "utf": "["}, {"txt": "]", "tex": "]", "utf": "]"})

class Op(object):

    opsBool = {True: 1, False: -1}
    ops = {0: 'X', 1: '|', -1: '&'}
    sym_tks = {0: DSS.getFmtKey("OP_X"), 1: DSS.getFmtKey("OP_OR"), -1: DSS.getFmtKey("OP_AND")}


    DSS.addKeySymbols("OP_X", "X")
    DSS.addKeySymbols("OP_OR", {"txt": " or ", "tex": " \\lor{} ", "utf": " "+SYM.SYM_OR+" "})
    DSS.addKeySymbols("OP_AND", {"txt": " and ", "tex": " \\land{} ", "utf": " "+SYM.SYM_AND+" "})

    def __init__(self, nval=0):
        if type(nval) == bool:
            self.val = self.opsBool[nval]
        elif nval in Op.ops:
            self.val = nval
        elif isinstance(nval, Op):
            self.val = nval.val
        else:
            raise Exception('Uninterpretable operator: %s !' % nval)

    def isSet(self):
        return self.val != 0

    def flip(self):
        self.val *= -1

    def copy(self):
        return Op(self.val)

    def other(self):
        return Op(-self.val)

    def __int__(self):
        return self.val

    def isOr(self):
        return self.val == 1

    def isAnd(self):
        return self.val == -1

    def __str__(self):
        return Op.ops[self.val]

    def disp(self):
        return str(self)
    
    def getDispTk(self, protect=False):
        return DSS.prepare(Op.sym_tks[self.val], protect=protect)
    
    def styledDisp(self, style="", protect=True, str_before="", str_after=""):
        return str_before + DSS.replaceStyleTks(self.getDispTk(protect), style) + str_after

    def compare(self, other):
        if type(other) is bool:
            return cmp_vals(self.val, self.opsBool[other])
        if other is None or not isinstance(other, Op):
            return 1
        return cmp_vals(self.val, other.val)

    def __eq__(self, other):
        return self.compare(other) == 0

    def toKey(self):
        return (self.val, )
    def __hash__(self):
        return hash(self.toKey())
    

class Neg(object):
    sym = ['', '! ']
    sym_tks = {False: DSS.getFmtKey("NONNEG"), True: DSS.getFmtKey("NEG")}
    DSS.addKeySymbols("NONNEG", {"txt": "", "tex": "", "utf": ""})
    DSS.addKeySymbols("NEG", {"txt": "not ", "tex": "\\neg ", "utf": SYM.SYM_NOT+" "})

    def __init__(self, nNeg=False):
        if nNeg == True or nNeg < 0:
            self.neg = -1
        else:
            self.neg = 1

    def copy(self):
        return Neg(self.neg)

    def boolVal(self):
        return self.neg < 0

    def flip(self):
        self.neg *= -1

    def __int__(self):
        return self.neg

    def toKey(self):
        return (self.neg, )
    def __hash__(self):
        return hash(self.toKey())
    
    def __str__(self):
        return Neg.sym[self.boolVal()]
    
    def disp(self):
        return str(self)

    def getDispTk(self, protect=False):
        return DSS.prepare(Neg.sym_tks[self.boolVal()], protect=protect)
    
    def styledDisp(self, style="", protect=True, str_before="", str_after=""):
        return str_before + DSS.replaceStyleTks(self.getDispTk(protect), style) + str_after


class Term(object):

    pattVName = VARIABLE_MARK+"%d"
    type_id = 0
    type_letter = '-'
    type_name = '-'
    
    def __init__(self, ncol, vX=None, vY=None):
        self.col = ncol

    def toKey(self):
        return (self.colId(), self.typeId(), -1, -1)
    def __hash__(self):
        return hash(self.toKey())

    
    def __eq__(self, other):
        return isinstance(other, Term) and self.toKey() == other.toKey()

    def __lt__(self, other):
        return isinstance(other, Term) and self.toKey() < other.toKey()

    def __le__(self, other):
        return isinstance(other, Term) and self.toKey() <= other.toKey()

    def __gt__(self, other):
        return not isinstance(other, Term) or self.toKey() > other.toKey()

    def __ge__(self, other):
        return not isinstance(other, Term) or self.toKey() >= other.toKey()

    def valRange(self):
        return [None, None]

    def specifies(self):
        return None

    def specValues(self):  # always a list
        return []

    def isAnon(self):
        return False

    def getAdjusted(self, bounds):
        tmp = self.copy()
        tmp.setRange(bounds)
        return tmp

    def setRange(self, newR):
        pass

    def getCanonical(self, is_neg=False):
        # term with negation and complement information
        return (self, is_neg, False)

    def getComplement(self):
        return None

    def isComplement(self, term):
        return False

#     def simple(self, neg):
#         return neg

    def copy(self):
        return Term(self.col)

    def colId(self):
        return self.col

    def replaceCol(self, cid):
        self.col = cid
    
    def typeId(self):
        return self.type_id

    def cmpCol(self, other):
        if other is None or not isinstance(other, Term):
            return 1
        return cmp_vals(self.col, other.col)

    def cmpType(self, other):
        if other is None or not isinstance(other, Term):
            return 1
        return cmp_vals(self.typeId(), other.typeId())

    def getDispTermTk(self, tk, protect=False):
        return DSS.getTermTk(tk, self.type_letter, protect=protect)
    
    def getDispCol(self, names=None, fmts=None):
        if type(names) == list and len(names) > 0:
            return getNameCol(self.col, names)
        else:
            return Term.pattVName % self.col
        
    def getDispColTk(self, names=None, fmts=None, protect=False):
        if type(names) == list and len(names) > 0:
            return DSS.prepare(getNameCol(self.col, names), "TEXT", protect=protect)
        else:
            return DSS.prepare(Term.pattVName % self.col, protect=protect)

    def getDispNeg(self, neg=None):
        if neg is None:
            neg = ""
        if type(neg) == bool:
            neg = Neg(neg)
        if type(neg) is Neg:
            return str(neg)
        return neg
    
    def getDispNegTk(self, neg=None, protect=False):
        if neg is None:
            neg = ""
        if type(neg) == bool:
            neg = Neg(neg)
        if type(neg) is Neg:
            return neg.getDispTk(protect)
        return neg

    def disp(self, neg=None, names=None, fmts=None):
        return self.getDispNeg(neg) + self.getDispCol(names, fmts)

    def __str__(self):
        return self.disp()
        
    def getDispTk(self, neg=None, names=None, fmts=None, protect=False):
        return DSS.prepare(self.getDispNegTk(neg) + self.getDispTermTk(self.getDispColTk(names, fmts)), protect=protect)
  
    def styledDisp(self, neg=None, names=None, fmts=None, style="", protect=True, str_before="", str_after=""):
        return str_before + DSS.replaceStyleTks(self.getDispTk(neg, names, fmts, protect), style) + str_after
    

class BoolTerm(Term):
    type_id = 1
    type_letter = 'B'
    type_name = 'Boolean'

    DSS.addTermSymbols(type_letter, {"txt": "", "tex": "", "utf": ""}, {"txt": "", "tex": "", "utf": ""})
    
    def valRange(self):
        return [True, True]

    def specifies(self):
        return True

    def specValues(self):  # always a list
        return [True]

    def copy(self):
        return BoolTerm(self.col)

    def toKey(self):
        return (self.colId(), self.typeId(), 0, 0)
    
    def truthEval(self, variableV):
        return variableV




class CatTermONE(Term):  # LEGACY
    type_id = 2
    type_letter = 'C'
    type_name = 'Categorical'
    
    sym_tks = {False: DSS.getFmtKey("CAT_IS"), True: DSS.getFmtKey("CAT_ISNOT")}
    DSS.addKeySymbols("CAT_IS", {"txt": " = ", "tex": " = ", "utf": " = "})
    DSS.addKeySymbols("CAT_ISNOT", {"txt": " != ", "tex": " \\neq{} ", "utf": " "+SYM.SYM_NEQ+" "})

    
    def __init__(self, ncol, ncat, vY=None):
        self.col = ncol
        self.cat = ncat

    def getCat(self):
        return self.cat

    def valRange(self):
        return [self.getCat(), self.getCat()]

    def specifies(self):
        return self.getCat()

    def specValues(self):  # always a list
        return [self.getCat()]

    def setRange(self, cat):
        self.cat = cat

    def copy(self):
        return CatTermONE(self.col, self.cat)

    def truthEval(self, variableV):
        return variableV == self.cat

    def getCatsStr(self, op_curl="{", cl_curl="}", sep=",", op_any="", cl_any=""):
        return op_any + "%s" % self.cat + cl_any
    # def getCatsBin(self):
    #     return 2**self.cat

    def hashCat(self):
        # return self.getCatsBin()
        return hash(self.getCatsStr())

    def toKey(self):
        return (self.colId(), self.typeId(), 1, self.hashCat())

        
    def disp(self, neg=None, names=None, fmts=None):
        return self.getDispNeg(neg) + self.getDispCol(names, fmts) + "=" + self.getCat()
    
    def getDispCatTk(self, neg=None, protect=False):
        negBVal = False
        if type(neg) == bool:
            negBVal = neg
        if type(neg) is Neg:
            negBVal = neg.boolVal()
        return self.sym_tks[negBVal] + DSS.prepare(self.getCat(), "TEXT", protect=protect)
    
    def getDispTk(self, neg=None, names=None, fmts=None, protect=False):
        return DSS.prepare(self.getDispTermTk(self.getDispColTk(names, fmts) + self.getDispCatTk(neg)), protect=protect)


class CatTerm(Term):
    type_id = 2
    type_letter = 'C'
    type_name = 'Categorical'
    
    # DSS.addTermSymbols(type_letter, {"txt": "[", "tex": "[", "utf": "["}, {"txt": "]", "tex": "]", "utf": "]"})

    sym_one_tks = {False: DSS.getFmtKey("CAT_IS"), True: DSS.getFmtKey("CAT_ISNOT")}
    # DSS.addKeySymbols("CAT_IS", {"txt": "=", "tex": "=", "utf": "="})
    # DSS.addKeySymbols("CAT_ISNOT", {"txt": "!=", "tex": "\\neq{}", "utf": SYM.SYM_NEQ})

    sym_multi_tks = {False: DSS.getFmtKey("CAT_IN"), True: DSS.getFmtKey("CAT_NOTIN")}
    DSS.addKeySymbols("CAT_IN", {"txt": " in ", "tex": " \\in ", "utf": " "+ SYM.SYM_EIN+" "})
    DSS.addKeySymbols("CAT_NOTIN", {"txt": " != ", "tex": " \\neq{} ", "utf": " "+SYM.SYM_NIN+" "})
    DSS.addKeySymbols("CAT_SEP", {"txt": ", ", "tex": ", ", "utf": ", "})

    @classmethod
    def is_multi_cat(tcl, cat):
        return type(cat) in [tuple, list, set]
    @classmethod
    def is_int_cat(tcl, cat):
        return (type(cat) is int) or numpy.issubdtype(type(cat), numpy.integer)
    
    def __init__(self, ncol, ncat, vY=None):
        self.col = ncol
        if self.is_multi_cat(ncat):
            self.cat = ncat
        else:
            self.cat = set([ncat])

    def getCat(self):
        return self.cat

    def nbCats(self):
        return len(self.cat)

    def valRange(self):
        return sorted(self.getCat())

    def specifies(self):
        return self.getCat()

    def specValues(self):  # always a list
        return sorted(self.getCat())

    def setRange(self, cat):
        if self.is_multi_cat(cat):
            self.cat = cat
        else:
            self.cat = set([cat])

    def copy(self):
        return CatTerm(self.col, self.cat)

    def truthEval(self, variableV):
        return variableV in self.cat

    def getCatsStr(self, op_curl="{", cl_curl="}", sep=",", op_any="", cl_any=""):
        if len(self.cat) == 1:
            return op_any + sep.join(["%s" % c for c in self.cat]) + cl_any
        else:
            return op_curl + op_any + sep.join(sorted(["%s" % c for c in self.cat])) + cl_any + cl_curl
    # def getCatsBin(self):
    #     bincat = 0
    #     for c in self.cat:
    #         bincat += 2**c
    #     return bincat

    def hashCat(self):
        # return self.getCatsBin()
        return hash(self.getCatsStr())

    def toKey(self):
        return (self.colId(), self.typeId(), self.nbCats(), self.hashCat())

    def disp(self, neg=None, names=None, fmts=None):
        return self.getDispNeg(neg) + self.getDispCol(names, fmts) + "=" + self.getCatsStr()
    
    def getDispCatTk(self, neg=None, protect=False):
        negBVal = False
        if type(neg) == bool:
            negBVal = neg
        if type(neg) is Neg:
            negBVal = neg.boolVal()
        if len(self.cat) == 1:
            cats_tk = DSS.prepare(sorted(["%s" % c for c in self.cat]), encl_each="TEXT", protect=protect)
            return self.sym_one_tks[negBVal] + cats_tk
        else:
            cats_tk = DSS.prepare(sorted(["%s" % c for c in self.cat]), encl_all="CURL", sep="CAT_SEP", encl_each="TEXT", protect=protect)
            return self.sym_multi_tks[negBVal] + cats_tk
    
    def getDispTk(self, neg=None, names=None, fmts=None, protect=False):
        return DSS.prepare(self.getDispTermTk(self.getDispColTk(names, fmts) + self.getDispCatTk(neg)), protect=protect)

        
class NumTerm(Term):
    type_id = 3
    type_letter = 'N'
    type_name = 'Numerical'
    
    bound_styles_trimm = ["tex"]
    bound_styles_adjprec = ["tex"]

    sym_leq_tk = DSS.getFmtKey("NUM_LEQ")
    DSS.addKeySymbols("NUM_LEQ", {"txt": "<", "tex": "\\leq{}", "utf": " "+ SYM.SYM_LEQ+" "})

    
    def __init__(self, ncol, nlowb, nupb):
        if numpy.isinf(nlowb) and numpy.isinf(nupb) or nlowb > nupb:
            raise Warning('Unbounded numerical term !')
        self.col = ncol
        self.lowb = nlowb
        self.upb = nupb

    def getCanonical(self, is_neg=False):
        # term with negation and complement information
        if numpy.isinf(self.lowb):
            return (NumTerm(self.col, self.upb, float("Inf")), not is_neg, True)
        return (self, is_neg, False)

    def getComplement(self):
        if numpy.isinf(self.lowb):
            return NumTerm(self.col, self.upb, float("Inf"))
        elif numpy.isinf(self.upb):
            return NumTerm(self.col, float("-Inf"), self.lowb)
        return None

    def isComplement(self, term):
        if term.type_id == self.type_id:
            return (numpy.isinf(self.lowb) and numpy.isinf(term.upb) and term.lowb == self.upb) or \
                   (numpy.isinf(term.lowb) and numpy.isinf(self.upb) and self.lowb == term.upb)
        else:
            return False
#     def simple(self, neg):
#         if neg:
#             if self.lowb == float('-Inf'):
#                 self.lowb = self.upb
#                 self.upb = float('-Inf')
#                 neg = False
#             elif self.upb == float('-Inf'):
#                 self.upb = self.lowb
#                 self.lowb = float('-Inf')
#                 neg = False
#         return neg

    def valRange(self):
        return [self.lowb, self.upb]

    def specifies(self):
        return [self.lowb, self.upb]

    def specValues(self):  # always a list
        return [self.lowb, self.upb]

    def setRange(self, bounds):
        if numpy.isinf(bounds[0]) and numpy.isinf(bounds[1]) or bounds[0] > bounds[1]:
            raise Warning('Unbounded numerical term !')
        self.lowb = bounds[0]
        self.upb = bounds[1]

    def copy(self):
        return NumTerm(self.col, self.lowb, self.upb)

    def getUpb(self):
        return self.upb

    def getLowb(self):
        return self.lowb

    def isUpbounded(self):
        return not numpy.isinf(self.upb)

    def isLowbounded(self):
        return not numpy.isinf(self.lowb)

    def truthEval(self, variableV):
        return self.lowb <= variableV and variableV <= self.upb

    def toKey(self):
        return (self.colId(), self.typeId(), self.lowb, self.upb)

    def prepareBoundValue(self, low=True, details={}, style=None):
        if low:
            notInf = self.lowb > float('-Inf')
            val = float(self.lowb)
        else:
            notInf = self.upb < float('Inf')
            val = float(self.upb)
        bd = ""
        if notInf:
            if details.get("as_time", False):
                bd = "@%s" % TimeTools.format_time(val, time_prec=details.get("fmt", {}).get("time_prec"))

            else:
                if details.get("adjprec", False):
                    bd = ('%'+details.get("fmt", {}).get("prec", "")+'f') % val
                else:
                    bd = "%s" % val
                    
                if details.get("trimm", False):
                    bd = bd.rstrip("0").rstrip(".")
        return bd, val, notInf
    
    def getDispBound(self, low=True, details={}):
        bd, val, notInf = self.prepareBoundValue(low, details)
        if notInf:
            if low:
                return bd +'<'
            return '<' + bd
        return ''
    
    def getDispBoundTk(self, low=True, details={}):
        bd, val, notInf = self.prepareBoundValue(low, details)
        if notInf:
            if low:
                return bd + self.sym_leq_tk
            return self.sym_leq_tk + bd
        return ''

    def disp(self, neg=None, names=None, fmts=None):
        details = {"as_time": TimeTools.asTimeVar(self.col, names, fmts)}
                   # "trimm": style in self.bound_styles_trimm,
                   # "adjprec": style in self.bound_styles_adjprec}
        if fmts is not None:
            details["fmt"] = getFmtCol(self.col, fmts)
               
        return self.getDispNeg(neg) + self.getDispBound(True, details) + self.getDispCol(names, fmts) + self.getDispBound(False, details)
    
    def getDispTk(self, neg=None, names=None, fmts=None, protect=False):
        details = {"as_time": TimeTools.asTimeVar(self.col, names, fmts)}
                   # "trimm": style in self.bound_styles_trimm,
                   # "adjprec": style in self.bound_styles_adjprec}
        if fmts is not None:
            details["fmt"] = getFmtCol(self.col, fmts)

        return DSS.prepare(self.getDispNegTk(neg) + self.getDispTermTk(self.getDispBoundTk(True, details) + self.getDispColTk(names, fmts) + self.getDispBoundTk(False, details)), protect=protect)    

# constant term, always evaluating to True if non-zero odd number of bangs, False otherwise
class ConstantTerm(Term):

    pattVName = "!"
    type_id = -1
    type_letter = '!'
    type_name = '!'

    def __init__(self, bangs):
        self.col = -1
        self.bangs = bangs

    def isAnon(self):
        return False

    def setTypeId(self, type_id):
        self.type_id = type_id

    def typeId(self):
        return self.type_id

    def copy(self):
        return ConstantTerm(self.bangs)

    def getAdjusted(self, bounds):
        tmp = self.copy()
        return tmp

    def hashBangs(self):
        return  hash(self.bangs)
    
    def toKey(self):
        return (-99.99, -1, self.hashBangs(), -1)

    def truthEval(self, variableV=None):
        if type(variableV) is bool:
            return variableV
        return len(self.bangs) > 0 and (len(self.bangs) % 2 == 0)

    def __str__(self):
        return self.bangs
    
    def disp(self, neg=None, names=None, fmts=None):
        return str(self)

    def getDispTk(self, neg=None, names=None, fmts=None, protect=False):
        return DSS.prepare(str(self), protect=protect)

class AnonTerm(Term):

    pattVName = "?"+VARIABLE_MARK+"%d"
    type_id = -1
    type_letter = '?'
    type_name = '?'

    types_classes = {}
    for c in Term.__subclasses__():
        types_classes[c.type_id] = c

    def __init__(self, ncol, type_id=None, vY=None):
        self.col = ncol
        if type_id is not None:
            self.type_id = type_id
        else:
            self.type_id = AnonTerm.type_id

    def isAnon(self):
        return True

    def isXpr(self):
        return False

    def setTypeId(self, type_id):
        self.type_id = type_id

    def typeId(self):
        return self.type_id

    def copy(self):
        return AnonTerm(self.col, self.type_id)

    def getAdjusted(self, bounds):
        tmp = None
        if self.typeId() in self.types_classes:
            tmp = self.types_classes[self.typeId()](self.colId(), 0, 1)
            tmp.setRange(bounds)
        return tmp

    def toKey(self):
        return (self.colId(), self.typeId(), hash("??"), -1)

    def truthEval(self, variableV):
        return False

    def getDispCol(self, names=None, fmts=None):
        if type(names) == list and len(names) > 0:
            return "?"+getNameCol(self.col, names)
        else:
            return AnonTerm.pattVName % self.col
        
    def getDispColTk(self, names=None, fmts=None, protect=False):
        if type(names) == list and len(names) > 0:
            return DSS.prepare("?"+getNameCol(self.col, names), "TEXT", protect=protect)
        else:
            return DSS.prepare(AnonTerm.pattVName % self.col, protect=protect)


class XprTerm(AnonTerm):

    pattVName = XPR_MARK+" %s"
    type_id = -1
    type_letter = 'X'
    type_name = 'Xpr'

    def __init__(self, xpr, val_col=None, name_col=None):
        self.val_col = val_col
        self.name_col = name_col
        self.xpr = xpr
        if val_col is not None:
            self.col = val_col.getId()
            self.type_id = val_col.typeId()
        else:
            self.col = -1
            self.type_id = XprTerm.type_id
            
    def isAnon(self):
        return True

    def isXpr(self):
        return True

    def copy(self):
        return XprTerm(self.xpr, self.val_col, self.name_col)

    def getName(self):
        return self.name_col

    def getCol(self):
        return self.val_col

    def getXpr(self):
        return self.xpr

    def getAdjusted(self, bounds):
        tmp = self.copy()
        return tmp

    def toKey(self):
        return (-99.99, self.typeId(), hash(self.xpr), -1)

    def truthEval(self, variableV):
        return False

    def __str__(self):
        tmp = XprTerm.pattVName % self.xpr
        if self.name_col is not None:
            tmp += " # %s" % self.name_col
        return tmp
    
    def disp(self, neg=None, names=None, fmts=None):
        return str(self)

    def getDispTk(self, neg=None, names=None, fmts=None, protect=False):
        return DSS.prepare(str(self), protect=protect)

    @classmethod
    def hasXpr(tcl, xxpr):
        return re.match(XPR_MARK, xxpr)

    @classmethod
    def parseXpr(tcl, xxpr, side, data):
        parts = xxpr.split("#")

        xx = parts[0].lstrip(XPR_MARK).strip()
        if len(xx) == 0:
            return None

        xname = None
        if len(parts) > 1:
            xname = re.sub(" ", "_", "#".join(parts[1:]).strip())

        xpr = re.sub(VARIABLE_MARK+"([0-9]+)", "row.getValue(side,\\1)", xx)
        vals = []
        for rid, row in enumerate(data.getRows()):
            try:
                v = eval(xpr, {}, {"row": row, "side": side, "rid": rid, "numpy": numpy})
            except:
                v = None
            vals.append(v)
        distinct_values = set(vals)
        col = None
        if len(distinct_values.difference([None])) > 0:
            if len(distinct_values.difference([True, False, None])) == 0:
                col = data.getColClassForName("Boolean").parseList(vals)
            else:
                try:
                    col = data.getColClassForName("Numerical").parseList(vals)
                except ValueError:
                    try:
                        col = data.getColClassForName("Categorical").parseList(vals)
                    except ValueError:
                        col = None
        if col is not None:
            return tcl(xx, col, xname)
        return None


class Literal(object):

    # types ordered for parsing
    termTypes = [{'class': NumTerm},
                 {'class': CatTerm},
                 {'class': BoolTerm}]

    @classmethod
    def fromString(tcl, part, names=None, ids_map=None):
        parser = RedQueryParser(parseinfo=False, variable_mark=VARIABLE_MARK)
        qs = QuerySemantics(names, ids_map)
        try:
            p = parser.parse(part, "literal", semantics=qs)
            if len(p) == 1:
                k, l = p.popitem()
                return l[0]
        except FailedParse as e:
            raise Exception("Failed parsing literal %s!\n\t%s" % (part, e))

    def __init__(self, nneg, nterm=None, ids_map=None):
        if isinstance(nterm, Term):
            self.term = nterm  # Already an Term instance
            self.neg = Neg(nneg)
        else:  # nneg is a string representing the literal, nterm is none or contains var names
            l = Literal.fromString(nneg, nterm, ids_map)
            if l is not None:
                self.term = l.term
                self.neg = l.neg
        self.key = None

    def copy(self):
        return Literal(self.getNeg().boolVal(), self.term.copy())

    def isAnon(self):
        return self.term.isAnon()

    def isXpr(self):
        return self.term.isAnon() and self.term.isXpr()

    def getAdjusted(self, bounds):
        return Literal(self.getNeg().boolVal(), self.term.getAdjusted(bounds))

    def valRange(self):
        if self.typeId() == 1 and self.isNeg():
            return [False, False]
        else:
            return self.getTerm().valRange()

    def specValues(self):
        return self.getTerm().specValues()

    def specifies(self):
        return self.getTerm().specifies()

    def __repr__(self):
        return "%s(\"%s\")" % (self.__class__.__name__, self.disp())

    def __str__(self):
        return self.disp()

    def disp(self, names=None, fmts=None):
        return self.getTerm().disp(self.getNeg(), names, fmts)

    def getDispTk(self, names=None, fmts=None, protect=False):
        return self.getTerm().getDispTk(self.getNeg(), names, fmts, protect)
    
    def styledDisp(self, names=None, fmts=None, style="", protect=True, str_before="", str_after=""):
        return self.getTerm().styledDisp(self.getNeg(), names, fmts, style, protect, str_before, str_after)
    
    def tInfo(self, names=None):
        return self.getTerm().tInfo(names)

    def toKey(self, side=None):
        if side is not None:
            return (side,)+self.toKey()
        if self.key is None:
            self.key = self.getTerm().toKey()+self.getNeg().toKey()
        return self.key
    def __hash__(self):
        return hash((self.getTerm().toKey(), self.getNeg().toKey()))

    def __eq__(self, other):
        return isinstance(other, Literal) and self.toKey() == other.toKey()

    def __ne__(self, other):
        return not isinstance(other, Literal) or self.toKey() != other.toKey()

    def __lt__(self, other):
        return isinstance(other, Literal) and self.toKey() < other.toKey()

    def __le__(self, other):
        return isinstance(other, Literal) and self.toKey() <= other.toKey()

    def __gt__(self, other):
        return not isinstance(other, Literal) or self.toKey() > other.toKey()

    def __ge__(self, other):
        return not isinstance(other, Literal) or self.toKey() >= other.toKey()

    def getTerm(self):
        return self.term

    def getNeg(self):
        return self.neg

    def getCanonical(self):
        return self.term.getCanonical(self.isNeg())

    def colId(self):
        return self.getTerm().colId()

    def replaceCol(self, cid):
        self.getTerm().replaceCol(cid)

    
    def typeId(self):
        return self.getTerm().typeId()

    def isNeg(self):
        return self.getNeg().boolVal()

    def setNeg(self, neg):
        self.key = None
        self.neg = Neg(neg)

    def flip(self):
        self.key = None
        self.neg.flip()

    def cmpFlip(self, term):
        if other is None or not isinstance(other, Literal):
            return 1
        elif cmp_vals(self.getTerm(), other.getTerm()) == 0:
            return 1-cmp_vals(self.getNeg(), other.getNeg())
        elif self.getTerm().isComplement(other.getTerm()) and self.getNeg() == other.getNeg():
            return 0
        else:
            return cmp_vals(self.getTerm(), other.getTerm())

    def truthEval(self, variableV):
        if self.isNeg():
            return not self.getTerm().truthEval(variableV)
        else:
            return self.getTerm().truthEval(variableV)


class QTree(object):

    branchN, branchY = (0, 1)

    def __init__(self, root_id=None, branches=None, fill=False, broken=False):
        self.root_id = root_id
        self.tree = {root_id: {"children": [[], []], "depth": 0}}
        self.leaves = set()
        self.supps = None
        self.Ys = None
        self.Xs = None
        self.max_depth = 0
        self.broken = broken
        if branches is not None:
            self.build(branches)
            if fill:
                self.fill()

    def copy(self):
        cp = QTree(root_id=self.root_id)
        cp.tree = copy.deepcopy(self.tree)
        cp.leaves = set(self.leaves)
        cp.supps = copy.deepcopy(self.supps)
        cp.Ys = copy.deepcopy(self.Ys)
        cp.Xs = copy.deepcopy(self.Xs)
        cp.max_depth = self.max_depth
        return cp

    def __str__(self):
        return self.dispNode(self.root_id)

    def isBroken(self):
        return self.broken

    def dispNode(self, nid):
        if self.isBroken():
            return "Broken Tree"
        strn = ""
        if "split" in self.tree[nid]:
            if self.supps is not None:
                sc = self.supps.get(nid, [])
                st = sum(sc)
                suppsn = " + ".join(map(str, self.supps.get(nid, []))) + (" = %d" % st)
            else:
                suppsn = ""
            if self.Ys is not None:
                yn = str(self.Ys.get(nid, ""))
            else:
                yn = ""
            strn += "%s'-- (%d) %d:N %s%s[%s]\t%s\n" % ("\t" * self.tree[nid]["depth"], nid,
                                                        self.tree[nid]["ynb"], self.tree[nid]["split"],
                                                        "\t" * self.max_depth, suppsn, yn)
        if "children" in self.tree[nid]:
            for ynb in [0, 1]:
                for node in self.tree[nid]["children"][ynb]:
                    strn += self.dispNode(node)
        else:
            if self.supps is not None:
                sc = [len(s) for s in self.supps.get((nid, "L"), [])]
                st = sum(sc)
                suppsn = " + ".join(map(str, self.supps.get(nid, []))) + (" = %d" % st)
            else:
                suppsn = ""
            if self.Ys is not None:
                yn = str(self.Ys.get(nid, "")) + "---" + str(self.Ys.get((nid, "L"), ""))
            else:
                yn = ""
            try:
                strn += "%s'-- (%s) %d:L %s%s[%s]\t%s\n" % ("\t" * self.tree[nid]["depth"], nid,
                                                            self.tree[nid]["ynb"], self.tree[nid]["leaf"],
                                                            "\t" * self.max_depth, suppsn, yn)
            except TypeError:
                pdb.set_trace()
                print(self.max_depth)
        return strn

    def hasNode(self, node):
        return node in self.tree

    def isRootNode(self, node):
        return node == self.root_id

    def isLeafNode(self, node):
        return node in self.tree and "leaf" in self.tree[node]

    def isLeafInNode(self, node):
        return node in self.tree and self.tree[node].get("leaf", None) >= 0

    def isLeafOutNode(self, node):
        return node in self.tree and self.tree[node].get("leaf", None) == -1

    def isSplitNode(self, node):
        return node in self.tree and "split" in self.tree[node]

    def isParentNode(self, node):
        return node in self.tree and "children" in self.tree[node]

    def isEmpty(self):
        return len(self.tree) == 1 and self.isLeafNode(self.root_id)

    def mkEmpty(self):
        self.tree = {self.root_id: {"leaf": -1, "depth": 0, "ynb": -1, "parent": -1}}
        self.leaves = set([self.root_id])

    def getMaxDepth(self):
        return self.max_depth

    def getLeaves(self):
        return self.leaves

    def getNodeSplit(self, node):
        return self.tree[node]["split"]

    def getNodeBranch(self, node):
        return self.tree[node]["ynb"]

    def getNodeLeaf(self, node):
        return self.tree[node]["leaf"]

    def getNodeParent(self, node):
        return self.tree[node]["parent"]

    def getNodeChildren(self, node, ynb):
        return self.tree[node]["children"][ynb]

    def nbNodeChildren(self, node, ynb):
        return len(self.tree[node]["children"][ynb])

    def trimNodeChildren(self, node, ynb):
        while len(self.tree[node]["children"][ynb]) > 0:
            c = self.tree[node]["children"][ynb].pop()
            if self.isParentNode(c):
                self.trimNodeChildren(c, 0)
                self.trimNodeChildren(c, 1)
            else:
                self.leaves.discard(c)

    def getNodeXY(self, node):
        if self.Ys is not None and self.Xs is not None:
            return (self.Xs[self.tree[node]["depth"]], self.Ys[node])
        return (None, None)

    def getBranchQuery(self, node):
        cn = node
        buk = []
        while self.getNodeParent(cn) is not None:
            prt = self.getNodeParent(cn)
            neg = cn in self.getNodeChildren(prt, self.branchN)
            tmp = self.getNodeSplit(prt)
            if neg and tmp.type_id == NumTerm.type_id and tmp.getComplement() is not None:
                buk.insert(0, Literal(not neg, tmp.getComplement()))
                ## print(neg, tmp, "=>", buk[0])
            else:
                buk.insert(0, Literal(neg, tmp))
                ## print(neg, tmp, "->", buk[0])
            cn = prt
        return buk

    def getQuery(self):
        buks = []
        for node in self.leaves:
            if self.isLeafInNode(node):
                buks.append((node, self.getBranchQuery(node)))
        buks.sort(key=lambda x: (self.getNodeLeaf(x[0]), x[0]))
        qu = Query()
        if len(buks) == 1:
            if len(buks[0][1]) == 1:
                qu.op = Op(0)
            else:
                qu.op = Op(-1)
            qu.buk = buks[0][1]
        else:
            qu.op = Op(1)
            qu.buk = []
            for x in buks:
                if len(x[1]) == 1:
                    qu.buk.append(x[1][0])
                else:
                    qu.buk.append(x[1])
        return qu

    def getSimpleQuery(self):
        cp = self.copy()
        cp.recSimply(cp.root_id)
        return cp.getQuery().toTree().getQuery()

    def recSimply(self, node):
        if self.isLeafNode(node):
            if self.isLeafInNode(node):
                return 1
            else:
                return -1

        else:
            pr = [None, None]
            for ynb in [0, 1]:
                if ynb == QTree.branchY or not self.isRootNode(node):
                    chlds = self.getNodeChildren(node, ynb)
                    if len(chlds) > 0:
                        tmppr = set([self.recSimply(c) for c in chlds])
                        pr[ynb] = max(tmppr)
                        if pr[ynb] == 1:
                            self.trimNodeChildren(node, ynb)
                            self.addLeafNode(node, ynb, 1)
            if pr[0] == pr[1]:
                return pr[0]
            else:
                return 0

    def getBottomX(self):
        if self.Xs is not None:
            return self.Xs[-1]
        return None

    def getNodeSupps(self, node):
        if self.supps is not None:
            return self.supps[node]
        return None

    def getNodeSuppSets(self, node):
        if self.supps is not None and self.isLeafNode(node):
            return self.supps[(node, "L")]
        return None

    def setNodeLeaf(self, node, bid):
        if node in self.tree and "leaf" in self.tree[node]:
            self.tree[node]["leaf"] = bid

    def setNodeLeaf(self, node, bid):
        if node in self.tree and "leaf" in self.tree[node]:
            self.tree[node]["leaf"] = bid

    def addSplitNode(self, pid, ynb, split):
        if pid in self.tree:
            tid = len(self.tree)
            self.tree[tid] = {"split": split, "children": [[], []], "parent": pid,
                              "depth": self.tree[pid]["depth"]+1, "ynb": ynb}
            self.tree[pid]["children"][ynb].append(tid)
            if self.tree[tid]["depth"] > self.max_depth:
                self.max_depth = self.tree[tid]["depth"]
            return tid

    def addLeafNode(self, pid, ynb, bid):
        if pid in self.tree:
            tid = len(self.tree)
            self.tree[tid] = {"leaf": bid, "parent": pid,
                              "depth": self.tree[pid]["depth"]+1, "ynb": ynb}
            self.tree[pid]["children"][ynb].append(tid)
            self.leaves.add(tid)
            if self.tree[tid]["depth"] > self.max_depth:
                self.max_depth = self.tree[tid]["depth"]
            return tid

    def build(self, branches):
        if len(branches) > 0:
            commons = {}
            for bi, branch in enumerate(branches):
                for li, l in enumerate(branch):
                    if l.getTerm() in commons:
                        key = l.getTerm()
                        cpm = False
                    else:
                        key = l.getTerm().getComplement()
                        cpm = True
                        if key is None or key not in commons:
                            key = l.getTerm()
                            cpm = False
                            commons[key] = [[], []]
                    # is it the yes or no branch?
                    if (not cpm and not l.isNeg()) or (cpm and l.isNeg()):
                        commons[key][self.branchY].append((bi, li))
                    else:
                        commons[key][self.branchN].append((bi, li))
            self.recTree(range(len(branches)), commons, pid=None, fynb=QTree.branchY)
        else:
            self.mkEmpty()

    def recTree(self, bids, commons, pid, fynb):
        mc = max([len(vs[0])+len(vs[1]) for vs in commons.values()])
        kks = [k for (k, vs) in commons.items() if len(vs[0])+len(vs[1]) == mc]
        # TODO choose the split
        # if len(kks) > 1:
        #     pdb.set_trace()
        # pdb.set_trace()
        pick = sorted(kks, key=lambda x: (commons[x], x.toKey()))[0]

        split_commons = [{}, {}, {}]
        to_ynbs = [[v[0] for v in commons[pick][0]], [v[0] for v in commons[pick][1]]]
        to_ynbs.append([bid for bid in bids if bid not in to_ynbs[0] and bid not in to_ynbs[1]])
        for k, vs in commons.items():
            vvs = [[[], []], [[], []], [[], []]]
            if k != pick:
                for yn_org in [0, 1]:
                    for v in vs[yn_org]:
                        if v[0] in to_ynbs[0]:
                            vvs[0][yn_org].append(v)
                        elif v[0] in to_ynbs[1]:
                            vvs[1][yn_org].append(v)
                        else:
                            vvs[2][yn_org].append(v)
            for ynb in [0, 1, 2]:
                if len(vvs[ynb][0])+len(vvs[ynb][1]) > 0:
                    split_commons[ynb][k] = vvs[ynb]

        tid = self.addSplitNode(pid, fynb, pick)
        for ynb in [0, 1]:
            if len(split_commons[ynb]) > 0:
                self.recTree(to_ynbs[ynb], split_commons[ynb], tid, ynb)
            elif len(to_ynbs[ynb]) > 0:
                if len(to_ynbs[ynb]) == 1:
                    bid = to_ynbs[ynb][0]
                else:
                    print("Not exactly one branch ending in %d: %s!" % (tid, to_ynbs[ynb]))
                    bid = None
                self.addLeafNode(tid, ynb, bid)

        if len(split_commons[2]) > 0:
            self.recTree(to_ynbs[2], split_commons[2], pid, fynb)
        elif len(to_ynbs[2]) > 0:
            print("Unexpected something")
            # pdb.set_trace()

    def fill(self):
        basic_nodes = list(self.tree.keys())
        for ni in basic_nodes:
            if "children" in self.tree[ni] and ni != self.root_id:
                for ynb in [0, 1]:
                    if len(self.tree[ni]["children"][ynb]) == 0:
                        self.addLeafNode(ni, ynb, -1)

    def computeSupps(self, side, data, subsets=None):
        self.supps = {}
        if subsets is None:
            subsets = [data.rows()]
        self.recSupps(side, data, self.root_id, subsets)
        # print("SUPPORT %d" % side)
        # print([(n, s) for (n, s) in self.supps.items() if type(n) is int])
        # pdb.set_trace()

    def recSupps(self, side, data, node, subsets):
        self.supps[node] = [len(s) for s in subsets]
        if self.isLeafNode(node):
            self.supps[(node, "L")] = subsets
        else:
            supps_node = [None, None, None]
            if self.isSplitNode(node):
                supp, miss = data.literalSuppMiss(side, self.tree[node]["split"])
                supps_node[QTree.branchY] = [supp & s for s in subsets]
                supps_node[QTree.branchN] = [(s - supp) - miss for s in subsets]
                supps_node[-1] = [s & miss for s in subsets]

            else:
                supps_node[QTree.branchY] = subsets
                supps_node[QTree.branchN] = [set() for s in subsets]
                supps_node[-1] = [set() for s in subsets]

            if self.isParentNode(node):
                for ynb in [0, 1]:
                    # cs = tree[node]["children"][ynb]
                    # if len(cs) == 0 and node is not None:
                    #     supps[(node, "X%d" % ynb)] = supps_node[ynb]

                    for child in self.getNodeChildren(node, ynb):
                        self.recSupps(side, data, child, supps_node[ynb])

    def positionTree(self, side, all_width=1.0, height_inter=[1., 2.]):
        if self.isEmpty():
            self.Xs = [-2*(0.5-side)*1.7]
            self.Ys = {self.root_id: .5*(height_inter[0]+height_inter[1])}
            return

        mdepth = self.getMaxDepth()
        width = all_width/(mdepth+1)
        self.Xs = [-2*(0.5-side)*(i+2)*width for i in range(mdepth+2)][::-1]
        leaves = []
        self.getLeavesYs(side, None, [0., 1.], leaves)
        leaves.sort(key=lambda x: x[-1])
        if len(leaves) < 2:
            width = 0
        else:
            width = (height_inter[1]-height_inter[0])/(len(leaves)-1)
        self.Ys = {}
        for li, leaf in enumerate(leaves):
            self.Ys[leaf[0]] = height_inter[0] + li*width
        self.setYs(side, None)

    def getLeavesYs(self, side, node, interval, leaves):
        if self.isLeafNode(node):
            leaves.append((node, (interval[1]+interval[0])/2.))
        elif self.isParentNode(node):
            eps = (interval[1]-interval[0])/10.
            mid = (interval[1]+interval[0])/2.
            for ynb, subint in [(QTree.branchN, (interval[0]+eps, mid)), (QTree.branchY, (mid, interval[1]-eps))]:
                if self.nbNodeChildren(node, ynb) > 0:
                    width = (subint[1]-subint[0])/self.nbNodeChildren(node, ynb)
                    for ci, child in enumerate(self.getNodeChildren(node, ynb)):
                        self.getLeavesYs(side, child, [subint[0]+ci*width, subint[0]+(ci+1)*width], leaves)

    def setYs(self, side, node):
        if self.isLeafNode(node):
            return self.Ys[node]

        elif self.isParentNode(node):
            ext_p = [[], []]
            for ynb in [0, 1]:
                for ci, child in enumerate(self.getNodeChildren(node, ynb)):
                    ext_p[ynb].append(self.setYs(side, child))
            if node is not None:
                self.Ys[node] = (numpy.max(ext_p[0]) + numpy.min(ext_p[1]))/2.01
            else:
                if len(ext_p[0])+len(ext_p[1]) > 0:
                    self.Ys[node] = numpy.mean(ext_p[0]+ext_p[1])
                else:
                    self.Ys[node] = 1.5
        return self.Ys[node]


class Query(object):
    diff_literals, diff_cols, diff_op, diff_balance, diff_length = range(1, 6)
    side = 0
    class_letter = "q"

    DSS.addKeySymbols("EMPTY_Q", {"txt": "[]", "tex": "", "utf": ""})
    
    # PROPS WHAT
    elem_letters = {"L": "Literals", "T": "Terms", "C": "Cols"}
    info_what = {"len": "len(self)", "depth": "self.max_depth()",
                 "containsOR": "self.usesOr()", "containsAND": "self.usesAnd()",
                 "containsAnon": "self.containsAnon()", "isTreeCompatible": "self.isTreeCompatible()"}
    for ei, el in elem_letters.items():
        info_what[ei+"set"] = "self.inv%s()" % el
        info_what[ei+"nb"] = "len(self.inv%s())" % el
    mtch_contain = "contains(?P<typ>["+"".join(list(elem_letters.keys()))+"])"

    Pwhat_match = "(" + "|".join(["query", mtch_contain]+list(info_what.keys())) + ")"
    @classmethod
    def hasPropWhat(tcl, what):
        return re.match(tcl.Pwhat_match, what) is not None

    @classmethod
    def pieceTogether(tcl, op, pieces, resort=True):
        if len(pieces) == 0:
            return tcl(0, [])
        elif len(pieces) == 1:
            return tcl(0, [pieces[0][1]])

        pp = dict(pieces)
        max_len = max([len(k) for k in pp.keys()])
        by_len = [[] for i in range(max_len+1)]
        for k, l in pp.items():
            by_len[len(k)].append(k)
        if len(by_len) > 2:
            pdb.set_trace()
        while len(by_len) > 1:
            tmp = by_len.pop()
            tmp.sort(key=lambda x: (x[-1] < 0, x[-1]))
            for t in tmp:
                if t[:-1] not in pp:
                    pp[t[:-1]] = []
                    by_len[-1].append(t[:-1])
                pp[t[:-1]].append(pp.pop(t))
        buk = pp[by_len[0][0]]
        q = tcl(op, buk)
        if resort:
            q.doSort()
        return q

    # PROPS WHICH
    Pwhich_match = "([0-9q]*)"
    # which is generic, search element for containment test, i.e. col id, var name
    @classmethod
    def testContainsWhich(tcl, which, eset, typ):
        if typ == "C":
            try:
                return int(which) in eset
            except TypeError:
                pass
        return False

    def __init__(self, OR=0, buk=None):
        self.op = Op(OR)
        if buk is not None:
            self.buk = buk
        else:
            self.buk = []

    def length(self, ex_anon=False):
        if len(self.buk) == 0:
            return 0
        if ex_anon:
            return recurse_numeric(self.buk, function=lambda x: int(isinstance(x, Literal) and not x.isAnon()))
        else:
            return recurse_numeric(self.buk, function=lambda x: int(isinstance(x, Literal)))

    def __len__(self):
        return self.length()

    def __hash__(self):
        if len(self) == 0:
            return 0
        return hash(tuple([self.op.toKey()]+[(i, e.toKey()) for (i,e) in sorted(self.indsLit())]))
        # return hash(self.op) + recurse_numeric(self.buk, function=lambda x, trace: hash(t)+sum(trace), args={"trace": []})

    def getProp(self, what, which=None, details={}):
        if what == "query":
            return self.prepareQuery(details)
        if what in self.info_what:
            return eval(self.info_what[what])
        tc = re.match(self.mtch_contain, what)
        if tc is not None and which is not None:
            typ = tc.group("typ")
            el = Query.elem_letters.get(typ)
            if el is not None:
                eset = eval("self.inv%s()" % el)
                return Query.testContainsWhich(which, eset, typ)
        return None

    def max_depth(self):  # Does the query involve some disjunction?
        if len(self) == 0:
            return 0
        return max(recurse_list(self.buk, function=lambda x, trace: len(trace), args={"trace": []}))

    def getOp(self):
        return self.op

    def usesOr(self):  # Does the query involve some disjunction?
        max_d = self.max_depth()
        return max_d > 1 or (len(self) > 1 and self.op.isOr())

    def usesAnd(self):  # Does the query involve some conjunction?
        max_d = self.max_depth()
        return max_d > 1 or (len(self) > 1 and not self.op.isOr())

    def getOuterOp(self):
        return self.op

    def opBuk(self, nb):  # get operator for bucket nb (need not exist yet).
        if nb % 2 == 0:  # even bucket: query operator, else other
            return self.op.copy()
        else:
            return self.op.other()

    def getBukElemAtR(self, path, buk=None, i=None):  # starting pos from end of path
        if i is None:
            i = len(path)-1
        if buk is None:
            buk = self.buk
        if path[i] < len(buk):
            if i == 0:
                return buk[path[i]]
            else:
                return self.getBukElemAt(path, buk[path[i]], i-1)
        return None

    def getBukElemAt(self, path, buk=None, i=None):  # starting pos from start of path
        if i is None:
            i = 0
        if buk is None:
            buk = self.buk
        if type(buk) is list:
            if len(path) == 0:
                return buk
            elif path[i] < len(buk):
                if i == len(path)-1:
                    return buk[path[i]]
                else:
                    return self.getBukElemAt(path, buk[path[i]], i+1)
        return None

    def setBukElemAtR(self, newE, path, buk=None, i=None):  # starting pos from end of path
        if i is None:
            i = len(path)-1
        if buk is None:
            buk = self.buk
        if len(path) == 0:
            return buk
        if path[i] < len(buk):
            if i == 0:
                tmp = buk[path[i]]
                buk[path[i]] = newE
                return tmp
            else:
                return self.setBukElemAt(newE, path, buk[path[i]], i-1)
        return None

    def setBukElemAt(self, newE, path, buk=None, i=None):  # starting pos from start of path
        if i is None:
            i = 0
        if buk is None:
            buk = self.buk
        if type(buk) is list and path[i] < len(buk):
            if i == len(path)-1:
                tmp = buk[path[i]]
                buk[path[i]] = newE
                return tmp
            else:
                return self.setBukElemAt(newE, path, buk[path[i]], i+1)
        return None

    def copy(self):
        c = Query()
        c.op = self.op.copy()
        c.buk = recurse_deep(self.buk, function=lambda x: x.copy())
        return c

    def push_negation(self):
        def evl(b, flip=False):
            if isinstance(b, Literal):
                if flip:
                    b.flip()
                return (False, b)
            else:
                now_flip = False
                neg = [bb for bb in b if isinstance(bb, Neg)]
                if len(neg) == 1:
                    b.remove(neg[0])
                    now_flip = True
                vs = []
                for bb in b:
                    sfliped, res = evl(bb, now_flip ^ flip)
                    if sfliped:
                        vs.extend(res)
                    else:
                        vs.append(res)
                return (now_flip, vs)
        if len(self) == 0:
            return
        sfliped, res = evl(self.buk, False)
        self.buk = res
        if sfliped:
            self.op.flip()
        # print(self)
        # pdb.set_trace()
        # print("-------")

    def flip(self):
        self.negate()

    def negate(self):
        if len(self) == 0:
            return
        neg = [bb for bb in self.buk if isinstance(bb, Neg)]
        if len(neg) == 1:
            self.buk.remove(neg[0])
        else:
            self.buk.insert(0, Neg(True))
        self.push_negation()
        # self.op.flip()
        # recurse_list(self.buk, function =lambda x: x.flip())

    def __eq__(self, other):
        return self.compare(other) == 0

    def compare(self, other):
        if other is None or not isinstance(other, Query):
            return 1
        try:
            if self.op.compare(other.op) == 0 and cmp_reclists(self.buk, other.buk) == 0:
                return 0
        except AttributeError:
            # Such error means the buckets are not identical...
            pass

        if len(self) < len(other):  # nb of literals in the query, shorter better
            return Query.diff_length
        elif len(self) == len(other):
            if len(self.buk) < len(other.buk):  # nb of buckets in the query, shorter better
                return Query.diff_balance
            elif len(self.buk) == len(other.buk):
                if self.op.compare(other.op) > 0:  # operator
                    return Query.diff_op
                elif self.op.compare(other.op) == 0:
                    if self.invCols() > other.invCols():  # literals in the query
                        return Query.diff_cols
                    elif self.invCols() == other.invCols():
                        return Query.diff_literals
                    else:
                        return -Query.diff_cols
                else:
                    return -Query.diff_op
            else:
                return -Query.diff_balance
        else:
            return -Query.diff_length

    @classmethod
    def comparePair(tcl, x0, x1, y0, y1):  # combined compare for pair
        if (x0.op.compare(y0.op) == 0 and x0.buk == y0.buk and x1.op.compare(y1.op) == 0 and x1.buk == y1.buk):
            return 0

        if len(x0) + len(x1) < len(y0) + len(y1):  # nb of terms in the query, shorter better
            return Query.diff_length

        elif len(x0) + len(x1) == len(y0) + len(y1):
            if len(x0.buk) + len(x1.buk) < len(y0.buk) + len(y1.buk):  # nb of sets of terms in the query, shorter better
                return Query.diff_balance
            elif len(x0.buk) + len(x1.buk) == len(y0.buk) + len(y1.buk):
                if max(len(x0), len(x1)) < max(len(y0), len(y1)):  # balance of the nb of terms in the query, more balanced is better
                    return Query.diff_balance
                elif max(len(x0), len(x1)) == max(len(y0), len(y1)):
                    if max(len(x0.buk), len(x1.buk)) < max(len(y0.buk), len(y1.buk)):  # balance of the nb of sets of terms in the query, more balanced is better
                        return Query.diff_balance

                    elif max(len(x0.buk), len(x1.buk)) == max(len(y0.buk), len(y1.buk)):
                        if x0.op.compare(y0.op) > 0:  # operator on the left
                            return Query.diff_op
                        elif x0.op.compare(y0.op) == 0:
                            if x1.op.compare(y1.op) > 0:  # operator on the right
                                return Query.diff_op
                            elif x1.op.compare(y1.op) == 0:
                                if x0.invCols() > y0.invCols():
                                    return Query.diff_cols
                                elif x0.invCols() == y0.invCols():
                                    if x1.invCols() > y1.invCols():
                                        return Query.diff_cols
                                    elif x1.invCols() == y1.invCols():
                                        return Query.diff_literals
                                return -Query.diff_cols
                        return -Query.diff_op
            return -Query.diff_balance
        return -Query.diff_length


    def replaceCols(self, cols_map={}):
        def replaceCol(lit, cols_map):
            if lit.colId() in cols_map:
                lit.replaceCol(cols_map[lit.colId()])
        recurse_list(self.buk, function=replaceCol, args={"cols_map": cols_map})


    def invCols(self, ex_anon=False):
        def getCol(lit, ex_anon):
            if (not ex_anon or not lit.isAnon()):
                return lit.colId()
        return set(recurse_list(self.buk, function=getCol, args={"ex_anon": ex_anon}))

    def invLiterals(self, ex_anon=False):
        def getLit(lit, ex_anon):
            if (not ex_anon or not lit.isAnon()):
                return lit
        return set(recurse_list(self.buk, function=getLit, args={"ex_anon": ex_anon}))

    def invTerms(self, ex_anon=False):
        def getTerm(lit, ex_anon):
            if (not ex_anon or not lit.isAnon()):
                return lit.getTerm()
        return set(recurse_list(self.buk, function=getTerm, args={"ex_anon": ex_anon}))

    def posLiterals(self, buk=None, op=None, preff=[], mdepth=None):
        if buk is None:
            buk = self.buk
            op = self.op.copy()
            mdepth = self.max_depth()
        if isinstance(buk, Literal):
            return [{"pth": tuple(preff), "mdepth": mdepth, "cid": buk.colId(), "op": op, "lit": buk}]
        elif isinstance(buk, Neg):
            return [{"pth": tuple(preff), "mdepth": mdepth, "cid": -1, "op": op, "lit": buk.boolVal()}]
        tmp = []
        for i in range(len(buk)):
            opp = op.copy()
            opp.flip()
            tmp.extend(self.posLiterals(buk[i], opp, [i]+preff, mdepth))
        return tmp

    def replace(self, depth_ind, replacement, init=None):
        if len(depth_ind) == 0:
            return replacement
        inr = depth_ind.pop()
        if init is None:
            init = self.buk
        src = init
        for i in depth_ind:
            src = src[i]
        tmp = [l for l in src]
        if replacement is None or (type(replacement) is list and len(replacement) == 0):
            tmp.pop(inr)
        else:
            tmp[inr] = replacement
        return self.replace(depth_ind, tmp, init=init)

    def indsLit(self, depth_ind=[], current_el=None, only_anon=False):
        if current_el is None and len(depth_ind) == 0:
            current_el = self.buk
        results = []
        if type(current_el) is list:
            for ni, next_el in enumerate(current_el):
                results.extend(self.indsLit(depth_ind+[ni], next_el, only_anon))
        elif isinstance(current_el, Literal) and (not only_anon or current_el.isAnon()):
            results.append((tuple(depth_ind), current_el))
        return results

    def minusInd(self, depth_ind):
        qq = Query(self.op.isOr(), self.replace(list(depth_ind), None))
        qq.unfold()
        return qq

    def minusInds(self, depth_inds):
        src = self.buk
        for ind in sorted(depth_inds, reverse=True):
            src = self.replace(list(ind), None, src)
        qq = Query(self.op.isOr(), src)
        qq.unfold()
        return qq

    def minusOneLiteral(self):
        xps = []
        for ind, lit in self.indsLit():
            q = self.minusInd(ind)
            xps.append((ind, q))
        return xps

    def minusOneLiteralSupps(self, side, data=None, restrict=None, sm_lits=None):
        xps = []
        if sm_lits is None:
            sm_lits = {}
        for ind, lit in self.indsLit():
            q = self.minusInd(ind)
            dt = {"ind": ind, "q": q, "lit": lit, "lk": lit.toKey(side)}
            if data is not None:
                dt["sm_q"] = q.recompute(side, data, restrict, sm_lits)
                dt["sm_lit"] = sm_lits.get(dt["lk"])
            xps.append(dt)
        return xps

    def minusAnon(self, keep_inds=None):
        q, dropped, drop_inds = self, [], []
        if self.containsAnon():
            inds = self.indsLit()
            for (ind, lit) in inds:
                if lit.isAnon() or (keep_inds is not None and ind in keep_inds):
                    drop_inds.append(ind)
                    dropped.append((ind, lit))
            if len(drop_inds) == len(inds):
                q = Query()
            elif len(drop_inds) > 0:
                q = self.minusInds(drop_inds)
        return q, dropped

    def minusAnonButOne(self):
        dropped, drop_inds = [], []
        if self.containsAnon():
            inds = self.indsLit()
            for (ind, lit) in inds:
                if lit.isAnon():
                    drop_inds.append(ind)
                    dropped.append((ind, lit))
        xps = []
        for i in range(len(drop_inds)):
            mq = self.minusInds(drop_inds[:i]+drop_inds[i+1:])
            minds = mq.indsLit(only_anon=True)
            if len(minds) == 1:
                xps.append((mq, dropped[i][0], dropped[i][1], minds[0][0]))
            else:
                pdb.set_trace()
        return xps

    def unfoldRec(self, buk):
        if isinstance(buk, Query):
            return buk.unfoldRec(buk.buk)

        if isinstance(buk, Literal):
            return buk, False
        tmp = []
        for bi, bb in enumerate(buk):
            tpb, fl = self.unfoldRec(bb)
            if fl is not None:
                if fl:
                    tmp.extend(tpb)
                else:
                    tmp.append(tpb)
        if type(tmp) is list:
            if len(tmp) == 1:
                return tmp[0], not isinstance(tmp[0], Literal)
            elif len(tmp) == 0:
                return tmp, None
        return tmp, False

    def unfold(self):
        new_buk, opflip = self.unfoldRec(self.buk)
        if isinstance(new_buk, Literal):
            self.buk = [new_buk]
        else:
            self.buk = new_buk
        if len(self.buk) < 2:
            self.op = Op()
        elif opflip:
            self.op.flip()

    def partsOCD(self, depth_ind):
        op = self.opBuk(len(depth_ind)-1)
        conj_part = []
        disj_part = []
        if len(depth_ind) > 1:
            up_to = len(depth_ind)-1
            if op.isOr():
                disj_part = list(self.getBukElemAt(depth_ind[:up_to]))
                disj_part.pop(depth_ind[up_to])
                up_to -= 1

            while up_to >= 0:
                c = list(self.getBukElemAt(depth_ind[:up_to]))
                x = c.pop(depth_ind[up_to])
                conj_part.extend(c)
                up_to -= 2

            qc = Query(False, conj_part)
            qd = Query(self.op.isOr(), self.replace(list(depth_ind[:-1]), disj_part))
        else:
            new_buk = list(self.buk)
            if len(depth_ind) > 0:
                new_buk.pop(depth_ind[0])
            if self.op.isOr():
                disj_part = new_buk
            elif self.op.isAnd():
                conj_part = new_buk
            qc = Query(self.op.isOr(), conj_part)
            qd = Query(self.op.isOr(), disj_part)
        return op, qc, qd

    def makeIndexesNew(self, format_str):
        if len(self) == 0:
            return ""
        return recurse_list(self.reorderedLits()[1], function=lambda term, trace: format_str % {'col': term.colId(), 'buk': ".".join(map(str, trace))}, args={"trace": []})

    def makeIndexes(self, format_str):
        if len(self) == 0:
            return ""
        return recurse_list(self.reorderedLits()[1], function=lambda term, trace: format_str % {'col': term.colId(), 'buk': len(trace)}, args={"trace": []})

    def isTreeCompatible(self):
        # a DNF with several conjunctions
        # or a single disjunction
        # or a single conjunction
        # or a single literal
        cp = self.copy()
        cp.push_negation()
        branches = cp.collectTreeBranches([])
        return branches is not None

    def collectTreeBranches(self, bl=None):
        if bl is None:
            bl = self.buk
        branches = []
        if len(self) == 0:
            pass
        elif self.max_depth() == 2 and self.op.isOr():
            for buk in bl:
                if type(buk) is list:
                    branches.append(list(buk))
                else:
                    branches.append([buk])
        elif self.max_depth() == 1 and self.op.isOr():
            for buk in bl:
                branches.append([buk])
        elif (self.max_depth() == 1 and not self.op.isOr()) or \
                (len(self.buk) == 1 and isinstance(self.buk[0], Literal)):
            branches.append(list(bl))
        else:
            return None
        return branches

    def toTree(self, fill=False):
        broken = False
        cp = self.copy()
        cp.push_negation()
        branches = cp.collectTreeBranches()
        if branches is None:
            print("Not a tree form!", self.disp(), self.buk)
            broken = True
            #raise Warning("Not tree form!")
        return QTree(branches=branches, fill=True, broken=broken)

    # return the truth value associated to a configuration

    def truthEval(self, config={}):
        def evl(b, op, config={}):
            if isinstance(b, Literal):
                return b.colId() in config and b.truthEval(config[b.colId()])
            else:
                vs = [evl(bb, op.other(), config) for bb in b]
                if op.isOr():
                    return numpy.sum(vs) > 0
                else:
                    return numpy.prod(vs) > 0
        if len(self) == 0:
            return True
        cp = self.copy()
        cp.push_negation()
        return evl(cp.buk, cp.op, config)

    # return the support associated to a query
    def recompute(self, side, data=None, restrict=None, sm_lits=None):
        def evl(b, op, side, data, restrict=None, sm_lits=None):
            if isinstance(b, Literal):
                if sm_lits is None:
                    sm = data.literalSuppMiss(side, b)
                    if restrict is None:
                        return sm
                    else:
                        return sm[0] & restrict, sm[1] & restrict
                else:
                    bk = b.toKey(side)
                    if bk not in sm_lits:
                        sm = data.literalSuppMiss(side, b)
                        if restrict is not None:
                            sm = (sm[0] & restrict, sm[1] & restrict)
                        sm_lits[bk] = sm
                    return sm_lits[bk]

            else:
                vs = [evl(bb, op.other(), side, data, restrict, sm_lits) for bb in b]
                return SParts.partsSuppMissMass(op.isOr(), vs)

        if len(self) == 0 or data == None or self.containsAnon():
            return (set(), set())
        else:
            if side == -1:
                # in case the query is a condition and the data does not have condition cols, prepare from geo is possible
                if not data.isConditional() and data.isGeospatial():
                    data.prepareGeoCond()
                if not data.isConditional():
                    # impossible to get condition
                    return (set(), set())

            cp = self.copy()
            cp.push_negation()
            return evl(cp.buk, cp.op, side, data, restrict, sm_lits)

    def proba(self, side, data=None, restrict=None):
        def evl(b, op, side, data, restrict=None):
            if isinstance(b, Literal):
                if restrict is None:
                    return len(data.supp(side, b))/float(data.nbRows())
                else:
                    return len(data.supp(side, b) & restrict)/float(len(restrict))
            else:
                vs = [evl(bb, op.other(), side, data, restrict) for bb in b]
                return SParts.updateProbaMass(vs, op.isOr())

        if data is None or self.containsAnon():
            pr = -1
        elif len(self) == 0:
            pr = 1
        else:
            cp = self.copy()
            cp.push_negation()
            pr = evl(cp.buk, cp.op, side, data, restrict)
        return pr

    def probaME(self, dbPr=None, side=None, epsilon=0):
        def evl(b, op, dbPr, side, epsilon):
            if isinstance(b, Literal):
                return dbPr.pointPrLiteral(side, literal, epsilon)
            else:
                vs = [evl(bb, op.other(), dbPr, side, epsilon) for bb in b]
                return SParts.updateProbaMass(vs, op.isOr())

        if dbPr is None:
            pr = -1
        elif len(self) == 0:
            pr = 1
        else:
            cp = self.copy()
            cp.push_negation()
            pr = evl(cp.buk, cp.op, dbPr, side, epsilon)
        return pr

    # RESORT TODO FOR DEBUGGING
    def extend(self, op, literal, resort=True):
        if len(self) == 0:
            self.append(literal)
        elif len(self) == 1:
            self.append(literal)
            if type(op) is bool:
                op = Op(op)
            self.op = op
        elif op == self.op:
            self.append(literal)
        else:
            self.op = self.op.other()
            self.buk = [self.buk, literal]
        if resort:
            self.doSort()

    def append(self, literal):
        self.buk.append(literal)

    def asDisLit(self):
        if not self.op.isOr():
            if self.op.isAnd():
                self.op.flip()
                self.buk = [self.buk]
            else:
                return self.buk[0]
        return self

    def appendBuk(self, buk, op=None, resort=False):
        bid = None
        if op is None:
            op = Op(1)
        if len(self) == 0:
            bid = len(self.buk)
            self.buk.extend(buk)
            self.op = op.other()
        elif len(self) == 1 and self.buk != buk:
            bid = 1
            self.buk = [self.buk, buk]
            self.op = op
        elif self.op == op:
            if buk not in self.buk:
                bid = len(self.buk)
                self.buk.append(buk)
        else:
            if self.buk != buk:
                bid = 1
                self.op = self.op.other()
                self.buk = [self.buk, buk]
        if resort:
            self.doSort()
            bid = None
        return bid

    def doSort(self):
        def soK(x):
            if type(x) is list:
                return -1
            else:
                return x.colId()
        self.buk.sort(key=lambda x: soK(x))

    def getSingleLiteral(self):
        if len(self.buk) == 1 and isinstance(self.buk[0], Literal):
            return self.buk[0]

    def listLiterals(self):
        def evl(b, lits):
            for bb in b:
                if isinstance(bb, Literal):
                    lits.append(bb)
                elif not isinstance(bb, Neg):
                    evl(bb, lits)
        lits = []
        if len(self) > 0:
            evl(self.buk, lits)
        return lits

    def listLiteralsDetails(self):
        def evl(b, lits, path):
            for bi, l in enumerate(b):
                if isinstance(l, Literal):
                    key, neg, cpm = l.getCanonical()
                    if key not in lits:
                        lits[key] = []
                    lits[key].append((tuple(path+[bi]), cpm, neg))
                    # if l.getTerm() in lits:
                    #     key = l.getTerm()
                    #     cpm = False
                    # else:
                    #     key = l.getTerm().getComplement()
                    #     cpm = True
                    #     if key is None or key not in lits:
                    #         key = l.getTerm()
                    #         cpm = False
                    #         lits[key] = []
                    # lits[key].append((tuple(path+[bi]), cpm, not l.isNeg()))
                elif not isinstance(l, Neg):
                    evl(l, lits, path+[bi])
        lits = {}
        path = []
        if len(self) > 0:
            evl(self.buk, lits, path)
        return lits

    def isBasis(self):
        if len(self) == 1:
            ll = self.listLiterals()
            return ll[0].isAnon()
        return False

    def typeId(self):
        if len(self) == 1:
            ll = self.listLiterals()
            return ll[0].typeId()
        return None

    def containsAnon(self):
        if len(self) > 0:
            return any([l.isAnon() for l in self.listLiterals()])
        return False

    def reorderedLits(self, b=None):
        if b is None:
            b = self.buk
        if isinstance(b, Literal):
            return ((b.colId(), b.isNeg(), b.getTerm().valRange()), b)
        elif isinstance(b, Neg):
            return ((-1, ), b)
        else:
            if len(b) == 0:
                return ()
            vs = [self.reorderedLits(bb) for bb in b]
            vs.sort(key=lambda x: x[0])
            return (vs[0][0], [v[1] for v in vs])

    def reorderLits(self):
        if len(self) > 0:
            self.buk = self.reorderedLits()[1]

    @classmethod
    def dispBuk(tcl, buk):
        if isinstance(buk, Literal) or isinstance(buk, Neg):
            return "%s" % buk
        else:
            return "[" + "; ".join([tcl.dispBuk(b) for b in buk])+"]"
        return "-"

    def disp(self, names=None, fmts=None):
        def evl(b, op, names, fmts=None, outer=False):
            if isinstance(b, Literal):
                return b.disp(names=names, fmts=fmts)
            if isinstance(b, Neg):
                return "!NEG!"
            else:
                vs = [evl(bb, op.other(), names, fmts) for bb in b]
                if len(vs) == 1:
                    return vs[0]
                else:
                    jstr = " " + op.disp() + " "
                    pref, suff = ("", "")
                    ### HERE SPACING
                    spc = " "
                    if "!NEG!" in vs:
                        vs.remove("!NEG!")
                        pref, suff = (Neg(True).disp()+ " ("+spc, spc+")")
                    elif not outer:
                        pref, suff = ("("+spc, spc+")")
                    return pref + jstr.join(vs) + suff

        if len(self) == 0:
            return "[]"
        else:
            str_tk = evl(self.buk, self.op, names, fmts, outer=True)
            return str_tk

    def __str__(self):
        return self.disp()
    
    def getDispTk(self, names=None, fmts=None, protect=False):
        def evl(b, op, names, fmts=None, protect=False, outer=False):
            if isinstance(b, Literal):
                return b.getDispTk(names=names, fmts=fmts) #, protect=protect)
            if isinstance(b, Neg):
                return "!NEG!"
            else:
                vs = [evl(bb, op.other(), names, fmts) for bb in b]
                if len(vs) == 1:
                    return DSS.prepare(vs[0], protect=protect and outer)
                else:
                    jstr = op.getDispTk()
                    if "!NEG!" in vs:
                        vs.remove("!NEG!")
                        front_neg = Neg(True).getDispTk()
                        str_tk = DSS.prepare(front_neg + DSS.prepare(vs, encl_all="BRCK", sep=jstr), protect=protect and outer)
                    elif not outer:
                        str_tk = DSS.prepare(vs, encl_all="BRCK", sep=jstr, protect=False)
                    else:
                        str_tk = DSS.prepare(vs, sep=jstr, protect=protect)
                    return str_tk


        if len(self) == 0:
            return DSS.getFmtKey("EMPTY_Q")
        else:
            str_tk = evl(self.buk, self.op, names, fmts, protect=protect, outer=True)
            # str_tk = DSS.cleanEnclosed(str_tk)
            return str_tk
            # old code to write the query justified in length lenField
            #### string.ljust(qstr, lenField)[:lenField]

    def styledDisp(self, names=None, fmts=None, style="", protect=True, str_before="", str_after=""):
        return str_before + DSS.replaceStyleTks(self.getDispTk(names, fmts, protect), style) + str_after

        
    def prepareQuery(self, details={}):
        style = details.get("style")
        side = details.get("side")
        names = None
        fmts = None
        if side is not None:
            if details.get("named", False) and "names" in details:
                names = details["names"][side]
            if "fmts" in details:
                fmts = details["fmts"][side]
        if style is None:
            return self.disp(names=names, fmts=fmts)
        else:
            return self.styledDisp(names=names, fmts=fmts, style=style)

    def algExp(self):
        def evl(b, op, tmap):
            if isinstance(b, Literal):
                if b.getTerm() not in tmap:
                    key = b.getTerm().getComplement()
                    if b.isNeg():
                        return "t[%s]" % tmap[key]
                    else:
                        return "(not t[%s])" % tmap[key]
                elif b.isNeg():
                    return "(not t[%s])" % tmap[b.getTerm()]
                else:
                    return "t[%s]" % tmap[b.getTerm()]
            if isinstance(b, Neg):
                return "!NEG!"
            else:
                vs = [evl(bb, op.other(), tmap) for bb in b]
                if len(vs) == 1:
                    return vs[0]
                else:
                    if op.isOr():
                        jstr = " or "
                    else:
                        jstr = " and "
                    if "!NEG!" in vs:
                        vs.remove("!NEG!")
                        pref = "( not ("
                        suff = "))"
                    else:
                        pref = "( "
                        suff = " )"
                    return pref + jstr.join(vs) + suff

        if len(self) == 0:
            return 'False', {}
        else:
            tmap = dict([(v, i) for (i, v) in enumerate(sorted(self.listLiteralsDetails().keys()))])
            vs = [evl(bb, self.op.other(), tmap) for bb in self.buk]
            if len(vs) == 1:
                return vs[0], tmap
            else:
                if self.op.isOr():
                    jstr = " or "
                else:
                    jstr = " and "
                if "!NEG!" in vs:
                    vs.remove("!NEG!")
                    pref = "( not ("
                    suff = "))"
                else:
                    pref = "( "
                    suff = " )"
                return pref + jstr.join(vs) + suff, tmap

    def truthTable(self):
        # return a map of literals to columns (tmap), and the combinations of truth values of the literal that evaluate to True
        def recTT(lstr, vlist, nbvar):
            if nbvar == 0:
                if eval(lstr, {}, {"t": vlist}) == 1:
                    return [vlist]
                else:
                    return []
            else:
                ####
                return recTT(lstr, [False]+vlist, nbvar-1)+recTT(lstr, [True]+vlist, nbvar-1)

        lstr, tmap = self.algExp()
        tb = recTT(lstr, [], len(tmap))
        return numpy.array(tb, dtype=numpy.int), tmap

    def algNormalized(self):
        tmp = Query()
        if len(self) > 0:
            tt, tmap = self.truthTable()
            stt = tt.copy()
            stt = simplerTT(stt)
            tlist = sorted(tmap.keys(), key=lambda x: tmap[x])
            branches = []
            for bi in range(stt.shape[0]):
                branches.append([Literal(1-stt[bi, ti], tlist[ti]) for ti in range(stt.shape[1]) if stt[bi, ti] != -1])
            if len(branches) > 0:
                if all([len(b) == 0 for b in branches]):  # branches contain only -1 -> always True
                    branches = [[Literal("!!")]]
                # for b in branches:
                #     print([t.disp() for t in b])
                # pdb.set_trace()
                qt = QTree(branches=branches)
                tmp = qt.getQuery()
                # ## TODO! Check correct!
                # if True:  # checking that all went fine
                #     tto, tmapo = tmp.truthTable()
                #     idsc = []
                #     dets = []
                #     for l in tmapo.keys():
                #         if l in tmap:
                #             idsc.append(tmap[l])
                #             dets.append((l, True))
                #         else:
                #             lc = l.getComplement()
                #             if lc is not None and lc in tmap:
                #                 idsc.append(lc)
                #                 dets.append((l, False))
                #     if len(idsc) != len(tmapo):
                #         print("--- OUPS! something went wrong when normalizing %s" % self)
                #         print("---", idsc, tmap, tmapo)
                #     pdb.set_trace()
                #     if all([x[-1] for x in dets]):  # no complement involved
                #         dropC = [i for (l, i) in tmap.items() if l not in tmapo]
                #         ctt = tt
                #         if len(dropC) > 0:
                #             keep_rows = numpy.all(tt[:, dropC] == 1, axis=1)
                #             ctt = tt[keep_rows, :][:, idsc]
                #         if numpy.sum(ctt != tto) > 0:
                #             print("--- OUPS! something went wrong when normalizing %s" % self)
                #             print("---", idsc, tmap, tmapo)

        if len(self) == 0 and len(tmp) == 0:
            return tmp, False
        else:
            tmp.reorderLits()
            comp = self.reorderedLits()[1]
            if comp == tmp.buk and self.op == tmp.op:
                # actually unchanged
                return tmp, False
            else:
                return tmp, True

    def isXpr(self):
        return len(self.buk) == 1 and isinstance(self.buk[0], Literal) and self.buk[0].isXpr()

    def getXprTerm(self):
        if self.isXpr():
            return self.buk[0].getTerm()

    @classmethod
    def hasXpr(tcl, part):
        return XprTerm.hasXpr(part)

    @classmethod
    def parseXpr(tcl, part, side, data):
        if XprTerm.hasXpr(part):
            xpr_term = XprTerm.parseXpr(part, side, data)
            if xpr_term is not None:
                return Query(buk=[Literal(False, xpr_term)])
        return None

    @classmethod
    def parse(tcl, part, names=None, ids_map=None):
        if len(part.strip()) == 0 or part.strip() == "[]":
            return Query()
        qs = QuerySemantics(names, ids_map)
        parser = RedQueryParser(parseinfo=False, variable_mark=VARIABLE_MARK)
        try:
            tmp = parser.parse(part, "query", semantics=qs)
        except FailedParse as e:
            tmp = Query()
            raise Exception("Failed parsing query %s!\n\t%s" % (part, e))
        return tmp

# GENERATE PARSER:
#     python -m grako -m RedQuery -o redquery_parser.py redquery.ebnf; sed -i 's/division, absolute_import, unicode_literals/division, unicode_literals/' redquery_parser.py (!!! REMOVE ABSOLUTE_IMPORT FROM GENERATED FILE)
# RUN:
#     python redquery_parser.py queries.txt QUERIES


class QuerySemantics(object):

    def __init__(self, names=None, ids_map=None):
        self.names = names
        self.ids_map = ids_map

    def query(self, ast):
        buk = []
        OR = 0
        if "conjunction" in ast:
            buk = ast["conjunction"]
            OR = False
        elif "disjunction" in ast:
            buk = ast["disjunction"]
            OR = True
        elif "literal" in ast:
            buk = list(ast["literal"].values())[0]
        if "mass_neg" in ast:
            buk.insert(0, Neg(True))
        return Query(OR, buk)

    def conjunction(self, ast):
        tmp = []
        for e in ast:
            if len(e) == 1:
                tmp.extend(e)
            else:
                tmp.append(e)
        return tmp

    def disjunction(self, ast):
        tmp = []
        for e in ast:
            if len(e) == 1:
                tmp.extend(e)
            else:
                tmp.append(e)
        return tmp

    def conj_item(self, ast):
        if "mass_neg" in ast:
            del ast["mass_neg"]
            return [Neg(True)]+list(ast.values())[0]
        return list(ast.values())[0]

    def disj_item(self, ast):
        if "mass_neg" in ast:
            del ast["mass_neg"]
            return [Neg(True)]+list(ast.values())[0]
        return list(ast.values())[0]

    def categorical_literal(self, ast):
        return [Literal(("neg" in ast) ^ ("cat_false" in ast.get("cat_test", {})),
                        CatTerm(self.parse_vname(ast.get("variable_name")), self.parse_categories(ast.get("categories"))))]

    def realvalued_literal(self, ast):
        if re.match("@", ast.get("lower_bound", "")):
            lower = TimeTools.parse_time(ast.get("lower_bound")[1:])
        else:
            lower = float(ast.get("lower_bound", "-inf"))
        if re.match("@", ast.get("upper_bound", "")):
            upper = TimeTools.parse_time(ast.get("upper_bound")[1:])
        else:
            upper = float(ast.get("upper_bound", "inf"))
        return [Literal("neg" in ast, NumTerm(self.parse_vname(ast.get("variable_name")), lower, upper))]

    def boolean_literal(self, ast):
        return [Literal("neg" in ast, BoolTerm(self.parse_vname(ast.get("variable_name"))))]

    def anonymous_literal(self, ast):
        return [Literal("neg" in ast, AnonTerm(self.parse_vname(ast.get("variable_name"))))]

    def constant_literal(self, ast):
        return [Literal(False, ConstantTerm(ast.get("constant_bangs")))]

    def variable_name(self, ast):
        return ast

    def categories(self, ast):
        return ast

    def parse_categories(self, ast):
        if "category" in ast:
            return set([ast["category"]])
        elif "catlist" in ast:
            return set(ast["catlist"])
        return set()

    def parse_vname(self, vname):
        tmp = re.match(VARIABLE_MARK+"(?P<id>\d+)$", vname)
        if tmp is not None:
            vv = int(tmp.group("id"))
            if self.ids_map is not None:
                return self.ids_map.get(vv, vv)
            return vv
        elif self.names is not None:
            if vname in self.names:
                return self.names.index(vname)
            # print("No match")
            raise Exception("No matching variable")
        else:
            print(vname)
            # pdb.set_trace()
            raise Exception("No variables names provided when needed!")


if __name__ == '__main__':
    # # items = [("disjunction", Op(False)), ("conjunction", Op(True)), ("nonnegation", Neg()), ("negation", Neg(True))]
    # # for lgd, item in items:
    # #     print("Next item: %s" % lgd)
    # #     print(item)
    # #     print(item.styledDisp(style="utf"))
    # #     print(item.styledDisp(style="tex"))
    # #     print(item.styledDisp(style="tex", str_before=">>", str_after="<<"))
    # #     print(item.styledDisp())
    # #     print(item.disp())


    # # items = [("BoolTerm", BoolTerm(1)), ("CatTermONE", CatTermONE(0, "C1")), ("CatTerm", CatTerm(0, ["C1", "C4"])), ("NumTerm", NumTerm(0, 3.523, 6.22))]
    # # items = [("CatTerm", CatTerm(0, ["C1", "C4"])), ("CatTerm", CatTerm(0, "C1"))]
    # # items = [("NumTerm", NumTerm(0, 3.523, 6.22)), ("NumTerm", NumTerm(0, 3.52300, float("Inf")))] , fmts={0: {"prec": ".8"}}
    # items = [("AnonTerm", AnonTerm(0)), ("XprTerm", XprTerm("3 * v4"))]
    # names=["A", "B", "C"]
    # for lgd, item in items:
    #     print("Next item: %s" % lgd)
    #     print(item)
    #     print(item.getDispTk(neg=True, names=names))
    #     print(item.getDispTk(names=names))
    #     print(item.getDispTk(names=names))
    #     print(item.getDispTk())
    #     print(item.disp())

    # exit()

    # # items = [("BoolTerm", BoolTerm(1)), ("CatTermONE", CatTermONE(0, "C1")), ("CatTerm", CatTerm(0, ["C1", "C4"])), ("NumTerm", NumTerm(0, 3.523, 6.22))]
    # # items = [("CatTerm", CatTerm(0, ["C1", "C4"])), ("CatTerm", CatTerm(0, "C1"))]
    # # items = [("NumTerm", NumTerm(0, 3.523, 6.22)), ("NumTerm", NumTerm(0, 3.52300, float("Inf")))] , fmts={0: {"prec": ".8"}}
    # items = [("AnonTerm", AnonTerm(0)), ("XprTerm", XprTerm("3 * v4"))]
    # names=["A", "B", "C"]
    # for lgd, item in items:
    #     print("Next item: %s" % lgd)
    #     print(item)
    #     print(item.styledDisp(neg=True, names=names, style="utf"))
    #     print(item.styledDisp(names=names, style="tex"))
    #     print(item.styledDisp(names=names, style="tex", str_before=">>", str_after="<<"))
    #     print(item.styledDisp())
    #     print(item.disp())
    # exit()


    from classData import Data
    from classQuery import QuerySemantics  # import for Literal instance type test
    import sys
    rep = "/home/egalbrun/short/rsmallZx/"
    data = Data([rep+"data_LHS.csv", rep+"data_RHS.csv", {}, "nan"], "csv")
    qsLHS = QuerySemantics(data.getNames(0))
    qsRHS = QuerySemantics(data.getNames(1))
    parser = RedQueryParser(parseinfo=False)

    with open(rep+"redescriptions.csv") as f:
        header = None
        for line in f:
            if header is None:
                header = line.strip().split("\t")
            elif len(line.strip().split("\t")) >= 3:
                resLHS = parser.parse(line.strip().split("\t")[1], "query", semantics=qsLHS)
                resRHS = parser.parse(line.strip().split("\t")[2], "query", semantics=qsRHS)
                print("----------")
                print(line.strip())
                print(resLHS.disp(names=data.getNames(0)))
                print(resRHS.disp(names=data.getNames(1)))
                
                # print("ORG   :", resLHS, "---", resRHS)
                # print(len(resLHS.recompute(0, data)[0]))
                # print(len(resRHS.recompute(1, data)[0]))

                # cp = resLHS.copy()
                # resLHS.push_negation()
                # print("COPY  :", cp)
                # print("PUSHED:", resLHS)
                # cp.negate()
                # print("NEG   :", cp)
                # print(resLHS.recompute(0, data))
                # print(resRHS.recompute(1, data))
                # pdb.set_trace()
                # print(len(resLHS))
                # print(resLHS.listLiterals())
                # tmp = resLHS.copy()
                # print(tmp)
                # tmp.negate()
                # print(tmp)
                # print(resLHS.styledDisp(names=data.getNames(0), style="utf"), "\t", resRHS.disp())
                # print(resLHS.makeIndexesNew('%(buk)s:%(col)i:'))
                # resLHS.reorderLits()
                # print(resLHS.styledDisp(names=data.getNames(0), style="utf"), "\t", resRHS.disp())
