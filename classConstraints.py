import re
import os.path
import numpy

try:
    import toolXML
    from classAction import ActionFactory
    from classSParts import cmp_lower
    from classPreferencesManager import PreferencesReader
    from csv_reader import getFp, FOLDS_PATT

except ModuleNotFoundError:
    from . import toolXML
    from .classAction import ActionFactory
    from .classSParts import cmp_lower
    from .classPreferencesManager import PreferencesReader
    from .csv_reader import getFp, FOLDS_PATT

import pdb


class Constraints(object):

    @classmethod
    def applyVarsMask(tcl, data, params):
        if data is not None:
            return data.applyDisableMasks(params.get("mask_vars_LHS"), params.get("mask_vars_RHS"), params.get("mask_rows"))

    special_cstrs = {}

    params_sets = {"data": set(["parts_type", "method_pval"])}

    def __init__(self, params, data=None, AR=None, filenames={}):
        self.deps = []
        self.folds = None
        if AR is None:
            AR = params.get("AR")
        self._pv = PreferencesReader.copyParams(params)

        self._pv.update(self.prepareValues(params))
        if filenames is not None:
            self._pv.update(filenames)
        self.resetAR(AR)
        self.resetDataDependent(data, params)

    ## SSETTS_PARAMS = set(["parts_type", "method_pval"])
    def reset(self, params, data=None, AR=-1, dtv=None):
        changed = {}
        self._pv.update(params)
        self._pv.update(self.prepareValues(params))

        if not isinstance(AR, ActionsRegistry):
            AR = params.get("AR", -1)
        changedAR = self.resetAR(AR)
        if changedAR:
            changed["actions_filters"] = True
        if data is not None and (dtv is None or len(self.params_sets["data"].intersection(dtv)) > 0):
            self.resetDataDetails(data)
            changed["reset_all"] = True
        self._pv.update(self.prepareValuesDataDependent(params))
        return changed

    def prepareValues(self, params, vals=None):
        def getNegBool(p):
            return p == "negative"

        def getOpsBool(p):
            return p == "disjunction"

        # preparing query types
        if vals is None:
            vvs = {}
        else:
            vvs = vals
        for side_id in [0, 1]:
            for type_id in [1, 2, 3]:
                kp = "neg_query_s%d_%d" % (side_id, type_id)
                if vals is None or kp in params:
                    vvs[kp] = [getNegBool(v) for v in params.get(kp, [])]

            kp = "ops_query_s%d" % side_id
            if vals is None or kp in params:
                vvs[kp] = [getOpsBool(v) for v in params.get(kp, [])]

        # preparing score coeffs
        if vals is None:
            vvs["score_coeffs"] = {}
        for k in ["impacc", "rel_impacc", "pval_red", "pval_query", "pval_fact"]:
            if vals is None or k in params:
                vvs["score_coeffs"][k] = params.get("score.%s" % k, 0)
        return vvs

    def prepareValuesDataDependent(self, params):
        # scaling support thresholds
        min_itm_c, min_itm_in, min_itm_out = self.scaleSuppParams(params.get("min_itm_c"), params.get("min_itm_in"), params.get("min_itm_out"))
        _, min_fin_in, min_fin_out = self.scaleSuppParams(-1, params.get("min_fin_in"), params.get("min_fin_out"))
        return {"min_itm_c": min_itm_c, "min_itm_in": min_itm_in, "min_itm_out": min_itm_out,
                "min_fin_in": min_fin_in, "min_fin_out": min_fin_out}

    def resetDataDetails(self, data):
        if data is not None:
            self.N = data.nbRows()
            data.resetSSetts(self.getCstr("parts_type"), self.getCstr("method_pval"))
            self.ssetts = data.getSSetts()
            self._pv["parts_type"] = self.ssetts.getTypeParts()
        else:
            self.N = -1
            self.ssetts = None

    def resetDataDependent(self, data, params):
        self.applyVarsMask(data, params)
        if data is not None:
            data.setCompat(self.getCstr("var_compat"))
        self.resetDataDetails(data)
        self._pv.update(self.prepareValuesDataDependent(params))

    def resetAR(self, AR=None):
        ar_fns = []
        for f in self._pv.get("actions_rdefs", "").split(";"):
            ff = f.strip()
            if len(ff) > 0:
                ar_fns.append(ff)

        if AR is None:
            self.AR = ActionsRegistry()
            changed = True
        elif isinstance(AR, ActionsRegistry):
            self.AR = AR
            ar_fns = []
            changed = True
        else:
            changed = False

        if len(ar_fns) > 0:
            self.AR.extend(ar_fns)
            changed = True
        return changed

    def setFolds(self, data):
        fcol = data.getColsByName(FOLDS_PATT)
        if len(fcol) == 1:
            self.folds = data.getFoldsStats(fcol[0][0], fcol[0][1])

    def scaleF(self, f):
        if f == -1 or f is None:
            return -1
        if f >= 1:
            return int(f)
        elif f >= 0 and f < 1 and self.N != 0:
            return int(round(f*self.N))
        return 0

    def scaleSuppParams(self, min_c, min_in=None, min_out=None):
        sc_min_c = self.scaleF(min_c)
        if min_in == -1:
            sc_min_in = sc_min_c
        else:
            sc_min_in = self.scaleF(min_in)
        if min_out == -1:
            sc_min_out = sc_min_in
        else:
            sc_min_out = self.scaleF(min_out)
        return (sc_min_c, sc_min_in, sc_min_out)

    def getSSetts(self):
        return self.ssetts

    def getCstrSimple(self, k):
        return self._pv[k]

    def getCstr(self, k, **kargs):
        if k in self.special_cstrs:
            return eval("self.%s" % self.special_cstrs[k])(**kargs)

        k_bak = k
        if "side" in kargs:
            k += "_s%d" % kargs["side"]
        if "type_id" in kargs:
            k += "_%d" % kargs["type_id"]

        if k in self._pv:
            return self._pv[k]
        else:
            return self._pv.get(k_bak, kargs.get("default"))

    # STATUS TEST TO KNOW WHAT IS ALLOWED
    @classmethod
    def expandDefStatus(tcl):
        return {"init": False, "other_contains_OR": False, "contains_OR": False, "other_type_id": 0, "type_id": 0, "cond": False}

    @classmethod
    def getExpandedStatus(tcl, status=0):
        if type(status) is dict:
            return status
        xpd = tcl.expandDefStatus()
        if cmp_lower(status, 0):
            xpd["init"] = True
        return xpd

    @classmethod
    def isStatusInitStage(tcl, status):
        if type(status) is dict:
            return status.get("init", False)
        return cmp_lower(status, 0)

    @classmethod
    def isStatusCond(tcl, status):
        if type(status) is dict:
            return status.get("cond", False)
        return False

    @classmethod
    def getStatusPair(tcl, col, side, fixTerm):
        status = tcl.expandDefStatus()
        status["init"] = True
        status["type_id"] = col.typeId()
        status["other_type_id"] = fixTerm.typeId()
        return status

    @classmethod
    def getStatusCond(tcl, pair=False):
        status = tcl.expandDefStatus()
        if pair:
            status["init"] = True
        status["cond"] = True
        return status

    @classmethod
    def getStatusRed(tcl, red=None, side=None, force_ops=None):
        if red is not None and side is not None:
            status = tcl.expandDefStatus()
            status["init"] = (red.length(side) == 0)
            status["force_ops"] = force_ops
            status["contains_OR"] = red.usesOr(side)
            status["other_contains_OR"] = red.usesOr(1-side)
            if red.length(1-side) == 1:
                status["other_type_id"] = red.query(1-side).invTerms().pop().typeId()
            if red.length(side) == 1:
                status["type_id"] = red.query(side).invTerms().pop().typeId()
            return status
        return 0

    # special constraints (not just lookup)
    def allw_ops(self, side, currentRStatus=0):
        if self.isStatusCond(currentRStatus):
            return [False]
        if self.isStatusInitStage(currentRStatus):
            return [True]
        else:
            xpd = self.getExpandedStatus(currentRStatus)
            if xpd.get("force_ops") is not None:
                return xpd.get("force_ops")
            tmp = self.getCstr("ops_query", side=side)
            if xpd["other_contains_OR"] and self._pv["single_side_or"]:
                tmp = [o for o in tmp if not o]
            return tmp
    special_cstrs["allw_ops"] = "allw_ops"

    def allw_negs(self, side, type_id, currentRStatus=0):
        if self.isStatusCond(currentRStatus):
            return [False]
        else:
            return self.getCstr("neg_query", side=side, type_id=type_id)
    special_cstrs["allw_negs"] = "allw_negs"

    def neg_query_init(self, side, currentRStatus=0):
        if self.isStatusInitStage(currentRStatus):
            xpd = self.getExpandedStatus(currentRStatus)
            if xpd["other_type_id"] > 0:
                return True in self.getCstr("neg_query", side=(1-side), type_id=xpd["other_type_id"])
        return False
    special_cstrs["neg_query_init"] = "neg_query_init"

    def getActionList(self, k, action_substitutions=[]):
        return self.AR.getActionListRecurse(k, action_substitutions=action_substitutions)

    def getActionsRegistry(self):
        return self.AR

 # FOLDS
    # def filter_folds(self, red):
    #    if self.folds is None:
    #        return False

    #    bcountI = numpy.bincount(self.folds["folds"][list(red.getSuppI())], minlength=self.folds["nb_folds"])
    #    bcountU = numpy.bincount(self.folds["folds"][list(red.getSuppU())], minlength=self.folds["nb_folds"])
    #    bcountU[bcountU == 0] = 1
    #    accs = bcountI/(1.*bcountU)
    #    print "--------------------"
    #    print red.disp()
    #    print accs
    #    if len(numpy.where(accs >= red.getAcc())[0]) > 1:
    #        return False
    #        bb = accs # bcount/self.folds["counts_folds"]
    #        # bpr = bcount/float(numpy.sum(bcount))
    #        # entropS = -numpy.sum(numpy.log(bpr)*bpr)
    #        bpr = bb/numpy.max(bb)
    #        score = numpy.sum(bpr)
    #        print score
    #        # entropM = -numpy.sum(numpy.log(bpr)*bpr)
    #        if score > 1.5:
    #            return False
    #    return True

    # Dependencies between variables (ex, single dataset)

    def setDeps(self, deps=[]):
        self.deps = deps

    def getDeps(self, col=None):
        if col is None:
            return self.deps
        else:
            return self.deps[col]

    def hasDeps(self):
        return len(self.deps) > 0

            
class ActionsRegistry:

    pref_dir = os.path.dirname(os.path.abspath(__file__))
    def_file_basic = pref_dir+"/actions_rdefs_basic.xml"
    default_def_files = [def_file_basic]

    def __init__(self, actions_fns=[], strict=False):
        if not strict:
            actions_fns = self.default_def_files + actions_fns
        self.actions_lists = {}
        self.actions_srcs = {}
        self.parsed_fns = []
        self.setupFDefsFiles(actions_fns)

    def setupFDefsFiles(self, actions_fns):
        for actions_fn in actions_fns:
            default = actions_fn in self.default_def_files
            try:
                fp, fcl = getFp(actions_fn)
                self.readActionsFile(fp, default)
                if fcl:
                    fp.close()
                    if not default:
                        self.parsed_fns.append(actions_fn)
                else:
                    self.parsed_fns.append("package")
            except IOError:
                print("Cannot read actions defs from file %s!" % actions_fn)
                
    def extend(self, actions_fns=[]):
        self.setupFDefsFiles(actions_fns)

    def setActionList(self, fk, flist, fcompact=None, default=False):
        self.actions_lists[fk] = flist
        if not default:
            self.actions_srcs[fk] = fcompact
        else:
            self.actions_srcs[fk] = None

    def delActionList(self, fk):
        if fk in self.actions_lists:
            del self.actions_lists[fk]

    def getActionsKeys(self, public_only=True):
        if public_only:
            return [k for k in self.actions_lists.keys() if not re.match("_", k)]
        return list(self.actions_lists.keys())

    def getActionsKeysSimple(self, public_only=True, patt=None):
        ks = []
        for k in self.getActionsKeys(public_only):
            l = self.getActionListRecurse(k, just_what=True)
            if len(l) == 1 and (patt is None or re.match(patt, l[0])):
                ks.append(k)
        return ks

    def getActionsKeysMulti(self, public_only=True, patt=None):
        ks = []
        for k in self.getActionsKeys(public_only):
            l = self.getActionListRecurse(k, just_what=True)
            if len(l) > 1 and (patt is None or re.match(patt, l[0])):
                ks.append(k)
        return ks

    def clearActions(self):
        self.actions_lists = {}

    def hasActionList(self, fk):
        return fk in self.actions_lists

    def getActionList(self, fk):
        priv = "_"+fk
        if fk not in self.actions_lists and priv in self.actions_lists:
            return self.actions_lists[priv]
        return self.actions_lists.get(fk, [])

    def getActionListRecurse(self, fk, action_substitutions=[], just_what=False):
        actions = []
        for action in self.getActionList(fk):
            if type(action) is dict and "actions" in action:
                actions.extend(self.getActionListRecurse(action["actions"], action_substitutions=action_substitutions, just_what=just_what))
            else:
                a_what = action.getWhat()
                for from_what, to_what in action_substitutions:
                    if action.isWhat(from_what):
                        if to_what is not None:
                            if just_what:
                                a_what = self.getSubstitutionWhat(to_what)
                            else:
                                action = action.doSubstitution(to_what)
                        else:
                            a_what = None
                            action = None
                if just_what:
                    if a_what is not None:
                        actions.append(a_what)
                else:
                    if action is not None:
                        actions.append(action)
        return actions


    def readActionsFile(self, actions_fp, default=False):
        doc = toolXML.parseXML(actions_fp)
        for n in toolXML.subTags(doc, "actionlist"):
            ak, alist = self.parseActionListNode(n)
            if ak is not None:
                self.setActionList(ak, alist, n.toxml(), default)
        
    def parseActionListNode(self, n):
        name = toolXML.getTagData(n, "name")
        actions = []
        for child in toolXML.children(n):
            if toolXML.isElementNode(child):
                if toolXML.tagName(child) == "action":
                    a = ActionFactory.initFromXml(child)
                    if a is not None:
                        actions.append(a)

                elif toolXML.tagName(child) == "list":
                    lname = toolXML.getTagData(child, "name")
                    if self.hasActionList(lname):
                        actions.append({"actions": lname})
        return name, actions
    
    def actionsToStr(self):
        srcs = [src for k, src in self.actions_srcs.items() if src is not None]
        xps = "\n".join(srcs)
        if len(xps) > 0:
            head = "<!-- Action definitions read from %s-->" % ";".join(self.parsed_fns)
            xps = "<?xml version=\"1.0\" encoding=\"utf-8\"?>\n<root>\n"+head+"\n"+xps+"\n</root>"
        return xps
