import argparse
import glob
import re
import tempfile
import os.path
import sys
try:
    import toolXML
except ModuleNotFoundError:
    from . import toolXML
import pdb

USAGE_DEF_HEADER = "Clired redescription mining"
USAGE_DEF_URL_HELP = "http://cs.uef.fi/siren/help/"

EXT_SIREN = ".siren"
# PATT_QALT = "suff(ix)?(?P<suff>.+)"
PATT_QSUFF = "\+(?P<suff>.+)"

conf_names = {"miner": 0, "inout": 0, "dataext": 1, "rnd": 1, "folds": 1}
conf_names_data = ["inout", "dataext"]


class TmpDirManager:
    TOKEN_TMP_DIR = "__TMP_DIR__"

    @classmethod
    def isTmpDirToken(tcl, v):
        return v == tcl.TOKEN_TMP_DIR

    @classmethod
    def getTmpDir(tcl, prefix):
        return tempfile.mkdtemp(prefix=prefix)

    @classmethod
    def inTmpDir(tcl, path):
        tmp_dir = tempfile.gettempdir()
        try:
            return os.path.commonpath([path, tmp_dir]) == tmp_dir
        except ValueError:
            return False


def fillNamesForConfDefs(conf_defs=None, wdata=True):
    if conf_defs is None:
        conf_defs = [k for (k, v) in conf_names.items() if v == 0]
    elif any([(type(c) is int) for c in conf_defs]):
        cd = set()
        for cc in conf_defs:
            if type(cc) is int:
                cd.update([k for (k, v) in conf_names.items() if v == cc])
            else:
                cd.add(cc)
        conf_defs = list(cd)
    if not wdata:
        conf_defs = [c for c in conf_defs if c not in conf_names_data]
    return conf_defs


class CHelpFormatter(argparse.HelpFormatter):

    def _get_help_string(self, action):
        return action.help

    def _format_action_invocation(self, action):
        if not action.option_strings:
            default = self._get_default_metavar_for_positional(action)
            metavar, = self._metavar_formatter(action, default)(1)
            return metavar

        else:
            parts = []

            # if the Optional doesn't take a value, format is:
            #    -s, --long
            if action.nargs == 0:
                if len(action.option_strings) == 1 and action.option_strings[0] == ("--no-%s" % action.dest):
                    parts.append("--[no-]%s" % action.dest)
                else:
                    parts.extend(action.option_strings)

            # if the Optional takes a value, format is:
            #    -s ARGS, --long ARGS
            else:
                default = self._get_default_metavar_for_optional(action)
                args_string = self._format_args(action, default)
                for option_string in action.option_strings:
                    parts.append('%s %s' % (option_string, args_string))

            return ', '.join(parts)

    def _format_args(self, action, default_metavar):
        get_metavar = self._metavar_formatter(action, default_metavar)
        if action.nargs is None:
            result = '%s' % get_metavar(1)
        elif action.nargs == argparse.OPTIONAL:
            result = '[%s]' % get_metavar(1)
        elif action.nargs == argparse.ZERO_OR_MORE:
            result = '[%s ...]' % get_metavar(1)
        elif action.nargs == argparse.ONE_OR_MORE:
            result = '%s [...]' % get_metavar(1)
        elif action.nargs == argparse.REMAINDER:
            result = '...'
        elif action.nargs == argparse.PARSER:
            result = '%s ...' % get_metavar(1)
        elif action.nargs == argparse.SUPPRESS:
            result = ''
        else:
            try:
                formats = ['%s' for _ in range(action.nargs)]
            except TypeError:
                raise ValueError("invalid nargs value") from None
            result = ' '.join(formats) % get_metavar(action.nargs)
        return result


class CAction(argparse.Action):

    suppress_sections = None
    only_sections = None
    only_core = False

    def prepareHelp(self, hlp):
        if self._section is not None:
            if CAction.only_core and self._core is None:
                return argparse.SUPPRESS
            if CAction.only_sections is not None and not any([re.match(s, self._section) is not None for s in CAction.only_sections]):
                return argparse.SUPPRESS
            if CAction.suppress_sections is not None and any([re.match(s, self._section) is not None for s in CAction.suppress_sections]):
                return argparse.SUPPRESS
        return self._help

    @property
    def help(self):
        return self.prepareHelp(self._help)

    @help.setter
    def help(self, hlp):
        if hlp == argparse.SUPPRESS:
            self._help = hlp
            self._section = None
            self._core = None
        else:
            tmp = re.match("(\[(?P<section>.*)\])?\s*(?P<core>CORE)?\s*(?P<hlp>.*)$", hlp)
            self._help = tmp.group("hlp")
            self._section = tmp.group("section")
            self._core = tmp.group("core")


class CStoreAction(argparse._StoreAction, CAction):
    def __init__(self,
                 option_strings,
                 dest,
                 nargs=None,
                 const=None,
                 default=None,
                 type=None,
                 choices=None,
                 required=False,
                 help=None,
                 metavar=None):
        argparse._StoreAction.__init__(self,
                                       option_strings=option_strings,
                                       dest=dest,
                                       nargs=nargs,
                                       const=const,
                                       default=default,
                                       type=type,
                                       choices=choices,
                                       required=required,
                                       help=help,
                                       metavar=metavar)


class CStoreConstAction(argparse._StoreConstAction, CAction):

    def __init__(self,
                 option_strings,
                 dest,
                 const,
                 default=None,
                 required=False,
                 help=None,
                 metavar=None):
        argparse._StoreConstAction.__init__(self,
                                            option_strings=option_strings,
                                            dest=dest,
                                            const=const,
                                            default=default,
                                            required=required,
                                            help=help,
                                            metavar=metavar)


class CStoreBoolAction(CStoreConstAction):
    

    def __call__(self, parser, namespace, values, option_string=None):
        if len(values) == 1 and BooleanCParameter.map_str.get(values[0].lower()) == False:
            setattr(namespace, self.dest, not self.const)
        else:
            setattr(namespace, self.dest, self.const)


class CStoreBoolTAction(CStoreBoolAction):
    
    def __init__(self,
                 option_strings,
                 dest,
                 default=True,
                 required=False,
                 help=None,
                 metavar=None):
        CStoreConstAction.__init__(self,
                              option_strings=option_strings,
                              dest=dest,
                              const=True,
                              default=default,
                              required=required,
                              help=help,
                              metavar=metavar)
        
class CStoreBoolFAction(CStoreBoolAction):

    def __init__(self,
                 option_strings,
                 dest,
                 default=False,
                 required=False,
                 help=None,
                 metavar=None):
        CStoreConstAction.__init__(self,
                              option_strings=option_strings,
                              dest=dest,
                              const=False,
                              default=default,
                              required=required,
                              help=help,
                              metavar=metavar)

class CFilenameAction(CAction):
    def __call__(self, parser, namespace, values, option_string=None):
        items = getattr(namespace, self.dest, None)
        if items is None:
            items = {}
        for filename in values:
            if "src_folder" not in items:
                items["src_folder"] = os.path.dirname(os.path.abspath(filename))
            p, ext = os.path.splitext(filename)
            if ext == EXT_SIREN:
                k = "pack_file"
            elif ext == ".conf" or ext == ".xml":
                k = "conf_file"
            elif ext == ".queries":
                k = "queries_file"
            elif ext == ".csv":
                k = "data_file"
                if re.search("queries", p) or len(items.get(k, [])) > 1:
                    k = "queries_file"
            elif ext == ".txt":
                k = "conf_file"
                if re.search("queries", p) or len(items.get(k, [])) > 0:
                    k = "queries_file"
            elif re.match(PATT_QSUFF, p):
                k = "queries_file"
            else:
                k = "unknown"
            if k not in items:
                items[k] = []
            items[k].append(filename)
        setattr(namespace, self.dest, items)


class CRangeAction(CAction):
    def __call__(self, parser, namespace, values, option_string=None):
        # if self.const[0] <= values <= self.const[1]:
        setattr(namespace, self.dest, values)


class COptionAction(CAction):
    def __call__(self, parser, namespace, values, option_string=None):
        v = values[0]  # nargs -> always exactly one value
        if v not in self.choices:
            v = self.choices[int(v)]
        setattr(namespace, self.dest, v)


class COptionsAction(CAction):
    def __call__(self, parser, namespace, values, option_string=None):
        items = getattr(namespace, self.dest, None)
        if items is None:
            items = set()
        for v in values:
            if v not in self.choices:
                v = self.choices[int(v)]
            items.add(v)
        setattr(namespace, self.dest, items)


class CParameter(object):
    type_id = "X"
    type_str = "X"
    value_types = {"text": str, "boolean": bool, "integer": int, "float": float, "color": str}

    def __init__(self,  name=None, label=None, default=None, value_type=None, legend=None):
        self._name = name
        self._label = label
        self._default = default
        self._legend = legend
        self._value_type = value_type
        self._pfr = None

    def hasOptions(self):
        return False

    def getCardinality(self):
        return "unique"

    def isCore(self):
        return re.match("CORE", self._legend) is not None

    def setPathFromRoot(self, pfr):
        self._prf = tuple(pfr)

    def getPathFromRoot(self, sep=":"):
        return self.getPath(sep=sep, from_pos=0)

    def getPath(self, sep=":", from_pos=1):
        if sep is None:
            return self._prf[from_pos:]
        return sep.join(self._prf[from_pos:])

    def getConfDef(self):
        return self._prf[0]

    def suppress(self, only_core=False, only_sections=None, suppress_sections=None, is_def=None, defaults=None):
        if only_core and not self.isCore():
            return True
        path = self.getPathFromRoot()
        if only_sections is not None and not any([re.match(s, path) is not None for s in only_sections]):
            return True
        if suppress_sections is not None and any([re.match(s, path) is not None for s in suppress_sections]):
            return True
        if is_def is None:
            return False

        if only_core:
            return False  # must be core to get here, keep regardless
        return is_def and not defaults

    def getId(self):
        return self._name

    def getValueType(self):
        return self._value_type

    def valToText(tcl, value):
        return str(value)

    def textToVal(self, txt):
        return toolXML.parseToType(txt, self.getValueType())

    def textToIndex(self, txt):
        return None

    def getTypeStr(self):
        return self.type_str

    def getTypeDets(self):
        return self.type_str

    def getName(self):
        return self._name

    def getLegend(self):
        return re.sub("CORE *", "", self._legend)

    def getArgLegend(self, sep=":"):
        if sep is not None:
            return "[%s] %s" % (self.getPathFromRoot(sep), self._legend)
        return self._legend

    def getInfo(self):
        return "%s (%s)" % (self.getLegend(), self.getTypeDets())

    def getLabel(self):
        return self._label

    def parseNode(self, node):
        self._name = toolXML.getTagData(node, "name")
        self._legend = toolXML.getTagData(node, "legend")
        if self._name is None:
            raise Exception("Name for param undefined!")
        self._label = toolXML.getTagData(node, "label")
        if self._label is None:
            raise Exception("Label for param %s undefined!" % self._name)
        tmp_vt = toolXML.getTagData(node, "value_type")
        if tmp_vt in self.value_types:
            self._value_type = self.value_types[tmp_vt]
        if self._value_type is None:
            raise Exception("Value type for param %s undefined!" % self._name)

    def __str__(self):
        return "Parameter (%s): %s" % (self.type_id, self._name)

    def getArgparseBase(self):
        return {"help": self.getArgLegend(), "default": argparse.SUPPRESS, "type": self.getValueType(), "action": CStoreAction}

    def getArgparseKargs(self):
        return [(["--"+self.getName()], self.getArgparseBase())]

    # DEFAULT
    def getDefaultValue(self):
        return self._default

    def isDefaultValue(self, value):
        if TmpDirManager.isTmpDirToken(self.getDefaultValue()):
            return TmpDirManager.inTmpDir(value) or value == self.getDefaultValue()
        return value == self.getDefaultValue()

    def getDefaultIndex(self):
        return None

    def isDefaultIndex(self, value):
        return value == self.getDefaultIndex()

    def getDefaultText(self):
        return self.valToText(self._default)

    def getDefaultStr(self):
        return self.getDefaultText()

    # PARAM
    def getParamValue(self, raw_value):
        return self.textToVal(raw_value)

    def getParamValueOrDefault(self, raw_value):
        v = self.getParamValue(raw_value)
        if v is None:
            return self.getDefaultValue()
        return v

    def getParamIndex(self, raw_value):
        return None

    def getParamText(self, value):
        return self.valToText(value)

    def getParamTextRaw(self, raw_value):
        v = self.getParamValue(raw_value)
        if v is not None:
            return self.valToText(v)


class OpenCParameter(CParameter):
    type_id = "open"
    type_str = "open text"

    def __init__(self, name=None, label=None, default=None, value_type=None, length=None, legend=None):
        CParameter.__init__(self, name, label, default, value_type, legend)
        self._length = length

    def parseNode(self, node):
        CParameter.parseNode(self, node)
        self._length = toolXML.getTagData(node, "length", int)
        et = node.getElementsByTagName("default")
        if len(et) > 0:
            self._default = toolXML.getValue(et[0], self._value_type)
        if self._default is None:
            raise Exception("Default value for param %s undefined!" % self._name)


class RangeCParameter(CParameter):
    type_id = "range"
    type_str = "range"

    def __init__(self, name=None, label=None, default=None, value_type=None, range_min=None, range_max=None, legend=None):
        CParameter.__init__(self, name, label, default, value_type, legend)
        self._range_min = range_min
        self._range_max = range_max

    def parseNode(self, node):
        CParameter.parseNode(self, node)
        self._range_min = toolXML.getTagData(node, "range_min", self._value_type)
        self._range_max = toolXML.getTagData(node, "range_max", self._value_type)
        et = node.getElementsByTagName("default")
        if len(et) > 0:
            self._default = toolXML.getValue(et[0], self._value_type)
        if self._default is None or self._range_min is None or self._range_max is None \
                or self._default < self._range_min or self._default > self._range_max:
            raise Exception("Default value for param %s not in range!" % self._name)

    def getTypeDets(self):
        return "%s in [%s, %s]" % (self.type_str, self._range_min, self._range_max)

    def getArgparseBase(self):
        kargs = CParameter.getArgparseBase(self)
        kargs["action"] = CRangeAction
        kargs["const"] = (self._range_min, self._range_max)
        return kargs

    def getParamValue(self, raw_value):
        tmp = self.textToVal(raw_value)
        if tmp is not None and tmp >= self._range_min and tmp <= self._range_max:
            return tmp
        return None


class SingleOptionsCParameter(CParameter):
    type_id = "single_options"
    type_str = "single option"

    def __init__(self, name=None, label=None, default=None, value_type=None, options=None, legend=None):
        CParameter.__init__(self, name, label, default, value_type, legend)
        self._options = options

    def hasOptions(self):
        return True

    def parseNode(self, node):
        CParameter.parseNode(self, node)
        et = node.getElementsByTagName("options")
        if len(et) > 0:
            self._options = toolXML.getValues(et[0], self._value_type)
        et = node.getElementsByTagName("default")
        if len(et) > 0:
            self._default = toolXML.getValue(et[0], int)
        if self._default is None or self._default < 0 or self._default >= len(self._options):
            raise Exception("Default value for param %s not among options!" % self._name)

    def getTypeDets(self):
        return "%s in {%s}" % (self.type_str, ", ".join(self.getOptionsText()))

    def getArgparseBase(self):
        kargs = CParameter.getArgparseBase(self)
        kargs["action"] = COptionAction
        kargs["nargs"] = 1
        kargs["choices"] = self.getOptions()
        return kargs

    def textToIndex(self, txt):
        try:
            tmp = toolXML.parseToType(txt, int)
            if tmp is not None and tmp >= 0 and tmp < len(self._options):
                return tmp
        except ValueError:
            tmp = None

    def getOptions(self):
        return self._options

    def getOptionsText(self):
        return [self.valToText(v) for v in self._options]

    # DEFAULT
    def getDefaultValue(self):
        return self._options[self._default]

    def getDefaultIndex(self):
        return self._default

    def getDefaultText(self):
        return self.valToText(self._options[self._default])

    # PARAM
    def getParamValue(self, raw_value):
        if type(raw_value) is int:
            if raw_value >= -1 and raw_value < len(self._options):
                return self._options[raw_value]
            return None
        try:
            tmp = self.textToVal(raw_value)
            if tmp is not None and tmp in self._options:
                return tmp
        except ValueError:
            idx = self.textToIndex(raw_value)
            if idx is not None:
                return self._options[idx]

    def getParamIndex(self, raw_value):
        if type(raw_value) is int:
            if raw_value >= -1 and raw_value < len(self._options):
                return raw_value
            return None
        try:
            tmp = self.textToVal(raw_value)
            if tmp is not None and tmp in self._options:
                return self._options.index(tmp)
        except ValueError:
            return self.textToIndex(raw_value)


class BooleanCParameter(SingleOptionsCParameter):
    type_id = "boolean"
    type_str = "yes/no"
    opts_data = [False, True]
    map_str = {"yes": True, "no": False,
               "true": True, "false": False,
               "t": True, "f": False,
               "1": True, "0": False}
    inv_str = {True: "yes", False: "no"}

    def __init__(self, name=None, label=None, default=None, value_type=bool, options=None, legend=None):
        CParameter.__init__(self, name, label, default, value_type, legend)
        self._options = self.opts_data

    def hasOptions(self):
        return False

    def parseNode(self, node):
        CParameter.parseNode(self, node)
        et = node.getElementsByTagName("default")
        if len(et) > 0:
            self._default = toolXML.getValue(et[0], bool)
        if self._default is None:
            raise Exception("Default value for param %s not among options!" % self._name)

    def getTypeDets(self):
        return self.type_str

    def getArgparseBase(self, positive=True):
        if positive:
            return {"help": argparse.SUPPRESS, "default": argparse.SUPPRESS, "action": CStoreBoolTAction, "dest": self.getName()}
        else:
            return {"help": self.getArgLegend(), "default": argparse.SUPPRESS, "action": CStoreBoolFAction, "dest": self.getName()}

    def getArgparseKargs(self):
        return [(["--"+self.getName()], self.getArgparseBase(True)), (["--no-"+self.getName()], self.getArgparseBase(False))]

    def valToText(self, value):
        return self.inv_str.get(value)


class MultipleOptionsCParameter(SingleOptionsCParameter):
    type_id = "multiple_options"
    type_str = "multiple options"

    def getCardinality(self):
        return "multiple"

    def hasOptions(self):
        return True

    def parseNode(self, node):
        CParameter.parseNode(self, node)
        et = node.getElementsByTagName("options")
        if len(et) > 0:
            self._options = toolXML.getValues(et[0], self._value_type)
        et = node.getElementsByTagName("default")
        if len(et) > 0:
            self._default = toolXML.getValues(et[0], int)
        if self._default is None or (len(self._default) > 0 and (min(self._default) < 0 or max(self._default) >= len(self._options))):
            raise Exception("Some default value for param %s not among options!" % self._name)

    def getArgparseBase(self):
        kargs = CParameter.getArgparseBase(self)
        kargs["action"] = COptionsAction
        kargs["nargs"] = "*"
        kargs["choices"] = self.getOptions()
        return kargs

    # DEFAULT
    def getDefaultValue(self):
        return [self._options[d] for d in self._default]

    def getDefaultData(self):
        return [self._options[i] for i in self._default]

    def getDefaultText(self):
        return [self.valToText(self._options[i]) for i in self._default]

    def getDefaultStr(self):
        return "{" + ", ".join(self.getDefaultText()) + "}"

    # PARAM

    def getParamValue(self, raw_value):
        vs = set()
        for v in raw_value:
            if type(v) is int:
                if v >= -1 and v < len(self._options):
                    vs.add(self._options[v])
            else:
                try:
                    tmp = self.textToVal(v)
                    if tmp is not None and tmp in self._options:
                        vs.add(tmp)
                except ValueError:
                    idx = self.textToIndex(v)
                    if idx is not None:
                        vs.add(self._options[idx])
        return sorted(vs)

    def getParamIndex(self, raw_value):
        vs = set()
        for v in raw_value:
            if type(v) is int:
                if v >= -1 and v < len(self._options):
                    vs.add(v)
            else:
                try:
                    tmp = self.textToVal(v)
                    if tmp is not None and tmp in self._options:
                        vs.add(self._options.index(tmp))
                except ValueError:
                    idx = self.textToIndex(v)
                    if idx is not None:
                        vs.add(idx)
        return sorted(vs)

    def getParamText(self, value):
        return set([self.valToText(v) for v in value])

    def getParamTextRaw(self, raw_value):
        vs = self.getParamValue(raw_value)
        return set([self.valToText(v) for v in vs])


class ColorCParameter(CParameter):
    type_id = "color_pick"
    type_str = "color #RRGGBB"
    match_p = "^#(?P<rr>[0-9A-Fa-f][0-9A-Fa-f])(?P<gg>[0-9A-Fa-f][0-9A-Fa-f])(?P<bb>[0-9A-Fa-f][0-9A-Fa-f])$"

    def valToText(self, tuple_value):
        return "#"+"".join([v.replace("x", "")[-2:] for v in map(hex, tuple_value)])

    def textToVal(self, txt_value):
        if txt_value is not None:
            g = re.match(self.match_p, txt_value)
            if g is not None:
                try:
                    return (int(g.group("rr"), 16), int(g.group("gg"), 16), int(g.group("bb"), 16))
                except:
                    raise Warning("Could not parse color %s!" % txt_value)
        return None

    def parseNode(self, node):
        CParameter.parseNode(self, node)
        et = node.getElementsByTagName("default")
        if len(et) > 0:
            tmp = toolXML.getValue(et[0], self._value_type)
            self._default = self.textToVal(tmp)

        if self._default is None:
            raise Exception("Default value for param %s not correct!" % self._name)

    def getArgparseBase(self):
        kargs = CParameter.getArgparseBase(self)
        kargs["type"] = self.textToVal
        return kargs


class PreferencesManager(object):
    parameter_types = {"open": OpenCParameter,
                       "range": RangeCParameter,
                       "single_options": SingleOptionsCParameter,
                       "boolean": BooleanCParameter,
                       "multiple_options": MultipleOptionsCParameter,
                       "color_pick": ColorCParameter}
    MTCH_ST = "^(?P<basis>[^0-9]*)((_s(?P<side>[01])_(?P<typ>[0-9]))|(_s(?P<oside>[01]))|(_(?P<otyp>[0-9])))$"

    @classmethod
    def getConfDefSections(tcl, conf_defs=None, wdata=True):
        if conf_defs is None:
            return None
        return ["%s_confdef.xml" % c if c in conf_names else os.path.basename(c) for c in fillNamesForConfDefs(conf_defs, wdata)]

    def __init__(self, filenames):
        self.subsections = []
        self.pdict = {}

        if type(filenames) == str:
            filenames = [filenames]
        for filename in filenames:
            if filename is not None:
                doc = toolXML.parseXML(filename)
                if doc is not None:
                    src = os.path.basename(filename)
                    params = self.processDom(doc.documentElement, secs=[src])
                    if type(params) == dict and len(params) == 1 and "subsections" in params:
                        for si in range(len(params["subsections"])):
                            params["subsections"][si]["filename"] = filename
                        self.subsections.extend(params["subsections"])

    def __str__(self):
        strd = "Preferences manager:\n"
        for sec in self.getTopSections():
            strd += self.dispSection(sec)
        return strd

    def getTopSections(self, conf_filter=None, wdata=True):
        only_sections = self.getConfDefSections(conf_filter, wdata)
        sections = []
        for t in self.subsections:
            path = ":".join(t.get("pfr"))
            if only_sections is None or any([re.match(s, path) is not None for s in only_sections]):
                sections.append(t)
        return sections

    def dispSection(self, parameters, level=0):
        strs = ("\t"*level)+("[%s]\n" % parameters.get("name", ""))
        for k in self.parameter_types.keys():
            if len(parameters[k]) > 0:
                strs += ("\t"*level)+("   * %s:\n" % k)
            for item_id in parameters[k]:
                item = self.pdict[item_id]
                tmp_str = str(item)
                tmp_str.replace("\n", "\n"+("\t"*(level+1)))
                strs += ("\t"*(level+1))+tmp_str+"\n"
        if len(parameters["subsections"]) > 0:
            strs += ("\t"*level)+("   * subsections:\n")
        for k in parameters["subsections"]:
            strs += self.dispSection(k, level+1)
        return strs

    def getItem(self, item_id):
        return self.pdict.get(item_id, None)

    def getItems(self):
        return self.pdict.items()

    def getNameSidesTypes(self, name):
        tmp = re.match(self.MTCH_ST, name)
        if tmp is not None:
            return (tmp.group("basis"),
                    int(tmp.group("side") or tmp.group("oside") or -1),
                    int(tmp.group("typ") or tmp.group("otyp") or -1))
        return (name, -1, -1)

    def getItemsSidesTypes(self):
        dd = {}
        for k in self.pdict.keys():
            tmp = self.getNameSidesTypes(k)
            if tmp[1] != -1 or tmp[2] != -1:
                dd[k] = tmp
        return dd

    def getDefaultValues(self, only_core=False, only_sections=None, suppress_sections=None):
        return dict([(item_id, item.getDefaultValue())
                     for (item_id, item) in self.pdict.items() if not item.suppress(only_core=only_core, only_sections=only_sections, suppress_sections=suppress_sections)])

    def getListOptions(self):
        return [item.getName()+"=" for (item_id, item) in self.pdict.items()]

    def processDom(self, current, secs=[]):
        parameters = None
        name = None
        if toolXML.isElementNode(current):
            if toolXML.tagName(current) in ["root", "section"]:
                parameters = {"subsections": []}
                if toolXML.tagName(current) == "section":
                    parameters["name"] = toolXML.getTagData(current, "name")
                    parameters["pfr"] = tuple(secs)
                    for k in self.parameter_types.keys():
                        parameters[k] = []
                    secs = list(secs + [parameters["name"]])
                for child in toolXML.children(current):
                    tmp = self.processDom(child, secs)
                    if tmp is not None:
                        if type(tmp) == dict:
                            parameters["subsections"].append(tmp)
                        elif tmp.type_id in parameters.keys():
                            tmp_id = tmp.getId()
                            if tmp_id in self.pdict:
                                raise Exception("Encountered two parameters with same id %s!" % tmp_id)
                            else:
                                self.pdict[tmp_id] = tmp
                                parameters[tmp.type_id].append(tmp_id)
            if toolXML.tagName(current) == "parameter":
                name = toolXML.getTagData(current, "name")
                parameter_type = toolXML.getTagData(current, "parameter_type")
                if parameter_type in self.parameter_types.keys():
                    parameters = self.parameter_types[parameter_type]()
                    parameters.parseNode(current)
                    parameters.setPathFromRoot(secs)
        if parameters is not None:
            return parameters

    def dispParametersRec(self, parameters, pv, level=0, sections=True, helps=False, defaults=False, only_core=False, only_sections=None, suppress_sections=None, xml=True):
        indents = ""
        strd, header, footer = ("", "", "")
        if xml:
            if sections:
                indents = "\t"*(level+1)
                header = ("\t"*level)+"<section>\n"+indents+("<name>%s</name>\n" % parameters.get("name", ""))
                footer = ("\t"*level)+"</section>\n"
        else:
            if sections:
                footer = ("#"*level) + ("# %s" % parameters.get("name", "")) + "\n"
            if not sections and not helps:
                indents = " --"

        for k in self.parameter_types.keys():
            for item_id in parameters[k]:
                item = self.getItem(item_id)
                if item is not None:
                    if pv is None or item_id not in pv:
                        vs = item.getDefaultText()
                        is_def = True
                    else:
                        vs = item.getParamText(pv[item_id])
                        is_def = item.isDefaultValue(pv[item_id])

                    if item.suppress(only_core=only_core, only_sections=only_sections, suppress_sections=suppress_sections, is_def=is_def, defaults=defaults):
                        item = None

                if item is not None:  # and ((core and item.isCore()) or (defaults or not is_def))
                    if xml:
                        strd += indents+"<parameter>\n"
                        strd += indents+"\t<name>" + item.getName() + "</name>\n"
                        if sections:
                            strd += indents+"\t<label>" + item.getLabel() + "</label>\n"
                        if helps:
                            strd += indents+"\t<info>" + item.getInfo() + "</info>\n"
                        if type(vs) is set or type(vs) is list:
                            for v in vs:
                                strd += indents+"\t<value>" + v + "</value>\n"
                        else:
                            strd += indents+"\t<value>" + vs + "</value>\n"
                        strd += indents+"</parameter>\n"
                    else:
                        if type(vs) is set or type(vs) is list:
                            strd += indents + item.getName() + " " + " ".join(vs)
                        elif isinstance(item, BooleanCParameter):
                            strd += indents + "no-"*(vs == "no") + item.getName()
                        else:
                            strd += indents + item.getName() + " " + vs

                        if sections or helps:
                            strd += " #"
                            if sections:
                                strd += " [%s]" % item.getLabel()
                            if helps:
                                strd += " " + item.getInfo()
                            strd += "\n"

        for k in parameters["subsections"]:
            strd += self.dispParametersRec(k, pv, level+1, sections, helps, defaults,
                                           only_core=only_core, only_sections=only_sections, suppress_sections=suppress_sections, xml=xml)
        if len(strd) > 0:
            return header + strd + footer
        else:
            return ""

    def dispParameters(self, pv=None, conf_filter=None, sections=True, helps=False, defaults=False, only_core=False, only_sections=None, suppress_sections=None, xml=True, wdata=True):
        if only_sections is None and suppress_sections is None:
            only_sections = self.getConfDefSections(conf_filter, wdata)
        if (pv is None or conf_filter is not None) and defaults:
            helps = True
            if not only_core:
                sections = True

        strd = ""
        for subsection in self.subsections:
            # print("---SUBSECTION:", subsection.get("name"))
            strd += self.dispParametersRec(subsection, pv, 0, sections, helps, defaults,
                                           only_core=only_core, only_sections=only_sections, suppress_sections=suppress_sections, xml=xml)
        if xml:
            if len(strd) == 0:
                strd = "\t<!-- Using only default parameters -->\n"
            return "<root>\n"+strd+"</root>"
        else:
            return strd

    def collectParametersRec(self, parameters, sections=True, only_core=False, only_sections=None, suppress_sections=None):
        items = []
        for k in self.parameter_types.keys():
            for item_id in parameters[k]:
                item = self.getItem(item_id)
                if not item.suppress(only_core=only_core, only_sections=only_sections, suppress_sections=suppress_sections):
                    items.append((item_id, item))
        if len(items) > 0 and sections:
            items.insert(0, ("[section]", tuple(parameters["pfr"] + (parameters["name"],))))

        for k in parameters["subsections"]:
            items.extend(self.collectParametersRec(k, sections,
                                                   only_core=only_core, only_sections=only_sections, suppress_sections=suppress_sections))
        return items

    def collectParameters(self, conf_filter=None, sections=True, only_core=False, only_sections=None, suppress_sections=None, wdata=True):
        if only_sections is None and suppress_sections is None:
            only_sections = self.getConfDefSections(conf_filter, wdata)
        items = []
        for subsection in self.subsections:
            items.extend(self.collectParametersRec(subsection, sections,
                                                   only_core=only_core, only_sections=only_sections, suppress_sections=suppress_sections))
        return items


class PreferencesReader(argparse.ArgumentParser):

    outer_arguments = ["template", "config", "pack_file", "conf_file"]
    @classmethod
    def copyParams(tcl, params):
        return dict(params.items())

    def __init__(self, pm, tasks=[], tasks_default=None, description=None, epilog=None):
        if description is None:
            description = USAGE_DEF_HEADER+"."
        if epilog is None:
            epilog = "\nFor more details see "+USAGE_DEF_URL_HELP+"formats.html and "+USAGE_DEF_URL_HELP+"preferences.html#mining-parameters"
        argparse.ArgumentParser.__init__(self, description=description, epilog=epilog, formatter_class=CHelpFormatter)
        self.pm = pm

        self.add_argument("filename", nargs="*", type=str, action=CFilenameAction, default=argparse.SUPPRESS, help="Preferences file (.conf/.xml), data file (.csv), queries file (.queries) or Siren package (.siren) containing the data, preferences, etc.")

        self.add_argument("--config", action="store_true", help="generate a default preferences file")
        self.add_argument("--template", action="store_true", help="generate a basic preferences template")
        self.add_argument("--debug", action="store_true", help="turn on debugging mode")

        choices_task = set(tasks)
        if tasks_default is not None:
            choices_task.add(tasks_default)
            default_task = tasks_default
        else:
            default_task = argparse.SUPPRESS
        if len(choices_task) > 0:
            self.add_argument("--task", type=str, choices=choices_task, default=default_task, help="task to perform")

        g = self
        for (item_id, item) in self.collectParameters():
            if item_id == "[section]":
                g = self.add_argument_group(":".join(item[1:]))  # drop file src
            else:
                for (args, kargs) in item.getArgparseKargs():
                    g.add_argument(*args, **kargs)

    def getManager(self):
        return self.pm

    def getTopSections(self, conf_filter=None):
        return self.getManager().getTopSections(conf_filter=conf_filter)

    def getDefaultValues(self, only_core=False, only_sections=None, suppress_sections=None):
        return self.getManager().getDefaultValues(only_core=only_core, only_sections=only_sections, suppress_sections=suppress_sections)

    def dispParameters(self, pv=None, conf_filter=None, sections=True, helps=False, defaults=False, only_core=False, only_sections=None, suppress_sections=None):
        return self.getManager().dispParameters(pv=pv, conf_filter=conf_filter, sections=sections, helps=helps, defaults=defaults,
                                                only_core=only_core, only_sections=only_sections, suppress_sections=suppress_sections)

    def collectParameters(self, conf_filter=None, sections=True, only_core=False, only_sections=None, suppress_sections=None):
        return self.getManager().collectParameters(conf_filter=conf_filter, sections=sections, only_core=only_core, only_sections=only_sections, suppress_sections=suppress_sections)

    def dispRST(self, conf_filter=None, sections=True, helps=True, defaults=True, only_core=False, only_sections=None, suppress_sections=None):
        rst = ""
        for (item_id, item) in self.collectParameters(conf_filter=conf_filter, sections=sections, only_core=only_core, only_sections=only_sections, suppress_sections=suppress_sections):
            if item_id == "[section]":
                rst += "\n.. rubric:: %s\n" % ": ".join(item[1:])
                # rst += "\n%s\n%s\n\n%s (%s)\n" % (item[-1], "-"*len(item[-1]), " ".join(item[1:-1]), item[0])
            else:
                rst += "\n.. cparam:: %s" % item.getName()
                if item.isCore():
                    rst += "\n    :core:"
                rst += "\n    :label: %s" % item.getLabel()
                rst += "\n    :conf_def: %s" % item.getConfDef()
                rst += "\n    :path: %s" % item.getPath(sep=", ")
                rst += "\n    :type: %s" % item.getTypeStr()
                if item.hasOptions():
                    rst += "\n    :choices: " + ", ".join(item.getOptionsText())
                if helps:
                    rst += "\n\n    * " + item.getLegend()
                    rst += "\n    * " + item.getTypeDets()
                if defaults:
                    rst += "\n    * default: " + item.getDefaultStr()
                if helps:
                    rst += "\n    * %s" % item.getConfDef()
                rst += "\n"
        return rst

    def format_help(self, only_core=True, only_sections=None, suppress_sections=None):
        CAction.only_core = only_core
        CAction.only_sections = only_sections
        CAction.suppress_sections = suppress_sections
        return argparse.ArgumentParser.format_help(self)

    def processHelp(self, params, conf_filter=None):
        if params.get("usage", False):
            return self.format_usage()
        elif params.get("template"):
            return self.getManager().dispParameters(pv=params, conf_filter=conf_filter, sections=False, helps=True, defaults=False, only_core=True)
        elif params.get("config"):
            return self.getManager().dispParameters(pv=None, conf_filter=conf_filter, sections=True, helps=True, defaults=True, only_core=False)

    def parse_arguments(self, arguments=[]):
        parsed, leftover = self.parse_known_args(arguments)
        return vars(parsed), leftover

    def apply_actions(self, action_tuples, namespace=None):
        if namespace is None:
            namespace = argparse.Namespace()

        for action, argument_strings, option_string in action_tuples:
            argument_values = self._get_values(action, argument_strings)
            if argument_values is not argparse.SUPPRESS:
                action(self, namespace, argument_values, option_string)

        return vars(namespace)

    def _check_value(self, action, value):
        # converted value must be one of the choices (if specified)
        if action.choices is not None and value not in action.choices:
            if value == "-":
                value = "-1"
            try:
                v = int(value)
            except ValueError:
                v = None
            if v is None or v < -1 or v >= len(action.choices):
                raise argparse.ArgumentError(action, "invalid choice: %s (choose from %s)" % (value, ", ".join(map(repr, action.choices))))
        elif action.const is not None and type(action.const) is tuple and not (action.const[0] <= value <= action.const[1]):
            raise argparse.ArgumentError(action, "invalid value: %s (choose from [%s, %s])" % (value, action.const[0], action.const[1]))

    def getPreferences(self, arguments=None, conf_file=None, pv=None):
        preferences_mod = set()
        if pv is None:
            pv = self.getManager().getDefaultValues()

        if arguments is not None and len(arguments) > 1:
            params_cmdline, leftover_args = self.parse_arguments(arguments[1:])
        else:
            params_cmdline, leftover_args = ({}, [])
        conf_files = []
        if conf_file is not None:
            conf_files.append(conf_file)
        if "conf_file" in params_cmdline.get("filename", {}):
            conf_files.extend(params_cmdline["filename"]["conf_file"])

        for conf_file in conf_files:
            acts_f = self.collectActsFromXMLFile(conf_file)
            if acts_f is None:
                acts_f = self.collectActsFromTxtFile(conf_file)
            if acts_f is None:
                print("%s is not a valid configuration file!" % conf_file)
            else:
                params_cfile = self.apply_actions(acts_f)
                preferences_mod.update(params_cfile.keys())
                pv.update(params_cfile)
        preferences_mod.update(params_cmdline.keys())
        preferences_mod.difference(self.outer_arguments)
        pv.update(params_cmdline)
        return pv, leftover_args, preferences_mod

    def collectActsFromTxtFile(self, conf_file):
        action_tuples = []
        if isinstance(conf_file, str):
            fp = open(conf_file)
            filepos = None
        else:
            fp = conf_file
            filepos = fp.tell()
        for line in fp:
            if not re.match("#", line):
                parts = line.strip().split()
                if len(parts) > 0:
                    option_string = "--"+parts[0]
                    if option_string in self._option_string_actions:
                        if len(parts) == 1:
                            action_tuples.append((self._option_string_actions[option_string], [], option_string))
                        else:
                            for p in parts[1:]:
                                action_tuples.append((self._option_string_actions[option_string], [p], option_string))
                    else:
                        raise Warning("Ignore (probably outdated) argument [%s]" % parts[0])
        if filepos is not None:
            fp.seek(filepos)
        else:
            fp.close()
        return action_tuples

    def collectActsFromXMLFile(self, conf_file):
        try:
            doc = toolXML.parseXML(conf_file)
        except Exception as instX:
            return None

        if doc is not None:
            action_tuples = []
            try:
                for current in doc.documentElement.getElementsByTagName("parameter"):
                    name = toolXML.getTagData(current, "name")
                    values = toolXML.getValues(current)

                    option_string = "--"+name
                    if option_string in self._option_string_actions:
                        if len(values) == 0:
                            action_tuples.append((self._option_string_actions[option_string], [], option_string))
                        else:
                            for p in values:
                                action_tuples.append((self._option_string_actions[option_string], [p], option_string))
            except Exception as instX:
                return None
            return action_tuples

    # def readParametersDict(self, params_dict, type_cast=False):
    #     pv = {}
    #     non_matched = {}
    #     tt = self.pm.getItemsSidesTypes()
    #     bb_matched = {}
    #     for tk, t in tt.items():
    #         if not t[0] in bb_matched:
    #             bb_matched[t[0]] = []
    #         bb_matched[t[0]].append((tk, t[1], t[2]))

    #     for name, values in params_dict.items():
    #         item = self.pm.getItem(name)
    #         if item is not None:
    #             self.prepareItemVal(pv, item, name, values, type_cast=type_cast)
    #         else:
    #             tmp = self.pm.getNameSidesTypes(name)
    #             if tmp[0] in bb_matched:
    #                 non_matched[name] = tmp
    #     if len(non_matched):
    #         fill_matched = {}
    #         kk = sorted(non_matched.keys())
    #         for k in kk:
    #             cbasis, cside, ctyp = non_matched[k]
    #             for elem in bb_matched[cbasis]:
    #                 if (elem[1] == cside or elem[1] == -1 or -1 == cside) and (elem[2] == ctyp or elem[2] == -1 or -1 == ctyp):
    #                     fill_matched[elem[0]] = k
    #         # print("fill_matched", fill_matched)
    #         for name, nn in fill_matched.items():
    #             values = params_dict[nn]
    #             item = self.pm.getItem(name)
    #             if item is not None:
    #                 self.prepareItemVal(pv, item, name, values)
    #     return pv


def getPathsForConfDefs(conf_defs=None, wdata=True):
    conf_defs = fillNamesForConfDefs(conf_defs, wdata)
    pref_dir = os.path.dirname(os.path.abspath(__file__))
    conf_paths = ["%s/%s_confdef.xml" % (pref_dir, c) if c in conf_names else c for c in conf_defs]
    return conf_paths


def getPreferencesManager(conf_defs=None):
    conf_paths = getPathsForConfDefs(conf_defs)
    return PreferencesManager(conf_paths)


def getPreferencesReader(conf_defs=None, tasks=[], tasks_default=None, description=None, epilog=None):
    return PreferencesReader(getPreferencesManager(conf_defs), tasks, tasks_default, description, epilog)


def fillParamsDocRST(paramdoc_file_in, conf_files, conf_patt, insert_patt="^\.\. ", before_block=None, after_block=None):
    if before_block is None:
        before_block = []
    if after_block is None:
        after_block = []

    pr = getPreferencesReader(conf_defs=conf_files)
    map_files = {}
    for filename in conf_files:
        map_files[os.path.basename(filename)] = filename

    content_lines = []
    current_block = []
    with open(paramdoc_file_in) as fp:
        for line in fp:
            tmp = re.search(insert_patt, line)
            if tmp is not None:
                content_lines.extend(current_block)
                current_block = []
            else:
                for match in re.finditer(conf_patt, line):
                    conf_filter = match.group("conf_filter")
                    if conf_filter in map_files:
                        blck = pr.dispRST(conf_filter=[conf_filter]).split("\n")
                        if len(blck) > 0:
                            current_block.extend(before_block)
                            current_block.extend(blck)
                            current_block.extend(after_block)
                        del map_files[conf_filter]
            content_lines.append(line.rstrip())

    if len(map_files) > 0:
        blck = pr.dispRST(conf_filter=list(map_files.keys())).split("\n")
        if len(blck) > 0:
            current_block.extend(before_block)
            current_block.extend(blck)
            current_block.extend(after_block)
    content_lines.extend(current_block)

    return content_lines


if __name__ == "__main__":

    pref_dir = os.path.dirname(os.path.abspath(__file__))
    conf_defs = glob.glob(pref_dir + "/*_confdef.xml")
    pr = getPreferencesReader(conf_defs)
    params, leftover_args, preferences_mod = pr.getPreferences(sys.argv)
    print(params)
