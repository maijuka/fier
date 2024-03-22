try:
    import toolXML
    from classSParts import cmp_fun_map

except ModuleNotFoundError:
    from . import toolXML
    from .classSParts import cmp_fun_map
import pdb

def flipValue(v):
    nv = v
    if type(v) is bool:
        v = not nv
    elif nv is not None:
        try:
            v = -nv
        except TypeError:
            try:
                v = reversed(nv)
            except TypeError:
                v = nv
    return v

def doComparison(vA, cop, vB):
    if cop in cmp_fun_map:
        return cmp_fun_map[cop](vA, vB)
    return False

def getItemVal(exp, item, other=None, constraints=None, details={}, swap=False):
    if other is not None and swap:
        v = other.getExpProp(exp, details)
        return v, "[%s:%s:]%s" % (other.getShortId(), exp, v)
    if item is not None:
        v = item.getExpProp(exp, details)
        return v, "[%s:%s:]%s" % (item.getShortId(), exp, v)

def getPairVal(exp, item, other=None, constraints=None, details={}, swap=False):
    if other is None:
        raise Exception("Pair comparison missing other item")
    if item is not None and other is not None:
        if swap:
            (item, other) = (other, item)
        v = item.getExpPairProp(other, exp, details)
        return v, "[%s:%s:%s]%s" % (item.getShortId(), exp, other.getShortId(), v)

class ACstrValue:
    def_atts = {"flip": False, "type": None}
    dynamic_types = {"item": getItemVal,
                     "pair": getPairVal}
    def_types_parity = 0
    types_parity = {"item": 1,
                    "pair": 2}

        
    def __init__(self, value, cstrs=None, atts=None, is_default=False, value_type=None):
        self.is_default = is_default
        if cstrs is None:
            self.cstrs = {}
        else:
            self.cstrs = cstrs
        if atts is None:
            self.atts = {}
        else:
            self.atts = atts
        if toolXML.isXMLNode(value):
            self.parseNode(value, value_type=value_type)
        else:
            self.value = value
        if value_type is not None and self.atts.get("type") is None:
            self.atts["type"] = value_type

    def hasConstraints(self):
        return len(self.cstrs) > 0
    def getSimpleCstr(self):
        if not self.isDynamicType() and len(self.cstrs) == 1:
            c = list(self.cstrs.keys())[0]
            if self.value == "{"+c+"}":
                return c
        return None
    def isSimpleCstr(self):
        return self.getSimpleCstr() is not None

    def hasUnfilledConstraints(self):
        return any([v is None for v in self.cstrs.values()])
    def isDefault(self):
        return self.is_default
    def getType(self):
        return self.atts.get("type")
    def getParity(self):
        return self.types_parity.get(self.getType(), self.def_types_parity)
    def isDynamicType(self):
        return self.getType() in self.dynamic_types
    def getDynamicFunction(self):
        return self.dynamic_types.get(self.getType())
    def getAtt(self, k):
        return self.atts.get(k)
            
    def __repr__(self):
        return "%s(%s, %s, %s, %s)" % (self.__class__.__name__, repr(self.value), repr(self.cstrs), repr(self.atts), self.is_default)
    
    def parseNode(self, xml_node, value_type=None):
        self.atts = toolXML.getAttsDict(xml_node, self.def_atts)
        if "type" in self.atts and self.atts["type"] not in self.dynamic_types:
            try:
                att_vtype = eval(self.atts["type"])
            except NameError:
                del self.atts["type"]
            else:                
                if value_type is None and type(att_vtype) is type:
                    self.atts["type"] = att_vtype
                    value_type = att_vtype
                    
        pieces = []
        for child in toolXML.children(xml_node):
            if toolXML.isElementNode(child):
                if toolXML.tagName(child) == "cstr":
                    cstr = toolXML.getChildrenText(child)
                    pieces.append("{"+cstr+"}")
                    self.cstrs[cstr] = None
            else:
                txt = toolXML.getNodeText(child)
                if txt is not None:
                    pieces.append(txt)
        self.value = "".join(pieces)
        if len(self.cstrs) == 0 and value_type is not None:
            self.value = toolXML.parseToType(self.value, value_type=value_type)

    def copy(self, constraints=None, override=True, cstrs=None):
        if cstrs is None:
            cstrs = dict(self.cstrs)
        if constraints is not None:
            self.fillCstrs(constraints=constraints, override=override, cstrs=cstrs)
        return self.__class__(self.value, cstrs, dict(self.atts), self.is_default)
    
    def fillCstrs(self, constraints, override=True, cstrs=None):
        if cstrs is None: ## fill in-place
            cstrs = self.cstrs
        for c in self.cstrs.keys():
            v = constraints.getCstr(c)
            if v != cstrs.get(c) and (cstrs.get(c) is None or override):
                cstrs[c] = v
        return cstrs
    def getFilledCstrs(self, constraints, override=True):
        cstrs = dict(self.cstrs)
        return self.fillCstrs(constraints=constraints, override=override, cstrs=cstrs)

    def getRawValue(self, constraints=None, override=True):
        if not self.isDynamicType() and not self.hasConstraints():
            return self.value
        c = self.getSimpleCstr()
        if c is not None: ## a simple cstr, not a complex expression, return the corresponding value
            if constraints is not None:
                v = constraints.getCstr(c)
                if v != self.cstrs[c] and (self.cstrs[c] is None or override):
                    return v
            return self.cstrs[c]
        return self.getExp(constraints=constraints, override=override)
    
    def getExp(self, constraints=None, override=True):
        if constraints is None:
            cstrs = self.cstrs
        else:
            cstrs = self.getFilledCstrs(constraints, override=override)
        return self.value.format(**cstrs)

    def evaluateAndTrack(self, item=None, other=None, constraints=None, swap=False):
        if self.isDynamicType():
            dfun = self.getDynamicFunction()
            exp = self.getExp(constraints)
            v, t = dfun(exp, item, other, constraints, swap=swap)
            if self.getAtt("flip"):
                v = flipValue(v)
                t = "NOT "+t
            return v, t
        else:
            v = self.getRawValue(constraints)
            if self.getAtt("flip"):
                v = flipValue(v)
            if v is not None:
                return v, "%s" % v
        return None, "??"

    def evaluate(self, item=None, other=None, constraints=None, swap=False):
        return self.evaluateAndTrack(item, other, constraints, swap)[0]
    
class AOperator(ACstrValue):
    def_operator = "="

    def __init__(self, value=None, cstrs={}, atts={}, is_default=False):
        ACstrValue.__init__(self, value or self.def_operator, cstrs={}, atts={}, is_default=is_default or (value is None))

class AValue(ACstrValue):
    pass

class AComparison:
    
    def __init__(self, valueA, operator=None, valueB=None):
        self.valueA = valueA
        if operator is None:
            self.operator = AOperator()
        else:
            self.operator = operator
        if valueB is None:
            self.valueB = valueA.copy()
        else:
            self.valueB = valueB

    def hasUnfilledConstraints(self):
        return self.valueA.hasUnfilledConstraints() or self.operator.hasUnfilledConstraints() or self.valueB.hasUnfilledConstraints()

    def getParity(self):
        return min(2, self.valueA.getParity()+self.valueB.getParity())
    
    def __repr__(self):
        return "%s(%s, %s, %s)" % (self.__class__.__name__, repr(self.valueA), repr(self.operator), repr(self.valueB))

    def copy(self, constraints=None, override=True, cstrs=None):
        return self.__class__(self.valueA.copy(constraints=constraints, override=override, cstrs=cstrs),
                              self.operator.copy(constraints=constraints, override=override, cstrs=cstrs),
                              self.valueB.copy(constraints=constraints, override=override, cstrs=cstrs))

    def evaluate(self, item=None, other=None, constraints=None, swap=False):
        vA = self.valueA.evaluate(item, other, constraints, swap=swap)
        vB = self.valueB.evaluate(item, other, constraints, swap=not swap)
        cop = self.operator.evaluate(constraints = constraints)
        return doComparison(vA, cop, vB)

    def evaluateAndTrack(self, item=None, other=None, constraints=None, swap=False):
        vA, tA = self.valueA.evaluateAndTrack(item, other, constraints, swap=swap)
        vB, tB = self.valueB.evaluateAndTrack(item, other, constraints, swap=not swap)
        cop, tcop = self.operator.evaluateAndTrack(constraints = constraints)
        v = doComparison(vA, cop, vB)
        return v, "%s%s?%s -> %s" % (tA, tcop, tB, v)
    
def parseValueNode(xml_node):
    return AValue(xml_node)

def parseComparisonNode(xml_node):
    kwargs = {}
    for child in toolXML.children(xml_node):
        if kwargs is not None and toolXML.isElementNode(child):
            if toolXML.tagName(child) == "value":
                v = AValue(child)
                if v is not None:
                    if "valueA" in kwargs:
                        if "valueB" in kwargs:
                            kwargs = None
                        else:
                            kwargs["valueB"] = v
                    else:
                        kwargs["valueA"] = v
            elif toolXML.tagName(child) == "operator":
                v = AOperator(child)
                if v is not None:
                    if "operator" in kwargs:
                        kwargs = None
                    else:
                        kwargs["operator"] = v
    if kwargs is not None and "valueA" in kwargs:
        return AComparison(**kwargs)


def all_subclasses(cls):
    return cls.__subclasses__() + [g for s in cls.__subclasses__() for g in all_subclasses(s)]    

class ActionFactory(object):
    
    top_class = None
    substitute_comp_classes = []
    
    @classmethod
    def getActionClass(tcl, what, top_class=None):
        if top_class is None:
            top_class = tcl.top_class
        for c in all_subclasses(top_class):
            if c.what == what:
                return c
        return None    

    @classmethod
    def initFromXml(tcl, xml_node=None):
        what = toolXML.getAttData(xml_node, "what")
        c = tcl.getActionClass(what)
        if c is not None:
            return c(xml_node=xml_node)
        return None

    @classmethod
    def getSubstituteActionClass(tcl, from_ac, to_what):
        for sc in tcl.substitute_comp_classes:
            if isinstance(from_ac, sc):
                return tcl.getActionClass(to_what, top_class=sc)
        return None

    
class Action(object):

    what = "?"
    def_args = {"reverse": False}
    block_types = {"value": parseValueNode,
                   "comparison": parseComparisonNode}
    action_parity = None
        
    
    def __init__(self, xml_node=None, args=None, blocks=None):
        if args is None:
            self.args = {}
        else:
            self.args = args
        if blocks is None:
            self.blocks = []
        else:
            self.blocks = blocks

        if xml_node is not None:
            if type(xml_node) is str:
                xxml_node = toolXML.parseStr(xml_node)
                self.parseActionNode(xxml_node)
            else:
                self.parseActionNode(xml_node)
        self.setSrc(xml_node)

    def getSrc(self):
        return self.src
    def setSrc(self, xml_node):
        self.src = None
        if type(xml_node) is str:
            self.src = xml_node
        elif toolXML.isXMLNode(xml_node):
            self.src = toolXML.getSrc(xml_node)

    def getWhat(self):
        return self.what
    def isWhat(self, what):
        return self.what == what

    def getDefArgsItems(self):
        return self.def_args.items()

    def getNbBlocks(self):
        return len(self.blocks)
        
    def parseBlockNode(self, xml_node):
        if toolXML.isXMLNode(xml_node) and toolXML.isElementNode(xml_node) and toolXML.tagName(xml_node) in self.block_types:
            b = self.block_types[toolXML.tagName(xml_node)](xml_node)
            if self.action_parity is not None and b.getParity() > self.action_parity:
                raise Warning("Too high parity for action %s! parity %s (%s)" % (self.getWhat(), b.getParity(), b))
            return b
        return None

    def parseActionNode(self, xml_node=None):
        if toolXML.isXMLNode(xml_node):
            for k, dv in self.getDefArgsItems():
                attv = toolXML.getAttData(xml_node, k, value_type=type(dv))
                if attv is not None:
                    self.args[k] = AValue(attv, value_type=type(dv))
                else:
                    tagn = toolXML.getTagNode(xml_node, k)
                    if tagn is not None:
                        self.args[k] = AValue(tagn, value_type=type(dv))
                    elif k not in self.args:
                        self.args[k] = AValue(dv, is_default=True, value_type=type(dv))
            for child in toolXML.children(xml_node):
                b = self.parseBlockNode(child)
                if b is not None:
                    self.blocks.append(b)

    def __repr__(self):
        return "%s(%s, %s)" % (self.__class__.__name__, repr(self.args), repr(self.blocks))
    
    def copy(self, constraints=None, override=True, to_class=None):
        blocks = [b.copy(constraints=constraints, override=override) for b in self.blocks]
        args = dict([(k, v.copy(constraints=constraints, override=override)) for (k,v) in self.args.items()])
        if to_class is None:
            to_class = self.__class__
        new = to_class(args=args, blocks=blocks)
        new.setSrc(self.getSrc())
        return new

    def fillCstrs(self, constraints, override=True, cstrs=None):
        for b in self.blocks:
            b.fillCstrs(constraints, override=override)
        for k in self.args.keys():
            self.args[k].fillCstrs(constraints, override=override)
    
    def getSubstitutionClass(self, to_what):
        return ActionFactory.getSubstituteActionClass(self, to_what)

    def getSubstitutionWhat(self, to_what):
        sc = ActionFactory.getSubstituteActionClass(self, to_what)
        if sc is not None:
            return sc.getWhat()
        return None
    def doSubstitution(self, to_what):
        to_class = self.getSubstitutionClass(to_what)
        if to_class is not None:
            return self.copy(to_class=to_class)
        return None

    def getArg(self, k, constraints=None, override=True):
        if k in self.args:
            return self.args[k].getRawValue(constraints=constraints, override=override)
        elif k in self.def_args:
            return self.def_args[k]

    ### WARNING "reverse" argument must be taken into account in the parent calling function of evaluate (e.g. when sorting)
    ### but is already taken into account in computeSatisfy
    def evaluateAndTrack(self, item=None, other=None, constraints=None, swap=False):
        return [b.evaluateAndTrack(item, other=other, constraints=constraints, swap=swap) for b in self.blocks]
    
    def evaluate(self, item=None, other=None, constraints=None, swap=False):
        return [b.evaluate(item, other=other, constraints=constraints, swap=swap) for b in self.blocks]

    def computeSatisfyAndTrack(self, item=None, other=None, constraints=None, swap=False):
        reverse = self.getArg("reverse", constraints)
        ts = []
        for bi, b in enumerate(self.blocks):
            v, t = b.evaluateAndTrack(item, other=other, constraints=constraints, swap=swap)
            if not v:
                return reverse, [(bi, t)]
            ts.append((bi, t))
        return not reverse, ts    

    def computeSatisfy(self, item=None, other=None, constraints=None, swap=False):
        return self.computeSatisfyAndTrack(item, other, constraints, swap)[0]
    
class ActionUnary(Action):
    action_parity = 1


class ActionBinary(Action):
    action_parity = 2
        
ActionFactory.top_class = Action
ActionFactory.substitute_comp_classes = [ActionUnary, ActionBinary]


# sed '
# s:actionlist\t\(.*\)$:</actionlist>\n<actionlist><name>\1</name>:;
# s:^list\t\(.*\)$:<list><name>\1</name></list>:;
# s:^\([^\t]*\)\t\(.*\)$:<action what=\"\1\">\2</action>:;
# s!\([a-zA-Z0-9_][a-zA-Z0-9_]*\)\=\([0-9a-zA-Z_:][0-9a-zA-Z_:]*\)[ \t]!<\1>\2</\1>!g;
# s!PAIR:\([^ <>=\t]*\)!<value type=\"pair\">\1</value>!g;
# s!ITEM:\([^ <>=\t]*\)!<value type=\"item\">\1</value>!g;
# s!\(CSTR:[^: <\t]*\)!<value>\1</value>!g;
# s%</value>\([<>=][<>=]*\)<value>%</value><operator><![CDATA[\1]]></operator><value>%g;
# s!><\([a-z]\)!>\n<\1!g;
# s!CSTR:\([^: <\t]*\)!<cstr>\1</cstr>!g;
# s!\t!\n!g' actions_rdefs_basic_alt.txt > actions_rdefs_basic_alt.xml
