import xml.dom.minidom
import re
import pdb

BOOL_MAP = {'true': True, 'false': False, 't': True, 'f': False, '0': False, '0.0': False, '1': True, 0:False, 1:True, None: None}

def getSrc(node):
    return node.toxml()

def parseStr(xml_str):
    return xml.dom.minidom.parseString(xml_str)

def parseXML(filename):
    try:
        doc = xml.dom.minidom.parse(filename)
    except Exception as inst:
        doc = None
        print("File %s could not be read! (%s)" % (filename, inst))
    return doc

def isXMLNode(node):
    return isinstance(node, xml.dom.Node)

def isElementNode(node):
    return node.nodeType == xml.dom.Node.ELEMENT_NODE

def tagName(node):
    return node.tagName

def children(node):
    return node.childNodes

def subTags(node, tag):
    return node.getElementsByTagName(tag)

def getNodeText(node):
    if node.nodeType == node.TEXT_NODE:
        return node.data
    if node.nodeType == node.CDATA_SECTION_NODE:
        return node.data.strip("\"")    

def getNodesText(nodelist):
    return "".join([getNodeText(node) or "" for node in nodelist])

def getChildrenText(element):
    return getNodesText(element.childNodes)

def parseToType(raw_value, value_type=None):
    if value_type is not None and type(raw_value) != value_type:
        try:
            if value_type is bool and raw_value in BOOL_MAP:
                return BOOL_MAP.get(raw_value)            
            raw_value = value_type(raw_value)
        except ValueError:
            raw_value = None
    return raw_value

def getTagNode(node, tag, check_uniq=True):
    tag_nodes = []
    for child in children(node):
        if isElementNode(child) and tagName(child) == tag:
            tag_nodes.append(child)
    if len(tag_nodes) == 1 or (len(tag_nodes) > 1 and not check_uniq):
        return tag_nodes[0]
    return None

def getTagData(node, tag, value_type=None, default_value=None, check_uniq=True):
    tag_node = getTagNode(node, tag, check_uniq=check_uniq)
    if tag_node is not None:
        return parseToType(getChildrenText(tag_node), value_type)
    return default_value

def getAttsDict(node, def_atts=None, set_defaults=False):    
    if set_defaults and def_atts is not None:
        atts = dict(def_atts)
    else:
        atts = {}
    if node.hasAttributes():
        for k, v in node.attributes.items():
            if def_atts is not None and k in def_atts and def_atts[k] is not None:
                v = parseToType(v, type(def_atts[k]))
            atts[k] = v
    return atts

def getAttData(node, att, value_type=None, default_value=None):
    if att in node.attributes:
        return parseToType(node.attributes[att].value, value_type)
    return default_value

def getValue(node, value_type=None):
    return getTagData(node, "value", value_type)

def getValues(node, value_type=None, tag_name="value"):
    values = []
    for valuen in node.getElementsByTagName(tag_name):
        tmp = parseToType(getChildrenText(valuen), value_type)
        if tmp is not None:
            values.append(tmp)
    return values

def parseValues(node, sep=",", strip=True, content_type=None):
    vs = []
    if re.match("values?$", tagName(node)):
        if tagName(node) == "values":
            vs = getChildrenText(node).split(node.attributes.get("sep", sep))
        else:
            vs = [getChildrenText(node)]
        if node.attributes.get("strip", strip):
            vs = [v.strip() for v in vs]
        ct = node.attributes.get("content_type", content_type)
        vs = [parseToType(v, ct) for v in vs if len(v) > 0]
    return vs
