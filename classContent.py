import sys
import os.path
import re
import multiprocessing
# from ctypes import c_bool, c_int

try:
    from toolICList import ICList
    from toolICDict import ICDict
    from classAction import Action, ActionUnary, ActionBinary # for BatchCollection
except ModuleNotFoundError:
    from .toolICList import ICList
    from .toolICDict import ICDict
    from .classAction import Action, ActionUnary, ActionBinary # for BatchCollection
    
import pdb

#####################################################
# Classes defined here, for managing content
# class Item(object):
#       basic content item
# class Container(Item, ICList):
#       basic container item
# class StoredContainer(Container):
#       List of element that can be of different origin (read from file, results from run, traces of history, etc.)
#       ContentTable assumes the displayed content, i.e. different lists in the interface, are subclasses of StoredContainer
#
# class ContentCollection(object):
#       Content items arranged in several container collections, with corresponding list of containers, allows sorting etc.
#       items are Item instances, collections are Container instances, list of collection is itself a Container instance
#       Data is a subclass of ContentCollection
#
# class TrackedContentCollection(ContentCollection):
#       adds tracking of changes to content collection, especially, tracking actions performed on items and their results
#
# class StoredContentCollection(TrackedContentCollection):
#       containers are StoredContainer instances
#
# class BatchCollection(TrackedContentCollection):
#       allows manipulating the lists through action lists, e.g. special sequences of sorting, filtering etc.
#       typically used for storing redescriptions
#       RCollection used by Miner is a subclass of BatchCollection (and Souvenirs)
#
# class StoredRCollection(StoredContentCollection, BatchCollection):
#       combines StoredContentCollection and BatchCollection, i.e. StoredContainer instances that can be manipulated through actions
#       used for storing redescriptions that are displayed in the interface
#       redescription are stored in a StoredRCollection by the DataWrapper

#  Example use of BatchCollection
#         bc = BatchCollection(lists = {"A": lA, "B": lB, ...})
#         # where lA, lB, are list of items to initialize the collections
#         # equivalent with lists = [("A", lA), ("B", lB), ...]
#         # while {0: l0, 1: l1, ...}, [(0, l0), (1, l1), ...] and [l0, l1, ...] are all equivalent
#
#         # same as initializing empty, adding lists and then adding items to them
#         bc = BatchCollection()
#         bc.newList(iid="A")
#         for itA in lA:
#             bc.addItem(itA, "B")
#         # or
#         bc.newList(iid="B")
#         bc.addItems(lB, "B")
#
#         # more
#         bc.clearList("A") # remove the list, not the items it contains
#         l = bc.getItems() # get all stored items
#         lA = bc.getItems("A") # get all items in a list
#         l = bc.getItems(iids) # get items for ids
#         liA = bc.getIidsList("A") # get ids of all items in a list
#         item = bc.getItem(iid) # get item for id
#         # merge items from "B" into "A" following "final" actions specifications
#         ids_sel = bc.selected(constraints.getActionList("final"), ids="A", new_ids="B", trg_lid="A", constraints=constraints)
#
#         # not using containers
#         bc = BatchCollection(items)
#         ids_sel = bc.selected(constraints.getActionList("final"), constraints=constraints)
#####################################################


def strToLower(value):
    try:
        return value.lower()
    except AttributeError:
        pass
    return value


class Item(object):

    class_letter = "#"
    next_uid = -1
    step_uid = 1
    mp_lock = False
    @classmethod
    def generateNextUid(tcl):
        if tcl.mp_lock:
            with tcl.next_uid.get_lock():
                tcl.next_uid.value += tcl.step_uid
                v = tcl.next_uid.value
            return v
        tcl.next_uid += tcl.step_uid
        return tcl.next_uid

    @classmethod
    def setUidGen(tcl, nv=-1, step=1, mp_lock=False):
        if tcl.mp_lock == mp_lock and step is None and nv is None:
            return
        if nv is None:
            if tcl.mp_lock:
                nv = tcl.next_uid.value
            else:
                nv = tcl.next_uid
        if step is None:
            step = tcl.step_uid
        tcl.mp_lock = mp_lock
        if tcl.mp_lock:
            tcl.next_uid = multiprocessing.Value('i', nv)
        else:
            tcl.next_uid = nv
        tcl.step_uid = step

    def __init__(self, iid=None):
        self.setUid(iid)

    def setUid(self, iid=None):
        if iid is None:
            iid = self.generateNextUid()
        self.uid = iid

    def getUid(self):
        return self.uid

    def getShortId(self, details={}):
        if details is not None and details.get("aim") == "sort":
            return self.getUid()
        return "%s%s" % (self.class_letter, self.getUid())

    def getFieldV(self, field="uid", details={}):
        if field == "uid":
            return self.getUid()
        tmp = None
        item = self
        if (type(field) is list or type(field) is tuple) and len(field) > 1:
            methode = eval(field[1])
            if len(field) > 2 and field[2] is not None:
                details.update(field[2])
        else:
            methode = eval(field)
        if callable(methode):
            try:
                tmp = methode(details)
            except IndexError:
                methode(details)
        else:
            tmp = methode
        if tmp is None and "replace_none" in details:
            return details["replace_none"]
        return tmp

    def __str__(self):
        return "Item %s" % self.getShortId()


class Container(Item, ICList):

    def __init__(self, iid=None, name=None, data=[], isChanged=False, content_handle=None):
        Item.__init__(self, iid)
        self.resetSortP()
        if name is None:
            self.makeName()
        else:
            self.name = name
        ICList.__init__(self, data, isChanged)
        self.setContentHandle(content_handle)

    def clear(self):
        self.resetSortP()
        self.reset()

    def resetSortP(self):
        self.sortP = (None, False)

    def setContentHandle(self, content_handle):
        self.content_handle = content_handle

    def getContentHandle(self):
        return self.content_handle

    def hasContentHandle(self):
        return self.content_handle is not None

    def getItemForIid(self, iid):
        if self.hasContentHandle():
            return self.content_handle.get(iid)
        return iid

    def getItemsForIids(self, iids=None):
        if iids is None:
            iids = self
        return [self.getItemForIid(iid) for iid in iids]

    def getItemFieldV(self, iid, field, details):
        if self.hasContentHandle():
            item = self.getItemForIid(iid)
            if item is not None:
                return item.getFieldV(field, details)
        return None

    def makeName(self):
        self.name = "L#%s" % self.getUid()

    def getName(self):
        return self.name

    def setName(self, text):
        self.name = text

    def getShortStr(self):
        return "(%d%s) %s" % (len(self), "*"*(1*self.isChanged), self.getName())

    def __str__(self):
        return "%s (%d)" % (self.getName(), len(self))

    def setIids(self, data):
        ICList.__init__(self, data, isChanged=True)

    def getPosIid(self, iid):
        if iid in self:
            return self.index(iid)

    def getSortInfo(self):
        return self.sortP

    def setSort(self, colS=None, direct=None):
        old_sortP = self.sortP
        if direct is not None:
            self.sortP = (colS, direct)
        elif colS is not None:
            if self.sortP[0] == colS:
                self.sortP = (self.sortP[0], not self.sortP[1])
            else:
                self.sortP = (colS, False)
        else:
            self.sortP = (None, False)
        if self.sortP != old_sortP and self.sortP[0] is not None:
            return True
        return False

    def updateSort(self, fields=[], details={}):
        default_cmp = None
        map_v = {}
        sdetails = {}
        sdetails.update(details)
        sdetails.update({"aim": "sort", "replace_none": None})
        field = self.sortP[0]
        if type(field) is int and field < len(fields):
            field = fields[field]
        if field is not None and self.hasContentHandle():
            vs = [(x, strToLower(self.getItemFieldV(x, field, sdetails))) for x in self]
            map_v = dict([(v[0], (v[1] is None, v[1] if v[1] is not None else v[0])) for v in vs])
        else:
            map_v = dict([(p, p) for p in self])
        self.sort(key=lambda x: map_v.get(x) or default_cmp, reverse=self.sortP[1])


class StoredContainer(Container):
    NFILE = "file"
    NRUN = "run"
    NMANU = "manual"
    NHIST = "history"
    NBUFF = "buffer"
    NPACK = "package"
    STORED_TYPES_NAMES = [NFILE, NRUN, NMANU, NHIST, NBUFF]
    stored_types_map = dict([(v, k) for (k, v) in enumerate(STORED_TYPES_NAMES)])

    @classmethod
    def getDefaultSrc(tcl):
        return (tcl.NMANU, None, 0)

    @classmethod
    def getHistSrc(tcl):
        return (tcl.NHIST, None, 0)

    @classmethod
    def getBufferSrc(tcl):
        return (tcl.NBUFF, None, 0)

    @classmethod
    def getPackSrc(tcl):
        return (tcl.NFILE, tcl.NPACK, 1)

    @classmethod
    def makeSrc(tcl, src=None):
        if src is None:
            return tcl.getDefaultSrc()
        elif src == tcl.NPACK:
            return tcl.getPackSrc()
        elif src == tcl.NHIST:
            return tcl.getHistSrc()
        elif src == tcl.NBUFF:
            return tcl.getBufferSrc()
        return src

    def __init__(self, iid=None, src=None, name=None, data=[], isChanged=False, content_handle=None):
        self.src = self.makeSrc(src)
        Container.__init__(self, iid, name, data, isChanged, content_handle)

    def makeName(self):
        Container.makeName(self)
        if self.src[0] == self.NRUN:
            self.name = "%s#%s" % (self.src[1][1].title(), self.src[1][0])
        elif self.isHist():
            self.name = self.NHIST
            self.uid = self.NHIST
        elif self.isBuffer():
            self.name = self.NBUFF
            self.uid = self.NBUFF
        elif self.src[0] == self.NFILE:
            if self.inPack():
                self.name = "::"
            else:
                self.name = "//"
            if self.src[1] is not None:
                if os.path.exists(self.src[1]):
                    self.name += os.path.basename(self.src[1])
                else:
                    self.name += self.src[1]

    def updateSrc(self, src):
        self.src = self.makeSrc(src)
        self.makeName()

    def matches(self, what={}):
        match = True
        if "namepIn" in what:
            match &= (re.match(what["namepIn"], self.getName()) is not None)
        if "namepOut" in what:
            match &= not (re.match(what["namepOut"], self.getName()) is not None)

        if "srcTypIn" in what:
            srcTyp = what["srcTypIn"]
            if type(srcTyp) is list:
                match &= (self.getSrcTyp() in srcTyp)
            elif type(srcTyp) is int:
                match &= (self.getSrcTyp() == STORED_TYPES_NAMES[srcTyp])
            else:
                match &= (self.getSrcTyp() == srcTyp)
        if "srcTypOut" in what:
            srcTyp = what["srcTypOut"]
            if type(srcTyp) is list:
                match &= not (self.getSrcTyp() in srcTyp)
            elif type(srcTyp) is int:
                match &= not (self.getSrcTyp() == STORED_TYPES_NAMES[srcTyp])
            else:
                match &= not (self.getSrcTyp() == srcTyp)
        return match

    def getSrc(self):
        return self.src

    def getSrcTyp(self):
        return self.src[0]

    def getSrcTypId(self):
        return self.stored_types_map[self.src[0]]

    def getSrcInfo(self):
        return self.src[1]

    def getSrcPath(self):
        if self.hasSrcPath():
            return self.src[1]
        else:
            return None

    def getPackPath(self):
        if self.hasSrcPath():
            return self.src[1]
        else:
            return None

    def hasSrcPath(self):
        return self.src[0] == self.NFILE and not type(self.src[1]) is int and os.path.exists(self.src[1])

    def inPack(self):
        return self.src[0] == self.NFILE and self.src[2] == 1

    def isHist(self):
        return self.src[0] == self.NHIST

    def isBuffer(self):
        return self.src[0] == self.NBUFF

    def setSrc(self, src_typ, src_info=None, inpck=False):
        self.src = (src_typ, src_info, inpck)
        self.makeName()


class ContentCollection(object):

    name_class = "Content"
    container_class = Container

    def nbTracks(self):
        return -1

    def __init__(self, items=None, lists=None):
        self.items = ICDict()
        self.containers = ICDict()
        self.clist = Container()
        self.map_iids = {}
        self.keep_lids = set()
        if items is not None:
            self.initItems(items)
        elif lists is not None:
            self.initLists(lists)

    def updateSort(self, lid, fields=[], details={}):
        if lid in self.containers:
            self.containers[lid].updateSort(fields, details)
            for pos, iid in enumerate(self.containers[lid]):
                self.map_iids[iid][lid] = pos

    def initItems(self, items):
        self.clear()
        if type(items) is list:
            self.items = ICDict(enumerate(items))
        elif type(items) is dict:
            self.items = ICDict(items.items())

    def initLists(self, lists):
        self.clear()
        iter_lists = []
        if type(lists) is list:
            if len(lists) > 0 and type(lists[0]) is tuple:
                iter_lists = lists
            else:
                iter_lists = enumerate(lists)
        elif type(lists) is dict:
            iter_lists = lists.items()
        for lid, lls in iter_lists:
            self.newList(iid=lid)
            for ll in lls:
                self.addItem(ll, lid)

    def getClassName(self):
        return "%s collection" % self.name_class

    def __str__(self):
        return "%s: %d items in %d/%d lists" % (self.getClassName(), self.nbItems(), self.nbDispLists(), self.nbTotLists())
    # hex(id(self))

    def dispDetails(self):
        print("Items: ", self.items)
        print("Containers: ", self.containers)
        print("CList: ", self.clist)
        print("Map: ", self.map_iids)

    def clear(self):
        self.items.reset()
        self.map_iids.clear()
        if len(self.keep_lids) == 0:
            self.containers.reset()
            self.clist.clear()
        else:
            lids = self.containers.keys()
            for lid in lids:
                if lid not in self.keep_lids:
                    try:
                        self.clist.remove(lid)
                    except ValueError:
                        pass
                    del self.containers[lid]

    def prune(self):
        mks = self.map_iids.keys()
        for mk in mks:
            if len(self.map_iids[mk]) == 0:
                del self.map_iids[mk]
        delete_iids = [iid for iid in self.items.keys() if iid not in self.map_iids]
        for iid in delete_iids:
            del self.items[iid]
        delete_lids = [lid for lid in self.containers.keys() if self.toPrune(lid)]
        for lid in delete_lids:
            del self.containers[lid]
            try:
                self.clist.remove(lid)
            except ValueError:
                pass

    def toPrune(self, lid):
        return self.getLen(lid) == 0 and lid not in self.keep_lids
    # HANDLING OF ITEMS
    #######################

    def getIidForLidPos(self, lid, pos=-1):
        if lid in self.containers:
            if pos >= -1 and pos < len(self.containers[lid]):
                return self.containers[lid][pos]

    def popIidForLidPos(self, lid, pos=-1):
        iid = None
        if lid in self.containers:
            if pos < -1 or pos >= len(self.containers[lid]):
                pos = None
            elif pos == -1:
                pos = len(self.containers[lid])-1
            if pos is not None:
                for pos_follow in range(pos+1, len(self.containers[lid])):
                    iid_follow = self.containers[lid][pos_follow]
                    self.map_iids[iid_follow][lid] -= 1
                iid = self.containers[lid].pop(pos)
                del self.map_iids[iid][lid]
        return iid

    def getPosForLidIid(self, lid, iid):
        return self.map_iids.get(iid, {}).get(lid)

    def getLidPosForIid(self, iid):
        return self.map_iids.get(iid, {}).items()

    def setIids(self, lid, iids):
        if lid in self.containers:
            for iid in self.containers[lid]:
                self.map_iids[iid].pop(lid)
            for pos, iid in enumerate(iids):
                if iid not in self.map_iids:
                    self.map_iids[iid] = {}
                self.map_iids[iid][lid] = pos
            self.containers[lid].setIids(iids)

    def insertIidAtLidPos(self, iid, trg_lid=None, trg_pos=-1):
        pos = None
        if trg_lid is not None and trg_lid in self.containers:
            if trg_lid in self.map_iids.get(iid, {}):
                return self.map_iids[iid][trg_lid]

            if trg_pos == -1 or trg_pos > len(self.containers[trg_lid]):
                pos = len(self.containers[trg_lid])
                self.containers[trg_lid].append(iid)
            else:
                pos = trg_pos
                for pos_follow in range(trg_pos, len(self.containers[trg_lid])):
                    iid_follow = self.containers[trg_lid][pos_follow]
                    self.map_iids[iid_follow][trg_lid] += 1
                self.containers[trg_lid].insert(trg_pos, iid)
        if pos is not None:
            if iid not in self.map_iids:
                self.map_iids[iid] = {}
            self.map_iids[iid][trg_lid] = pos
        return pos

    def substituteItem(self, iid, item, backup=False):
        old_item = None
        if self.hasIid(iid):
            old_item = self.getItem(iid)
            if old_item is not None and item is not None and old_item.compare(item) != 0:
                old_item.setUid()
                item.setUid(iid)
                self.items[iid] = item
                for lid in self.map_iids.get(iid, {}).keys():
                    self.containers[lid].isChanged = True
            else:
                old_item = None
        return old_item

    def addItem(self, item, trg_lid=None, trg_pos=-1):
        if isinstance(item, Item):
            iid = item.getUid()
            # if iid not in self.items:
            self.items[iid] = item
        else:
            iid = item
        return self.insertIidAtLidPos(iid, trg_lid, trg_pos)

    def addItems(self, items, trg_lid=None, trg_pos=-1):
        return [self.addItem(item, trg_lid, trg_pos) for item in items]

    def deleteItem(self, item):
        iid = None
        removed = None
        if isinstance(item, Item):
            iid = item.getUid()
        else:
            iid = item
        if iid is not None:
            if iid in self.map_iids:
                lidsPos = self.map_iids[iid].items()
                for lid, pos in lidsPos:
                    for pos_follow in range(pos+1, len(self.containers[lid])):
                        iid_follow = self.containers[lid][pos_follow]
                        self.map_iids[iid_follow][lid] -= 1
                    self.containers[lid].pop(pos)
                del self.map_iids[iid]
            removed = self.items.pop(iid)
        return removed

    def copyIid(self, src_lid, src_pos=None, iid=None, trg_lid=None, trg_pos=-1):
        if src_pos is None and iid is None:
            return None
        if iid is not None:
            pos = getPosForLidIid(src_lid, iid)
            if src_pos is None:
                src_pos = pos
            elif src_pos != pos:
                return None
        iid = getIidForLidPos(src_lid, src_pos)
        if iid is not None:
            item = self.items.get(iid)
            if item is not None:
                item_cp = item.copy()
                return self.addItem(item_cp, trg_lid, trg_pos)

    def moveIid(self, src_lid, src_pos=None, iid=None, trg_lid=None, trg_pos=-1):
        if src_pos is None and iid is None:
            return None
        if iid is not None:
            pos = self.getPosForLidIid(src_lid, iid)
            if src_pos is None:
                src_pos = pos
            elif src_pos != pos:
                return None
        iid = self.popIidForLidPos(src_lid, src_pos)
        if iid is not None:
            return self.insertIidAtLidPos(iid, trg_lid, trg_pos)

    def copyIids(self, src_lid, src_pos=None, iids=None, trg_lid=None, trg_pos=-1):
        if src_pos is None and iids is None:
            return None
        if iids is not None:
            src_p = []
            for iid in iids:
                pos = self.getPosForLidIid(src_lid, iid)
                if pos is not None:
                    src_p.append(pos)
            if src_pos is None:
                src_pos = src_p
            elif src_pos != src_p:
                return None
        popped = []
        for sp in sorted(src_pos, reverse=True):
            popped.append(self.getIidForLidPos(src_lid, sp))
        inserted = []
        if trg_pos == -1:
            trg_pos = len(self.containers[trg_lid])
        offset = 0
        for iid in popped:
            if iid is not None:
                item = self.items.get(iid)
                if item is not None:
                    item_cp = item.copy()
                    inserted.append(self.addItem(item_cp, trg_lid, trg_pos+offset))
                    offset += 1
        return inserted

    def moveIids(self, src_lid, src_pos=None, iids=None, trg_lid=None, trg_pos=-1):
        if src_pos is None and iids is None:
            return None
        if iids is not None:
            src_p = []
            for iid in iids:
                pos = self.getPosForLidIid(src_lid, iid)
                if pos is not None:
                    src_p.append(pos)
            if src_pos is None:
                src_pos = src_p
            elif src_pos != src_p:
                return None
        popped = []
        for sp in sorted(src_pos, reverse=True):
            popped.append(self.popIidForLidPos(src_lid, sp))
        inserted = []
        if trg_pos == -1:
            trg_pos = len(self.containers[trg_lid])
        offset = 0
        for iid in popped:
            if iid is not None:
                inserted.append(self.insertIidAtLidPos(iid, trg_lid, trg_pos+offset))
                offset += 1
        return inserted

    # HANDLING OF LISTS
    #######################
    def getLidForPos(self, pos=-1):
        if pos >= -1 and pos < len(self.clist):
            return self.clist[pos]

    def popLidForPos(self, pos=-1):
        if pos >= -1 and pos < len(self.clist):
            return self.clist.pop(pos)

    def getPosForLid(self, lid):
        if lid in self.clist:
            return self.clist.index(lid)

    def insertLidAtPos(self, lid, trg_pos=-1):
        pos = None
        if trg_pos is not None:
            if trg_pos == -1 or trg_pos > len(self.clist):
                pos = len(self.clist)
                self.clist.append(lid)
            else:
                pos = trg_pos
                self.clist.insert(trg_pos, lid)
        return pos

    def addList(self, lst, trg_pos=-1):
        lid = lst.getUid()
        self.containers[lid] = lst
        self.addFromList(lst)
        self.insertLidAtPos(lid, trg_pos)
        return lid

    def deleteList(self, lst):
        lid = None
        removed = None
        if isinstance(lst, Item):
            lid = lst.getUid()
        else:
            lid = lst
        if lid is not None:
            pos = self.getPosForLid(lid)
            if pos is not None:
                self.popLidForPos(pos)
                self.clearList(lid)
                removed = self.containers.pop(lid)
        return removed

    def addFromList(self, lst):
        ll = None
        if isinstance(lst, Item):
            ll = lst
        elif lst in self.containers:
            ll = self.containers[lst]
        if isinstance(ll, Container) and len(ll) > 0:
            lid = ll.getUid()
            for pos, iid in enumerate(ll):
                if iid not in self.map_iids:
                    self.map_iids[iid] = {}
                self.map_iids[iid][lid] = pos

    def clearList(self, lst):
        ll = None
        if isinstance(lst, Item):
            ll = lst
        elif lst in self.containers:
            ll = self.containers[lst]
        if isinstance(ll, Container):
            lid = ll.getUid()
            for iid in ll:
                del self.map_iids[iid][lid]
            ll.reset()

    def getItem(self, iid):
        return self.items.get(iid)

    def getIids(self):
        return list(self.items.keys())

    def hasIid(self, iid):
        return iid in self.items

    def getComplementIids(self, iids):
        return list(set(self.items.keys()) - set(iids))

    def getIntersectIids(self, iids):
        return list(set(self.items.keys()) & set(iids))

    def nbItems(self):
        return len(self.items)

    def nbTotLists(self):
        return len(self.containers)

    def nbDispLists(self):
        return len(self.clist)

    def getIidsList(self, lid=None):
        if lid in self.containers:
            return list(self.containers[lid])
        return []

    def getIidsListAbove(self, lid, iid):
        iids = []
        if lid in self.containers:
            for iidn in self.containers[lid]:
                if iidn == iid:
                    return iids
                iids.append(iidn)
        return iids

    def getIidsListBelow(self, lid, iid):
        iids = []
        if lid in self.containers:
            found = False
            for iidn in self.containers[lid]:
                if iidn == iid:
                    found = True
                if found:
                    iids.append(iidn)
        return iids

    def getList(self, lid=None):
        if lid is None:
            return self.clist
        elif lid in self.containers:
            return self.containers[lid]

    def newList(self, **args):
        lst = self.container_class(content_handle=self.items, **args)
        return self.addList(lst, trg_pos=-1)

    def getOrdLids(self):
        return list(self.clist)

    def getLids(self):
        return self.containers.keys()

    def hasLid(self, lid):
        return lid in self.containers

    def showingLid(self, lid):
        return lid in self.clist

    def getLen(self, lid=None):
        if lid is None:
            return len(self.clist)
        elif lid in self.containers:
            return len(self.containers[lid])
        return 0

    def getItems(self, lid=None):
        if lid is None:
            return self.items.values()
        try:
            if lid in self.containers:
                return [self.getItem(iid) for iid in self.containers[lid]]
        except TypeError:  # not hashable, maybe a list of ids?..
            pass
        try:
            return [self.getItem(iid) for iid in lid]  # lid is assumed to be a list of iids
        except TypeError:  # not iterable, no idea what is being asked for...
            pass
        return []

    def lenBuffer(self):
        return 0

    def getNonEmpyLids(self, check_changed=False, inc_hist=False, inc_buffer=False):
        return []

    def getItemFieldV(self, iid, field, details):
        item = self.getItem(iid)
        if item is not None:
            return item.getFieldV(field, details)
        return None

    def getContentInfo(self, lids=None, iids=None):
        info = {"path": None, "is_buffer": False, "is_hist": False, "in_pack": False}
        lid = None
        if iids is not None:
            lid = None
        if lids is not None:
            if type(lids) is int:
                lid = lids
                lids = [lid]
            elif type(lids) is list and len(lids) == 1:
                lid = lids[0]
            if lid is not None:
                ll = self.getList(lid)
                if isinstance(ll, StoredContainer):
                    info.update({"path": ll.getSrcPath(), "is_buffer": ll.isBuffer(), "is_hist": ll.isHist(), "in_pack": ll.inPack()})
            for llid in lids:
                iids.extend(self.getIidsList(llid))
        elif iids is None:  # both are None, list everything
            iids = self.getIids()
        info.update({"lid": lid, "lids": lids, "iids": iids, "nb": len(iids)})
        return info


class TrackedContentCollection(ContentCollection):

    name_class = "Tracked content"
    to_str_fields = ["do", "action", "src", "trg", "rationales"]

    def __init__(self, items=None, lists=None, tracks_switch="OFF"):
        self.stored_tracks = []
        self.resetLatestTrackI()
        self.setTracksSwitch(tracks_switch=tracks_switch)
        ContentCollection.__init__(self, items, lists)

    def setTracksSwitch(self, tracks_switch="OFF"):
        self.tracks_switch = tracks_switch
        
    def getTracksSwitch(self):
        return self.tracks_switch

    def resetLatestTrackI(self, ti=None):
        if ti is None:
            ti = len(self.stored_tracks)
        self.last_ti = ti
        return ti
    def getLatestTrackI(self):
        return self.last_ti
    
    def getLatestTracksRange(self, update=True):
        ti = self.getLatestTrackI()
        if update:
            self.resetLatestTrackI()
        return (ti, len(self.stored_tracks))
    
    def getLatestTracks(self, update=True):
        tf, tt = self.getLatestTracksRange(update=update)
        return self.stored_tracks[tf:tt]
    
    def getTracks(self):
        return self.stored_tracks

    def importTracks(self, tracks, source=None):
        for t in tracks:
            if "@" not in t:
                t["@"] = source
                self.stored_tracks.append(t)
            else:
                print(">>> Check tracks!!")
                # pdb.set_trace()

    def clearTracks(self):
        del self.stored_tracks[:]

    def keepTrack(self, track, mute_tracking=None):
        if self.tracks_switch == "OFF" or mute_tracking:
            return False
        return True
            
    def addTrack(self, track, mute_tracking=None):
        if self.keepTrack(track, mute_tracking):
            self.stored_tracks.append(track)

    def tracksToStr(self, trange=None, line_break="\n"):
        if trange is None:
            tf, tt = (0, len(self.stored_tracks))
        elif type(trange) is bool:
            tf, tt = self.getLatestTracksRange(update=trange)
        elif type(trange) is int:
            tf, tt = (trange, len(self.stored_tracks))
        else:
            tf, tt = trange

        xps = "---- NO TRACKS"
        if tt > tf:            
            xps = "---- TRACKS:"
            for ti in range(tf, tt):
                t = self.stored_tracks[ti]
                at = ""
                if "@" in t:
                    at = " {@%s}" % t["@"]
                xps += "%s%d%s\t%s" % (line_break, ti, at, "\t".join(["%s: %s" % (k, t[k]) for k in self.to_str_fields if k in t]))
        return xps

    def nbTracks(self):
        return len(self.stored_tracks)

    def __str__(self):
        return "%s, %d tracks" % (ContentCollection.__str__(self), self.nbTracks())

    def dispDetails(self):
        ContentCollection.dispDetails(self)
        print("Tracks: ", self.nbTracks())

    def clear(self):
        self.clearTracks()
        ContentCollection.clear(self)


class StoredContentCollection(TrackedContentCollection):

    name_class = "Stored content"
    container_class = StoredContainer

    def __init__(self):
        TrackedContentCollection.__init__(self)
        self.map_src_ids = {}
        # prepare buffer list
        lstBuff = self.container_class(src=StoredContainer.NBUFF)
        self.addList(lstBuff, trg_pos=None)
        # prepare history list
        lstHist = self.container_class(src=StoredContainer.NHIST)
        self.addList(lstHist, trg_pos=-1)
        self.special_uids = {StoredContainer.NBUFF: lstBuff.getUid(), StoredContainer.NHIST: lstHist.getUid()}
        self.keep_lids = set(self.special_uids.values())
        self.buff_pos = []

    def clear(self):
        ContentCollection.clear(self)
        srcs = self.map_src_ids.items()
        for src, lid in srcs:
            if not self.hasLid(lid):
                del self.map_src_ids[src]

    def getUidForSrc(self, src):
        return self.map_src_ids.get(src)

    def addList(self, lst, trg_pos=-1):
        ContentCollection.addList(self, lst, trg_pos)
        lid = lst.getUid()
        self.map_src_ids[lst.getSrc()] = lid
        return lid

    def getBufferLid(self):
        return self.special_uids[StoredContainer.NBUFF]

    def cutIid(self, src_lid, src_pos=None, iid=None):
        self.buff_pos.append([-1])
        return self.moveIid(src_lid, src_pos, iid, trg_lid=self.getBufferLid(), trg_pos=-1)

    def pasteIid(self, trg_lid, trg_pos=-1):
        src_pos = None
        if len(self.buff_pos) > 0:
            if len(self.buff_pos[-1]) == 1:
                src_pos = self.buff_pos.pop()[0]
            elif len(self.buff_pos[-1]) > 1:
                src_pos = self.buff_pos[-1].pop()
        if src_pos is not None:
            return self.moveIid(src_lid=self.getBufferLid(), src_pos=src_pos, trg_lid=trg_lid, trg_pos=trg_pos)

    def copyToBufferIid(self, src_lid, src_pos=None, iid=None):
        self.buff_pos.append([-1])
        return self.copyIid(src_lid, src_pos, iid, trg_lid=self.getBufferLid(), trg_pos=-1)

    def cutIids(self, src_lid, src_pos=None, iids=None):
        pos = self.moveIids(src_lid, src_pos, iids, trg_lid=self.getBufferLid(), trg_pos=-1)
        self.buff_pos.append(pos)
        return pos

    def pasteIids(self, trg_lid, trg_pos=-1):
        src_pos = None
        if len(self.buff_pos) > 0:
            src_pos = self.buff_pos.pop()
        if src_pos is not None:
            return self.moveIids(src_lid=self.getBufferLid(), src_pos=src_pos, trg_lid=trg_lid, trg_pos=trg_pos)

    def copyToBufferIids(self, src_lid, src_pos=None, iids=None):
        pos = self.copyIids(src_lid, src_pos, iids, trg_lid=self.getBufferLid(), trg_pos=-1)
        self.buff_pos.append(pos)
        return pos

    def lenBuffer(self):
        return len(self.containers.get(self.getBufferLid(), []))

    def clearBuffer(self):
        self.buff_pos = []
        return self.clearList(self.getBufferLid())

    def getHistLid(self):
        return self.special_uids[StoredContainer.NHIST]

    def addItemToHist(self, item):
        self.addItem(item, self.getHistLid())

    def copyToHistIid(self, src_lid, src_pos=None, iid=None):
        if src_lid != self.getHistLid():
            return self.copyIid(src_lid, src_pos=None, iid=None, trg_lid=self.getHistLid(), trg_pos=-1)

    def lenHist(self):
        return len(self.containers.get(self.getHistLid(), []))

    def clearHist(self):
        return self.clearList(self.getHistLid())

    def insertLidAtPos(self, lid, trg_pos=-1):
        pos = None
        if trg_pos is not None:
            if trg_pos == -1 or trg_pos > len(self.clist):
                if len(self.clist) != 0 and self.clist[-1] == self.getHistLid() and lid != self.getHistLid():
                    pos = len(self.clist)-1
                    self.clist.insert(pos, lid)
                else:
                    pos = len(self.clist)
                    self.clist.append(lid)
            else:
                pos = trg_pos
                self.clist.insert(trg_pos, lid)
        return pos

    def substituteItem(self, iid, item, backup=False):
        old_item = ContentCollection.substituteItem(self, iid, item, backup)
        if old_item is not None and backup:
            self.addItemToHist(old_item)
        return old_item

    def appendItemSrc(self, item, trg_src):
        trg_lid = self.getUidForSrc(trg_src)
        if trg_lid is None:
            trg_lid = self.newList(src=trg_src)
        return self.addItem(item, trg_lid)

    def updateSrc(self, from_src, to_src):
        lid = self.getUidForSrc(from_src)
        if lid is not None:
            self.containers[lid].updateSrc(to_src)
            del self.map_src_ids[from_src]
            self.map_src_ids[to_src] = lid

    def isListToPack(self, lid):
        if lid is not None:
            ll = self.getList(lid)
            if ll is not None:
                return ll.inPack()
        return False

    def addDelListToPack(self, lid):
        src, new_src = (None, None)
        if lid is not None:
            ll = self.getList(lid)
            if ll is not None:
                src = ll.getSrc()
                if ll.inPack():
                    if os.path.exists(src[1]):
                        new_src = ('file', src[1], 0)
                    else:
                        new_src = ('manual', None, 0)
                else:
                    if src[0] == "file":  # and os.path.exists(src[1]):
                        new_src = ('file', src[1], 1)
                    else:
                        new_src = ('file', "%s.csv" % ll.name, 1)
        if new_src is not None:
            self.updateSrc(src, new_src)
        return new_src

    def getNonEmpyLids(self, check_changed=False, inc_hist=False, inc_buffer=False):
        changed = []
        for lid in self.getLids():
            if (inc_hist or lid != self.getHistLid()) and (inc_buffer or lid != self.getBufferLid()):
                if len(self.getList(lid)) and (not check_changed or self.getList(lid).isChanged):
                    changed.append(lid)
        return changed


class BatchCollection(TrackedContentCollection):

    name_class = "Batch"
    actions_methods = {}

    def prepareIids(self, ids=None):
        iids = []
        if ids is None:
            iids = self.getIids()
        elif type(ids) is list:
            iids = [iid for iid in ids if self.hasIid(iid)]
        elif self.hasLid(ids):
            iids = self.getIidsList(ids)
        # if complement:
        #     iids = self.getComplementIids(iids)
        return iids
    
    def applyFunction(self, action, i, data=None, constraints=None, mute_tracking=None, changes=False):
        item = self.getItem(i)
        xps_item, modi = action.applyFunction(item, data, constraints)
        if modi and changes:
            self.items._isChanged = True
        if modi:
            track = {"do": "apply", "trg": [i], "action": action.getWhat()}
            self.addTrack(track, mute_tracking)
        return xps_item, modi

    def computeSatisfy(self, action, i, otheri=None, constraints=None, mute_tracking=None):
        item = self.getItem(i)
        other = None
        if otheri is not None:
            other = self.getItem(otheri)
        out, ts = action.computeSatisfyAndTrack(item, other=other, constraints=constraints)            
        track = {"do": "satisfy", "trg": [i], "out": out, "rationales": ts, "action": action.getWhat()}
        if otheri is not None:
            track["src"] = [otheri]            
        self.addTrack(track, mute_tracking)
        return out
    
    ### corresponding Action subclass
    class ActionApplyGeneric(ActionUnary):
        def_args = {"function": "identity", "reverse": False}
    
        def applyFunction(self, item, data=None, constraints=None):
            funct = self.getArg("function", constraints)
            xps_item, modi = None, False
            if callable(funct):
                if data is not None:
                    xps = funct(item, data)
                else:
                    xps = funct(item)
            if type(funct) == str and len(funct) > 0:
                if funct == "identity":
                    xps = item
                try:
                    xps = eval("item."+funct)
                except TypeError:
                    xps = None
    
            if callable(xps):
                f = xps
                if data is not None:
                    xps = f(data)
                else:
                    xps = f()
    
            if isinstance(xps, Item):
                xps_item, modi = xps, True
            elif xps is not None:
                try:
                    if len(xps) >= 2 and isinstance(xps[0], Item) and type(xps[1]) is bool:
                        xps_item, modi = xps[0], xps[1]
                except (ValueError, IndexError):
                    xps_item, modi = xps, False
            return xps_item, modi

    
    def applyBulkIds(self, action, ids, new_ids=None, data=None, constraints=None, mute_tracking=None, changes=False):
        xps_ids = []
        for i in ids:
            xps_ids.append(i)
            xps_item, modi = self.applyFunction(action, i, data=data, constraints=constraints, mute_tracking=True, changes=changes)
            if xps_item is not None and modi:
                self.addItem(xps_item)
                xps_ids.append(xps_item.getUid())
        track = {"do": "apply", "trg": list(xps_ids), "action": action.getWhat()}
        self.addTrack(track, mute_tracking)                
        return xps_ids
    actions_methods["applyBulk"] = applyBulkIds

    ### corresponding Action subclass
    class ActionApplyBulk(ActionApplyGeneric):
        what = "applyBulk"

    def applyIds(self, action, ids, new_ids=None, data=None, constraints=None, mute_tracking=None, changes=False):
        org_ids = ids
        if action.getNbBlocks() == 0:
            if action.getArg("reverse", constraints):
                ids = self.getComplementIids(ids)
        else:
            ids = [i for i in ids if self.computeSatisfy(action, i, constraints=constraints, mute_tracking=True)]
        for i in ids:
            self.applyFunction(action, i, data=data, constraints=constraints, mute_tracking=True, changes=changes)
        track = {"do": "apply", "trg": list(ids), "action": action.getWhat()}
        self.addTrack(track, mute_tracking)
        return ids
    actions_methods["apply"] = applyIds
    
    ### corresponding Action subclass
    class ActionApply(ActionApplyGeneric):
        what = "apply"

    def sortIds(self, action, ids, new_ids=None, data=None, constraints=None, mute_tracking=None):
        # sort
        if len(ids) > 1:
            track = {"do": "sort", "src": list(ids), "action": action.getWhat()}
            if action.getNbBlocks() > 0:
                map_vids = dict([(iid, tuple(action.evaluate(self.getItem(iid), constraints=constraints))) for iid in ids])
                ids.sort(key=lambda x: map_vids[x], reverse=action.getArg("reverse", constraints))
            elif action.getArg("reverse", constraints):
                ids.reverse()
            track["trg"] = list(ids)
            self.addTrack(track, mute_tracking)
        return ids
    actions_methods["sort"] = sortIds

    ### corresponding Action subclass
    class ActionSort(ActionUnary):
        what = "sort"
        def_args = {"reverse": False}

        
    def filterSingleIds(self, action, ids, new_ids=None, data=None, constraints=None, mute_tracking=None):
        # tests ids if they are not new (new_ids = None means everything is new), so old one are kept regardless
        track = {"do": "filter", "src": list(ids), "action": action.getWhat()}
        if new_ids is not None:
            track["src"] = new_ids
        fids = [i for i in ids if (new_ids is not None and i not in new_ids) or not self.computeSatisfy(action, i, constraints=constraints, mute_tracking=True)]
        track["trg"] = list(fids)
        self.addTrack(track, mute_tracking)                
        return fids
    actions_methods["filterSingle"] = filterSingleIds

    class ActionFilterUnary(ActionUnary):
        def_args = {"reverse": False}
    
    ### corresponding Action subclass
    class ActionFilterSingle(ActionFilterUnary):
        what = "filterSingle"
        
           
    def filterLastIds(self, action, ids, new_ids=None, data=None, constraints=None, mute_tracking=None):
        # check last element against all previous in list
        max_filter = action.getArg("max", constraints)
        if action.getNbBlocks() > 0 and (len(ids) >= max_filter+2):
            filter_out = []
            ts_out = []
            posC = len(ids)-1
            posP = 0
            while len(filter_out) <= max_filter and posP < posC:
                v, ts = action.computeSatisfyAndTrack(self.getItem(ids[posC]), other=self.getItem(ids[posP]), constraints=constraints)
                if v:
                    filter_out.append(ids[posP])
                    ts_out.append(ts)
                posP += 1
            if len(filter_out) > max_filter:
                track = {"do": "filter", "trg": [ids[posC]], "src": filter_out, "out": True, "rationales": ts, "action": action.getWhat()}
                self.addTrack(track, mute_tracking)
                return True
        return False
    actions_methods["filterLast"] = filterLastIds

    class ActionFilterBinary(ActionBinary):
        def_args = {"reverse": False, "max": 0}
    
    ### corresponding Action subclass
    class ActionFilterLast(ActionFilterBinary):
        what = "filterLast"

    def filterToFirstIds(self, action, ids, new_ids=None, data=None, constraints=None, mute_tracking=None):
        # filter list checking elements against first in list
        if action.getNbBlocks() > 0 and len(ids) >= 2:
            # start from where the condition might be broken
            posC = 0
            posP = 1
            while posP < len(ids):
                v, ts = action.computeSatisfyAndTrack(self.getItem(ids[posP]), other=self.getItem(ids[posC]), constraints=constraints)
                if v:
                    track = {"do": "filter", "trg": [ids[posP]], "src": [ids[posC]], "out": True, "rationales": ts, "action": action.getWhat()}
                    self.addTrack(track, mute_tracking)
                    ids.pop(posP)
                else:
                    posP += 1
        return ids
    actions_methods["filterToFirst"] = filterToFirstIds

    ### corresponding Action subclass
    class ActionFilterToFirst(ActionFilterBinary):
        what = "filterToFirst"
    
    def filterPairsIds(self, action, ids, new_ids=None, data=None, constraints=None, mute_tracking=None):
        # filter list checking elements against previous in list
        before_ids = list(ids)
        max_filter = action.getArg("max", constraints)
        if action.getNbBlocks() > 0 and (len(ids) >= max_filter+2):
            # start from where the condition might be broken
            posC = max_filter+1
            kept_new = (new_ids is None or ids[0] in new_ids)
            while posC < len(ids):
                if (kept_new or ids[posC] in new_ids):
                    if self.filterLastIds(action, ids[:posC+1], constraints=constraints, mute_tracking=mute_tracking):
                        ids.pop(posC)
                    else:
                        kept_new = True
                        posC += 1
                else:
                    posC += 1
        return ids
    actions_methods["filterPairs"] = filterPairsIds

    ### corresponding Action subclass
    class ActionFilterPairs(ActionFilterBinary):
        what = "filterPairs"
    
    def cutIds(self, action, ids, new_ids=None, data=None, constraints=None, mute_tracking=None):
        max_cut = action.getArg("max", constraints)
        cut = abs(max_cut)
        cut_org = abs(max_cut)
        direction_cut = action.getArg("direction", constraints)
        if cut is not None and cut < len(ids):
            ts_out = []
            if action.getNbBlocks() > 0 and direction_cut != 0:
                stop = cut <= 0 or cut >= len(ids)
                while not stop:
                    v, ts = action.computeSatisfyAndTrack(self.getItem(ids[cut-1]), other=self.getItem(ids[cut]), constraints=constraints)
                    ts_out.append(ts)
                    if v:
                        cut += direction_cut
                        stop = cut <= 0 or cut >= len(ids)
                    else:
                        stop = True
            if max_cut < 0:
                track = {"do": "cut", "trg": ids[:cut], "src": ids[cut:cut_org+1], "out": (cut, cut_org), "rationales": ts_out, "action": action.getWhat()}
                del ids[:cut]
            else:
                track = {"do": "cut", "trg": ids[cut:], "src": ids[cut_org-1:cut], "out": (cut_org, cut), "rationales": ts_out, "action": action.getWhat()}
                del ids[cut:]
            self.addTrack(track, mute_tracking)
        return ids
    actions_methods["cut"] = cutIds
    
    ### corresponding Action subclass
    class ActionCut(ActionFilterBinary):
        what = "cut"
        def_args = {"reverse": False, "max": 0, "direction": 0}

    def doAction(self, action, ids=None, complement=False, new_ids=None, data=None, constraints=None, mute_tracking=None):
        if len(ids) > 0:
            if action.getWhat() in self.actions_methods:
                return self.actions_methods[action.getWhat()](self, action, ids, new_ids=new_ids, data=data, constraints=constraints, mute_tracking=mute_tracking)
            else:
                raise Exception('Oups action method does not exist (%s)!' % action.getWhat())
        return []

    def doActions(self, actions, ids=None, complement=False, new_ids=None, data=None, constraints=None, mute_tracking=None):
        al_name = None
        ids = self.prepareIids(ids)
        nb_src_ids = len(ids)
        if new_ids is not None:
            new_ids = self.prepareIids(new_ids)
            ids.extend(new_ids)
        before_ids = list(ids)
        
        if type(actions) is str and constraints is not None:
            al_name = actions
            actions = constraints.getActionList(al_name)
            if len(actions) > 0:
                track = {"do": "START", "action": al_name, "src": before_ids[:nb_src_ids], "trg": before_ids[nb_src_ids:]}
                self.addTrack(track, mute_tracking)
        
        for action in actions:
            ids = self.doAction(action, ids, complement, new_ids, data, constraints, mute_tracking=mute_tracking)
        if complement:
            ids = (set(before_ids)-set(ids))

        if al_name is not None and len(actions) > 0:
            track = {"do": "END", "action": al_name, "trg": list(ids)}
            self.addTrack(track, mute_tracking)            
        return ids

    def selected(self, actions=[], ids=None, complement=False, new_ids=None, trg_lid=None, new_only=False, data=None, constraints=None, mute_tracking=None):
        ids_org = ids
        # applies a sequence of action to select a sequence of elements from the batch
        ids = self.doActions(actions, ids, complement, new_ids, data, constraints, mute_tracking=mute_tracking)
        if new_only and new_ids is not None:
            ids = [i for i in ids if i in new_ids]
        if trg_lid is not None:
            if self.hasLid(trg_lid) and trg_lid == ids_org:
                self.clearList(trg_lid)
            for iid in ids:
                self.addItem(iid, trg_lid)
        return ids
    
    def selectedItems(self, actions=[], ids=None, complement=False, new_ids=None, data=None, constraints=None, mute_tracking=None):
        return [self.getItem(iid) for iid in self.selected(actions, ids, complement, new_ids, data=data, constraints=constraints, mute_tracking=mute_tracking)]


class StoredRCollection(StoredContentCollection, BatchCollection):
    name_class = "StoredR content"

    # def __init__(self):
    #     StoredContentCollection.__init__(self)


def unitTest():
    items = [Item() for i in range(4)]
    content = dict([(item.getUid(), item) for item in items])
    x = StoredContainer()
    x.extend(content.keys())
    print("List:", x)
    print("List sort info:", x.getSortInfo())
    print("List elements:")
    for iid in x:
        print(iid, content[iid])

    print("List set sort k")
    x.setSort("k", direct=True)
    print("List sort info:", x.getSortInfo())
    print("List elements:")
    for iid in x:
        print(iid, content[iid])
    print("Locate first element:", x.getPosIid(items[0].getUid()))
    print("List shortStr:", x.getShortStr())


if __name__ == '__main__':
    unitTest()
