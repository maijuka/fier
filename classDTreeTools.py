import numpy
from sklearn import tree

try:
    from classQuery import  *
except ModuleNotFoundError:
    from .classQuery import  *

import pdb

class DTree:
    class_name = "DTree"

    @classmethod
    def isNumericalFeat(tcl, feat):
        try:
            v = int(feat)
            return True
        except (ValueError, TypeError):
            return False
    @classmethod
    def isSpecialFeat(tcl, feat):
        return any([isinstance(feat, c) for c in [Literal, Term, Query]])
    @classmethod
    def isEmptyFeat(tcl, feat):
        return feat is None
    @classmethod
    def testValue(tcl, value):
        return value != 0

    
    def isEmpty(self):
        return self.getNodeCount() == 0
    def hasFeature(self):
        return any([v is not None for v in self.getFeature()])

    def getNodeIds(self):
        return range(self.getNodeCount())
    def __str__(self):
        # tstr = "%s %s (%d)" % (self.class_name, self.getTid(), self.getNodeCount())
        tstr = "%s (%d)" % (self.class_name, self.getNodeCount())
        for i in self.getNodeIds():
            if self.isLeaf(i):
                tstr += "\n\t#%d: %s\t[%d+%d=%d]" % (i, self.getClass(i), self.getLabelCount(i, 0), self.getLabelCount(i, 1), self.getLabelCount(i, 0) + self.getLabelCount(i, 1))
            else:
                tstr += "\n\t#%d %s=%s <%s,%s>\t[%d+%d=%d]" % (i, self.getFeature(i), self.getThreshold(i), self.getChildrenLeft(i), self.getChildrenRight(i), self.getLabelCount(i, 0), self.getLabelCount(i, 1), self.getLabelCount(i, 0)+ self.getLabelCount(i, 1))
        return tstr
    # def setTid(self, tid=None):
    #     self.tid = tid
    # def getTid(self):
    #     return self.tid
    def getTree(self):
        return self
    def isValidId(self, node_id):
        return node_id > -1 and node_id < self.getNodeCount()
    def isLeaf(self, node_id):
        return self.isValidId(node_id) and self.getTree().tclass[node_id] != -1

    def getNodeCount(self):
        return self.getTree().node_count
    def getChildrenLeft(self, node_id=None):
        if node_id is None:
            return [self.getChildrenLeft(i) for i in self.getNodeIds()]
        return self.getTree().children_left[node_id]
    def getChildrenRight(self, node_id=None):
        if node_id is None:
            return [self.getChildrenRight(i) for i in self.getNodeIds()]
        return self.getTree().children_right[node_id]
    def findValue(self, k, values):
        if k in values:
            return values.index(k)
        return None
    def getParent(self, node_id=None):
        if node_id is None:
            return [self.getParent(i) for i in self.getNodeIds()]
        idx = self.findValue(node_id, self.getTree().children_left)
        if idx is None:
            idx = self.findValue(node_id, self.getTree().children_right)
            if idx is None:
                return None
            return (1, idx)
        return (0, idx)
    def getDepth(self, node_id=None):
        if node_id is None:
            return max([self.getDepth(i)  for i in self.getNodeIds() if self.isLeaf(i)])
        pidx = self.getParent(node_id)
        if pidx is None:
            return 0
        else:
            return self.getDepth(pidx[1]) + 1

    def getFeature(self, node_id=None):
        if node_id is None:
            return [self.getFeature(i) for i in self.getNodeIds()]
        return self.getTree().feature[node_id]
    def getFeatures(self):
        return set([v for v in self.getTree().feature if v > -1])

    def getThreshold(self, node_id=None):
        if node_id is None:
            return [self.getThreshold(i) for i in self.getNodeIds()]
        return self.getTree().threshold[node_id]
    def getClass(self, node_id=None): ### return the class decision for leaf node
        if node_id is None:
            return [self.getClass(i) for i in self.getNodeIds()]
        return self.getTree().tclass[node_id]
    def getLabelCount(self, node_id=None, label=0):
        if node_id is None:
            return [self.getLabelCount(i, label) for i in self.getNodeIds()]
        if self.getTree().counts[node_id] is None or label not in [0,1]:
            return -1
        return self.getTree().counts[node_id][label]
    def getMajLabel(self, node_id=None): ### return majority training label among instances in node
        if node_id is None:
            return [self.getMajLabel(i) for i in self.getNodeIds()]
        if self.getLabelCount(node_id, 0) == self.getLabelCount(node_id, 1):
            return -1
        return int(self.getLabelCount(node_id, 0) < self.getLabelCount(node_id, 1))
    def getNbInstances(self):
        return self.nb_instances
    def getSuppSet(self, node_id=None):
        if node_id is None:
            return self.supp_sets
        return self.supp_sets[node_id]
    def getSuppVect(self, node_id=0):
        if node_id is None:
            return [self.getSuppVect(i) for i in self.getNodeIds()]
        if self.isValidId(node_id):            
            X = numpy.zeros(self.getNbInstances(), dtype=int)
            X[list(self.getSuppSet(node_id))] = 1
            return X
        return -numpy.ones(self.getNbInstances(), dtype=int)

    def resetInstances(self, n=0):
        self.nb_instances = n
        self.counts = [None for i in self.getNodeIds()]
        self.supp_sets = [None for i in self.getNodeIds()]
    def setCounts(self, node_id, counts):
        self.counts[node_id] = counts
    def setSuppSet(self, node_id, sset):
        self.supp_sets[node_id] = sset
    def applyMaskSuppSets(self, mask):
        map_mask = {}
        if type(mask) is list:
            map_mask = dict(enumerate(mask))            
        for node_id in self.getNodeIds():
            self.supp_sets[node_id] = set([map_mask.get(i, i) for i in self.supp_sets[node_id]])

        
    def computeSupp(self, data, node_id=None, ids=None, supp_vect=None, target=None):
        if node_id is None:
            node_id = 0
        if ids is None:
            ids = list(range(data.shape[0]))
        if supp_vect is None:
            supp_vect = numpy.zeros(data.shape[0], dtype=bool)
        if not self.isValidId(node_id):
            return supp_vect
        if target is not None:
            np = numpy.sum(target[ids])
            nn = len(ids) - np
            self.setCounts(node_id, [nn, np])
            self.setSuppSet(node_id, set(ids))
            
        if self.isLeaf(node_id):
            supp_vect[ids] = self.getClass(node_id)
        else:
            feat = self.getFeature(node_id)
            thres = self.getThreshold(node_id)
            lc = self.getChildrenLeft(node_id)
            rc = self.getChildrenRight(node_id)
            if self.isEmptyFeat(feat) or self.isSpecialFeat(feat):
                lids = list(set(range(data.shape[0])).difference(ids))
                rids = list(ids)
            else:
                if thres is None: pdb.set_trace()
                lids = [ids[x] for x in numpy.where(data[ids, feat] <= thres)[0]]
                rids = [ids[x] for x in numpy.where(data[ids, feat] > thres)[0]]
            # print("---", node_id, feat, len(lids), len(rids))
            self.computeSupp(data, node_id=lc, ids=lids, supp_vect=supp_vect, target=target)
            self.computeSupp(data, node_id=rc, ids=rids, supp_vect=supp_vect, target=target)        
        return supp_vect
    def getSupportVect(self, data, ids=None):
        return self.computeSupp(data, ids=ids)

    def collectLeaves(self, only_class=-1):
        return [i for i in self.getNodeIds() if ((only_class == -1 and self.getClass(i) != -1) or (only_class != -1 and self.getClass(i) == only_class))]
    def collectBranches(self, node_id=None, suffix=None, branches=None, only_class=-1, only_branch=-1):
        if branches is None:
            branches = []
        if suffix is None:
            suffix = []
        if node_id is None:
            if self.getNodeCount() > 0:
                node_id = 0
            else:
                return []
        if self.isLeaf(node_id):
            if only_class == -1 or self.getClass(node_id) == only_class:
                branches.append([(None, node_id)]+suffix)
        else:
            if only_branch == -1 or only_branch == 0: 
                self.collectBranches(node_id=self.getChildrenLeft(node_id),
                                     suffix=[(0, node_id)]+suffix, branches=branches, only_class=only_class, only_branch=only_branch)
            if only_branch == -1 or only_branch == 1: 
                self.collectBranches(node_id=self.getChildrenRight(node_id),
                                     suffix=[(1, node_id)]+suffix, branches=branches, only_class=only_class, only_branch=only_branch)
        return branches

    def exportNodes(self):
        return (list(self.getFeature()), list(self.getThreshold()), list(self.getClass()),
                list(self.getChildrenLeft()), list(self.getChildrenRight()))

    def copy(self, in_data=None, in_target=None):
        return None
    
    def pruneDeadBranches(self, in_data=None, in_target=None):
        leaf_pairs_matched = [(pc, lc, rc) for (pc, lc, rc) in zip(self.getNodeIds(), self.getChildrenLeft(), self.getChildrenRight()) if self.isLeaf(lc) and self.isLeaf(rc) and self.getClass(lc) == self.getClass(rc)]
        if len(leaf_pairs_matched) > 0:
            feats, thres, clss, chls, chrs = self.exportNodes()
            pop_nodes = []
            for (pc, lc, rc) in sorted(leaf_pairs_matched):
                ppc = self.getParent(pc)
                if ppc is None: ### no parent, root -> empty tree
                    return DTree([], in_data, in_target)
                else:
                    pop_nodes.extend((lc, rc))
                    feats[pc], thres[pc], chls[pc], chrs[pc] = (-1, -1, -1, -1)
                    clss[pc] = self.getClass(rc)
            if len(pop_nodes) > 0:
                keep_pos = [i for i in self.getNodeIds() if i not in pop_nodes]
                map_nds = dict([(v,k) for (k,v) in enumerate(keep_pos)])
                feats = [feats[i] for i in keep_pos]
                thres = [thres[i] for i in keep_pos]
                clss = [clss[i] for i in keep_pos]
                chls = [map_nds.get(chls[i],-1) for i in keep_pos]
                chrs = [map_nds.get(chrs[i],-1) for i in keep_pos]
                xx = DTree((feats, thres, clss, chls, chrs), in_data, in_target)
                return xx.pruneDeadBranches(in_data, in_target)
        return self
        
    
    def __init__(self, nodes=[], in_data=None, in_target=None):
        if type(nodes) is dict:
            self.initBasic(supp_pos=nodes.get("supp_pos"), feat=nodes.get("feat"), n=nodes.get("n", 0))
        elif len(nodes) == 5  and not type(nodes[0]) is dict:
            self.initCopy(nodes)
        else:
            self.initFromNodes(nodes)
        if not self.isEmpty() and in_data is not None and in_target is not None:
            X = self.computeSupp(data=in_data, target=in_target)
                        
    def initCopy(self, nodes):
        if len(nodes) == 5 and not type(nodes[0]) is dict:
            self.node_count = len(nodes[0])
            self.resetInstances(n=0)
            self.feature = nodes[0]
            self.threshold = nodes[1]
            self.tclass = nodes[2]
            self.children_left = nodes[3]
            self.children_right = nodes[4]            
        else:
            self.initFromNodes(nodes)
            
    def initFromNodes(self, nodes=[]):
        self.node_count = len(nodes)
        self.resetInstances(n=0)       
        self.feature = [-1 for i in self.getNodeIds()]
        self.threshold = [None for i in self.getNodeIds()]
        self.tclass = [-1 for i in self.getNodeIds()]
        self.children_left = [-1 for i in self.getNodeIds()]
        self.children_right = [-1 for i in self.getNodeIds()]
        for ni, node in enumerate(nodes):
            self.tclass[ni] = node.get("class", -1)
            if "attribute" in node:
                # if self.feature[node["parent"]] != -1 and self.feature[node["parent"]] != node["attribute"]: pdb.set_trace()            
                self.feature[node["parent"]] = node["attribute"]
                if self.testValue(node["value"]):
                    self.threshold[node["parent"]] = node["value"] - 0.5 #
                    self.children_right[node["parent"]] = ni
                    # self.counts[node["parent"]][0] = node["size"]
                else:
                    self.threshold[node["parent"]] = node["value"] + 0.5 #
                    self.children_left[node["parent"]] = ni
                    # self.counts[node["parent"]][1] = node["size"]

    def initBasic(self, supp_pos=None, feat=None, n=0):
        if supp_pos is not None or n > 0 or feat is not None:
            self.node_count = 3
            self.feature = [feat, -1, -1]
            self.tclass = [-1, 0, 1]
            self.threshold = [None, None, None]
            self.children_left = [1, -1, -1]
            self.children_right = [2, -1, -1]
            if supp_pos is None:
                self.supp_sets = [set(range(n)), set(), set(range(n))]
            elif type(supp_pos) is set:
                self.supp_sets = [set(range(n)), set(range(n)).difference(supp_pos), set(supp_pos)]
            else:
                self.supp_sets = [None, set(numpy.where(supp_pos==0)[0]), set(numpy.where(supp_pos==1)[0])]
                self.supp_sets[0] = self.supp_sets[1].union(self.supp_sets[2])
            nn, np = len(self.supp_sets[1]), len(self.supp_sets[2])
            self.counts = [[nn, np], [nn, 0], [0, np]]
            self.nb_instances = nn+np
        else:
            self.initFromNodes()

    def copy(self, in_data=None, in_target=None):
        return DTree(self.exportNodes(), in_data, in_target)

    
class DTreeToolsSKL:
    @classmethod
    def getNodes(tcl, skl_dtc):
        dtree = skl_dtc.tree_
        feats, thres = list(dtree.feature), list(dtree.threshold)
        chls, chrs = list(dtree.children_left), list(dtree.children_right)
        clss = []
        for node_id in range(len(dtree.children_left)):
            if dtree.children_left[node_id] == tree._tree.TREE_LEAF: ### node is a leaf
                clss.append(int(dtree.value[node_id][0][0] < dtree.value[node_id][0][1]))
            else:
                clss.append(-1)
        return feats, thres, clss, chls, chrs
    
    @classmethod
    def fitTree(tcl, in_data, in_target, in_depth, in_min_bucket, split_criterion="gini", random_state=0, logger=None):
        if logger is not None:
            logger.printL(10, "FITTING TREE\tnb input data size = %s, nb positive labels = %d, depth = %d, min node size = %d" % (in_data.shape, sum(in_target), in_depth, in_min_bucket), "X")
        if len(set(in_target)) == 1:
            if logger is not None:
                logger.printL(10, "All target labels identical, stopping.", "X")
            return DTree([], in_data, in_target)
        skl_dtc = tree.DecisionTreeClassifier(criterion=split_criterion, max_depth=in_depth, min_samples_leaf=in_min_bucket, random_state=0).fit(in_data, in_target)
        if skl_dtc.tree_.value.shape[2] == 1: ### did not split
            if logger is not None:
                logger.printL(10, "Did not split", "X")
            skl_dtc.tree_.value = numpy.array([[[numpy.sum(in_target==0), numpy.sum(in_target==1)]]])        
        dt = DTree(tcl.getNodes(skl_dtc), in_data, in_target)
        dtt = dt.pruneDeadBranches(in_data, in_target)
        if logger is not None:
            logger.printL(10, "Output tree:\t %s" % dtt, "X")
        return dtt
