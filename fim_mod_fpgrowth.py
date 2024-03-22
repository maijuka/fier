import sys
import datetime
import numpy

import pdb


def calc_jacc(main_pattern, group_patterns):
    return main_pattern[1]/sum([g[1]-main_pattern[1] for g in group_patterns]+[main_pattern[1]])


class PatternStore(dict):

    def add(self, main_pattern, group_patterns=None):
        if group_patterns is None:
            self[len(self)] = main_pattern
        else:
            self[len(self)] = (main_pattern, group_patterns)


class FPNode(object):

    # index -1 is used exclusively for tree roots
    NODE_ID = -1

    @classmethod
    def get_next_id(tcl):
        tcl.NODE_ID += 1
        return tcl.NODE_ID

    def __init__(self, value, count=1):
        self.id = self.get_next_id()
        self.value = value
        self.count = count

    def __str__(self):
        return "%s (%d) [n%s]" % (self.value, self.count, self.id)

    def get_id(self):
        return self.id

    def get_value(self):
        return self.value

    def get_count(self):
        return self.count

    def add_count(self, count=1):
        self.count += count


class FPTree(object):

    def __str__(self):
        return "FPTree (%s:%d) %s" % (self.get_suffix(), self.get_suffix_count(), self.get_mining_order())

    def str_long(self):
        xs = "FPTree (%s:%d) %s" % (self.get_suffix(), self.get_suffix_count(), self.get_mining_order())
        return xs + self.recurse_str()

    def recurse_str(self, node_id=-1, level=0):
        xs = ""
        if node_id in self.nodes:
            xs += "\n" + level*"  " + str(self.nodes[node_id]) + "*"*(len(self.children.get(node_id, [])) == 0)
        for value, child_id in self.children.get(node_id, {}).items():
            xs += self.recurse_str(child_id, level+1)
        return xs

    # FPTree are meant only to be read, except for the mining order, hence
    def copy(self, item_counts=None, mining_order=None, local_support_threshold=None):
        if item_counts is None:
            item_counts = self.item_counts
        if mining_order is None:
            mining_order = self.mining_order
        if local_support_threshold is None:
            local_support_threshold = self.local_support_threshold
        data = {"nodes": self.nodes,
                "children": self.children,
                "parent": self.parent,
                "links": self.links,
                "suffix_count": self.suffix_count,
                "local_support_threshold": local_support_threshold}
        return FPTree(data, item_counts, mining_order, suffix=self.suffix)

    def __init__(self, data, item_counts, mining_order=[], transaction_counts=None, suffix=[]):

        self.local_support_threshold = None
        self.suffix = suffix
        self.item_counts = item_counts
        self.mining_order = mining_order
        if type(data) is dict:
            self.local_support_threshold = data["local_support_threshold"]
            self.nodes = data["nodes"]
            self.children = data["children"]
            self.parent = data["parent"]
            self.links = data["links"]
            self.suffix_count = data["suffix_count"]

        else:
            if transaction_counts is None:
                order_map = dict([(v, i) for (i, v) in enumerate(mining_order)])
            else:
                order_map = None

            self.nodes = {}
            self.children = {}
            self.parent = {}
            self.links = {}

            current_id = -1  # root
            self.children[current_id] = {}
            if type(data) is list:
                self.suffix_count = self.insert_transactions(data, item_counts, mining_order, transaction_counts, order_map, current_id)
            else:
                self.suffix_count = self.insert_matrix(data, item_counts, mining_order, transaction_counts, order_map, current_id)

    def insert_matrix(self, mat, item_counts, mining_order, transaction_counts=None, order_map=None, current_id=-1):
        for ti in range(mat.shape[0]):
            sorted_items = list(numpy.where(mat[ti, :])[0])
            if order_map is not None:
                sorted_items.sort(key=lambda x: order_map[x], reverse=True)
            if transaction_counts is None:
                count = 1
            else:
                count = transaction_counts[ti]
            if len(sorted_items) > 0:
                self.insert_transaction(sorted_items, current_id, count)
        return mat.shape[0]

    def insert_transactions(self, transactions, item_counts, mining_order, transaction_counts=None, order_map=None, current_id=-1):
        suffix_count = 0
        for ti, transaction in enumerate(transactions):
            sorted_items = [x for x in transaction if x in item_counts]
            if order_map is not None:
                sorted_items.sort(key=lambda x: order_map[x], reverse=True)
            if transaction_counts is None:
                count = 1
            else:
                count = transaction_counts[ti]
            suffix_count += count
            if len(sorted_items) > 0:
                self.insert_transaction(sorted_items, current_id, count)
        return suffix_count

    def insert_transaction(self, items, current_id, count=1):
        child_id = self.get_child(current_id, items[0])
        if child_id is not None:
            self.nodes[child_id].add_count(count)
        else:
            child_id = self.add_child(current_id, items[0], count)

        if len(items[1:]) > 0:
            self.insert_transaction(items[1:], child_id, count)

    def get_single_leaf_id(self, node_id=-1):
        """
        If there is a single path in the tree,
        return the id of the leaf at the end, which is -1 (the root) if the tree does not contain any nodes, else return None.
        """
        if node_id == -1 and not self.has_frequent_item():
            return node_id
        nb_children = len(self.children[node_id])
        if nb_children > 1:
            return None
        elif nb_children == 0:
            return node_id
        else:
            return self.get_single_leaf_id(list(self.children[node_id].values())[0])

    def get_node_count(self, node_id):
        return self.nodes[node_id].get_count()

    def get_node_value(self, node_id):
        return self.nodes[node_id].get_value()

    def get_path_to_root(self, node_id=-1):
        path = []
        while node_id != -1:  # while not reached the root
            path.append((node_id, self.get_node_value(node_id), self.get_node_count(node_id)))
            node_id = self.parent[node_id]
        return path

    def get_itemset_to_root(self, node_id=-1):
        items = []
        while node_id != -1:  # while not reached the root
            items.append(self.get_node_value(node_id))
            node_id = self.parent[node_id]
        return items

    def get_all_leaves_id(self):
        return [k for k in self.nodes.keys() if len(self.children.get(k, {})) == 0]

    def get_all_trans_counts(self):
        # also intermediate nodes with counts higher than the sum of their children's counts
        trans = []
        counts = []
        for (node_id, node) in self.nodes.items():
            count = node.get_count() - sum([self.get_node_count(ch) for (v, ch) in self.children.get(node_id, {}).items()])
            if count > 0:
                trn = []
                while node_id != -1:  # while not reached the root
                    trn.append(self.get_node_value(node_id))
                    node_id = self.parent[node_id]
                trans.append(trn[::-1])
                counts.append(count)
        return trans, counts

    def has_frequent_item(self, item=None):
        if item is None:
            return len(self.item_counts) > 0
        return item in self.item_counts

    def get_item_counts(self):
        return self.item_counts

    def get_mining_order(self):
        return self.mining_order

    def get_local_support_threshold(self):
        return self.local_support_threshold

    def set_local_support_threshold(self, local_support_threshold, inplace=True):
        if self.local_support_threshold != local_support_threshold:
            item_counts = dict([(k, v) for (k, v) in self.item_counts.items() if v >= local_support_threshold])
            mining_order = [k for k in self.mining_order if k in item_counts]
            if inplace:
                self.local_support_threshold = local_support_threshold
                self.mining_order = mining_order
                self.item_counts = item_counts
            else:
                return self.copy(item_counts=item_counts, mining_order=mining_order, local_support_threshold=local_support_threshold)
        return self

    def filter_items(self, selected_items=None, remove=True, inplace=True):
        if selected_items is not None:
            if remove:
                mining_order = [k for k in self.mining_order if k not in selected_items]
                item_counts = dict([(k, v) for (k, v) in self.item_counts.items() if k not in selected_items])
            else:
                mining_order = [k for k in self.mining_order if k in selected_items]
                item_counts = dict([(k, v) for (k, v) in self.item_counts.items() if k in selected_items])

            if len(item_counts) != len(self.item_counts) or len(mining_order) != len(self.mining_order):
                if inplace:
                    self.mining_order = mining_order
                    self.item_counts = item_counts
                else:
                    return self.copy(item_counts=item_counts, mining_order=mining_order)
        return self

    def get_suffix_count(self):
        return self.suffix_count

    def get_suffix_len(self):
        return len(self.suffix)

    def get_suffix(self):
        return self.suffix

    def has_frequent_suffix(self, support_threshold=None):
        if self.get_local_support_threshold() is not None and self.get_local_support_threshold() > support_threshold:
            return len(self.suffix) > 0 and (support_threshold is None or self.suffix_count >= self.get_local_support_threshold())
        return len(self.suffix) > 0 and (support_threshold is None or self.suffix_count >= support_threshold)

    def get_suffixed_itemset(self, items=tuple()):
        # return tuple(sorted(self.suffix+list(items)))
        return items+tuple(self.suffix)

    def get_links(self, item):
        return self.links.get(item, [])

    def get_child(self, node_id, value):
        if value in self.children[node_id]:
            return self.children[node_id][value]
        return None

    def add_child(self, node_id, value, count=1):
        child = FPNode(value, count)
        child_id = child.get_id()
        self.nodes[child_id] = child
        self.children[child_id] = {}
        self.children[node_id][value] = child_id
        self.parent[child_id] = node_id
        if value not in self.links:
            self.links[value] = []
        self.links[value].append(child_id)
        return child_id


class FPGrowthMiner:

    def key_freq(x, map_itog, item_counts):
        return item_counts.get(x, 0)

    def key_grpfreq(x, map_itog, item_counts):
        return (map_itog.get(x, [-1])[0], item_counts.get(x, 0))

    sorting_def = "freq-asc"
    sorting_dets = {"freq-asc": {"fnct": key_freq, "reverse": False},
                    "freq-desc": {"fnct": key_freq, "reverse": True},
                    "grp-freq-asc": {"fnct": key_grpfreq, "reverse": False},
                    "grp-freq-desc": {"fnct": key_grpfreq, "reverse": True}
                    }

    @classmethod
    def get_sorting_choices(tcl):
        return tcl.sorting_dets.keys()

    @classmethod
    def get_sorting_default(tcl):
        return tcl.sorting_def

    @classmethod
    def sort_items(tcl, items, map_itog, item_counts, sorting_name):
        # what's the best order?
        if sorting_name in tcl.sorting_dets:
            fnct = tcl.sorting_dets[sorting_name]["fnct"]
            rvr = tcl.sorting_dets[sorting_name]["reverse"]
        else:
            fnct = tcl.sorting_dets[tcl.sorting_def]["fnct"]
            rvr = tcl.sorting_dets[tcl.sorting_def]["reverse"]
        return sorted(items, key=lambda x: fnct(x, map_itog, item_counts), reverse=rvr)

    bthres_names = ["min_supp", "max_supp", "min_len", "max_len"]
    bthres_dets = {"min_supp": {"shrt": "s", "help": "minimum support threshold", "default": 0.1, "total": "nb_rows"},
                   "max_supp": {"shrt": "t", "help": "maximum support threshold", "total": "nb_rows"},
                   "min_len": {"shrt": "n", "help": "minimum itemset length", "total": "nb_cols"},
                   "max_len": {"shrt": "m", "help": "maximum itemset length", "total": "nb_cols"}
                   }

    gthres_names = []
    gthres_dets = {}

    @classmethod
    def get_bthres_names(tcl):
        return tcl.bthres_names

    @classmethod
    def get_bthres_details(tcl, thres=None):
        if thres is None:
            return [(t, tcl.bthres_dets[t]) for t in tcl.get_bthres_names()]
        return tcl.bthres_dets.get(t)

    @classmethod
    def get_bthres_totn(tcl, thres, group=False):
        if thres in tcl.bthres_dets and "total" in tcl.bthres_dets[thres]:
            return ("group_"*group)+tcl.bthres_dets[thres]["total"]

    @classmethod
    def get_gthres_names(tcl):
        return tcl.gthres_names

    @classmethod
    def get_gthres_details(tcl, thres=None):
        if thres is None:
            return [(t, tcl.gthres_dets[t]) for t in tcl.get_gthres_names()]
        return tcl.gthres_dets.get(t)

    @classmethod
    def get_gthres_fnct(tcl, thres):
        if thres in tcl.gthres_dets and "fnct" in tcl.gthres_dets[thres]:
            return tcl.gthres_dets[thres]["fnct"]

    def get_nb_cols(self, gi=None):
        if gi is None:
            return self.data_shape["nb_cols"]
        else:
            return self.data_shape["group_nb_cols"][gi]

    def get_nb_rows(self, gi=None):
        return self.data_shape["nb_rows"]

    @classmethod
    def prepare_thresholds(tcl, params, nb_groups, data_shape=None, split=False):
        prep_thres = {}
        # scale thresholds depending on data size
        for p in tcl.get_bthres_names():
            if p in params:
                if data_shape is not None and tcl.get_bthres_totn(p) is not None and 0 < params[p] < 1:
                    prep_thres[p] = int(params[p]*data_shape[tcl.get_bthres_totn(p)])
                else:
                    prep_thres[p] = int(params[p])
            gp = "group_"+p
            if split and gp in params:
                vs = []
                if type(params[gp]) is list or type(params[gp]) is tuple:
                    ovs = params[gp]
                else:
                    ovs = params[gp].split(",")
                if len(ovs) == 1 or len(ovs) == nb_groups:
                    for vi, v in enumerate(ovs):
                        try:
                            v = float(v)
                            if data_shape is not None and tcl.get_bthres_totn(p) is not None and 0 < v < 1:
                                if tcl.get_bthres_totn(p, group=True) in data_shape:
                                    vs.append(int(v*data_shape[tcl.get_bthres_totn(p, group=True)][vi]))
                                else:
                                    vs.append(int(v*data_shape[tcl.get_bthres_totn(p)]))
                            else:
                                vs.append(int(v))
                        except (ValueError, TypeError):
                            vs.append(None)
                    if len(vs) == nb_groups:
                        prep_thres[gp] = vs
                    elif len(vs) == 1:
                        prep_thres[gp] = vs*nb_groups
        # Global thresholds, no groups, no scaling
        for p in tcl.get_gthres_names():
            if p in params:
                prep_thres[p] = params[p]
        return prep_thres

    def prepare_exclusions(self, params, map_itog=False, items_groups=None):
        param_x = params.get("item_exclusions")
        x_within_group = False
        if param_x is None:
            return {}, x_within_group
        elif type(param_x) is dict:  # ready dict given as param
            item_exclusions = param_x
        else:  # assume exclusions given as a string "i:j,k;l:j..."
            item_exclusions = {}
            for p in param_x.split(";"):
                if ":" in p:
                    itsA, itsB = p.split(":")
                    itsB = self.prep_items(itsB.split(","))
                    for itA in self.prep_items(itsA.split(",")):
                        if itA not in item_exclusions:
                            item_exclusions[itA] = set(itsB)
                        else:
                            item_exclusions[itA].update(itsB)
                else:
                    its = self.prep_items(p.split(","))
                    for ii, it in enumerate(its):
                        if it not in item_exclusions:
                            item_exclusions[it] = set(its)
                        else:
                            item_exclusions[it].update(its)
        # check whether exclusions are within groups, if groups are not disjoint no need to check, they are handled as cross-group
        if self.groups_disjoint:
            for it, X in item_exclusions.items():
                gis = self.map_itog.get(it, [])
                if len(gis) > 1:
                    return item_exclusions, True
                if not X.issubset(self.items_groups[gis[0]]):
                    return item_exclusions, True
        return item_exclusions, x_within_group

    def __init__(self, params, items_groups, data, split=False, log=None):
        self.log = log
        self.items_groups = items_groups
        self.map_itog, self.groups_disjoint = self.map_items_to_groups(items_groups)
        self.num_items = type(list(self.map_itog.keys())[0]) is int
        self.item_exclusions, self.x_within_group = self.prepare_exclusions(params)
        self.data_shape = {"nb_rows": len(data), "group_nb_cols": [len(g) for g in items_groups]}
        self.data_shape["nb_cols"] = sum(self.data_shape["group_nb_cols"])
        self.thresholds = self.prepare_thresholds(params, self.get_nb_groups(), self.data_shape, split)
        self.parameters = dict([(k, v) for (k, v) in params.items() if k not in self.thresholds])

    def are_items_numerical(self):
        return self.num_items

    def to_log(self, level, message):
        if type(self.log) is int:
            if level < self.log:
                print(message)
        elif self.log is not None:
            self.log.printL(level, message)

    def prep_items(self, its):
        if self.num_items:
            try:
                return [int(v) for v in its]
            except ValueError:
                raise Warning("Integer values are expected for the items!", its)
        return its

    def map_items_to_groups(self, items_groups):
        map_itog = {}
        groups_disjoint = False
        for gi, grp in enumerate(items_groups):
            for i in grp:
                if i not in map_itog:
                    map_itog[i] = [gi]
                elif gi != map_itog[i][-1]:
                    map_itog[i].append(gi)
                    groups_disjoint = True
        return map_itog, groups_disjoint

    def are_groups_disjoint(self):
        return self.groups_disjoint

    def excludes_within_group(self):
        return self.x_within_group

    def get_exclusions(self, item=None):
        if item is not None:
            return self.item_exclusions.get(item)
        return self.item_exclusions

    def get_itog(self, item=None):
        if item is not None:
            return self.map_itog.get(item, [])
        return self.map_itog

    def get_nb_groups(self):
        return len(self.items_groups)

    def get_item_group(self, gi=None):
        if gi is None:
            return self.items_groups
        elif gi < len(self.items_groups):
            return self.items_groups[gi]
        return []

    def get_parameter(self, pname, default=None):
        return self.parameters.get(pname, default)

    def get_threshold(self, pname, gi=None, default=None):
        if gi is None:
            return self.thresholds.get(pname, default)

        else:
            totn = self.get_bthres_totn(pname)
            pv = self.thresholds.get(pname)
            gpv = None
            if "group_"+pname in self.thresholds and self.thresholds["group_"+pname][gi] is not None:
                gpv = self.thresholds["group_"+pname][gi]

            # if main threshold is stricter, use instead
            if totn == "nb_rows":  # about support, value of group is larger than main
                if pv is not None and pname[:4] == "min_" and (gpv is None or gpv < pv):
                    return pv

            elif totn == "nb_cols":  # about length, value of group is smaller than main
                if pv is not None and pname[:4] == "max_" and (gpv is None or gpv > pv):
                    return pv

            if gpv is None:
                return default
            return gpv

    def check_pattern_bthres(self, pattern, pname, thv=None):
        # pattern -> (itemset, count)
        if thv is None:  # no threshold, satisfied by default
            return True
        val = None
        tt = self.get_bthres_totn(pname)
        if tt == "nb_rows":
            val = pattern[1]
        elif tt == "nb_cols":
            val = len(pattern[0])
        if val is None:  # no value, satisfied by default
            return True

        if pname[:4] == "min_":
            return thv <= val
        if pname[:4] == "max_":
            return val <= thv
        return False

    def check_pattern_thresholds(self, pattern, bthress=None):
        # pattern -> (itemset, count)
        if bthress is None:
            bthress = self.get_bthres_names()
        for thn in bthress:
            if not self.check_pattern_bthres(pattern, thn, self.get_threshold(thn)):
                # print("FILTER OUT", thn, pattern)
                return False
        return True

    def mk_tree(self, data, transaction_counts=None, mining_order=None, suffix=[], gi=None, prev_d=None, initial_filter=False):
        item_counts = self.collect_frequent_items(data, transaction_counts, gi, prev_d, initial_filter)
        if mining_order is None:
            mining_order = self.sort_items(item_counts.keys(), self.get_itog(), item_counts, self.get_parameter("sorting"))
        else:
            mining_order = [x for x in mining_order if x in item_counts]
        return FPTree(data, item_counts, mining_order, transaction_counts, suffix)

    def get_conditional_tree(self, fp_tree, item, gi=None, last_round=False):
        tot_counts = sum([fp_tree.get_node_count(node_id) for node_id in fp_tree.get_links(item)])
        if last_round or tot_counts < self.get_threshold("min_supp", gi):
            return FPTree([[]], {}, [], [tot_counts], [item]+fp_tree.get_suffix())

        conditional_tree_trans = []
        conditional_tree_counts = []

        for node_id in fp_tree.get_links(item):

            path = fp_tree.get_path_to_root(node_id)
            nids, items, counts = zip(*path)
            conditional_tree_trans.append(items[:0:-1])
            conditional_tree_counts.append(counts[0])

        return self.mk_tree(conditional_tree_trans, conditional_tree_counts, fp_tree.get_mining_order(),
                            suffix=[item]+fp_tree.get_suffix(), gi=gi, prev_d=fp_tree.get_item_counts())

    @classmethod
    def collect_items_counts(tcl, data, transaction_counts=None):
        if type(data) is list:
            item_counts = {}
            for ti, transaction in enumerate(data):
                c = 1 if transaction_counts is None else transaction_counts[ti]
                for item in transaction:
                    if item in item_counts:
                        item_counts[item] += c
                    else:
                        item_counts[item] = c
        else:
            item_counts = dict([(i, int(c)) for (i, c) in enumerate(data.sum(0))])
        return item_counts

    def collect_frequent_items(self, data, transaction_counts=None, gi=None, prev_d=None, initial_filter=False):
        item_counts = self.collect_items_counts(data, transaction_counts)
        min_supp = self.get_threshold("min_supp", gi)
        if min_supp is not None or prev_counts is not None:
            item_counts = dict([(k, v) for (k, v) in item_counts.items()
                                if (prev_d is None or k in prev_d) and (min_supp is None or min_supp <= item_counts[k])])

        if initial_filter and gi is None:  # we are in a main tree
            # has the data been prefiltered by min_supp? if not, do it
            for ggi, grp in enumerate(self.get_item_group()):
                gmin_supp = self.get_threshold("min_supp", ggi)
                if gmin_supp is not None and gmin_supp != min_supp:
                    for i in grp:
                        if i in item_counts and item_counts[i] < gmin_supp:
                            del item_counts[i]
        return item_counts

    def check_add(self, pattern, patterns):
        if self.check_pattern_thresholds(pattern):
            patterns.add(pattern)

    def mine_patterns(self, fp_tree, patterns=None):
        if patterns is None:
            patterns = PatternStore()
        if fp_tree.has_frequent_suffix(self.get_threshold("min_supp")):
            self.check_add((fp_tree.get_suffixed_itemset(), fp_tree.get_suffix_count()), patterns)

        if fp_tree.get_suffix_len() >= self.get_threshold("max_len", default=self.get_nb_cols()):
            return patterns

        single_leaf_id = fp_tree.get_single_leaf_id()
        if single_leaf_id is None:
            return self.mine_sub_trees(fp_tree, patterns)
        elif single_leaf_id != -1:
            return self.generate_patterns_path(fp_tree, single_leaf_id, patterns)
        return patterns

    def generate_patterns_path(self, fp_tree, leaf_id, patterns={}):
        path = fp_tree.get_path_to_root(leaf_id)
        nids, items, counts = zip(*path)

        for i, c in enumerate(counts):
            self.check_add((fp_tree.get_suffixed_itemset(items[i:]), c), patterns)
        return patterns

    def mine_sub_trees(self, fp_tree, patterns={}):
        for item in fp_tree.get_mining_order():
            cond_fp_tree = self.get_conditional_tree(fp_tree, item)
            self.mine_patterns(cond_fp_tree, patterns)
        return patterns

    def find_frequent_patterns(self, data, patterns=None):
        fp_tree = self.mk_tree(data)
        self.to_log(1, "Made tree, mining...")
        return self.mine_patterns(fp_tree, patterns)


class FPGrowthMinerGroups(FPGrowthMiner):

    gthres_names = ["min_jacc"]
    gthres_dets = {"min_jacc": {"shrt": "j", "help": "minimum jaccard coefficient", "default": 0.1, "fnct": calc_jacc}
                   }

    def check_pattern_gthres(self, main_pattern, group_patterns, pname, thv=None):
        # pattern -> (itemset, count)
        if thv is None:  # no threshold, satisfied by default
            return True
        val = None
        # Global threshold, with function
        tf = self.get_gthres_fnct(pname)
        if tf is not None:
            val = tf(main_pattern, group_patterns)
        if val is None:  # no value, satisfied by default
            return True

        if pname[:4] == "min_":
            return thv <= val
        if pname[:4] == "max_":
            return val <= thv
        return False

    def check_pattern_thresholds(self, main_pattern, group_patterns=None, bthress=None, gthress=None, gi=None):
        if bthress is None:
            bthress = self.get_bthres_names()
        if gthress is None:
            gthress = self.get_gthres_names()

        if gi is not None:
            for thn in bthress:
                if not self.check_pattern_bthres(main_pattern, thn, self.get_threshold(thn, gi)):
                    # print("FILTER OUT", gi, thn, main_pattern)
                    return False

        else:
            for thn in bthress:
                if not self.check_pattern_bthres(main_pattern, thn, self.get_threshold(thn)):
                    # print("FILTER OUT", "MAIN", thn, main_pattern, group_patterns)
                    return False
            if group_patterns is not None:
                for thn in bthress:
                    for gi, gpatt in enumerate(group_patterns):
                        if not self.check_pattern_bthres(gpatt, thn, self.get_threshold(thn, gi)):
                            # print("FILTER OUT", gi, thn, main_pattern, group_patterns)
                            return False

                for thn in gthress:
                    if not self.check_pattern_gthres(main_pattern, group_patterns, thn, self.get_threshold(thn)):
                        # print("FILTER OUT", thn, main_pattern, group_patterns)
                        return False
        return True

    def __init__(self, params, items_groups, data=None, split=True, log=None):
        FPGrowthMiner.__init__(self, params, items_groups, data, split=split, log=log)

    def mk_trees(self, data, transaction_counts=None, mining_order=None, suffix=[]):
        item_counts = self.collect_frequent_items(data, transaction_counts, prev_d=self.get_itog(), initial_filter=self.get_parameter("initial_filter", True))
        if mining_order is None:
            mining_order = self.sort_items(item_counts.keys(), self.get_itog(), item_counts, self.get_parameter("sorting"))
        else:
            mining_order = [x for x in mining_order if x in item_counts]
        main_tree = FPTree(data, item_counts, mining_order, transaction_counts, suffix)
        trans, counts = main_tree.get_all_trans_counts()
        group_trees = []
        for grp in self.get_item_group():
            gitem_counts = dict([(k, item_counts[k]) for k in grp if k in item_counts])
            # gmining_order = [x for x in mining_order if x in gitem_counts]
            group_trees.append(FPTree(trans, gitem_counts, [], counts, suffix))
        return main_tree, group_trees

    def check_add(self, main_pattern, group_patterns, patterns):
        if self.check_pattern_thresholds(main_pattern, group_patterns):
            patterns.add(main_pattern, group_patterns)

    def mine_patterns(self, fp_tree, group_trees, patterns=None):
        if patterns is None:
            patterns = PatternStore()

        all_have_freq_suff = fp_tree.has_frequent_suffix(self.get_threshold("min_supp"))
        for gi, g_tree in enumerate(group_trees):
            if not g_tree.has_frequent_suffix(self.get_threshold("min_supp", gi)):
                all_have_freq_suff = False
                if not g_tree.has_frequent_item():
                    # some tree is dead (there is just the root, nothing frequent), no need to go further
                    return patterns

        if all_have_freq_suff:
            main_pattern = (fp_tree.get_suffixed_itemset(), fp_tree.get_suffix_count())
            group_patterns = [(g_tree.get_suffixed_itemset(), g_tree.get_suffix_count()) for g_tree in group_trees]
            self.check_add(main_pattern, group_patterns, patterns)

        if not fp_tree.has_frequent_item() or fp_tree.get_suffix_len() >= self.get_threshold("max_len", default=self.get_nb_cols()):
            return patterns

        return self.mine_sub_trees(fp_tree, group_trees, patterns)

    def generate_patterns_path(self, fp_tree, group_trees, leaf_id, gleaf_ids, patterns={}):
        support_counts = {}
        for gi, g_tree in enumerate(group_trees):
            support_counts[gi] = {}
            if g_tree.has_frequent_suffix(self.get_threshold("min_supp", gi)):
                support_counts[gi][None] = (g_tree.get_suffixed_itemset(), g_tree.get_suffix_count())
            if gleaf_ids[gi] != -1:
                gpath = g_tree.get_path_to_root(gleaf_ids[gi])
                gnids, gitems, gcounts = zip(*gpath)
                for i, c in enumerate(gcounts):
                    support_counts[gi][gitems[i]] = (g_tree.get_suffixed_itemset(gitems[i:]), c)
        latest_item = [None for g in group_trees]

        path = fp_tree.get_path_to_root(leaf_id)
        nids, items, counts = zip(*path)

        for i, c in enumerate(counts):
            main_pattern = (fp_tree.get_suffixed_itemset(items[i:]), c)
            group_patterns = []
            for gi in range(len(group_trees)):
                if items[i] in support_counts[gi]:
                    latest_item[gi] = items[i]
                if latest_item[gi] not in support_counts[gi]:
                    raise Warning("either suffix or item should be there...")
                group_patterns.append(support_counts[gi][latest_item[gi]])
            self.check_add(main_pattern, group_patterns, patterns)
        return patterns

    def mine_sub_trees(self, fp_tree, group_trees, patterns={}):
        max_len = self.get_threshold("max_len", default=self.get_nb_cols())
        remaining_ext = max_len - fp_tree.get_suffix_len()
        remaining_gext = [self.get_threshold("max_len", gi, default=self.get_nb_cols(gi))-g_tree.get_suffix_len() for gi, g_tree in enumerate(group_trees)]

        nb_items = len(fp_tree.get_mining_order())
        prj_items = fp_tree.get_mining_order()
        if self.get_parameter("symmetric_data", False) and fp_tree.get_suffix_len() == 0:
            prj_items = [i for i in prj_items if self.get_itog(i) == [0]]

        for ii, item in enumerate(prj_items):
            level = 2 + 3*len(fp_tree.get_suffix())
            if level < 10:
                self.to_log(level, "Projecting on item %s [%d/%d] (%d) %s" % (item, ii+1, nb_items, fp_tree.get_suffix_count(), fp_tree.get_suffix()))

            for gi in self.get_itog(item):
                if group_trees[gi].has_frequent_item(item):
                    if self.get_parameter("pruning", True):
                        self.mine_sub_trees_item(item, gi, fp_tree, group_trees, remaining_ext, remaining_gext, patterns)
                    else:
                        self.mine_sub_trees_item_no_pruning(item, gi, fp_tree, group_trees, remaining_ext, remaining_gext, patterns)

        return patterns

    def mine_sub_trees_item(self, item, gi, fp_tree, group_trees, remaining_ext, remaining_gext, patterns={}):

        # items to be removed from possible extensions
        remove_items = None
        if remaining_gext[gi] == 0:
            raise Warning("Should not get here?! remaining_gext[gi] == 0", item, gi)
            return patterns
        elif remaining_gext[gi] == 1:
            cond_g_tree = self.get_conditional_tree(group_trees[gi], item, gi, last_round=True)
            cpatt = (cond_g_tree.get_suffixed_itemset(), cond_g_tree.get_suffix_count())
            if not self.check_pattern_thresholds(cpatt, gi=gi):
                return patterns
            remove_items = self.get_item_group(gi)
        else:
            cond_g_tree = self.get_conditional_tree(group_trees[gi], item, gi)
            remove_items = self.get_exclusions(item)
            cond_g_tree.filter_items(remove_items)

        if remaining_ext == 0:
            raise Warning("Should not get here! remaining_ext == 0", item)
            return patterns
        elif remaining_ext == 1:
            cond_fp_tree = self.get_conditional_tree(fp_tree, item, last_round=True)
            cpatt = (cond_fp_tree.get_suffixed_itemset(), cond_fp_tree.get_suffix_count())
            if not self.check_pattern_thresholds(cpatt, gi=None):
                return patterns
        else:
            cond_fp_tree = self.get_conditional_tree(fp_tree, item)
            if self.excludes_within_group():
                cond_fp_tree.filter_items(remove_items)
        if not cond_fp_tree.has_frequent_suffix(self.get_threshold("min_supp")):
            return patterns

        cond_group_trees = list(group_trees)
        cond_group_trees[gi] = cond_g_tree
        if not self.excludes_within_group():
            remove_items = self.get_exclusions(item)
            for ggi in range(len(cond_group_trees)):
                if ggi != gi:
                    cond_group_trees[ggi] = cond_group_trees[ggi].filter_items(remove_items, inplace=False)

        # check accuracy
        min_jacc = self.get_threshold("min_jacc")
        if min_jacc is not None and cond_fp_tree.has_frequent_item():
            supp_bounds = sorted([(min_jacc*ctree.get_suffix_count(), ggi) for (ggi, ctree) in enumerate(cond_group_trees) if not ctree.has_frequent_item()])
            if len(supp_bounds) > 0:
                supp_bound, gi_bound = supp_bounds[-1]

                if supp_bound > self.get_threshold("max_supp", default=supp_bound+1):
                    return patterns
                if supp_bound > self.get_threshold("min_supp"):
                    if supp_bound > cond_fp_tree.get_suffix_count():
                        return patterns
                    else:
                        cond_fp_tree.set_local_support_threshold(supp_bound)

                for (ggi, ctree) in enumerate(cond_group_trees):
                    if ggi != gi_bound and ctree.has_frequent_item():
                        if supp_bound > self.get_threshold("min_supp", ggi):
                            if supp_bound > ctree.get_suffix_count():
                                return patterns
                        else:
                            cond_group_trees[ggi] = ctree.set_local_support_threshold(supp_bound, inplace=False)

        return self.mine_patterns(cond_fp_tree, cond_group_trees, patterns)

    def mine_sub_trees_item_no_pruning(self, item, gi, fp_tree, group_trees, remaining_ext, remaining_gext, patterns={}):

        # items to be removed from possible extensions
        remove_items = None
        cond_g_tree = self.get_conditional_tree(group_trees[gi], item, gi)
        remove_items = self.get_exclusions(item)
        cond_g_tree.filter_items(remove_items)

        cond_fp_tree = self.get_conditional_tree(fp_tree, item)
        if self.excludes_within_group():
            cond_fp_tree.filter_items(remove_items)
        if not cond_fp_tree.has_frequent_suffix(self.get_threshold("min_supp")):
            return patterns

        cond_group_trees = list(group_trees)
        cond_group_trees[gi] = cond_g_tree
        if not self.excludes_within_group():
            remove_items = self.get_exclusions(item)
            for ggi in range(len(cond_group_trees)):
                if ggi != gi:
                    cond_group_trees[ggi] = cond_group_trees[ggi].filter_items(remove_items, inplace=False)

        return self.mine_patterns(cond_fp_tree, cond_group_trees, patterns)

    def find_frequent_patterns(self, data, patterns=None):
        main_tree, group_trees = self.mk_trees(data)
        self.to_log(1, "Made trees, mining... %s" % main_tree)
        return self.mine_patterns(main_tree, group_trees, patterns)


if __name__ == "__main__":

    # python fim_mod_fpgrowth.py -s 50 --no-split datasets/abalone/abalone[L].dat datasets/abalone/abalone[R].dat
    # python fim_mod_fpgrowth.py -s 50 -T 4,5 -N 3 -M 0.1,- --nb_show 10  datasets/abalone/abalone[L].dat datasets/abalone/abalone[R].dat
    # python fim_mod_fpgrowth.py -j 0.9 -s 50 -M 2,2 --no-pruning --nb_show -1 datasets/abalone/abalone[L].dat datasets/abalone/abalone[R].dat
    # python fim_mod_fpgrowth.py -j 0.7 -s 50 -M 2,2 --sorting freq-desc --canon_show --nb_show -1 datasets/abalone/abalone[L].dat datasets/abalone/abalone[R].dat > out_50-.7_freq-descZ.txt
    import argparse

    parser = argparse.ArgumentParser(description='Mine itemset patterns')
    parser.add_argument("-o", "--output", type=str, help="folder to store the output", default="out.txt")
    for lng, dets in FPGrowthMinerGroups.get_bthres_details():
        parser.add_argument("-"+dets["shrt"], "--"+lng,
                            type=dets.get("type", float),
                            default=dets.get("default", argparse.SUPPRESS),
                            help=dets["help"])
        parser.add_argument("-"+dets["shrt"].upper(), "--group_"+lng,
                            type=str, default=argparse.SUPPRESS,
                            help=dets["help"]+"s for the different groups, comma separated")
    for lng, dets in FPGrowthMinerGroups.get_gthres_details():
        parser.add_argument("-"+dets["shrt"], "--"+lng,
                            type=dets.get("type", float),
                            default=dets.get("default", argparse.SUPPRESS),
                            help=dets["help"])

    parser.add_argument("--sorting", type=str, choices=FPGrowthMinerGroups.get_sorting_choices(), default=FPGrowthMinerGroups.get_sorting_default(), help="criterion to sort the items")

    parser.add_argument("--bm", action="store_true", help="the data come as a binary matrix")
    parser.add_argument("--fbm", action="store_true", help="force reading the data as a binary matrix")

    parser.add_argument("--log", type=int, help="level of logging", default=argparse.SUPPRESS)
    parser.add_argument("--nb_show", type=int, help="number of results to display in summary", default=argparse.SUPPRESS)
    parser.add_argument("--canon_show", action="store_true", help="display in canonical order")

    group = parser.add_mutually_exclusive_group()
    group.add_argument("--split", action="store_true", help="split items in groups", default=argparse.SUPPRESS)
    group.add_argument("--no-split", action="store_false", dest="split", help=argparse.SUPPRESS, default=argparse.SUPPRESS)

    group = parser.add_mutually_exclusive_group()
    group.add_argument("--pruning", action="store_true", help="try to prune candidates during enumeration", default=argparse.SUPPRESS)
    group.add_argument("--no-pruning", action="store_false", dest="pruning", help=argparse.SUPPRESS, default=argparse.SUPPRESS)

    parser.add_argument('input', metavar='I', type=str, nargs='+', help='names of the files containing the data')

    params = vars(parser.parse_args())

    if params["bm"]:
        # TODO to be implemented
        pass
    else:
        transactions = []
        items_groups = []
        for gi, fname in enumerate(params["input"]):
            items = set()
            with open(fname) as fp:
                offset = sum([len(i) for i in items_groups])
                for li, line in enumerate(fp):
                    if params["fbm"]:
                        trn = [offset+int(v) for v in line.strip().split()]
                    else:
                        trn = ["%s.%s" % (gi, v) for v in line.strip().split()]
                    if li == len(transactions):
                        transactions.append(trn)
                    elif li < len(transactions):
                        transactions[li].extend(trn)
                    items.update(trn)
                items_groups.append(sorted(items))

        if params["fbm"]:
            mat = numpy.zeros((len(transactions), sum([len(i) for i in items_groups])))
            for ti, t in enumerate(transactions):
                mat[ti, t] = 1
            item_counts = mat.sum(0)
            # keep_cis = numpy.where(item_counts >= support_threshold)[0]
            keep_cis = numpy.where(item_counts >= 0)[0]
            sids = sorted(keep_cis, key=lambda x: (item_counts[x], x))
            mat = mat[:, sids]
            map_items = dict([(v, i) for (i, v) in enumerate(sids)])
            mapped_items_groups = [[map_items[i] for i in gis if i in map_items] for gis in items_groups]
            items_str = {}
            for gi, gis in enumerate(items_groups):
                for i in gis:
                    if i in map_items:
                        items_str[map_items[i]] = "%s.%s" % (gi, i)

            data = mat
            data_items_groups = mapped_items_groups
        else:
            data = transactions
            data_items_groups = items_groups
            items_str = {}

    split = params.get("split", True) and len(data_items_groups) > 1

    tic = datetime.datetime.now()
    if split:
        fp_miner = FPGrowthMinerGroups(params, data_items_groups, data, log=params.get("log"))
    else:
        fp_miner = FPGrowthMiner(params, data_items_groups, data)
    patts = fp_miner.find_frequent_patterns(data)
    elps = datetime.datetime.now()-tic

    if "nb_show" in params:
        lines = []
        for i, (m, g) in sorted(patts.items()):
            if split:
                add_measures = ["%.4f" % calc_jacc(m, g)]
                tshow = [m]+g
            else:
                add_measures = []
                tshow = [(m, g)]
            if params.get("canon_show"):
                lines.append("\t".join([" ".join(["%s" % items_str.get(x, x) for x in sorted(gg[0])]) + " (%d)" % gg[1] for gg in tshow]+add_measures))
            else:
                lines.append("\t".join([" ".join(["%s" % items_str.get(x, x) for x in gg[0]]) + " (%d)" % gg[1] for gg in tshow]+add_measures))

        if params.get("canon_show"):
            lines.sort()
        nb_show = params["nb_show"]
        if 0 <= nb_show < len(lines):
            lines = lines[:nb_show]
        print("\n".join(lines))

    print("--- nb_trans=%s\tnb_patt=%d\ttime=%s" % (len(data), len(patts), elps))
