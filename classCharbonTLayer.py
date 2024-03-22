import numpy

try:
    from classCharbon import CharbonTree
    from classDTreeTools import DTree
    from classDTreeTools import DTreeToolsSKL as DTT
    # from classDTreeToolsDP import DTreeToolsDP as DTT
    from classQuery import *
    from classRedescription import *
except ModuleNotFoundError:
    from .classCharbon import CharbonTree
    from .classDTreeTools import DTree
    from .classDTreeTools import DTreeToolsSKL as DTT
    # from .classDTreeToolsDP import DTreeToolsDP as DTT
    from .classQuery import *
    from .classRedescription import *

import pdb


def isBootstrapTree(tree):
    return tree.isEmpty() or tree.isEmptyFeat(tree.getFeature(0))


def isSpecialRootTree(tree):
    return tree.isSpecialFeat(tree.getFeature(0))


def init_tree(data, side, more={}, cols_info=None, tid=None):
    parent_tree = {"id": tid,
                   "candidates": list(range(data[side].shape[1])),
                   "involv": []}
    if cols_info is not None:
        vid = None
        invol_narrow = more.get("involved", [])
        if len(invol_narrow) == 1:
            vid = invol_narrow[0]
        ttm = [cols_info[side][c][1] for c in invol_narrow]
        invol = [kk for (kk, vv) in cols_info[side].items() if vv[1] in ttm]

        parent_tree["involv"] = invol
        parent_tree["tree"] = DTree({"supp_pos": more["target"], "feat": more.get("src")})
        parent_tree["support"] = more["target"]
    else:
        parent_tree["tree"] = DTree({"n": data[side].shape[0]})
        parent_tree["support"] = parent_tree["tree"].getSuppVect()
    return parent_tree


def initialize_treepile(data, side_ini, more={}, cols_info=None):
    trees_pile = [[[]], [[]]]
    trees_store = {}

    PID = 0
    trees_pile[1-side_ini][-1].append(PID)
    trees_store[PID] = init_tree(data, 1-side_ini, tid=PID)

    PID += 1
    trees_pile[side_ini][-1].append(PID)
    trees_store[PID] = init_tree(data, side_ini, more, cols_info, tid=PID)

    PID += 1
    return trees_pile, trees_store, PID


def get_trees_pair(data, trees_pile, trees_store, side_ini, max_var, min_bucket, split_criterion="gini", PID=0, singleD=False, cols_info=None):
    current_side = 1-side_ini
    while len(trees_pile[0][-1]) > 0 or len(trees_pile[1][-1]) > 0:
        if len(trees_pile[current_side][-1]) > 0:
            # try extending ...
            trees_pile[current_side].append([])
            # ... if max depth not exhausted #### account for dummy tree on non init side when counting depth
            if (len(trees_pile[current_side]) - 1*(current_side != side_ini)) <= max_var[current_side]:

                opp_match_level = -1
                if len(trees_pile[1-current_side][-1]) == 0:
                    opp_match_level = -2
                target = numpy.sum([trees_store[tid]["support"] for tid in trees_pile[1-current_side][opp_match_level]], axis=0)
                for gpid in trees_pile[current_side][-2]:
                    gp_tree = trees_store[gpid]
                    candidates = [v for v in gp_tree["candidates"] if v not in gp_tree.get("involv", [])]
                    if singleD:
                        for ggid in trees_pile[1-current_side][opp_match_level]:
                            for vv in trees_store[ggid]["involv"]:
                                try:
                                    candidates.remove(vv)
                                except ValueError:
                                    pass

                    leaves, dt = [], None
                    if len(candidates) > 0:
                        leaves = gp_tree["tree"].collectLeaves()
                        dt = data[current_side][:, candidates]
                    for leaf in leaves:
                        mask = list(gp_tree["tree"].getSuppSet(leaf))
                        # print("SUB TARGET", current_side, sum(target[mask]), candidates)
                        if sum(target[mask]) > min_bucket and (len(mask)-sum(target[mask])) > min_bucket:
                            tree_rpart = DTT.fitTree(dt[mask, :], target[mask], in_depth=1,
                                                     in_min_bucket=min_bucket, split_criterion=split_criterion, random_state=0)
                            # print("SUB TREE", tree_rpart)
                            if not tree_rpart.isEmpty():
                                tree_rpart.applyMaskSuppSets(mask)
                                # CHECK SUPPORT
                                # if numpy.sum(split_tree["over_supp"][mask] != tree_rpart.getSupportVect(dt[mask,:])) > 0:
                                vrs = tree_rpart.getFeatures()
                                support = tree_rpart.computeSupp(dt, ids=mask)
                                ninvolved = [candidates[c] for c in vrs]
                                # if cols_info is None:
                                #     # ncandidates = [vvi for vvi in candidates if vvi not in vrs]
                                #     ninvolved = list(vrs)
                                # else:
                                #     ttm = [cols_info[current_side][c][1] for c in vrs]
                                #     # ncandidates = [vvi for vvi in candidates if cols_info[current_side][vvi][1] not in ttm]
                                #     ninvolved = [vvi for (vvi, vv) in cols_info[current_side].items() if vv[1] in ttm]
                                split_tree = {"id": PID, "tree": tree_rpart, "candidates": candidates,
                                              "branch": (gp_tree["id"], leaf),
                                              "support": support, "mask": mask,
                                              "involv": ninvolved}
                                trees_pile[current_side][-1].append(PID)
                                trees_store[PID] = split_tree
                                PID += 1

        current_side = 1-current_side
    return trees_pile, trees_store, PID


def generate_steps_seq(supps, side_ini):
    accs = []
    depths = []
    iSides = [0, 0]
    side_increment = side_ini
    if len(supps[1-side_increment]) == 1:
        iSides[1-side_increment] = -1
    while iSides[0] < len(supps[0]) and iSides[1] < len(supps[1]):
        accs.append(CharbonTree.getJacc(supps[0][iSides[0]], supps[1][iSides[1]]))
        depths.append(tuple(iSides))
        iSides[side_increment] += 1

        if iSides[1-side_increment] != -1:
            if iSides[side_increment]+1 >= len(supps[side_increment]):
                # other side has not reached end yet, keep going
                iSides[side_increment] = -1
                # else will stop
            side_increment = 1 - side_increment
    return accs, depths


def prune_trees(trees_pile, trees_store, side_ini, min_impr):
    levels = [[], []]
    supps = [[], []]
    counts = [[], []]
    for side in [0, 1]:
        if len(trees_pile[side][0]) == 1 and isBootstrapTree(trees_store[trees_pile[side][0][0]]["tree"]):
            trees_pile[side].pop(0)
        while len(trees_pile[side]) > 0 and len(trees_pile[side][-1]) == 0:
            trees_pile[side].pop()

        for tl in trees_pile[side]:
            levels[side].append(numpy.sum([trees_store[tid]["support"] for tid in tl], axis=0))
            if len(supps[side]) > 0:
                Z = numpy.zeros(levels[side][-1].shape)
                for tid in tl:
                    if "mask" in trees_store[tid]:
                        Z[trees_store[tid]["mask"]] += 1
                supps[side].append(supps[side][-1].copy())
                supps[side][-1][Z > 0] = levels[side][-1][Z > 0]
                if numpy.sum(Z > 1) > 0 or numpy.sum(levels[side][-1][Z == 0]) > 0:
                    raise Exception("Something is wrong with support masking in the layeredtrees construction !")
                    # pdb.set_trace()
                    # print(numpy.sum(Z > 1))
                # assert(numpy.sum(Z > 1) == 0)
                # assert(numpy.sum(levels[side][-1][Z==0]) == 0)
            else:
                supps[side].append(levels[side][-1].copy())
            counts[side].append(supps[side][-1].sum())

    accs, depths = generate_steps_seq(supps, side_ini)
    if len(depths) > 0 and depths[-1] not in [(-1, len(supps[1])-1), (len(supps[0])-1, -1)]:
        raise Exception("Something went wrong when pairing depth in layeredtrees !")
        # pdb.set_trace()
        # accs, depths = generate_steps_seq(supps, side_ini)

    ji = CharbonTree.pickStep(accs, min_impr)

    if ji is not None:
        dLHS, dRHS = depths[ji]
        if dLHS > -1 and (dLHS+1) < len(trees_pile[0]):
            del trees_pile[0][dLHS+1:]
            # print("PRUNED LHS")
        if dRHS > -1 and (dRHS+1) < len(trees_pile[1]):
            del trees_pile[1][dRHS+1:]
            # print("PRUNED RHS")
    # print(trees_pile)
    # for pid, t in trees_store.items():
    #     print("\n--- TREE %s %s (%s)" % (pid, t.get("branch"), t["support"].support()))
    #     print(t["tree"])
    # pdb.set_trace()
    # new_suppv = numpy.sum([trees_store[tid]["support"] for tid in trees_pile[side_ini][-1]], axis=0)


def graft_trees(trees_store, pile, in_data=None):
    list_nodes = {}
    tids = []
    ids_init = None
    for ti, t in enumerate(pile):
        if ti == 0 and len(t) == 1:
            if isBootstrapTree(trees_store[t[0]]["tree"]):
                continue  # Empty tree to bootstrap, drop
            elif isSpecialRootTree(trees_store[t[0]]["tree"]):
                ids_init = list(numpy.where(trees_store[t[0]]["support"])[0])
            tids.extend(t)
        else:
            tids.extend(t)
    graft_points = {}
    for tid in tids:
        list_nodes[tid] = list(trees_store[tid]["tree"].getNodeIds())
        if trees_store[tid].get("branch", [None])[0] in tids:
            list_nodes[trees_store[tid]["branch"][0]][trees_store[tid]["branch"][1]] = None
            graft_points[trees_store[tid]["branch"]] = (tid, 0)
    map_nodes = {}
    for tid in tids:
        for node_id, v in enumerate(list_nodes[tid]):
            if v is not None:
                map_nodes[(tid, node_id)] = len(map_nodes)
    feats, thres, clss, chls, chrs = ([], [], [], [], [])
    for tid in tids:
        ct = trees_store[tid]["tree"]
        # print("%d\t%s" % (tid, ct))
        for (i, v) in enumerate(list_nodes[tid]):
            if v is not None:
                if trees_store[tid].get("candidates") is None or ct.isEmptyFeat(ct.getFeature(i)) or ct.isSpecialFeat(ct.getFeature(i)):
                    feats.append(ct.getFeature(i))
                else:
                    feats.append(trees_store[tid]["candidates"][ct.getFeature(i)])
                thres.append(ct.getThreshold(i))
                clss.append(ct.getClass(i))
                for cid, chs in [(ct.getChildrenLeft(i), chls), (ct.getChildrenRight(i), chrs)]:
                    if (tid, cid) in graft_points:
                        chs.append(map_nodes.get(graft_points[(tid, cid)], -1))
                    else:
                        chs.append(map_nodes.get((tid, cid), -1))
        # print("\tFeat: %s\n * Thres: %s\n * Clss: %s\n * ChL: %s\n * ChR: %s" % (feats, thres, clss, chls, chrs))
    dtc = DTree((feats, thres, clss, chls, chrs))
    if dtc.isEmpty():
        return (None, None)
    suppv = dtc.getSupportVect(in_data, ids=ids_init)
    return dtc, suppv


class CharbonTLayer(CharbonTree):

    name = "TreeLayer"

    def getTreeCandidates(self, side, data, more, in_data, cols_info):
        # print("\n\n>>> getTreeCandidates", side, more["side"], more["involved"], more["src"], more["target"].sum())
        trees_pile, trees_store, PID = initialize_treepile(in_data, side, more, cols_info=cols_info)
        trees_pile, trees_store, PID = get_trees_pair(in_data, trees_pile, trees_store, side,
                                                      max_var=[self.constraints.getCstr("max_var", side=0), self.constraints.getCstr("max_var", side=1)],
                                                      min_bucket=self.constraints.getCstr("min_itm_c"),
                                                      split_criterion=self.constraints.getCstr("split_criterion"),
                                                      PID=PID, singleD=data.isSingleD(), cols_info=cols_info)

        prune_trees(trees_pile, trees_store, side, self.constraints.getCstr("min_impr"))
        # print("PILE", trees_pile)
        # for ti, tree in trees_store.items():
        #     print("------ %s %s" % (ti, tree["tree"]))
        dtc0, suppv0 = graft_trees(trees_store, trees_pile[0], in_data[0])
        dtc1, suppv1 = graft_trees(trees_store, trees_pile[1], in_data[1])
        if dtc0 is not None and dtc1 is not None:
            results = {0: {"tree": dtc0, "support": suppv0}, 1: {"tree": dtc1, "support": suppv1}}
            redex = self.get_redescription(results, data, cols_info)
            # print(">>> GOT REDESCRIPTION", redex)
            # print(">>> RED SUPP COUNTS\t", redex.getLenL(), redex.getLenR())
            return redex
        return None
