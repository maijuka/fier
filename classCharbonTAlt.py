try:
    from classCharbon import CharbonTree
    from classDTreeTools import DTreeToolsSKL as DTT
    # from classDTreeToolsDP import DTreeToolsDP as DTT
    from classQuery import *
    from classRedescription import Redescription
except ModuleNotFoundError:
    from .classCharbon import CharbonTree
    from .classDTreeTools import DTreeToolsSKL as DTT
    # from .classDTreeToolsDP import DTreeToolsDP as DTT
    from .classQuery import *
    from .classRedescription import Redescription

import pdb


class CharbonTCW(CharbonTree):
    name = "TreeCartWheel"

    def prepareVars(self, prev_result, current_side, singleD, in_data, cols_info):
        if prev_result.get("src") is not None and singleD:
            vrs = prev_result["src"].invCols()
            if cols_info is None:
                ncandidates = [vvi for vvi in range(in_data[current_side].shape[1]) if vvi not in vrs]
            else:
                ncandidates = [vvi for vvi in range(in_data[current_side].shape[1]) if cols_info[current_side][vvi][1] not in vrs]
            feed_data = in_data[current_side][:, ncandidates]

        elif prev_result.get("tree") is not None and singleD:
            vrs = [prev_result["candidates"][x] for x in prev_result["tree"].getFeatures()]
            if cols_info is None:
                ncandidates = [vvi for vvi in range(in_data[current_side].shape[1]) if vvi not in vrs]
            else:
                ttm = [cols_info[current_side][c][1] for c in vrs]
                ncandidates = [vvi for vvi in range(in_data[current_side].shape[1]) if cols_info[current_side][vvi][1] not in ttm]
            feed_data = in_data[current_side][:, ncandidates]
        else:
            ncandidates = None
            feed_data = in_data[current_side]
        return feed_data, ncandidates

    def getTreeCandidates(self, side, data, more, in_data, cols_info):
        if side is None:
            jj0, results0 = self.getSplit(0, in_data, more["target"], data.isSingleD(), cols_info, more.get("src"))
            jj1, results1 = self.getSplit(1, in_data, more["target"], data.isSingleD(), cols_info, more.get("src"))
            if jj0 > jj1:
                jj, results = (jj0, results0)
            else:
                jj, results = (jj1, results1)
        else:
            jj, results = self.getSplit(side, in_data, more["target"], data.isSingleD(), cols_info, more.get("src"))

        if results is not None and results[0].get("tree") is not None and results[1].get("tree") is not None:
            redex = self.get_redescription(results, data, cols_info)
            return redex
        return None

    def splitting_with_depth(self, in_data, in_target, in_depth, in_min_bucket,  split_criterion="gini"):
        dtc = DTT.fitTree(in_data, in_target, in_depth, in_min_bucket, split_criterion, random_state=0, logger=self.logger)
        suppv = dtc.getSupportVect(in_data)
        if sum(suppv) < in_min_bucket or len(suppv)-sum(suppv) < in_min_bucket:
            return None, None
        return dtc, suppv

    def getSplit(self, side, in_data, target, singleD=False, cols_info=None, src=None):
        best = (0, None)
        results = {0: {}, 1: {}}
        current_tids = [0, 1]
        current_side = 1-side
        if sum(target) >= self.constraints.getCstr("min_itm_c") and len(target)-sum(target) >= self.constraints.getCstr("min_itm_c"):
            results[current_tids[side]]["support"] = target
            results[current_tids[side]]["src"] = src
            rounds = 0
        else:
            rounds = -1

        while rounds < self.constraints.getCstr("max_rounds") and rounds >= 0:
            rounds += 1
            feed_data, ncandidates = self.prepareVars(results[current_tids[1-current_side]], current_side, singleD, in_data, cols_info)
            dtc, suppv = self.splitting_with_depth(feed_data, results[current_tids[1-current_side]].get("support"),
                                                   self.constraints.getCstr("max_var", side=current_side), self.constraints.getCstr("min_itm_c"),
                                                   split_criterion=self.constraints.getCstr("split_criterion"))

            if dtc is None or (results[current_tids[current_side]].get("tree") is not None and results[current_tids[1-current_side]].get("tree") is not None
                               and results[current_tids[current_side]].get("support") is not None and all(results[current_tids[current_side]].get("support") == suppv)):
                # nothing found or no change
                rounds = -1
            else:
                current_dt = {"support": suppv, "tree": dtc, "candidates": ncandidates}
                current_tids[current_side] = len(results)
                results[current_tids[current_side]] = current_dt
                current_side = 1-current_side
                jj = CharbonTree.getJacc(results[current_tids[0]].get("support"), results[current_tids[1]].get("support"))
                if jj > best[0] and results[current_tids[current_side]].get("tree") is not None:
                    best = (jj, (current_tids[0], current_tids[1]))
        if best[1] is not None:
            return (best[0], {0: results[best[1][0]], 1: results[best[1][1]]})
        return best


class CharbonTSprit(CharbonTCW):
    name = "TreeSprit"

    def getSplit(self, side, in_data, target, singleD=False, cols_info=None, src=None):
        best = (0, None)
        results = {0: {}, 1: {}}
        current_tids = [0, 1]
        current_side = 1-side
        if sum(target) >= self.constraints.getCstr("min_itm_c") and len(target)-sum(target) >= self.constraints.getCstr("min_itm_c"):
            results[current_tids[side]]["support"] = target
            results[current_tids[side]]["src"] = src
            rounds = 0
        else:
            rounds = -1

        depth = [2, 2]
        while rounds >= 0 and (depth[0] <= self.constraints.getCstr("max_var", side=0) or depth[1] <= self.constraints.getCstr("max_var", side=1)):
            rounds += 1
            feed_data, ncandidates = self.prepareVars(results[current_tids[1-current_side]], current_side, singleD, in_data, cols_info)
            dtc, suppv = self.splitting_with_depth(feed_data, results[current_tids[1-current_side]].get("support"),
                                                   depth[current_side], self.constraints.getCstr("min_itm_c"),
                                                   split_criterion=self.constraints.getCstr("split_criterion"))

            if dtc is None or (results[current_tids[current_side]].get("tree") is not None and results[current_tids[1-current_side]].get("tree") is not None
                               and results[current_tids[current_side]].get("support") is not None and all(results[current_tids[current_side]].get("support") == suppv)):
                # nothing found or no change
                rounds = -1
            else:
                depth[current_side] += 1
                current_dt = {"support": suppv, "tree": dtc, "candidates": ncandidates}
                current_tids[current_side] = len(results)
                results[current_tids[current_side]] = current_dt
                current_side = 1-current_side
                jj = CharbonTree.getJacc(results[current_tids[0]].get("support"), results[current_tids[1]].get("support"))
                if jj > best[0] and results[current_tids[current_side]].get("tree") is not None:
                    best = (jj, (current_tids[0], current_tids[1]))
        if best[1] is not None:
            return (best[0], {0: results[best[1][0]], 1: results[best[1][1]]})
        return best


class CharbonTSplit(CharbonTCW):

    name = "TreeSplit"

    def getTreeCandidates(self, side, data, more, in_data, cols_info):
        # print("\n\n>>> getTreeCandidates", side, more["side"], more["involved"], more["src"], more["target"].sum())
        results = self.getSplit(side, in_data, more["target"], data.isSingleD(), cols_info, more.get("src"))
        return self.extractRed(results, data, cols_info)

    def extractRed(self, results, data, cols_info):
        r = None
        accs = [0]
        others = [None]
        for ri in range(1, len(results)):
            ori = ri-1  # compute accuracy with
            if ri % 2 == 0 and not self.constraints.getCstr("splittrees_update_target"):
                ori = 1
            if results[ri].get("tree") is not None and results[ri].get("support") is not None and results[ri-1] is not None and results[ori].get("support") is not None:
                accs.append(CharbonTree.getJacc(results[ri]["support"], results[ori]["support"]))
                others.append(ori)
            else:
                accs.append(0)
                others.append(None)

        ji = CharbonTree.pickStep(accs, self.constraints.getCstr("min_impr"))
        if ji is not None and others[ji] is not None:
            r = self.get_redescription({results[others[ji]]["side"]: results[others[ji]],
                                        results[ji]["side"]: results[ji]}, data, cols_info)

        # r_fin = None
        # if len(results) > 1 and results[-1] is not None and results[-2] is not None:
        #     r_fin = self.get_redescription({results[-1]["side"]: results[-1],
        #                                     results[-2]["side"]: results[-2]}, data, cols_info)

        # print("\nStats (%d/%d)\t%s" % (ji, len(accs)-1, accs))
        # print(r)
        # if r is None or (ji < (len(results)-1) and not r.sameRectangles(r_fin)):
        #     print("INSTEAD OF")
        #     print(r_fin)
        return r

    def splitting_with_depth_both(self, side, in_data, in_target, singleD=False, cols_info=None, in_depth=1, in_min_bucket=0, split_criterion="gini", results=[]):
        trg = in_target
        change = [-1, -1]
        for ci, current_side in enumerate([1-side, side]):
            dtc, suppv = (None, None)
            side_depth = min(self.constraints.getCstr("max_var", side=current_side), in_depth)
            feed_data, ncandidates = self.prepareVars(results[-1], current_side, singleD, in_data, cols_info)
            dtc, suppv = self.splitting_with_depth(feed_data, trg, side_depth, in_min_bucket, split_criterion)
            if dtc is None:  # split failed, fall back
                if ci == 1:  # complement the split obtained
                    results.append(results[-2])
                # print(">>> SPLIT FAILED")
                return change
            else:
                # no previous results to compare to
                if results[-2] is None or results[-2].get("support") is None:
                    change[current_side] = 1
                ### change in support
                elif any(results[-2].get("support") != suppv):
                    change[current_side] = 1
                else:
                    change[current_side] = 0
                # print(">>> GOT TREE\tside=%d depth=%d change=%d supp_sum=%d" % (current_side, side_depth, change[current_side], suppv.sum()))
            results.append({"support": suppv, "tree": dtc, "candidates": ncandidates, "side": current_side})
            trg = suppv
        return change

    def getSplit(self, side, in_data, target, singleD=False, cols_info=None, src=None):
        depth = self.constraints.getCstr("splittrees_init_depth")
        trg = target
        results = [None, {"support": trg, "tree": None, "candidates": None, "side": side, "src": src}]

        if sum(target) >= self.constraints.getCstr("min_itm_c") and len(target)-sum(target) >= self.constraints.getCstr("min_itm_c"):
            while depth > 0 and depth <= max(self.constraints.getCstr("max_var", side=0), self.constraints.getCstr("max_var", side=1)):
                change = self.splitting_with_depth_both(side if self.constraints.getCstr("splittrees_dyn_side") else 1,
                                                        in_data, trg, singleD, cols_info,
                                                        in_depth=depth, in_min_bucket=self.constraints.getCstr("min_itm_c"),
                                                        split_criterion=self.constraints.getCstr("split_criterion"),
                                                        results=results)
                if -1 not in change and (not 0 in change or (not self.constraints.getCstr("splittrees_both_sides_change") and 1 in change)):
                    if self.constraints.getCstr("splittrees_update_target"):
                        trg = results[-1].get("support")
                    else:  # source for the actual target
                        results[-1]["src"] = results[1]["src"]
                    depth += 1
                elif depth == 1:
                    depth = 2  # try depth 2 anyway, keeping target
                else:
                    depth = 0
        return results
        # # Check if we have both vectors (split was successful on the left and right matrix)
        # if results[0].get("tree") is None or results[1].get("tree") is None:
        #     if prev_results is not None:
        #         # Check if left tree was able to split
        #         if results[0].get("tree") is None:
        #             results[0] = prev_results[0]
        #         if results[1].get("tree") is None:
        #             results[1] = prev_results[1]
        #     depth = 0
        # else:
        #     # Either have no previous results or have successful splits and have to check whether trees changed
        #     if prev_results is None or \
        #         ((prev_results[0].get("support") is None or any(prev_results[0].get("support") != results[0].get("support"))) and \
        #          (prev_results[1].get("support") is None or any(prev_results[1].get("support") != results[1].get("support")))):
        #         prev_results = results
        #         depth += 1
        #     else:
        #         depth = 0
