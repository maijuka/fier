import numpy
try:

    from classData import Data
    from classCharbon import CharbonXaust
    from classQuery import *
    from classRedescription import Redescription
    from fim_mod_fpgrowth import FPGrowthMiner, FPGrowthMinerGroups
except ModuleNotFoundError:
    from .classData import Data
    from .classCharbon import CharbonXaust
    from .classQuery import *
    from .classRedescription import Redescription
    from .fim_mod_fpgrowth import FPGrowthMiner, FPGrowthMinerGroups
import pdb


class FIStore:
    def __init__(self, lits_dets, log=None, no_groups=False):
        self.log = log
        self.lits_dets = lits_dets
        self.patterns = []
        self.no_groups = no_groups

    def to_log(self, level, message):
        if type(self.log) is int:
            if level < self.log:
                print(message)
        elif self.log is not None:
            self.log.printL(level, message)

    def __len__(self):
        return len(self.patterns)

    def add(self, main_pattern, group_patterns=None):
        if self.no_groups:
            lsupp_inter, lsupp_LHS, lsupp_RHS = (0, main_pattern[1], 0)
            q_LHS, q_RHS = (main_pattern[0], ())
        else:
            lsupp_inter, lsupp_LHS, lsupp_RHS = (main_pattern[1], group_patterns[0][1], group_patterns[1][1])
            q_LHS, q_RHS = (group_patterns[0][0], group_patterns[1][0])
        acc = lsupp_inter/(lsupp_LHS+lsupp_RHS-lsupp_inter)
        self.patterns.append((q_LHS, q_RHS, acc, lsupp_LHS-lsupp_inter, lsupp_RHS-lsupp_inter, lsupp_inter))
        if len(self.patterns) % 100 == 0:
            self.to_log(3, "Collected %d patterns" % len(self.patterns))
        elif len(self.patterns) % 10 == 0:
            self.to_log(8, "Collected %d patterns" % len(self.patterns))

    def sort(self):
        self.patterns.sort(key=lambda x: (x[2], x[5], x[3], x[4]), reverse=True)

    def mkQueries(self, q_LHS, q_RHS):
        lits0 = [self.lits_dets[q][1].copy() for q in q_LHS]
        lits1 = [self.lits_dets[q][1].copy() for q in q_RHS]
        return (Query(False, lits0), Query(False, lits1))

    def getPatterns(self):
        return self.patterns

    def getPattern(self, i):
        return self.patterns[i]

    def getQueryPair(self, i):
        return self.mkQueries(self.patterns[i][0], self.patterns[i][1])

    def getRed(self, i, data):
        return Redescription.fromQueriesPair(self.getQueryPair(i), data, copyQ=False)

    def getReds(self, data):
        return [self.getRed(i, data) for i in range(len(self))]


class CharbonXFIM(CharbonXaust):

    name = "XaustFIM"

    def isIterative(self):
        return False

    def computeExpansions(self, data, initial_terms, logger=None):
        sorted_terms = sorted(initial_terms, key=lambda c: (c.getSuppLen(), c.getSide(), c.getCid()))
        no_groups = False
        lits_dets = [(c.getSide(), c.getLit()) for c in sorted_terms]
        if data.isSymmetric():
            if self.constraints.getCstr("min_fin_acc") > -1:
                lits_dets = []
                for c in sorted_terms:
                    lits_dets.extend([(c.getSide(), c.getLit()), (1-c.getSide(), c.getLit())])
            else:
                no_groups = True

        from_col = {}
        items_groups = [[], []]
        for i, (side, lit) in enumerate(lits_dets):
            items_groups[side].append(i)
            if data.isSingleD():
                k = (None, lit.colId())
            else:
                k = (side, lit.colId())
            if k not in from_col:
                from_col[k] = []
            from_col[k].append(i)

        item_exclusions = {}
        for k, its in from_col.items():
            if len(its) > 1:
                for i in its:
                    if i not in item_exclusions:
                        item_exclusions[i] = set(its)
                    else:
                        item_exclusions[i].update(its)

        bmat = data.getLitsToBinaryMatrix(lits_dets)

        cs = FIStore(lits_dets, log=logger, no_groups=no_groups)
        params = {}
        params["single_data"] = data.isSingleD()
        params["symmetric_data"] = data.isSymmetric()
        if self.constraints.getCstr("min_itm_in") > -1:
            params["min_supp"] = self.constraints.getCstr("min_itm_in")
        if self.constraints.getCstr("min_itm_out") > -1:
            params["group_max_supp"] = [data.nbRows() - self.constraints.getCstr("min_itm_out")]*2
        if self.constraints.getCstr("max_var_s0") > -1 or self.constraints.getCstr("max_var_s1") > -1:
            params["group_max_len"] = [self.constraints.getCstr("max_var_s0") if self.constraints.getCstr("max_var_s0") > -1 else None,
                                       self.constraints.getCstr("max_var_s1") if self.constraints.getCstr("max_var_s1") > -1 else None]
        if self.constraints.getCstr("max_var_s0") > -1 and self.constraints.getCstr("max_var_s1") > -1:
            params["max_len"] = self.constraints.getCstr("max_var_s0") + self.constraints.getCstr("max_var_s1")
        if self.constraints.getCstr("min_fin_acc") > -1:
            params["min_jacc"] = self.constraints.getCstr("min_fin_acc")
        if len(item_exclusions) > 0:
            params["item_exclusions"] = item_exclusions
        params["initial_filter"] = False

        if no_groups:
            fp_miner = FPGrowthMiner(params, items_groups, bmat, log=logger)
        else:
            fp_miner = FPGrowthMinerGroups(params, items_groups, bmat, log=logger)
        patts = fp_miner.find_frequent_patterns(bmat, patterns=cs)
        return patts.getReds(data)
