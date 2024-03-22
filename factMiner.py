import re
# import numpy

try:
    from classRedescription import Redescription
    from classMiner import Miner
    from classMultiprocMiner import MultiprocMiner
    from classMinerLSH import MinerLSH

except ModuleNotFoundError:
    from .classRedescription import Redescription
    from .classMiner import Miner
    from .classMultiprocMiner import MultiprocMiner
    from .classMinerLSH import MinerLSH
    
import pdb

#######################################################################
# MINER INSTANCIATION
#######################################################################


def instMiner(data, params, logger=None, mid=None, qin=None, cust_params={}, filenames={}):
    if params["nb_processes"] > 1:
        Redescription.setUidGen(nv=None, step=None, mp_lock=True)
        return MultiprocMiner(data, params, logger, mid, qin, cust_params, filenames)
    else:
        if params.get("lsh_extensions", False):
            return MinerLSH(data, params, logger, mid, qin, cust_params, filenames)
        return Miner(data, params, logger, mid, qin, cust_params, filenames)


class StatsMiner:
    def __init__(self, data, params, logger=None, cust_params={}):
        self.data = data
        self.params = params
        self.logger = logger
        self.cust_params = cust_params

    def run_stats(self):
        rp = Redescription.getRP()
        modifiers = rp.getModifiersForData(self.data)
        modifiers["wfolds"] = True
        list_fields = rp.getListFields("basic", modifiers)
        exp_dict = rp.getExpDict(list_fields)
        stats_fields = [f for f in list_fields if not re.search("query", f) and not re.search("rid", f)]

        folds_info = self.data.getFoldsInfo()
        stored_folds_ids = sorted(folds_info["fold_ids"].keys(), key=lambda x: folds_info["fold_ids"][x])
        summaries = {}
        nbfolds = len(stored_folds_ids)
        for kfold in range(nbfolds):
            ids = {"learn": [], "test": []}
            for splt in range(nbfolds):
                if splt == kfold:
                    ids["test"].append(stored_folds_ids[splt])
                else:
                    ids["learn"].append(stored_folds_ids[splt])
            self.data.assignLT(ids["learn"], ids["test"])
            miner = instMiner(self.data, self.params, self.logger)
            try:
                miner.full_run()
            except KeyboardInterrupt:
                self.logger.printL(1, "Interrupted!", "status")

            stats = []
            reds = []
            for red in miner.rcollect.getItems("F"):
                red.recompute(self.data)  # to set the rsets
                evals_dict = rp.compEVals(red, exp_dict, details={})
                stats.append([evals_dict[f] for f in stats_fields])
                reds.append(red)
            summaries[kfold] = {"reds": reds, "stats": stats}

        reds_map = []
        stack_stats = []
        all_stats = {}
        for k, rr in summaries.items():
            all_stats[k] = rr["stats"]
            stack_stats.extend(rr["stats"])
            for ri, r in enumerate(rr["reds"]):
                found = None
                i = 0
                while i < len(reds_map) and found is None:
                    if reds_map[i]["red"].compare(r) == 0:
                        found = i
                    i += 1
                if found is None:
                    reds_map.append({"red": r, "pos": [(k, ri)]})
                else:
                    reds_map[found]["pos"].append((k, ri))
        reds_map.sort(key=lambda x: x["red"].getAcc(), reverse=True)
        reds_list = []
        for rr in reds_map:
            rr["red"].setTrack(rr["pos"])
            reds_list.append(rr["red"])
        all_stats[-1] = stack_stats
        return reds_list, all_stats, summaries, list_fields, stats_fields
