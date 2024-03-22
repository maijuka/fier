#!/usr/bin/python

import os.path
import sys
import re
import datetime
import numpy
import tempfile

try:
    from toolLog import Log
    from classData import Data
    from classQuery import Query
    from classRedescription import Redescription
    from classConstraints import Constraints
    from classPackage import IOTools
    from classContent import BatchCollection
    from factMiner import instMiner, StatsMiner
    from classRndFactory import RndFactory
except ModuleNotFoundError:
    from .toolLog import Log
    from .classData import Data
    from .classQuery import Query
    from .classRedescription import Redescription
    from .classConstraints import Constraints
    from .classPackage import IOTools
    from .classContent import BatchCollection
    from .factMiner import instMiner, StatsMiner
    from .classRndFactory import RndFactory

import pdb

# Run mining
# python exec_clired.py ~/short/rsmallZx/preferences_alpha.xml
# format results
# to filename generated from original
# python exec_clired.py --task printout --queries_alternate suff_named.tex ~/short/rsmallZx/preferences_alpha.xml
# or, equivalently
# python exec_clired.py --task printout ~/short/rsmallZx/preferences_alpha.xml suff_named.tex
# to different filename
# python exec_clired.py --task printout ~/short/rsmallZx/preferences_alpha.xml queries-v1.txt


def run(kw, loaded):
    params, data, logger, filenames = (loaded["params"], loaded["data"], loaded["logger"], loaded["filenames"])
    if "preferences_reader" in loaded and "params" in loaded:
        log_params = loaded["preferences_reader"].getManager().dispParameters(pv=loaded["params"], conf_filter=loaded["params"]["conf_defs"], sections=False, helps=False, defaults=False, only_core=False, xml=False, wdata=False)
        logger.printL(1, "Run with parameters:"+log_params, "log")

    miner = instMiner(data, params, logger, filenames=filenames)
    try:
        miner.full_run()
    except KeyboardInterrupt:
        # miner.initial_pairs.saveToFile()
        logger.printL(1, 'Stopped...', "log")

    mined_items = miner.rcollect.getItems("F")
    if len(mined_items) > 0:
        IOTools.outputResults(filenames, mined_items, data)
    logger.clockTac(0, None)


def run_filter(kw, loaded):
    params, data, logger, filenames = (loaded["params"], loaded["data"], loaded["logger"], loaded["filenames"])
    reds, srcs_reds, all_queries_src = IOTools.collectLoadedReds(loaded)
    constraints = Constraints(params, data, filenames=filenames)
    bc = BatchCollection(reds, tracks_switch=constraints.getCstr("tracks_switch"))

    actions = constraints.getActionList("postprocess")
    ids = bc.selected(actions, constraints=constraints)
    if constraints.getCstr("tracks_switch") != 'OFF':
        if logger is None:
            logger = Log(verbosity=params["verbosity"], output=filenames["logfile"])
        logger.printL(1, bc.tracksToStr(line_break="\n\t"), "log", "P")
    IOTools.outputResults(filenames, bc.getItems(ids), data)

    
def run_printout(kw, loaded):
    params, data, logger, filenames = (loaded["params"], loaded["data"], loaded["logger"], loaded["filenames"])
    reds, srcs_reds, all_queries_src = IOTools.collectLoadedReds(loaded)
    constraints = Constraints(params, data, filenames=filenames)
    # IOTools.saveAsPackage("/home/egalbrun/short/Zxx.siren", data, preferences=params, preferences_reader=loaded["preferences_reader"], reds=loaded.get("reds"), AR=constraints.getActionsRegistry())
    # data.saveExtensions(details={"dir": "/home/egalbrun/Desktop/"})
    IOTools.outputResults(filenames, reds, data)
    return loaded, reds, srcs_reds, all_queries_src


def run_test(kw, loaded):  # TODO
    params, data, logger, filenames = (loaded["params"], loaded["data"], loaded["logger"], loaded["filenames"])
    reds, srcs_reds, all_queries_src = IOTools.collectLoadedReds(loaded)

    for i, red in enumerate(reds):
        qt0, c0 = red.query(0).algNormalized()
        print(">>\t", red.query(0), "\tvs\t", qt0)
        qt1, c1 = red.query(1).algNormalized()
        print(">>\t", red.query(1), "\tvs\t", qt1)

        
def run_expand(kw, loaded):
    params, data, logger, filenames = (loaded["params"], loaded["data"], loaded["logger"], loaded["filenames"])
    reds, srcs_reds, all_queries_src = IOTools.collectLoadedReds(loaded)

    miner = instMiner(data, params, logger, filenames=filenames)
    collect_reds = []
    try:
        rcollect = miner.part_run({"reds": reds, "task": kw})  # , "side": 1})
        collect_reds.extend(rcollect.getItems("P"))
    except KeyboardInterrupt:
        # miner.initial_pairs.saveToFile()
        logger.printL(1, 'Stopped...', "log")
    IOTools.outputResults(filenames, collect_reds, data)
    logger.clockTac(0, None)


def run_folds(kw, loaded):
    params, data, logger, filenames = (loaded["params"], loaded["data"], loaded["logger"], loaded["filenames"])

    if "package" in filenames:
        head_p, tail_p = os.path.split(filenames["package"])
        root_p, ext_p = os.path.splitext(tail_p)
        filenames["basis"] = os.path.join(head_p, root_p)

    suff = None
    fold_cands = dict([(data.col(c[0], c[1]).getName(), c) for c in data.findCandsFolds(strict=False)])
    if len(params.get("folds_col", "")) > 0 and params.get("folds_col") in fold_cands:
        fci = fold_cands[params.get("folds_col")]
        logger.printL(2, "Using existing cross-fold subsets from side %s col %s" % fci, "log")
        sss = data.extractFolds(fci[0], fci[1])
        nb_folds = len(sss)
        suff = data.col(fci[0], fci[1]).getName()
    else:
        logger.printL(2, "Generating cross-fold subsets...", "log")
        sss = data.getFold(nbsubs=params["nb_folds"], coo_dim=params["coo_dim"], grain=params["grain"])
        suff = data.addFoldsCol()
        if suff is None:
            raise Warning("Creation of subsets failed!")
        nb_folds = params["nb_folds"]
        flds_pckgf = filenames["basis"] + ("_fold-%d:%s_empty.siren" % (nb_folds, suff))
        IOTools.saveAsPackage(flds_pckgf, data, preferences=params, preferences_reader=loaded["preferences_reader"])

    if suff is not None:
        print("SIDS", suff, sorted(data.getFoldsInfo()["fold_ids"].items(), key=lambda x: x[1]))
        print(data)
        flds_pckgf = filenames["basis"] + ("_fold-%d:%s.siren" % (nb_folds, suff))
        flds_statf = filenames["basis"] + ("_fold-%d:%s.txt" % (nb_folds, suff))

        stM = StatsMiner(data, params, logger)
        reds_list, all_stats, summaries, list_fields, stats_fields = stM.run_stats()

        rp = Redescription.getRP()
        flds_fk = filenames["basis"] + ("_fold-%d:%s-kall.txt" % (nb_folds, suff))
        with open(flds_fk, "w") as f:
            f.write(rp.printRedList(reds_list, fields=list_fields+["track"]))

        for fk, dt in summaries.items():
            flds_fk = filenames["basis"] + ("_fold-%d:%s-k%d.txt" % (nb_folds, suff, fk))
            with open(flds_fk, "w") as f:
                f.write(rp.printRedList(dt["reds"], fields=list_fields+["track"]))

        nbreds = numpy.array([len(ll) for (li, ll) in all_stats.items() if li > -1])
        tot = numpy.array(all_stats[-1])
        if nbreds.sum() > 0:
            summary_mat = numpy.hstack([numpy.vstack([tot.min(axis=0), tot.max(axis=0), tot.mean(axis=0), tot.std(axis=0)]), numpy.array([[nbreds.min()], [nbreds.max()], [nbreds.mean()], [nbreds.std()]])])

            info_plus = "\nrows:min\tmax\tmean\tstd\tnb_folds:%d" % (len(all_stats)-1)
            numpy.savetxt(flds_statf, summary_mat, fmt="%f", delimiter="\t", header="\t".join(stats_fields+["nb reds"])+info_plus)
            # IOTools.saveAsPackage(flds_pckgf, data, preferences=params, preferences_reader=loaded["preferences_reader"], reds=reds_list)
        else:
            with open(flds_statf, "w") as fo:
                fo.write("No redescriptions found")
        # for red in reds_list:
        #     print(red.disp())


def run_rnd(kw, loaded):
    params, data, logger, filenames = (loaded["params"], loaded["data"], loaded["logger"], loaded["filenames"])

    select_red = None
    if len(params.get("select_red", "")) > 0:
        select_red = params["select_red"]
    prec_all = None
    if params.get("agg_prec", -1) >= 0:
        prec_all = params["agg_prec"]
    count_vname = params.get("count_vname", "COUNTS")

    rf = RndFactory(org_data=data)
    with_traits = False
    if "traits_data" in filenames:
        traits_data = Data([filenames["traits_data"], None]+filenames["add_info"], filenames["style_data"])
        rf.setTraits(traits_data)
        with_traits = True

    if params.get("rnd_seed", -1) >= 0:
        rf.setSeed(params["rnd_seed"])

    stop = False
    for rnd_meth in params["rnd_meth"]:
        nb_copies = params["rnd_series_size"]
        if rnd_meth == "none":
            nb_copies = 1

        for i in range(nb_copies):
            sub_filenames = dict(filenames)
            suff = "_%s-%d" % (rnd_meth, i)
            sub_filenames["basis"] += suff
            for k in ["queries", "queries_named", "support"]:
                if k in sub_filenames:
                    parts = sub_filenames[k].split(".")
                    parts[-2] += suff
                    sub_filenames[k] = ".".join(parts)

            Dsub, sids, back, store = rf.makeupRndData(rnd_meth=rnd_meth, with_traits=with_traits, count_vname=count_vname, select_red=select_red, prec_all=prec_all)
            logger.printL(2, "STARTING Random series %s %d" % (rnd_meth, i), "log")
            logger.printL(2, Dsub, "log")

            miner = instMiner(Dsub, params, logger, filenames=sub_filenames)
            try:
                miner.full_run()
            except KeyboardInterrupt:
                # miner.initial_pairs.saveToFile()
                logger.printL(1, 'Stopped...', "log")
                stop = True

            IOTools.outputResults(sub_filenames, miner.final, Dsub)
            logger.clockTac(0, None)
            if stop:
                exit()


############################################
TASKS_METH = {"mine": run, "filter": run_filter, "printout": run_printout, "test": run_test,
              "expand": run_expand, "improve": run_expand, "folds": run_folds, "rnd": run_rnd}
TASKS_DEFAULT = "mine"
TASKS_LOAD = {"mine": {"log": False},  # first is the default one
              "filter": {"with_log": False, "queries_basic_dest": "?"},
              "test": {"with_log": False, "queries_basic_dest": "?"},
              "printout": {"with_log": False, "queries_basic_dest": "in", "conf_defs": ["dataext"]},
              "rnd": {"conf_defs": ["dataext", "rnd"]},
              "folds": {"conf_defs": ["folds"]}}
CONF_DEFS = [0, "dataext", "rnd", "folds"]
CONF_FILTER = [0]

# wrappers for IOTools functions, with task params


def load_all(arguments, conf_defs=None, tasks=None, tasks_default=None, tasks_load=None, conf_filter=None):
    if conf_defs is None:
        conf_defs = CONF_DEFS
    if tasks is None:
        tasks = TASKS_METH.keys()
    if tasks_default is None:
        tasks_default = TASKS_DEFAULT
    if tasks_load is None:
        tasks_load = TASKS_LOAD
    if conf_filter is None:
        conf_filter = CONF_FILTER
    return IOTools.loadAll(arguments, conf_defs, tasks, tasks_default, tasks_load, conf_filter)

# determine which task to run, load data etc. if not provided, and run


def do_task(sargs, conf_defs=CONF_DEFS, tasks_meths=TASKS_METH, tasks_default=None, tasks_load=None, conf_filter=None):
    success, loaded = load_all(sargs, conf_defs, tasks_meths.keys(), tasks_default, tasks_load, conf_filter)
    if not success:
        print(loaded)  # print help message
        sys.exit(2)
    if "filenames" in loaded:
        print(loaded["filenames"])
    task = loaded["params"]["task"]
    loaded["params"]["conf_defs"] = conf_defs
    print("Running %s" % task)
    return tasks_meths[task](task, loaded)


if __name__ == "__main__":
    do_task(sys.argv)
