import sys
import pdb
import random
import datetime
import copy
from io import IOBase


def get_datetime_now():
    return datetime.datetime.now()
#     return 0  # datetime.datetime.now()


class Log(object):

    verbs = {"std": {"*": 1, "progress": 0, "result": 0, "error": 0, "tracks": 0},
             "comm": {"log": 1, "error": 1, "status": 1, "time": 1, "progress": 2, "result": 1, "tracks": 1},
             "inter": {"log": 1, "error": 0, "status": 1, "time": 0, "progress": 2, "result": 1, "tracks": 1},
             "quiet": {"*": 1, "error": 1, "status": 0, "result": 0, "progress": 0},
             "error": {"error": 1},
             "shut": {}}

    @classmethod
    def get_verbs(tcl, k="std", verb=None):
        v = tcl.verbs.get(k, {})
        if type(verb) is dict:
            v.update(verb)
        elif verb is not None:
            if "*" in v:
                v["*"] = verb
            if "log" in v:
                v["log"] = verb
        return v

    def __init__(self, k="std", verbosity=1, output='-', method_comm=None):
        self.tics = {None: get_datetime_now()}
        self.progress_ss = {"current": 0, "total": 0}
        self.out = []
        self.oqu = []
        self.verbosity = -1
        self.addOut(k, verbosity, output, method_comm)

    # FOR PICKLING !!
    def __getstate__(self):
        tmp = {}
        for k, v in self.__dict__.items():
            if k == 'out':
                tmp[k] = []
            else:
                tmp[k] = v
        return tmp

    def append(self, destination, message, type_message=None, source=None):
        destination.append(message)

    # THE CLOCK PART
    def getTic(self, id, name=None):
        if name is None:
            return self.tics[None]
        elif (id, name) in self.tics:
            return self.tics[(id, name)]
        else:
            return None

    def setTic(self, id, name):
        self.tics[(id, name)] = get_datetime_now()
        return self.tics[(id, name)]

    def setTac(self, id, name=None):
        if name is None:
            return (self.tics[None], get_datetime_now())
        elif (id, name) in self.tics:
            return (self.tics.pop((id, name)), get_datetime_now())

    def getTac(self, id, name):
        if name is None:
            return (self.tics[None], get_datetime_now())
        elif (id, name) in self.tics:
            return (self.tics[(id, name)], get_datetime_now())

    def clockTic(self, id, name=None, details={}):
        tic = self.setTic(id, name)
        if name is None:
            name = "\t"
        mess = "Start %s\t((at %s))" % (name, tic)
        if details is not None and len(details) > 0:
            mess += ("\t%s" % details)
        self.printL(1, mess, "time", id)

    def clockTac(self, id, name=None, details=""):
        tic, tac = self.getTac(id, name)
        if name is None:
            name = "\t"
        mess = "End %s\t((at %s, elapsed %s))" % (name, tac, tac-tic)
        if details is not None:
            mess += ("\t%s" % details)
        self.printL(1, mess, "time", id)
        return tac-tic
    # END CLOCK

    # THE PROGRESS PART
    def initProgressFull(self, constraints, explore_list=None, nbCols=[0, 0], level=-1, id=None):
        lim_max_init = constraints.getCstr("max_inits")
        if explore_list is not None:
            self.progress_ss["pairs_gen"] = sum([p[-1] for p in explore_list])
            if lim_max_init > -1:
                nb_max_init = min(len(explore_list), lim_max_init)
            else:
                nb_max_init = len(explore_list)
        else:
            self.progress_ss["pairs_gen"] = 0
            nb_max_init = max(0, lim_max_init)
        self.progress_ss["cand_var"] = 1
        self.progress_ss["cand_side"] = [nbCols[0] * self.progress_ss["cand_var"],
                                         nbCols[1] * self.progress_ss["cand_var"]]
        self.progress_ss["generation"] = constraints.getCstr("batch_cap")*sum(self.progress_ss["cand_side"])
        self.progress_ss["expansion"] = (constraints.getCstr("max_var", side=0)+constraints.getCstr("max_var", side=0)-2)*2*self.progress_ss["generation"]
        self.progress_ss["total"] = self.progress_ss["pairs_gen"] + nb_max_init * self.progress_ss["expansion"]
        self.progress_ss["current"] = 0
        if level > -1:
            self.printL(level, self.getProgress(), 'progress', id)

    def initProgressPart(self, constraints, reds, nbCols=[0, 0], level=-1, id=None):
        self.progress_ss["cand_var"] = 1
        self.progress_ss["cand_side"] = [nbCols[0] * self.progress_ss["cand_var"],
                                         nbCols[1] * self.progress_ss["cand_var"]]
        self.progress_ss["generation"] = constraints.getCstr("batch_cap") * sum(self.progress_ss["cand_side"])
        self.progress_ss["expansion"] = (constraints.getCstr("max_var", side=0) -
                                         min([constraints.getCstr("max_var", side=0)] + [len(r.queries[0]) for r in reds]) +
                                         constraints.getCstr("max_var", side=1) -
                                         min([constraints.getCstr("max_var", side=1)] + [len(r.queries[1]) for r in reds])) \
            * self.progress_ss["generation"]
        self.progress_ss["total"] = self.progress_ss["expansion"]
        self.progress_ss["current"] = 0
        if level > -1:
            self.printL(level, self.getProgress(), 'progress', id)

    def updateProgress(self, details={}, level=-1, id=None):
        if details is not None and len(details) > 0:
            if "pload" in details:
                self.progress_ss["current"] += details["pload"]
            elif len(details) == 1:
                if details["rcount"] > 0:
                    self.progress_ss["current"] += self.progress_ss["expansion"]

        if level > -1:
            self.printL(level, self.getProgress(), 'progress', id)

    def sendCompleted(self, id):
        self.printL(1, None, 'progress', id)

    def getProgress(self):
        return (self.progress_ss["total"], self.progress_ss["current"])
    # END PROGRESS

    def disp(self):
        tmp = "LOGGER"
        for out in self.out:
            tmp += "\n\t* %s -> %s" % (out["verbosity"],  out["destination"])
        for out in self.oqu:
            tmp += "\n\t* %s -> %s" % (out["verbosity"],  out["destination"])
        return tmp

    def resetOut(self):
        self.out = []
        self.oqu = []
        self.verbosity = -1

    def capVerbosity(self, v=0):
        for o in self.oqu+self.out:
            if type(o["verbosity"]) is int and o["verbosity"] > v:
                o["verbosity"] = v
            elif type(o["verbosity"]) is dict:
                ks = o["verbosity"].keys()
                for k in ks:
                    if o["verbosity"][k] > v:
                        o["verbosity"][k] = v

    def addOut(self,  k="std", verbosity=None, output='-', method_comm=None):
        if type(output) == str:
            if output in ['-', "stdout"]:
                tmp_dest = sys.stdout
            elif output == 'stderr':
                tmp_dest = sys.stderr
            else:
                try:
                    tmp_dest = open(output, 'w')
                except IOError:
                    return
        else:
            tmp_dest = output

        verbs = self.get_verbs(k, verbosity)
        max_v = max(verbs.values())
        if max_v > self.verbosity:
            self.verbosity = max_v

        # OK ADD OUTPUT
        if output is not None and type(output) is not str:
            self.oqu.append({"verbosity": verbs, "destination": tmp_dest, "method": method_comm})
        else:
            self.out.append({"verbosity": verbs, "destination": tmp_dest, "method": method_comm})
        return len(self.out)+len(self.oqu)-1

    def usesOutMethods(self):
        for out in self.out+self.oqu:
            if out["method"] is not None:
                return True
        return False

    def printL(self, level, message, type_message="*", source=None):
        for out in self.out+self.oqu:
            if (type_message in out["verbosity"] and level <= out["verbosity"][type_message]) or \
               (type_message not in out["verbosity"] and "*" in out["verbosity"] and level <= out["verbosity"]["*"]):

                if isinstance(out["destination"], IOBase):
                    if type_message == "*":
                        header = ""
                    else:
                        header = type_message
                    if source is None:
                        header += ""
                    else:
                        header += "@%s" % source
                    if len(header) > 0:
                        header = "[[%-10s]]\t" % header
                    out["destination"].write("%s%s\n" % (header, message))
                    out["destination"].flush()
                else:
                    out["method"](out["destination"], message,  type_message, source)

    def logResults(self, rcollect, lid, pid):
        if rcollect.getTracksSwitch() == "IN": ## log latest tracks?
            self.printL(1, rcollect.tracksToStr(True, line_break="\n\t"), "log", pid)

        if rcollect.getLen(lid) > 0:
            self.printL(1, {lid: rcollect.getItems(lid)}, 'result', pid)
            # self.printL(1, rcollect.getLatestTracks(), 'tracks', pid)
            if rcollect.getLen(lid) > 0:
                self.printL(1, "%d redescriptions [%s]" % (rcollect.getLen(lid), lid), 'status', pid)
                for red in rcollect.getItems(lid):
                    self.printL(10, "--- %s" % red, source=pid)
            else:
                self.printL(1, "No redescription [%s]" % lid, 'status', pid)
                
    def logAllTracksEnd(self, rcollect, lid, pid):
        if rcollect.getTracksSwitch() != "OFF":
            self.printL(1, "=========================================", "log", pid)
            if rcollect.getTracksSwitch() == "ON":
                self.printL(1, rcollect.tracksToStr(line_break="\n\t"), "log", pid)
            self.printL(1, "----- REDS:", "log", pid)
            fids = rcollect.getIidsList(lid)
            for i in rcollect.getItems():
                prefix = "+" if i.getUid() in fids else "-"
                self.printL(1, "%s %s" % (prefix, i), "log", pid)
