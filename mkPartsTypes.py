import numpy
import random
import re
import pdb

status_labels = {True: "x", False: "o", None: "m"}
io_labels = {True: "into", False: "out", None: "imiss", "T": "tot"}

pairs_ord = dict([(v, k) for (k, v) in enumerate([(True, False), (False, True), (True, True), (False, False),
                                                  (True, None), (None, True), (False, None), (None, False), (None, None)])])
io_ord = dict([(v, k) for (k, v) in enumerate(["T", True, False, None])])


def land(A, B):
    if A is None and B == False:
        return B
    return A and B


def lor(A, B):
    if A is None and B == False:
        return None
    return A or B


def rSings(l):
    return ", ".join(["E%s%s" % (status_labels[s[0]], status_labels[s[1]]) for s in sorted(l, key=lambda x: pairs_ord[x])])


def rPairs(l):
    return ", ".join(["(%s, E%s%s)" % (io_labels[p[-1]], status_labels[p[0]], status_labels[p[1]])
                      for p in sorted(l, key=lambda x: (io_ord[x[-1]], pairs_ord[(x[0], x[1])]))])


def parse_xom(s):
    parts = set()
    for match in re.finditer("E[xom][xom]", s):
        e = match.group()
        parts.add((e[1] == "x" if e[1] != "m" else None, e[2] == "x" if e[2] != "m" else None))
    return parts


def group_xom_iom(s, nb_states=3):  # nb iom is number of different part types into, out, imiss
    tots, parts = (set(), set())
    unixoms = {}
    for ss in s:
        if ss[:2] not in unixoms:
            unixoms[ss[:2]] = set()
        unixoms[ss[:2]].add(ss[2])
    for xom, ioms in unixoms.items():
        if len(ioms) == nb_states:  # all parts present
            tots.add((xom[0], xom[1], "T"))
        else:
            parts.update([(xom[0], xom[1], i) for i in ioms])
    return tots, parts


def mkPartsType(defI, defU, defL, defO=None, states=None):
    if states is None:
        states = [True, False, None]

    set_N = set()
    for A in states:
        for B in states:
            set_N.add((A, B))

    org_sets = {}
    org_sets["I"] = set_N.intersection(parse_xom(defI))
    org_sets["U"] = set_N.intersection(parse_xom(defU))
    org_sets["L"] = set_N.intersection(parse_xom(defL))
    if defO is None:
        org_sets["O"] = set_N.difference(org_sets["U"])
    else:
        org_sets["O"] = set_N.intersection(parse_xom(defO))

    org_sets["R"] = set([(x[1], x[0]) for x in org_sets["L"]])

    org_sets["D"] = org_sets["U"].difference(org_sets["I"])

    results = {False: dict([((A, B), []) for (A, B) in set_N]),
               True: dict([((A, B), []) for (A, B) in set_N])}
    for (A, B) in set_N:
        for X in states:
            results[False][(land(A, X), B)].append((A, B, X))
            results[True][(lor(A, X), B)].append((A, B, X))

    new_sets = {}
    for s in ["I", "U", "D", "O", "L", "R"]:
        for op in [False, True]:
            new_sets[(s, op)] = group_xom_iom(set().union(*[results[op][p] for p in org_sets[s]]), nb_states=len(states))
    new_sets[("C", False)] = (set(), set([(xl, xr, i) for (xl, xr, i) in new_sets[("O", False)][1] if (xl, xr) in org_sets["D"]]))
    new_sets[("C", True)] = (set(), set([(xl, xr, i) for (xl, xr, i) in new_sets[("I", True)][1] if (xl, xr) in org_sets["D"]]))

    dsp = ""
    for name, ss in [("inter", "I"), ("uncovered", "O"), ("diff", "D"), ("suppL", "L"), ("suppR", "R")]:
        dsp += "\n%sself.IDS_%s = [%s]" % (8*" ", name, rSings(org_sets[ss]))
    dsp += "\n"
    for name, ss, fv in [("fixnum", "I", 0), ("varnum", "I", 1), ("fixden", "U", 0), ("varden", "U", 1)]:
        dsp += "\n%sself.IDS_%s = [[%s],\n%s[%s]]" % (8*" ", name, rPairs(new_sets[(ss, False)][fv]), 27*" ", rPairs(new_sets[(ss, True)][fv]))
    for name, ss in [("out", "O"), ("cont", "C"), ("nsupp", "L")]:
        dsp += "\n%sself.IDS_%s = [[%s],\n%s[%s]]" % (8*" ", name, rPairs(list(new_sets[(ss, False)][0])+list(new_sets[(ss, False)][1])),
                                                      27*" ", rPairs(list(new_sets[(ss, True)][0])+list(new_sets[(ss, True)][1])))
    return dsp


if __name__ == "__main__":

    # set_N = parse_xom("Exo, Eox, Exx, Eoo, Exm, Emx, Eom, Emo, Emm")
    # print("# Truth table\n# A  B   A and B   A or B\n" + (23*"-"))
    # for (A, B) in set_N:
    #     print("# %s  %s%7s%9s" % (status_labels[A], status_labels[B], status_labels[land(A, B)], status_labels[lor(B, A)]))

    # states = [True, False, None]
    # assigns = {(False, 0): [], (True, 0): [], (False, 1): [], (True, 1): []}
    # for A in states:
    #     for B in states:
    #         for X in states:
    #             assigns[(False, 0)].append("((%s, E%s%s), E%s%s)" % (io_labels[X], status_labels[A], status_labels[B], status_labels[land(A, X)], status_labels[B]))
    #             assigns[(True, 0)].append("((%s, E%s%s), E%s%s)" % (io_labels[X], status_labels[A], status_labels[B], status_labels[lor(A, X)], status_labels[B]))
    #             assigns[(False, 1)].append("((%s, E%s%s), E%s%s)" % (io_labels[X], status_labels[A], status_labels[B], status_labels[A], status_labels[land(B, X)]))
    #             assigns[(True, 1)].append("((%s, E%s%s), E%s%s)" % (io_labels[X], status_labels[A], status_labels[B], status_labels[A], status_labels[lor(B, X)]))

    # for k in [(False, 0), (True, 0), (False, 1), (True, 1)]:
    #     print("assigns[%s] = [%s]" % (k, ", ".join(assigns[k])))

    typs = [("none", {"defI": "|Exx|", "defU": "|Exx|+|Exo|+|Eox|", "defL": "|Exo|+|Exx|", "states": [True, False]}),
            ("rejective",   {"defI": "|Exx|",
                             "defU": "|Exx|+|Exo|+|Eox|",
                             "defL": "|Exo|+|Exx|", "defO": "|Eoo|"}),
            ("optimistic",  {"defI": "|Exx|+|Exm|+|Emx|+|Emm|",
                             "defU": "|Exo|+|Eox|+|Exx|+|Exm|+|Emx|+|Emm|",
                             "defL": "|Exo|+|Exx|+|Exm|+|Emx|+|Emm|"}),
            ("pessimistic", {"defI": "|Exx|",
                             "defU": "|Exo|+|Eox|+|Exx|+|Exm|+|Emx|+|Eom|+|Emo|+|Emm|",
                             "defL": "|Exo|+|Exx|+|Exm|+|Emo|"}),  # Emm can't be in both support sets, nor in either one...
            ("positive",    {"defI": "|Exx|+|Exm|+|Emx|+|Emm|",
                             "defU": "|Exo|+|Eox|+|Exx|+|Exm|+|Emx|+|Eom|+|Emo|+|Emm|",
                             "defL": "|Exo|+|Exx|+|Exm|+|Emo|+|Emx|+|Emm|"}),
            ("negative",    {"defI": "|Exx|",
                             "defU": "|Exo|+|Eox|+|Exx|+|Exm|+|Emx|",
                             "defL": "|Exo|+|Exx|+|Exm|"})]

    for typ, args in typs:
        dsp = mkPartsType(**args)
        print("\n\n    def init_part_ids_%s(self):%s" % (typ, dsp))
