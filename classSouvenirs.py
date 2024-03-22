try:
    from classQuery import *
except ModuleNotFoundError:
    from .classQuery import *

import pdb


class Souvenirs(object):

    format_index_pref = ':%(lengthL)i:%(lengthR)i:%(side)i:%(op)i:'
    format_index_suff = '%(buk)i:%(col)i:'

    def __init__(self, nAvailableMo=None, nAmnesic=False):
        self.indexes = {}
        if nAvailableMo is None:
            self.active = False
            self.availableMo = [set(), set()]
            self.amnesic = True
        else:
            self.active = True
            self.availableMo = nAvailableMo
            self.amnesic = nAmnesic

    def isActive(self):
        return self.active

    def cutOffSide(self, side):
        self.availableMo[side] = None

    def wasCutOff(self, side):
        return self.availableMo[side] is None

    def getAvailableS(self):
        return range(len(self.availableMo))

    def copyAvailableCols(self, side=None):
        if side is None:
            return [self.copyAvailableCols(side) for side in self.getAvailableS()]
        else:
            if self.availableMo[side] is None:
                return None
            else:
                return set(self.availableMo[side])

    def getAvailableCols(self, side=None):
        if side is None:
            return [self.getAvailableCols(side) for side in self.getAvailableS()]
        else:
            if self.availableMo[side] is None:
                return set()
            else:
                return self.availableMo[side]

    def getNbAvailableCols(self, side=None):
        if side is None:
            return [self.getNbAvailableCols(side) for side in self.getAvailableS()]
        else:
            if self.availableMo[side] is None:
                return 0
            else:
                return len(self.availableMo[side])

    def getCommonAvailableCols(self, sides=[0, 1]):
        common = None
        for side in sides:
            if self.availableMo[side] is not None:
                if common is None:
                    common = set(self.availableMo[side])
                else:
                    common.intersection_update(self.availableMo[side])
        if common is None:
            return set()
        return common

    def isAmnesic(self):
        return self.amnesic

    def __str__(self):
        if self.isAmnesic():
            x = 'amnesic'
        else:
            x = '%i indexes' % len(self.indexes)
        return x+", Availables: " + " + ".join(["%i" % i for i in self.getNbAvailableCols()])

    def add(self, red):
        if len(red) > 2:
            ix = red.getUid()
            for indx in self.makeOwnIndexes(red):
                if indx in self.indexes:
                    self.indexes[indx].add(ix)
                else:
                    self.indexes[indx] = set([ix])

    def makeIndexes(self, red, lengthL, lengthR, side, op):
        return red.queries[side].makeIndexes((self.format_index_pref % {'lengthL': lengthL, 'lengthR': lengthR, 'side': side, 'op': op}) + self.format_index_suff)

    def makeOwnIndexes(self, red):
        indexes = []
        for side in [0, 1]:
            indexes.extend(self.makeIndexes(red, red.length(0), red.length(1), side, red.queries[side].opBuk(0)))
        return indexes

    def extOneStep(self, red, side):
        if red.containsAnon():
            return red.invColsSide(side, ex_anon=True)
        cols_ext = red.invColsSide(side)
        if not self.isAmnesic():
            other_side = 1-side
            lengthL = red.length(0)+(1-side)
            lengthR = red.length(1)+side

            queries_ids_other_side = self.lookForQueries(self.makeIndexes(red, lengthL, lengthR, other_side, red.queries[other_side].opBuk(0)))
            if len(queries_ids_other_side) > 0:
                #                pdb.set_trace()
                if red.length(side) == 1:
                    queries_ids = self.lookForQueries(self.makeIndexes(red, lengthL, lengthR, side, Op(True)))
                    queries_ids |= self.lookForQueries(self.makeIndexes(red, lengthL, lengthR, side, Op(False)))
                else:
                    queries_ids = self.lookForQueries(self.makeIndexes(red, lengthL, lengthR, side, red.queries[side].opBuk(0)))

                queries_ids &= queries_ids_other_side
                if len(queries_ids) > 0:

                    # print('EXTENSIONS-%i------%s <=> %s -------' % (side, red.queries[0].dispIds(),red.queries[1].dispIds()))
                    #                     for i in queries_ids:
                    #                         print(self.queriesList[i][0].dispIds() + '<=>' + self.queriesList[i][1].dispIds())
                    #                     print('------------------------')
                    cols_ext = self.colsExtending(queries_ids, side)  # includes already used cols
#             if len(cols_ext) > red.queries[side].length() :
#                 print('EXCLUDED-%i-------------' % side)
#                 print(cols_ext)
#                 print('------------------------')

        return cols_ext

    def initialAvailable(self, initialPairRed):
        return [self.getAvailableCols(0) - self.extOneStep(initialPairRed, 0),
                self.getAvailableCols(1) - self.extOneStep(initialPairRed, 1)]

    def colsExtending(self, queries_ids, side):
        cols = set()
        for idr in queries_ids:
            cols |= self.getItem(idr).invColsSide(side)
        return cols

    def lookForQueries(self, indexes_p):
        if len(indexes_p) > 0 and not self.isAmnesic():
            query_ids = set([-1])
            id_inds = 0
            while id_inds < len(indexes_p) and indexes_p[id_inds] in self.indexes and len(query_ids) > 0:
                if id_inds == 0:
                    query_ids = set(self.indexes[indexes_p[id_inds]])
                else:
                    query_ids &= self.indexes[indexes_p[id_inds]]
                id_inds += 1
            if id_inds != len(indexes_p):
                query_ids = set()
            return query_ids
        else:
            return set()

    def updateSeen(self, redList):
        for red in redList:
            if not self.isAmnesic():
                self.add(red)

    def hasSeen(self, red):
        lengthL = red.length(0)
        lengthR = red.length(1)

        indices = self.makeIndexes(red, lengthL, lengthR, 0, red.query(0).opBuk(0))
        indices += self.makeIndexes(red, lengthL, lengthR, 1, red.query(1).opBuk(0))

        matches = list(self.lookForQueries(indices))
        seen = False
        mi = 0
        while not seen and mi < len(matches):
            if self.getItem(matches[mi]).bothSidesIdentical(red):
                seen = True
            else:
                mi += 1

        if not seen and not self.isAmnesic():
            self.add(red)
        return seen
