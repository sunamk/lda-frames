#!/usr/bin/python
#
# Copyright (C) 2012 Jiri Materna <xmaterna@fi.muni.cz>
# Licensed under the GNU GPLv3 - http://www.gnu.org/licenses/gpl-3.0.html
#

"""
Dictionary class defined in this file serves for representing dictionary translating
predicates and realisations to numbers and vice versa.
"""

class Dictionary:

    def __init__(self):
        self.p2id = {}
        self.id2p = []
        self.r2id = {}
        self.id2r = []
        self.relations = []

    def addPredicate(self, pred):
        if pred not in self.p2id:
            self.id2p.append(pred)
            self.p2id[pred] = len(self.id2p)

    def addRealisation(self, real):
        if real not in self.r2id and not real == "-":
            self.id2r.append(real)
            self.r2id[real] = len(self.id2r)

    def pred2id(self, pred):
        try:
            return self.p2id[pred]
        except KeyError:
            return False

    def id2pred(self, i):
        try:
            return self.id2p[i-1]
        except IndexError:
            return False

    def real2id(self, real):
        if real == "-": return 0
        try:
            return self.r2id[real]
        except KeyError:
            return False

    def id2real(self, i):
        if i == 0: return "-"
        try:
            return self.id2r[i-1]
        except IndexError:
            return False

    def setRelations(self, rels):
        self.relations = rels

    def getRelations(self):
        return self.relations

