#
# Copyright (C) 2013 Jiri Materna <xmaterna@fi.muni.cz>
# Licensed under the GNU GPLv3 - http://www.gnu.org/licenses/gpl-3.0.html
#

"""
Library for accessing semantic frames stored in a database.
"""

import shelve
import anydbm
import numpy
from scipy.spatial.distance import euclidean


class LDAF():
    
    def __init__(self, path):
        try:
            self.frames = shelve.open(path + "frames.db", flag='r')
            self.frameDist = shelve.open(path + "framedist.db", flag='r')
            self.realDist = shelve.open(path + "realdist.db", flag='r')
            self._SQRT2 = numpy.sqrt(2)
        except anydbm.error, msg:
            print "Cannot open database files."
            raise

    def getRelations(self):
        return self.frames["__description__"]

    def getPredicates(self):
        return sorted(self.frameDist.keys())


    def getFrame(self, frameid):
        return self.frames[frameid]

    def getFrameDist(self, predicate, limit=None):
        try:
           if limit==None:
               return self.frameDist[predicate]
           else:
               return self.frameDist[predicate][:limit]
        except KeyError:
            return None
    
    def getRealDist(self, realization, limit=None):
       if limit==None:
           return self.realDist[realization]
       else:
           return self.realDist[realization][:limit]
        
    def getF(self):
        return self.frames["__F__"]

    def getR(self):
        return self.frames["__R__"]

    def getHellingerDistance(self, v1, v2):
        return  euclidean(numpy.sqrt(v1),  numpy.sqrt(v2)) / self._SQRT2
        
