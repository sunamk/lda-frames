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
import nltk
from scipy.spatial.distance import euclidean


class LDAF():
    
    def __init__(self, path):
        try:
            self.frames = shelve.open(path + "frames.db", flag='r')
            self.frameDist = shelve.open(path + "framedist.db", flag='r')
            self._SQRT2 = numpy.sqrt(2)
    
            realHist = shelve.open(path + "realdist.db", flag='r')
            self.realDist = {}
            for role, hist in realHist.iteritems():
                freqdist = nltk.FreqDist(hist)
                freqsum = freqdist.N()
                relativeFreqs = map(lambda x:(x[0], x[1]*1.0/freqsum), freqdist.items())
                smoothedFreqs = nltk.probability.LaplaceProbDist(freqdist)
                self.realDist[role] = (relativeFreqs, smoothedFreqs)
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
    
    def getRealDist(self, role, limit=None):
       if limit==None:
           return self.realDist[role][0]
       else:
           return self.realDist[role][0][:limit]

    def getRealProb(self, realization, role):
        return self.realDist[str(role)][1].prob(realization)

    def getFrameProbs(self, predicate, realizations):
        frameProbs = []
        frameDist = self.getFrameDist(predicate)
        for (f, p) in frameDist:
            frame = self.getFrame(str(f))
            if len(frame) != len(realizations):
                raise Exception("The number of slots is inconsistent.")
            for s, r in enumerate(frame):
                p *= self.getRealProb(realizations[s], r)
            frameProbs.append((f, p)) 
        return sorted(frameProbs, key=lambda x: -x[1])
        
        
    def getF(self):
        return self.frames["__F__"]

    def getR(self):
        return self.frames["__R__"]

    def getHellingerDistance(self, v1, v2):
        return  euclidean(numpy.sqrt(v1),  numpy.sqrt(v2)) / self._SQRT2
        
