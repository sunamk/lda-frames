#
# Copyright (C) 2012 Jiri Materna <xmaterna@fi.muni.cz>
# Licensed under the GNU GPLv3 - http://www.gnu.org/licenses/gpl-3.0.html
#

"""
Library for accessing semantic frames stored in a database.
"""

import shelve
import anydbm


class LDAF():
    
    def __init__(self, path):
        try:
            self.frames = shelve.open(path + "frames.db")
            self.frameDist = shelve.open(path + "framedist.db")
            self.realDist = shelve.open(path + "realdist.db")
        except anydbm.error, msg:
            print "Cannot open database files."
            print msg

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
        
        
