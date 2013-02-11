#!/usr/bin/python
#
# Copyright (C) 2013 Jiri Materna <xmaterna@fi.muni.cz>
# Licensed under the GNU GPLv3 - http://www.gnu.org/licenses/gpl-3.0.html
#

"""
This script takes as input frame distributions and generates a similarity vector
for each lexical unit stored in a database.

USAGE: ./generate_similarities.py path
    path                  Path to the data files.
    -h, --help            Print this help.
"""

import sys
import getopt
import anydbm
import numpy
import libldaf
import shelve
import cPickle

if __name__ == "__main__":
    try:
        opts, args = getopt.getopt(sys.argv[1:], "h", ["help"])
    except getopt.error, msg:
        print msg
        print "for help use --help"
        sys.exit(1)
    for o, a in opts:
        if o in ("-h", "--help"):
            print __doc__
            sys.exit(0)
    if len(args) != 1:
        print "Wrong number of parameters."
        print __doc__
        sys.exit(1)

    path = args[0]
    if not path.endswith("/"): path += "/"
    try:
        ldaf = libldaf.LDAF(path)
        dist_database = shelve.open(path+"framesimilarities.db", flag='n',
            protocol=cPickle.HIGHEST_PROTOCOL)
        print "Creating frame distribution vectors."
        predicates = ldaf.getPredicates()
        vectors = []
        for pred in predicates:
            vec = numpy.zeros(ldaf.getF())
            distribution = ldaf.getFrameDist(pred)
            for f in distribution:
                vec[f[0]-1] = f[1]
            vectors.append(vec)
        print "Computing predicate distances."
        for i, pred1 in enumerate(predicates):
            if (i % 10 == 0 and i != 0):
                print "Progress at predicate no. %d/%d." % (i, len(predicates))
            similarPredicates = []
            for j, pred2 in enumerate(predicates):
                if i == j:
                    continue
                distance = ldaf.getHellingerDistance(vectors[i], vectors[j])
                if distance < 1.0:
                    similarPredicates.append((pred2, 1-distance))
            dist_database[pred1] = sorted(similarPredicates, key=lambda k: -k[1])
        dist_database.close()
        print "Finished."
                    

    except anydbm.error, msg:
        sys.exit(2)

