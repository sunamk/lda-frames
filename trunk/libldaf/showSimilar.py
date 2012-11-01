#!/usr/bin/python
#
# Copyright (C) 2012 Jiri Materna <xmaterna@fi.muni.cz>
# Licensed under the GNU GPLv3 - http://www.gnu.org/licenses/gpl-3.0.html
#

"""
This script prints out similar lexical units based on the frame distribution similarity.
If the predicate is not given, the script prints out all supported predicates.


USAGE: ./showSimilar.py path [predicate]
    path                  Path to database files.
    predicate             Predicate the similar predicates will be shown for.
    -n, --number          Number of similar predicates to be shown (all by default).
    -s, --similarity      Minimum similarity for predicate to be shown (unlimited by default).
    -h, --help            Print this help.
"""

import sys
import getopt
import shelve
import anydbm

def printSimilarities(db, query, number, similarity):
    try:
        for i, pred  in enumerate(db[query]):
            if (i == number and i != 0) or pred[1] < similarity:
                break
            print "%s\t%s\t%f" % (query, pred[0], pred[1])
    except KeyError:
        print "No such predicate in the database (%s)." % query
        sys.exit(3)

if __name__ == "__main__":
    number = 0
    similarity = 0

    try:
        opts, args = getopt.getopt(sys.argv[1:], "hn:s:", ["help", "number=", "similarity"])
    except getopt.error, msg:
        print msg
        print "for help use --help"
        sys.exit(1)
    for o, a in opts:
        if o in ("-h", "--help"):
            print __doc__
            sys.exit(0)
        if o in ("-n", "--number"):
            number = int(a)
        if o in ("-s", "--similarity"):
            similarity = float(a)

    if len(args) < 1 or  len(args) > 2:
        print "Wrong number of parameters."
        print __doc__
        sys.exit(1)
    query = None
    if len(args) == 2:
        query = args[1]

    path  = args[0]
    if not path.endswith("/"): path += "/"
    try:
        similarities = shelve.open(path + "framesimilarities.db", flag='r')

        if query != None:
            printSimilarities(similarities, query, number, similarity)
        else:    
            for pred in similarities:
                printSimilarities(similarities, pred, number, similarity)


        similarities.close()
    except anydbm.error, msg:
        print "Cannot open database files."
        sys.exit(2)

