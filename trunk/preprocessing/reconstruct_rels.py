#!/usr/bin/python
#
# Copyright (C) 2013 Jiri Materna <xmaterna@fi.muni.cz>
# Licensed under the GNU GPLv3 - http://www.gnu.org/licenses/gpl-3.0.html


"""
This script reconstructs grammatical relations from relations.dat and relations.dict files
and prints them to stdout.

USAGE: ./reconstruct_rels.py path
    path            Path to input files.
    -h, --help            Print this help.
"""

import sys
import getopt
import cPickle
from ldafdict import Dictionary

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
    dat = open(path + "relations.dat")
    dictionary = cPickle.load(open(path + "relations.dict"))
    print ";" + "\t".join(dictionary.getRelations())
    u = 1
    for line in dat.xreadlines():
        predicate = dictionary.id2pred(u)
        for realisation in line.strip().split("\t"):
            print predicate + "\t" + "\t".join(map(lambda x: dictionary.id2real(int(x)),
                    realisation.split(" ")))
        u += 1


