#!/usr/bin/python
#
# Copyright (C) 2013 Jiri Materna <xmaterna@fi.muni.cz>
# Licensed under the GNU GPLv3 - http://www.gnu.org/licenses/gpl-3.0.html
#

"""
This script takes as input a grammatical relations file, and generates the input data file
for the sampler (where grammatical realisations are represented by numbers) and a pickle 
with dictionaries, translating realisations and predicates to numbers and vice versa.

USAGE: ./generate_samplerinput.py input.rel path
    input.rel             Input file with the grammatical relation data.
    path                  Output path.
    -h, --help            Print this help.
"""

import sys
import getopt
import cPickle
import random
from ldafdict import Dictionary


if __name__ == "__main__":
    try:
        opts, args = getopt.getopt(sys.argv[1:], "ho:d:", ["help"])
    except getopt.error, msg:
        print msg
        print "for help use --help"
        sys.exit(1)
    for o, a in opts:
        if o in ("-h", "--help"):
            print __doc__
            sys.exit(0)

    if len(args) != 2:
        print "Wrong number of parameters."
        print __doc__
        sys.exit(1)
    input_file = args[0]
    path = args[1]
    if not path.endswith("/"): path += "/"

    output_file = open(path + "relations.dat", "w")
    dict_file = open(path + "relations.dict", "w")

    dictionary = Dictionary()
    data = []

    progress = 0
    f = open(input_file)
    relations = f.readline()[1:].strip().split("\t")
    dictionary.setRelations(relations)
    for line in f.xreadlines():
        progress += 1
        if progress % 100000 == 0:
            print "Progress at line no. #%d." % progress
        line = line.strip()
        if line == "": continue
        items = line.split("\t")
        predicate = items[0]
        reals = items[1:]
        dictionary.addPredicate(predicate)
        pid = dictionary.pred2id(predicate)
        if pid > len(data): data.append([])
        for r in reals: dictionary.addRealisation(r)
        slots = map(lambda r: dictionary.real2id(r), reals)
        if len(slots) != len(relations):
            print "Wrong number of slots in input data: '%s'." % line
            sys.exit(2)
        data[pid-1].append(slots)

    #shuffle lines
    random.seed(None)
    permutation = zip(xrange(len(data)), random.sample(xrange(len(data)), len(data)))
    for swap in permutation:
        tmp = data[swap[0]]
        data[swap[0]] = data[swap[1]]
        data[swap[1]] = tmp
    dictionary.permutatePredicates(permutation)
    print "Writing data."
    for p in data:
        output_file.write("\t".join([" ".join(map(str, r)) for r in p]) + "\n")

    output_file.close()
    cPickle.dump(dictionary, dict_file, cPickle.HIGHEST_PROTOCOL)
    dict_file.close()
    f.close()
