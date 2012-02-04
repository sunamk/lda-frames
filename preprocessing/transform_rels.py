#!/usr/bin/python
#
# Copyright (C) 2012 Jiri Materna <xmaterna@fi.muni.cz>
# Licensed under the GNU GPLv3 - http://www.gnu.org/licenses/gpl-3.0.html


"""
This script reads a grammatical relation file and prints several statistics or generates
a new file using some filters.

USAGE: ./transform_rels.py file.rel ACTION
    file.rel              Input file with the grammatical relation data.
    ACTION                Using this parametr you can specify the desired action. 
                          It takes one of the following values [info, transform, 
                          predicate_hist, realisation_hist].

                          info -- print general statistics, transform -- transform the input
                          file using specified filters and print the result to stdout, 
                          predicate_hist -- print histogram of all predicates from the input 
                          file, realisation_hist -- print histogram of all argument 
                          realisations.

    -h, --help            Print this help.
    -i, --ignore_empty    Ignore lines containing empty realisation (- sign) at any position.
    -p, --pred_min p      Minimum frequency of predicates. Default is 1.
    -P, --pred_max P      Maximum frequency of predicates. Default is unlimited.
    -r, --real_min r      Minimum frequency of argument realisations. Default is 1.
    -R, --real_max R      Maximum frequency of argument realisations. Default is unlimited.   
    

"""

import sys
import getopt

pred_min = 1
pred_max = -1
real_min = 1
real_max = -1
ignore_empty = False
pred_blacklist_file = False
real_blacklist_file = False
pred_blacklist = set([])
real_blacklist = set([])

try:
    opts, args = getopt.getopt(sys.argv[1:], "hip:P:r:R:b:", ["help", "ignore_empty=", "pred_min=",
        "pred_max=", "real_min=", "real_max=", "pred_blacklist=", "real_blacklist="])
except getopt.error, msg:
    print msg
    print "for help use --help"
    sys.exit(1)

for o, a in opts:
    if o in ("-h", "--help"):
        print __doc__
        sys.exit(0)
    if o in ("-i", "--ignore_empty"):
        ignore_empty = True
    if o in ("-p", "--pred_min"):
        pred_min = int(a)
    if o in ("-P", "--pred_max"):
        pred_max = int(a)
    if o in ("-r", "--real_min"):
        real_min = int(a)
    if o in ("-R", "--real_max"):
        real_max = int(a)
    if o == "--pred_blacklist":
        pred_blacklist_file = a
    if o == "--real_blacklist":
        real_blacklist_file = a

if len(args) != 2:
    print "Wrong number of parameters."
    print __doc__
    sys.exit(1)

input_file = args[0]
action = args[1]
if not action in ["info", "transform", "predicate_hist", "realisation_hist"]:
    print "Unknown parameter '%s'." % action
    print __doc__
    sys.exit(1)

f = open(input_file)
rel_names = f.readline()[1:].strip().split("\t")

if (pred_blacklist_file):
    b = open(pred_blacklist_file)
    for w in b.xreadlines():
        w = w.strip()
        if w != "": pred_blacklist.add(w)
    b.close()

if (real_blacklist_file):
    b = open(real_blacklist_file)
    for w in b.xreadlines():
        w = w.strip()
        if w != "": real_blacklist.add(w)
    b.close()

if action == "info":
    print "Grammatical relations (%d): %s." % (len(rel_names), ", ".join(rel_names))

predicates = {}
realisations = {}
lines = 0

for line in f.xreadlines():
    line = line.strip()
    items = line.split("\t")
    if len(items) == 0:
        continue
    predicate = items[0]
    predicates[predicate] = predicates.get(predicate, 0) + 1
    rels = items[1:]
    for rel in rels:
        realisations[rel] = realisations.get(rel, 0) + 1
    lines += 1
f.close()

if pred_max == -1:
    pred_max = max(predicates.values())

if real_max == -1:
    real_max =  max(realisations.values())


if action == "info":
    print "Number of lines:  %d" % lines
    print "Number of unique predicates: %d." % len(predicates)
    print "Number of unique grammatical relation realisations: %d." % len(realisations)
    print "---"
    frequent_preds = len([w for (w,i) in predicates.items() if i >= pred_min and \
            i <= pred_max and not w in pred_blacklist])
    print "Number of predicates satisfying the given constraints %d (%d %%)." % \
            (frequent_preds, int(frequent_preds*100.0/len(predicates)))
    frequent_reals = len([w for (w,i) in realisations.items() if i >= real_min and \
            i <= real_max and not w in real_blacklist])
    print "Number of realisations satisfying the given constraints %d (%d %%)." % \
            (frequent_reals, int(frequent_reals*100.0/len(realisations)))

if action == "predicate_hist":
    hist = sorted(predicates.items(), key = lambda k: -k[1])
    for p in hist:
        if p in pred_blacklist:
            continue
        print "%d\t%s" % (p[1], p[0])
    sys.exit(0)

if action == "realisation_hist":
    hist = sorted(realisations.items(), key = lambda k: -k[1])
    for p in hist:
        if p in real_blacklist:
            continue
        print "%d\t%s" % (p[1], p[0])
    sys.exit(0)


f = open(input_file)
rel_names = f.readline()[1:].strip().split("\t")
constrained_lines = 0
for line in f.xreadlines():
    line = line.strip()
    items = line.split("\t")
    if len(items) == 0:
        continue
    predicate = items[0]
    if predicates[predicate] < pred_min or predicates[predicate] > pred_max or \
            predicate in pred_blacklist:
        continue
    rels = items[1:]
    skip = False
    for rel in rels:
        if realisations[rel] < real_min or realisations[rel] > real_max or \
                (rel == "-" and ignore_empty) or rel in real_blacklist:
            skip = True
            break
    if not skip:
        if action == "transform": print line
        constrained_lines += 1
            
if action == "info":
    print "Number of lines satisfying the given constrains: %d (%d %%)." % \
        (constrained_lines, int(constrained_lines*100.0/lines))

f.close()

