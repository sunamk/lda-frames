#!/usr/bin/python
#
# Copyright (C) 2013 Jiri Materna <xmaterna@fi.muni.cz>
# Licensed under the GNU GPLv3 - http://www.gnu.org/licenses/gpl-3.0.html

"""Prepare the output of Stanford Parser for evaluating SemaEval.

USAGE: ./semeval_prepare.py input_file verb > result
"""

import sys
import methods
import re

if len(sys.argv) != 3:
    sys.stderr.write(__doc__)
    sys.exit(1)

method = getattr(methods, sys.argv[2])
progress = 0

f = open(sys.argv[1])
head = f.readline().strip()
if not head == ";ALL PARSER STRUCTURES":
    sys.stderr.write("Input data is not valid.\n")
    sys.exit(1)

predicate = re.sub(r'.*/', '', sys.argv[1])
predicate = re.sub(r'\..*', '', predicate)

for l in f.xreadlines():
    progress += 1
    if progress % 10000 == 0:
        sys.stderr.write("Progress at line #%d.\n" % progress)
    l = l.strip()
    results = method(eval(l))
    p = ''
    for r in results:
        a = "\t".join(r)
        if r[0] == predicate:
            p = a
    if p == '':
        print predicate + "\t-\t-\t-\t-"
    else:
        print p
sys.exit(0)
