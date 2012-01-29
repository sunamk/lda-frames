#!/usr/bin/python
#
# Copyright (C) 2012 Jiri Materna <xmaterna@fi.muni.cz>
# Licensed under the GNU GPLv3 - http://www.gnu.org/licenses/gpl-3.0.html


"""
This script reads grammatical relation file from stdint and prints several statistics.

USAGE: cat ./rel2stats.py RELS

"""

import sys

if len(sys.argv) != 2:
    sys.stderr.write(__doc__)
    sys.exit(1)

f = open(sys.argv[1])
rel_names = f.readline()[1:].split("\t")
print "Grammatical relations (%d): %s" % (len(rel_names), ", ".join(rel_names))

for line in f.xreadlines():
    line = line.strip()
f.close()
