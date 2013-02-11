#!/usr/bin/python
#
# Copyright (C) 2013 Jiri Materna <xmaterna@fi.muni.cz>
# Licensed under the GNU GPLv3 - http://www.gnu.org/licenses/gpl-3.0.html


"""Transform resulting structures of Stanford parser to a specified relation format.

USAGE: ./parser_structs2rels.py parse_structs.prs subjobj >  subjobj.rel

"""

import sys
import methods

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
print method(None, desc=True)
    
for l in f.xreadlines():
    progress += 1
    if progress % 10000 == 0:
        sys.stderr.write("Progress at line #%d.\n" % progress)
    l = l.strip()
    try:
        results = method(eval(l))
    except:
        continue
    for r in results:
        a = "\t".join(r)
        print a
sys.exit(0)

