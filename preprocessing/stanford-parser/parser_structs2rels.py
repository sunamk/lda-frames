#!/usr/bin/python

"""Transform resulting structures of Stanford parser to a specified relation format.

USAGE: cat parse_structs.txt | parser_structs2rels.py subjobj >  subjobj.rel
"""

import sys
import methods

if len(sys.argv) != 2:
    print __doc__
    sys.exit(1)

method = getattr(methods, sys.argv[1])
progress = 0
for l in sys.stdin.xreadlines():
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
        print a.encode("utf-8")
sys.exit(0)

