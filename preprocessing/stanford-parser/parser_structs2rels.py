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
for l in sys.stdin.xreadlines():
    l = l.strip()
    results = method(eval(l))
    for r in results:
        a = "\t".join(r)
        print a.encode("utf-8")
sys.exit(0)

