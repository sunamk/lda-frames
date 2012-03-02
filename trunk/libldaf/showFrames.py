#!/usr/bin/python
#
# Copyright (C) 2012 Jiri Materna <xmaterna@fi.muni.cz>
# Licensed under the GNU GPLv3 - http://www.gnu.org/licenses/gpl-3.0.html
#

"""
This script shows f best frames along with r best role realizations for the predicate
given as parametr. If the predicate is not given, the script prints out all supported
predicates.

USAGE: ./showFrames path [predicate]
    path                  Path to database files.
    predicate             Predicate the frames will be shown for.
    -h, --help            Print this help.
    -f, --frames          Number of best frames to be shown.
    -r, --realizations    Number of best role realizations to be shown.
"""

import sys
import getopt
import libldaf

if __name__ == "__main__":
    frames = 3
    realizations = 10
    
    try:
        opts, args = getopt.getopt(sys.argv[1:], "hf:r:", ["help", "frames=", "roles="])
    except getopt.error, msg:
        print msg
        print "for help use --help"
        sys.exit(1)
    for o, a in opts:
        if o in ("-h", "--help"):
            print __doc__
            sys.exit(0)
        if o in ("-f", "--frames"):
            frames = int(a)
        if o in ("-r", "--realizations"):
            realizations = int(a)

    if len(args) < 1 or  len(args) > 2:
        print "Wrong number of parameters."
        print __doc__
        sys.exit(1)

    path  = args[0]
    if not path.endswith("/"): path += "/"
    if len(args) == 1:
        ldaf = libldaf.LDAF(path)
        print "\n".join(ldaf.getPredicates())
        sys.exit(0)

    predicate = args[1]

    ldaf = libldaf.LDAF(path)
    frameDist = ldaf.getFrameDist(predicate, frames)
    if frameDist == None:
        print "No such predicate in the database (%s)." % predicate
        sys.exit(2)

    print "Grammatical relations: " + ", ".join(ldaf.getRelations())
    for fp in frameDist:
        frame = ldaf.getFrame(str(fp[0]))
        print fp, frame
        roles = [ldaf.getRealDist(s, realizations) for s in frame]
        for i in range(realizations):
            for j in range(len(frame)):
                sys.stdout.write(str(roles[j][i]) + "".join([" " for k in range(50-len(str(roles[j][i]))) ]))
            print
        print
    

