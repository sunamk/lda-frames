#!/usr/bin/env python
#
# Copyright (C) 2014 Jiri Materna <xmaterna@fi.muni.cz>
# Licensed under the GNU GPLv3 - http://www.gnu.org/licenses/gpl-3.0.html
#

"""
This script exports LDA-Frames data to text format.

USAGE: ./exportData path [output_path]
    path                  Path to database files.
    output_path           Output path. If not specified, the files are stored in the input directory.
    -h, --help            Print this help.
"""

import sys
import getopt
import libldaf

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

    if len(args) < 1 or  len(args) > 2:
        print "Wrong number of parameters."
        print __doc__
        sys.exit(1)

    path  = args[0]
    output_path = path
    if not path.endswith("/"): path += "/"
    if not output_path.endswith("/"): output_path += "/"
    if len(args) == 2:
        output_path = args[1]
        if not output_path.endswith("/"): output_path += "/"

    ldaf = libldaf.LDAF(path)
    frame_set = set([])
    role_set = set([])

    verbs_f = open(output_path + "verbs.txt", "w")
    for v in ldaf.getPredicates():
        for f in ldaf.getFrameDist(v):
            verbs_f.write("%s\t%d\t%f\n" % (v, f[0], f[1]))
            frame_set.add(f[0])
    verbs_f.close()

    
    frames_f = open(output_path + "frames.txt", "w")
    for f in frame_set:
        fr = ldaf.getFrame(str(f))
        frames_f.write(str(f) + "\t" + "\t".join(fr) + "\n")
        for r in fr:
            role_set.add(int(r))
    frames_f.close()

    roles_f = open(output_path + "roles.txt", "w")
    for r in role_set:
        for real in ldaf.getRealDist(str(r)):
            roles_f.write("%d\t%s\t%f\n" % (r, real[0], real[1]))
    roles_f.close()


