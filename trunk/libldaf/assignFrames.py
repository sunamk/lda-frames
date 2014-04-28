#!/usr/bin/python

#
# Copyright (C) 2014 Jiri Materna <xmaterna@fi.muni.cz>
# Licensed under the GNU GPLv3 - http://www.gnu.org/licenses/gpl-3.0.html
#

"""
This scipt assigns the most likely frames to corpus realizations.

USAGE: ./assign_frames.py path input
    path                  Path to the data files.
    input                 Input file with corpus realizations.
    -h, --help            Print this help.
    -e, --eval            Evaluation file
"""


import sys
import getopt
import libldaf

if __name__ == "__main__":
    evaluate = False
    evalFileName = None
    mapping = {}
    
    try:
        opts, args = getopt.getopt(sys.argv[1:], "he:", ["help", "eval="])
    except getopt.error, msg:
        print msg
        print "for help use --help"
        sys.exit(1)
    for o, a in opts:
        if o in ("-h", "--help"):
            print __doc__
            sys.exit(0)
        if o in ("-e", "--eval"):
            evaluate = True
            evalFileName = a
    
    if len(args) != 2:
        print "Wrong number of parameters."
        print __doc__
        sys.exit(1)

    path  = args[0]
    if not path.endswith("/"): path += "/"
    ldaf = libldaf.LDAF(path, 0.01)

    inputFileName = args[1]
    inputFile = open(inputFileName)
    if evaluate:
        evalFile = open(evalFileName)
        for l in evalFile.xreadlines():
            l = l.strip().split("\t")
            if not l[0] in mapping:
                mapping[l[0]] = {l[1]:set([l[2]])}
            try:
                mapping[l[0]][l[1]].add(l[2])
            except KeyError:
                mapping[l[0]][l[1]] = set([l[2]])
        evalFile.close()

    correct = 0
    wrong = 0
    skipped = 0
    allRealizations = 0
    for l in inputFile.xreadlines():
        l = l.strip().split("\t")
        pred = l[0]
        reals = l[1:len(l)-1]
        frame = l[len(l)-1]
        allRealizations += 1
        if evaluate:
            s = "s"
            try:
                if str(ldaf.getFrameProbs(pred, reals)[0][0]) in mapping[pred][frame]:
                    correct += 1
                    s = "o"
                else:
                    s = "x"
                    wrong += 1
            except KeyError:
                skipped += 1
            print pred, s, reals, frame, ldaf.getFrameProbs(pred, reals)[:4]
                
        else:
            print ldaf.getFrameProbs(pred, reals)[:4]
    
    inputFile.close()

    if evaluate:
        sys.stderr.write("Accuracy: %f\n" %(correct*1.0/allRealizations))
        sys.stderr.write("Precision: %f\n" %(correct*1.0/(correct + wrong)))
        sys.stderr.write("Recall: %f\n" %(correct*1.0/(correct + skipped)))
