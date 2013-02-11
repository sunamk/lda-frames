#!/usr/bin/python
#
# Copyright (C) 2013 Jiri Materna <xmaterna@fi.muni.cz>
# Licensed under the GNU GPLv3 - http://www.gnu.org/licenses/gpl-3.0.html
#

"""
This script evaluates quality of frames infered for data generated using 'randomData.py'.

USAGE: ./evaluateRandom.py original_filename path
"""

import sys
from itertools import izip

class Evaluator:

    def parseOriginal(self, ofile):
        original = []
        for line in open(ofile).xreadlines():
            predicate = []
            for frame in line.strip().split("\t"):
                predicate.append(frame.split(" "))
            original.append(predicate)
        return original

    def parseFrames(self, ffile):
        frames = []
        for line in open(ffile).xreadlines():
            frames.append(map(int, line.strip().split(" ")))
        return frames
    
    def parseRoles(self, rfile):
        roles = []
        for line in open(rfile).xreadlines():
            roles.append(map(int, line.strip().split(" ")))
        return roles

    def bestMapping(self, original, frames, roles): 
        mapping = {}
        freqs = {}
        for orig, fr in izip(original, frames):
            for frame, frame_id in izip(orig, fr):
                for role, role_id in izip(frame, roles[frame_id-1]):
                    if not role_id in freqs:
                        freqs[role_id] = {}
                    freqs[role_id][role] = freqs[role_id].get(role, 0) + 1
        for role_id, mapp in freqs.iteritems():
            role = None
            freq = 0
            for r, f in mapp.iteritems():
                if f >= freq:
                    role = r
                    freq = f
            mapping[role_id] = role
        return mapping

    def computeScore(self, original, frames, roles, mapping):
        all_freq = 0
        correct_freq = 0
        for pred, (orig, fr) in enumerate(izip(original, frames)):
            for pos, (frame, frame_id) in enumerate(izip(orig, fr)):
                for slot, (role, role_id) in  enumerate(izip(frame, roles[frame_id-1])):
                    all_freq += 1
                    if role == mapping[role_id]:
                        correct_freq += 1
                    else:
                        sys.stderr.write("Incorrectly classified %d as %s at "
                                          "predicate %d, position %d, slot %d.\n" % (
                            role_id, role, pred+1, pos+1, slot+1))
        sys.stderr.write("Correctly assigned sematic roles: %d/%d.\n" % (correct_freq, all_freq))
        return float(correct_freq)/all_freq        

if __name__ == "__main__":
        
    if len(sys.argv) != 3:
        sys.stderr.write(
"""Incorrect number of parameters.

USAGE: ./evaluateRandom.py original_filename path
""")
        sys.exit(1)
    orig_file_name = sys.argv[1]
    path = sys.argv[2]
    if not path.endswith("/"): path += "/"
        
    evaluator = Evaluator()
    original = evaluator.parseOriginal(orig_file_name)
    frames = evaluator.parseFrames(path + "frames.smpl")
    roles = evaluator.parseRoles(path + "roles.smpl")
    mapping = evaluator.bestMapping(original, frames, roles)
    print evaluator.computeScore(original, frames, roles, mapping)
