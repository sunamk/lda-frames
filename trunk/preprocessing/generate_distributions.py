#!/usr/bin/python
#
# Copyright (C) 2012 Jiri Materna <xmaterna@fi.muni.cz>
# Licensed under the GNU GPLv3 - http://www.gnu.org/licenses/gpl-3.0.html
#

"""
This script takes as input sampled files and generates distributions over frames
and roles.

USAGE: ./generate_distributions.py path
    path                  Path to the data files.
    -h, --help            Print this help.
"""

import sys
import getopt
import cPickle
import shelve
from ldafdict import Dictionary
from itertools import izip


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
    if len(args) != 1:
        print "Wrong number of parameters."
        print __doc__
        sys.exit(1)

    path = args[0]
    if not path.endswith("/"): path += "/"
    dict_file_name = path + "relations.dict"
    data_file_name = path + "relations.dat"
    frames_file_name = path + "frames.smpl"
    roles_file_name = path + "roles.smpl"

    dict_file = None
    data_file = None
    frames_file = None
    roles_file = None
    try:
        dict_file = open(dict_file_name)
        data_file = open(data_file_name)
        frames_file = open(frames_file_name)
        roles_file = open(roles_file_name) 
    except IOError, msg:
        print msg
        sys.exit(1)

    dictionary = cPickle.load(dict_file)
    dict_file.close()

    roles = []
    
    print "Loading roles..."
    for l in roles_file.xreadlines():
        roles.append(l.strip().split(" "))
    roles_file.close()

    
    roles_hist = {}
    print "Generating frame distributions..."
    frame_database = shelve.open(path+"framedist.db", flag='n', 
            protocol=cPickle.HIGHEST_PROTOCOL)

    for u, (data_line, frames_line) in enumerate(izip(data_file.xreadlines(), 
                frames_file.xreadlines())):
        data = data_line.strip().split("\t")
        frames = map(int, frames_line.strip().split(" "))
        s = len(frames)
        frames_hist = {}
        
        for rtuple, frame in izip(data, frames):
            rtuple = map(lambda x: dictionary.id2real(int(x)), rtuple.split(" "))
            frames_hist[frame] = frames_hist.get(frame, 0) + 1
            for role, realisation in izip(roles[frame-1], rtuple):
                if not role in roles_hist:
                    roles_hist[role] = {}
                roles_hist[role][realisation] = roles_hist[role].get(realisation, 0) + 1
        
        dist = []
        for fr,f in frames_hist.iteritems():
            dist.append((fr,f*1.0/s))
        dist = sorted(dist, key=lambda x: -x[1])
        frame_database[dictionary.id2pred(u+1)] = dist


    frame_database.close()    
    data_file.close()
    frames_file.close()

    print "Generating realisation distributions..."
    role_database = shelve.open(path+"realdist.db", flag='n',
             protocol=cPickle.HIGHEST_PROTOCOL)
    for role, hist in roles_hist.iteritems():
        s = 0
        dist = []
        for f in hist.itervalues():
            s+= f
        for w,f in hist.iteritems():
            dist.append((w,f*1.0/s))
        dist = sorted(dist, key=lambda x: -x[1])
        role_database[role] = dist
    role_database.close()

    frames_database = shelve.open(path+"frames.db", flag='n',
            protocol=cPickle.HIGHEST_PROTOCOL)
    for i,f in enumerate(roles):
        frames_database[str(i+1)] = map(str, roles[i])
        frames_database["__description__"] = dictionary.getRelations()
    frames_database.close()
        
        

