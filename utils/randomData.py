#!/usr/bin/python
#
# Copyright (C) 2013 Jiri Materna <xmaterna@fi.muni.cz>
# Licensed under the GNU GPLv3 - http://www.gnu.org/licenses/gpl-3.0.html
#

"""
This script generates random trainnig data using a predefined set of semantic 
roles and frames.

USAGE: ./randomData.py filename [seed]
"""

import sys
import random

ROLES = {"Animal": ["dog", "cat", "mouse", "fish", "chicken", "jenny"],
         "Food": ["cake", "fish", "lunch", "chicken", "dinner"],
         "Institution": ["school", "state", "university", "police"],
         "Person": ["people", "man", "woman", "john", "john", "jenny"],
         "Vehicle": ["car", "bus", "motorcycle", "truck", "tractor"],
         "-": ["-"]
        }

PREDICATES = {"boil": [("Person", "Food")],
              "buy": [("Person", "Food"),("Person", "Animal")],
              "catch": [("Person", "Animal")],
              "crash": [("Person", "Vehicle")],
              "cook": [("Person", "Food"), ("Institution", "Food")],
              "drive": [("Person", "Vehicle")],
              "eat": [("Person", "Food"), ("Animal", "Food")],
              "kill":[("Person", "Animal"), ("Person", "Person")],
              "love": [("Person", "Person")],
              "paint": [("Person", "-"), ("Person", "Person")],
              "pay": [("Institution", "Person")],
              "produce": [("Institution", "Food")],
              "repair": [("Person", "Vehicle")],
              "sell": [("Institution", "Food"), ("Person", "Vehicle")],
              "sue": [("Person", "Person"), ("Person", "Institution")],
              "teach": [("Person", "Person"), ("Institution", "Person")],
              "warn":[("Person", "Person"), ("Institution", "Person")],
              "work":[("Person", "-")],
              "work_for": [("Person", "Institution"), ("Person", "Person")]
             }

seed = None
if len(sys.argv) > 3 or len(sys.argv) < 2:
    sys.stderr.write(
"""Incorrect number of parameters.

USAGE: ./randomData.py filename [seed]
""")
    sys.exit(1)
filename = sys.argv[1]
if len(sys.argv) == 3:
    seed = int(sys.argv[2])
random.seed(seed)

fo = open(filename, "w")
ft = open(filename + ".chck", "w")

fo.write(";SUBJECT\tOBJECT\n")
for pred, frames in PREDICATES.iteritems():
    check = []
    for i in xrange(1, random.randint(10, 30)):
        frame = random.choice(frames)
        subj = random.choice(ROLES[frame[0]])
        obj = random.choice(ROLES[frame[1]])
        fo.write("%s\t%s\t%s\n" % (pred, subj, obj))
        check.append(frame[0] + " " + frame[1])
    ft.write("\t".join(check) + "\n")
        
fo.close()
ft.close()
