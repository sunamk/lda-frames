#!/usr/bin/python
#
# Copyright (C) 2014 Jiri Materna <xmaterna@fi.muni.cz>
# Licensed under the GNU GPLv3 - http://www.gnu.org/licenses/gpl-3.0.html

"""
This script transforms the CLT data format into LDA-Frames input file.

USAGE: ./clt2ldaf input
    input                 Input file.
    -h, --help            Print this help.
"""

import sys
import getopt

def group(w):
    pronouns = set(["I", "you", "he", "she", "it", "we", "you", "they",
                    "me", "him", "her", "us", "them"])
    if w in pronouns:
        return "PRONOUN"
    try:
        float(w)
        return "NUMBER"
    except:
        pass
    return w
       
    

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

    input_file_name = sys.argv[1]
    input_file = None
    try:
        input_file = open(input_file_name)
    except:
        sys.stderr.write("Cannot open input file.\n")

    print ";subject\tfirst_object\tsecond_object\tphrase"
    for line in input_file:
        items = line.strip().split("\t")
        if len(items) == 18:
            items.append("--") 
        if len(items) != 19:
            sys.stderr.write("Wrong number of columns.\n")
            sys.exit(2) 

        #group pronouns and numbers
        items = map(group, items)
        
        indexes = [5, 2, 7, 11] 
        new_items = [items[i] for i in indexes]

        #concatenate phrases with prepositions
        if items[15] != "--" and items[18] != "--":
            new_items.append(items[18] + "_" + items[15])
        else:
            new_items.append(items[15])
        
        print "\t".join(new_items)
