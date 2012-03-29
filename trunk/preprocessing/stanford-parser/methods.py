#
# Copyright (C) 2012 Jiri Materna <xmaterna@fi.muni.cz>
# Licensed under the GNU GPLv3 - http://www.gnu.org/licenses/gpl-3.0.html


"""Definitions of supported methods for selecting frame realizations"""
import sys

def get_id(word):
    """Returns the numeric suffix from the parsed recognized words"""
    return int(word[word.rindex("-")+1:])

def parserstruct(parseResult, desc=False):
    """Prints just the parser result."""

    if desc == True:
        return ";ALL PARSER STRUCTURES"

    res = str(parseResult) + "\n"
    return res

def normalize(word):
    """Normalizes words."""
    if word == "CLAUSE":
        return word
    word = word.lower()
    if word == "i" or word == "we" or word == "you" or word == "he" or word == "she" or\
            word == "me" or word == "him" or word == "her" or word == "us":
        return "person"
    if word == "they" or word == "it" or word == "them" or word =="this" or word=="those" \
            or word == "these":
        return "PRONOUN"
    try:
        float(word.replace(",", "").replace("-", "").replace(" ", ""))
        return "NUMBER"
    except ValueError:
        pass
    return word

def subjobj(parseResult, desc=False):
    """Generates SUBJECT -- OBJECT grammatical relations."""
    if desc == True:
        return ";SUBJECT\tOBJECT"

    result = []
    for sentence in parseResult:
        lemmas = None
        deps = None
        try:
            lemmas = sentence["lemmas"]
            deps = sentence["deps"]
        except TypeError:
            sys.stderr.write("Bad input. Ignoring line.\n")
            continue
        except KeyError:
            sys.stderr.write("Bad input. Ignoring line.\n")
            continue
        rels = {}
        for dep in deps:
            l2 = ""
            try:
                l2 = lemmas[get_id(dep[2]) - 1]
            except ValueError:
                continue
            except IndexError:
                continue

            if  dep[0] == "nsubj" or dep[0] == "nsubjpass":
                try:
                    rels[dep[1]][0] = l2
                except KeyError:
                    rels[dep[1]] = ["-", "-"]
                    rels[dep[1]][0] = l2
            
            elif dep[0] == "csubj" or dep[0] == "csubjpass":
                try:
                    rels[dep[1]][0] = "CLAUSE"
                except KeyError:
                    rels[dep[1]] = ["-", "-"]
                    rels[dep[1]][0] = "CLAUSE"
            elif dep[0] == "dobj":
                try:
                    rels[dep[1]][1] = l2
                except KeyError:
                    rels[dep[1]] = ["-", "-"]
                    rels[dep[1]][1] = l2
        for k, v in rels.iteritems():
            try:
                entry = (normalize(lemmas[get_id(k) - 1]), normalize(v[0]), 
                    normalize(v[1]))
                if "" in entry:
                    continue
                result.append(entry)
            except ValueError:
                continue
            except IndexError:
                continue
    return result
        
                    
