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

def normalize(word, tag):
    """Normalizes words."""
    if word == "CLAUSE":
        return word
    word = word.lower()
    if tag=="SYM":
        return "SYMBOL"
    if tag=="NUM":
        return NUMBER
    try:
        float(word.replace(",", "").replace("-", "").replace(" ", ""))
        return "NUMBER"
    except ValueError:
        pass
    return word

def is_verb(tag):
    if tag=="VB" or tag=="VBD" or tag=="VBG" or \
             tag=="VBN" or tag=="VBP" or tag=="VBZ":
        return True
    else:
        return False

def is_adverb(tag):
    if tag=="RB" or tag=="RBR" or tag=="RBS":
        return True
    else:
        return False

def is_noun(tag):
    if tag=="NN" or tag=="NNS" or tag=="NNP" or tag=="NNPS":
        return True
    else:
        return False

def is_pronoun(tag):
    if tag=="PRP":
        return True
    else:
        return False

def is_symbol(tag):
    if tag=="SYM":
        return True
    else:
        return False

def is_number(tag):
    if tag=="NUM":
        return True
    else:
        return False

def is_adjective(tag):
    if tag=="JJ" or tag=="JJR" or tag=="JJS":
        return True
    else:
        return False


def verb(parseResult, desc=False):
    """Generates SUBJECT, ACC_OBJECT, DAT_OBJECT, ADV_MODIFIER grammatical relations."""
    if desc == True:
        return ";SUBJECT\tACC_OBJECT\tDAT_OBJECT\tADV_MODIFIER"

    result = []
    for sentence in parseResult:
        lemmas = None
        deps = None
        tags = None
        try:
            lemmas = sentence["lemmas"]
            deps = sentence["deps"]
            tags = sentence["tags"]
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
                tag = tags[get_id(dep[2]) - 1]
                l2 = normalize(lemmas[get_id(dep[2]) - 1], tag)
            except ValueError:
                continue
            except IndexError:
                continue

            if  dep[0] == "nsubj" and (is_noun(tag) or is_number(tag) or \
                    is_pronoun(tag) or is_symbol(tag)):
                try:
                    rels[dep[1]][0] = l2
                except KeyError:
                    rels[dep[1]] = ["-", "-", "-", "-"]
                    rels[dep[1]][0] = l2
            
            elif dep[0] == "csubj" or dep[0] == "csubjpass":
                try:
                    rels[dep[1]][0] = "CLAUSE"
                except KeyError:
                    rels[dep[1]] = ["-", "-", "-", "-"]
                    rels[dep[1]][0] = "CLAUSE"

            elif dep[0] == "dobj" and (is_noun(tag) or is_number(tag) or \
                    is_pronoun(tag) or is_symbol(tag)):
                try:
                    rels[dep[1]][1] = l2
                except KeyError:
                    rels[dep[1]] = ["-", "-", "-", "-"]
                    rels[dep[1]][1] = l2
            elif dep[0] == "iobj" and (is_noun(tag) or is_number(tag) or \
                    is_pronoun(tag) or is_symbol(tag)):
                try:
                    rels[dep[1]][2] = l2
                except KeyError:
                    rels[dep[1]] = ["-", "-", "-", "-"]
                    rels[dep[1]][2] = l2
            elif dep[0] == "advmod" and is_adverb(tag):
                try:
                    rels[dep[1]][3] = l2
                except KeyError:
                    rels[dep[1]] = ["-", "-", "-", "-"]
                    rels[dep[1]][3] = l2

        for k, v in rels.iteritems():
            try:
                entry = (normalize(lemmas[get_id(k) - 1], tags[get_id(k) - 1]), 
                    v[0], v[1],v[2], v[3])
                if "" in entry or not is_verb(tags[get_id(k) - 1]):
                    continue
                result.append(entry)
            except ValueError:
                continue
            except IndexError:
                continue
    return result
        
                    
