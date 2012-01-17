"""Definitions of supported methods for selecting frame realizations"""


def get_id(word):
    """Returns the numeric suffix from the parsed recognized words"""
    return int(word[word.rindex("-")+1:])

def parserstruct(parseResult):
    """Prints just the parser result."""
    res = str(parseResult) + "\n"
    return res

def normalize(word):
    """Normalizes words."""
    if word == "CLAUSE":
        return word
    result = word.lower()
    try:
        float(result.replace(",", "").replace("-", "").replace(" ", ""))
        return "NUMBER"
    except ValueError:
        return result


def subjobj(parseResult):
    """Generates SUBJECT -- OBJECT grammatical relations."""
    result = []
    for sentence in parseResult:
        lemmas = sentence["lemmas"]
        rels = {}
        for dep in sentence["deps"]:
            l2 = ""
            try:
                l2 = lemmas[get_id(dep[2]) - 1]
            except ValueError:
                continue

            if dep[0] == "agent" or dep[0] == "nsubj" or dep[0] == "nsubjpass" or \
                    dep[0] == "xsubj":
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
            elif dep[0] == "dobj" or dep[0] == "xcomp":
                try:
                    rels[dep[1]][1] = l2
                except KeyError:
                    rels[dep[1]] = ["-", "-"]
                    rels[dep[1]][1] = l2
            elif dep[0] == "ccomp":
                try:
                    rels[dep[1]][1] = "CLAUSE"
                except KeyError:
                    rels[dep[1]] = ["-", "-"]
                    rels[dep[1]][1] = "CLAUSE"
        for k, v in rels.iteritems():
            try:
                result.append((normalize(lemmas[get_id(k) - 1]), normalize(v[0]), 
                    normalize(v[1])))
            except ValueError:
                continue
    return result
        
                    
