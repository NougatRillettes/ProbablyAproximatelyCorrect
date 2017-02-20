import fileinput
import re
import random
import functools
import operator
from collections import Counter

class Indexer(dict):
    def __missing__(self,key):
        n = len(self)
        self[key] = n
        return n
cre = re.compile(r'([\w\*]+)')
cre2 = re.compile(r'(?P<reac>.*)->(?P<prod>[^(]*)(\((?P<prop>.*)\))?')
cre3 = re.compile(r'(?:(.*)\*)?(.*)')

indexer = Indexer()

f = fileinput.input()

reactions = []

def specList(s):
    items = cre.findall(s)
    res = []
    for i in items:
        finalParse = cre3.match(i).groups()
        res.append((int(finalParse[0] or 1),indexer[finalParse[1]]))
    return res

for l in f:
    reaction = {}
    mdic = cre2.match(l).groupdict()
    reaction['reactants'] = specList(mdic['reac'])
    reaction['products'] = specList(mdic['prod'])
    reaction['propensity'] = float(mdic['prop'] or 1)
    reactions.append(reaction)

print(indexer)

state = [0]*len(indexer)

def canFire(l):
    for (coef,spec) in l:
        if state[spec] < coef:
            return False
    return True

influences = []
for _ in indexer:
    influences.append(Counter())
    influences.append(Counter())

for _ in range(10000):
    doable = [r for r in reactions if canFire(r['reactants'])]
    if not doable:
        break
    weights = [r['propensity']*functools.reduce(operator.mul,[state[spec]**c for (c,spec) in r['reactants']],1) for r in doable]
    chosen = random.choices(doable,weights=weights)[0]
    for (coef,spec) in chosen['products']:
        influences[2*spec][tuple([i for (i,x) in enumerate(state) if x > 0])] += 1
    for (coef,spec) in chosen['reactants']:
        influences[2*spec+1][tuple([i for (i,x) in enumerate(state) if x > 0])] += 1
    for (coef,spec) in chosen['products']:
        state[spec] += coef
    for (coef,spec) in chosen['reactants']:
        state[spec] -= coef

print(influences)
