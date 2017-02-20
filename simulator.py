import fileinput
import re
import random
import functools
import operator
from collections import Counter
import math
import json
from sys import argv

class Indexer(dict):
    def __missing__(self,key):
        n = len(self)
        self[key] = n
        return n
cre = re.compile(r'([\w\*]+)')
cre2 = re.compile(r'(?P<reac>.*)->(?P<prod>[^(]*)(\((?P<prop>.*)\))?')
cre3 = re.compile(r'(?:(.*)\*)?(.*)')

indexer = Indexer()

f = ""
with open(argv[1]) as foo:
    f = foo.read()

reactions = []

def specList(s):
    items = cre.findall(s)
    res = []
    for i in items:
        finalParse = cre3.match(i).groups()
        res.append((int(finalParse[0] or 1),indexer[finalParse[1]]))
    return res

for l in f.strip().split('\n'):
    reaction = {}
    mdic = cre2.match(l).groupdict()
    reaction['reactants'] = specList(mdic['reac'])
    reaction['products'] = specList(mdic['prod'])
    reaction['propensity'] = float(mdic['prop'] or 1)
    reactions.append(reaction)


state = [0]*len(indexer)

def canFire(l):
    for (coef,spec) in l:
        if state[spec] < coef:
            return False
    return True

def bigL(h,k):
    return int(2*h*((2*len(indexer))**(k+1)+math.log(h))+1)

influences = []
for _ in indexer:
    influences.append(set())
    influences.append(set())

h = int(argv[2])
k = int(argv[3])
for _ in range(bigL(h,k)):
    doable = [r for r in reactions if canFire(r['reactants'])]
    if not doable:
        break
    weights = [r['propensity']*functools.reduce(operator.mul,[state[spec]**c for (c,spec) in r['reactants']],1) for r in doable]
    chosen = random.choices(doable,weights=weights)[0]
    for (coef,spec) in chosen['products']:
        influences[2*spec] |= {tuple([i for (i,x) in enumerate(state) if x > 0])}
    for (coef,spec) in chosen['reactants']:
        influences[2*spec+1] |= {tuple([i for (i,x) in enumerate(state) if x > 0])}
    for (coef,spec) in chosen['products']:
        state[spec] += coef
    for (coef,spec) in chosen['reactants']:
        state[spec] -= coef

result = {'indexer': indexer, 'influences': [list(x) for x in influences]}

print(json.dumps(result))
