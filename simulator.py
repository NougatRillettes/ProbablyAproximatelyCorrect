import fileinput
import re
import random
import functools
import operator
from collections import Counter
import math
import json
from sys import argv,stderr

class Indexer(dict):
    def __missing__(self,key):
        n = len(self)
        self[key] = n
        return n
cre = re.compile(r'([\w\*]+)')
cre2 = re.compile(r'(?P<reac>[^/]*)(/(?P<inhib>.*))?->(?P<prod>[^(]*)(\((?P<prop>.*)\))?')
cre3 = re.compile(r'(?:(.*)\*)?(.*)')

indexer = Indexer()

f = ""
with open(argv[1]) as foo:
    f = foo.read()

reactions = []

def specList(s):
    if not s:
        return []
    items = cre.findall(s)
    res = []
    for i in items:
        finalParse = cre3.match(i).groups()
        res.append((int(finalParse[0] or 1),indexer[finalParse[1]]))
    return res

lines = f.strip().split('\n')

for l in lines[1:]:
    reaction = {}
    mdic = cre2.match(l).groupdict()
    reaction['reactants'] = specList(mdic['reac'])
    reaction['products'] = specList(mdic['prod'])
    reaction['inhibitors'] = specList(mdic['inhib'])
    reaction['propensity'] = float(mdic['prop'] or 1)
    reactions.append(reaction)


state = [0]*len(indexer)
dicState = json.loads(lines[0])
for (k,v) in dicState.items():
    state[indexer[k]] = v

def canFire(r):
    for (coef,spec) in r['reactants']:
        if state[spec] < coef:
            return False
    for (_,spec) in r['inhibitors']:
        if state[spec] > 0:
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
loop_end = bigL(h,k)
for loop_n in range(loop_end):
    if loop_n % (2)**14 == 0:
        print(state,file=stderr)
        print(100*loop_n/loop_end,file=stderr)
    doable = [r for r in reactions if canFire(r)]
    if not doable:
        break
    weights = [r['propensity']*functools.reduce(operator.mul,[state[spec]**c for (c,spec) in r['reactants']],1) for r in doable]
    chosen = random.choices(doable,weights=weights)[0]
    tmp = [0]*len(indexer)
    for (coef,spec) in chosen['products']:
        tmp[spec] += coef
    for (coef,spec) in chosen['reactants']:
        tmp[spec] -= coef
    for (spec,aux) in enumerate(tmp):
        if aux > 0 and state[spec] == 0:
        #Warning, outputs reverse of teh influence
            influences[2*spec] |= {tuple([2*i + (1 if x > 0 else 0) for (i,x) in enumerate(state)])}
        elif state[spec] > 0 and -aux >= state[spec]:
            influences[2*spec+1] |= {tuple([2*i + (1 if x > 0 else 0) for (i,x) in enumerate(state)])}
    for (coef,spec) in chosen['products']:
        state[spec] += coef
    for (coef,spec) in chosen['reactants']:
        state[spec] -= coef

result = {'indexer': indexer, 'influences': [list(x) for x in influences]}

print(json.dumps(result))
