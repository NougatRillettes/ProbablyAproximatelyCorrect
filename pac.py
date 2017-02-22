from copy import copy
from collections import deque
import json
from sys import argv
from itertools import combinations

class Partition:
    def __init__(self):
        self.marked = False
        self.elems = []
        self.code = 0
        self.children = []
        self.parents = []

    def __repr__(self):
        return "(" + repr(self.elems) + " | " + repr(self.code) + " | " + repr(self.marked) +  ")"

def generate_k_partitions(k,n):
    table = {0:Partition()}
    current_gen = [table[0]]
    for _ in range(k):
        new_gen = []
        for x in current_gen:
            for i in range (2*n):
                if 2**i & x.code == 0:
                    code = x.code | (2**i)
                    if not code in table:
                        table[code] = Partition()
                        table[code].elems = copy(x.elems)
                        table[code].elems.append(i)
                        table[code].code = code
                        new_gen.append(table[code])
                    table[code].parents.append(x)
                    x.children.append(table[code])
            current_gen = new_gen
    return table

def code_of(l):
    return sum([2**i for i in l])

def mark_all_greater(code,table):
    stack = deque([table[code]])
    while stack:
        x = stack.pop()
        if not x.marked:
            x.marked = True
            stack.extend(x.parents)

def sample(l,table,k):
    samples = combinations(l,k) if len(l) >= k else [l]
    for s in samples:
        code = code_of(s)
        mark_all_greater(code,table)

def conjonc(table,n):
    res = []
    stack = deque(table[0].children)
    while stack:
        x = stack.popleft()
        stack.extend(x.children)
        if not x.marked:
            tauto = False
            for i in range(n):
                extract1 = bool((4**i)*2 & x.code)
                extract2 = bool((4**i) & x.code)
                tauto = extract2 and extract1
                #print(x.elems,extract1,extract2,tauto)
                if tauto:
                    break
            if not tauto:
                res.append(x.elems)
            sub_stack = deque(x.children)
            while sub_stack:
                y = sub_stack.pop()
                if not y.marked:
                    y.marked = True
                    sub_stack.extend(y.children)
            x.marked = True
    return res

dic = None
with open(argv[1]) as f:
    dic = json.loads(f.read())

indexer = dic['indexer']
influences = dic['influences']
n = len(indexer)
revIndexer = [None]*n
for (k,v) in indexer.items():
    revIndexer[v] = k

def varName(i):
    if i % 2 == 0:
        return revIndexer[i//2]
    else:
        return '!'+revIndexer[i//2]
mask = 0
for _ in range(n):
    mask = 4*mask +1

k = int(argv[2])
for i in range(n):
    table = generate_k_partitions(k,n)
    #print(table)
    for s in influences[2*i]:
        sample(s,table,k)
    print("{}+ : {}".format(revIndexer[i],[[varName(x) for x in c] for c in conjonc(table,n)]))
    table = generate_k_partitions(k,n)
    for s in influences[2*i+1]:
        sample(s,table,k)
    print("{}- : {}".format(revIndexer[i],[[varName(x) for x in c] for c in conjonc(table,n)]))
