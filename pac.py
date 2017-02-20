from copy import copy
from collections import deque

class Partition:
    def __init__(self):
        self.marked = False
        self.elems = []
        self.code = 0
        self.children = []
        self.parents = []

    def __repr__(self):
        return "(" + repr(self.elems) + " | " + repr(self.code) + " | " + repr(self.marked) +  ")"

table = {0:Partition()}

def generate_k_partitions(k,n):
    current_gen = [table[0]]
    for _ in range(k):
        new_gen = []
        for x in current_gen:
            for i in range (n):
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

def mark_all_greater(code):
    stack = deque([table[code]])
    while stack:
        x = stack.pop()
        if not x.marked:
            x.marked = True
            stack.extend(x.parents)

def sample(l,n):
    code = 2**(n)-1 ^ code_of(l)
    mark_all_greater(code)

def conjonc():
    res = []
    stack = deque(table[0].children)
    while stack:
        x = stack.popleft()
        stack.extend(x.children)
        if not x.marked:
            res.append(x.elems)
            sub_stack = deque(x.children)
            while sub_stack:
                y = sub_stack.pop()
                if not y.marked:
                    y.marked = True
                    sub_stack.extend(y.children)
            x.marked = True
    return res

generate_k_partitions(3,3)
sample([0,1],3)
sample([0,2],3)
