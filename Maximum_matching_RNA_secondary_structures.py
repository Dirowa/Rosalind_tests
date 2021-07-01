string = 'UAGUCAAAGGCAACUGUUACGUAAAAUCCUCAGACACUAAAAAAAGAAUUAACUCGAAAUUCUCGCGUGGGACAUUCCGUGCA'

from math import perm
from math import factorial

#count nucleotides
seq={}
for N in string:
    if N not in seq.keys():
        seq[N] =1
    else:
        seq[N] += 1

# calculate max possibilities by using permutations
# perm is calculated with max options , howmany active
AU = perm(max(seq["A"],seq["U"]),min(seq["A"],seq["U"]))
GC = perm(max(seq["G"],seq["C"]),min(seq["G"],seq["C"]))
print(GC*AU)


