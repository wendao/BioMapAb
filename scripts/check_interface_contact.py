import sys
import numpy as np
from collections import defaultdict

contact_cutoff = 8.0

chain_set = defaultdict(list)

for f in sys.argv[1:]:
    lines = open(f, 'r').readlines()
    name = f.split('.')[0]
    id = name.split('_')[-1]
    for l in lines:
        if l[:6] == "ENDMDL": break
        if l[:4] != "ATOM": continue
        if l[13:15] != "CB": continue
        CH = l[21:22]
        x = float(l[30:38])
        y = float(l[38:46])
        z = float(l[46:54])
        chain_set[id].append(np.array([x,y,z]))

def count_contact( i, j ):
    sum_contact = 0
    for ri in chain_set[i]:
        for rj in chain_set[j]:
            d = np.linalg.norm(ri-rj)
            if d < contact_cutoff:
                sum_contact += 1
                if d < 0.0001:
                    return 0
    return sum_contact

chains = chain_set.keys()
#chains.sort()

nc = len(chains)
for ci in xrange(nc-1):
    for cj in xrange(ci+1, nc):
        ch_i = chains[ci]
        ch_j = chains[cj]
        counts = count_contact( ch_i, ch_j )
        if counts > 4:
            print ch_i, ch_j, counts
