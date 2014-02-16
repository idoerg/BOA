

import os
import sys
from collections import *

"""
Bins bacteriocins by reading frame
"""
def overlaps(L):
    overlapSame = 0 #Overlap in same reading frame
    overlapDiff = 0 #Overlap in different reading frame
    for i in range(len(L)):
        for j in range(0,i):
            st1,en1 = L[i]
            st2,en2 = L[j]
            if ((st1<=st2 and st2<=en1) or
                (st1>=st2 and st1<=en2) or
                (st1>=st2 and en1<=en2) or
                (st2>=st1 and en2<=en1)):
                if st1%3==st2%3:
                    overlapSame+=1
                else:
                    overlapDiff+=1
    return overlapSame,overlapDiff

coords = defaultdict( list )# a list of bacteriocin coordinates for each species

for ln in sys.stdin:
    ln = ln.rstrip()
    toks = ln.split('\t')
    organism,st,en = toks[1],int(toks[2]),int(toks[3])
    coords[organism].append( (st,en) )

for k,v in coords.iteritems():
    overs = overlaps(v)
    print k,'\t','\t'.join(map(str,overs))
