import os
import sys

for ln in sys.stdin:
    ln = ln.rstrip()
    toks = ln.split('\t')
    organism = toks[1]
    seq = toks[8]
    print ">%s\n%s"%(organism,seq)
