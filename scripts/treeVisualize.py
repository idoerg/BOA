import re
import os
import sys
import site
from collections import Counter

accReg = re.compile("[A-Z]+_(\d+)")

cnts = Counter()

for ln in sys.stdin:
    ln = ln.rstrip()
    toks = ln.split('\t')
    assert len(toks)==7
    bacID,accession,st,end,strand,gtype,seq = toks
    st,end = int(st), int(end)
    seq = seq.replace('-','')
    accNum = accReg.findall(accession)[0]
    if accNum in cnts: #calculate average length
        num,denom = cnts[accNum]
        num+=len(seq)
        denom +=1
        cnts[accNum] = (num,denom)
    else:
        cnts[accNum] = ( len(seq), 1)

for accum,frac in cnts.items():
    num,denom = frac
    average = num/denom
    print "%d\t%d"%(int(accum),average)


