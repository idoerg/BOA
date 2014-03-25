"""
Takes blast output and reformats it to prepare for iTOL
"""

import os,sys
import argparse
import re
from collections import Counter
from accessionMap import AccessionGG



accReg = re.compile("([A-Z]+_\d+.\d)")
def bacteriocinLengths(blastFile,lengthsOut,accTable):
    lenCnts = Counter()
    accTable = AccessionGG(args.accession_table)#maps accession ids to green gene ids    
    with open(blastFile,'r') as handle:
        for ln in handle:
            ln = ln.rstrip()
            toks = ln.split('\t')
            bacID,refID,st,end,start,gtype,seq = toks
            st,end = int(st),int(end)
            seq = seq.replace('-','')
            accession = accReg.findall(refID)[0]
            gg = accTable.lookupGenbank(accession)
            if gg==None: continue
            #Print out average lengths
            if gg in lenCnts: #calculate average length
                num,denom = lenCnts[gg]
                num+=len(seq)
                denom +=1
                lenCnts[gg] = (num,denom)
    lenHandle = open(lengthsOut,'w')
    for accum,frac in lenCnts.items():
        num,denom = frac
        average = num/denom
        lenHandle.write("%d\t%d\n"%(int(accum),average))
    lenHandle.close()

def bacteriocinCounts(blastFile,countsOut,accTable):
    bacCnts = Counter()
    with open(blastFile,'r') as handle:
        for ln in handle:
            ln = ln.rstrip()
            toks = ln.split('\t')
            bacID,refID,st,end,start,gtype,seq = toks
            st,end = int(st),int(end)
            seq = seq.replace('-','')
            accession = accReg.findall(refID)[0]
            gg = accTable.lookupGenbank(accession)
            if gg==None: continue
            #Print out average lengths
            if gg in bacCnts: #calculate average length
                bacCnts[gg] +=1
    bacHandle = open(countsOut,'w')
    for accum,count in bacCnts.items():
        bacHandle.write("%d\t%d\n"%(int(accum),count))
    bacHandle.close()

def anchorgeneCounts(anchorGeneFile,countsOut,accTable):
    pass

if __name__=="__main__":
    parser = argparse.ArgumentParser(description=\
        'Format Blast output for iTOL visualization ')
    parser.add_argument(\
        '--accession-table', type=str, required=True,
        help='A table for converting accession IDs into green gene ids')
    parser.add_argument(\
        '--blast-input', type=str, required=True,
        help='Blast tab output from bacteriocin.py')
    parser.add_argument(\
        '--anchor-genes', type=str, required=False,
        help='Anchor genes from genbank files')
    parser.add_argument(\
        '--lengths', type=str, required=True,
        help='Average bacteriocin length for each species')
    parser.add_argument(\
        '--counts', type=str, required=True,
        help='Number of bacteriocin for each species')    
    args = parser.parse_args()
    accTable = AccessionGG(args.accession_table)#maps accession ids to green gene ids    
    bacteriocinLengths(args.blast_input,args.lengths,accTable)
    bacteriocinCounts(args.blast_input,args.counts,accTable)
    
