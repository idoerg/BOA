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
            #gg = accTable.lookupGenbank(accession)
            gg = accession
            if gg==None: continue
            if gg in lenCnts: #calculate average length
                num,denom = lenCnts[gg]
                num+=len(seq)
                denom +=1
                lenCnts[gg] = (num,denom)
            else:
                lenCnts[gg] = (len(seq),1)
    lenHandle = open(lengthsOut,'w')
    for accum,frac in lenCnts.items():
        num,denom = frac
        average = num/denom
        lenHandle.write("%s\t%d\n"%(accum,average))
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
            #gg = accTable.lookupGenbank(accession)
            gg = accession
            if gg==None: continue
            bacCnts[gg] +=1
    bacHandle = open(countsOut,'w')
    for accum,count in bacCnts.items():
        bacHandle.write("%s\t%d\n"%(accum,count))
    bacHandle.close()

def anchorgeneCounts(anchorGeneFile,countsOut,accTable):
    anchorCnts = Counter()
    with open(anchorGeneFile,'r') as handle:
        for ln in handle:
            ln = ln.rstrip()
            toks = ln.split('\t')
            refID = toks[1]
            accession = accReg.findall(refID)[0]
            #gg = accTable.lookupGenbank(accession)
            gg = accession
            if gg==None: continue
            anchorCnts[gg] +=1
    anchorHandle = open(countsOut,'w')
    for accum,count in anchorCnts.items():
        anchorHandle.write("%s\t%d\n"%(accum,count))
    anchorHandle.close()
            
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
        '--lengths', type=str, required=False,
        help='Average bacteriocin length for each species')
    parser.add_argument(\
        '--bacteriocin-counts', type=str, required=False,
        help='Number of bacteriocin for each species')    
    parser.add_argument(\
        '--anchor-gene-counts', type=str, required=False,
        help='Number of bacteriocin for each species')    
    args = parser.parse_args()
    accTable = AccessionGG(args.accession_table)#maps accession ids to green gene ids    
    bacteriocinLengths(args.blast_input,args.lengths,accTable)
    bacteriocinCounts(args.blast_input,args.bacteriocin_counts,accTable)
    anchorgeneCounts(args.anchor_genes,args.anchor_gene_counts,accTable)
    
