"""
Run validation tests

1) Build database/interval tree from all.ptt
2) Query HMMER entries with interval trees 
   to determine where overlaps are
3) Place HMMER hits into 3 categories 
   ( intergene, hypothetical protein, annotated )
"""
import os,site,sys
from collections import *
from bx.intervals import *
from annotated_genes import AnnotationTree
import re
base_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
for directory_name in os.listdir(base_path):
    site.addsitedir(os.path.join(base_path, directory_name))
    
"""Load predicted operons """
def loadOperons(fname):
    genes=[]
    with open(fname,'r') as handle:
        for ln in handle:
            if ln[0]=='-': continue
            ln = ln.rstrip()
            toks = ln.split('|')
            acc,clrname,full_evalue,hmm_st,hmm_end,env_st,env_end,description=toks
            full_evalue = float(full_evalue)
            hmm_st,hmm_end,env_st,env_end = map(int,[hmm_st,hmm_end,env_st,env_end])
            genes.append((acc,clrname,full_evalue,hmm_st,hmm_end,env_st,env_end,description))
    return genes

""" Group all context genes into the three categories """
def categorize(operon_file,ptt_file,outfile):
    outhandle = open(outfile,'w')
    genes = loadOperons(operon_file)
    tree = AnnotationTree(ptt_file)
    tree.build()
    for gene in genes:
        acc,clrname,full_evalue,hmm_st,hmm_end,env_st,env_end,description=gene
        env_st,env_end = map(int,[env_st,env_end])
        toks = re.split("\s+",description)
        org = ' '.join(toks[:3])
        function = tree.find(org,env_st*3,env_end*3)
        if len(function)>1 or len(function)==0:
            outhandle.write('|'.join(
                            map(str,
                                [acc,clrname,full_evalue,
                                 hmm_st,hmm_end,env_st,env_end,
                                 description,"intergene"])
                            )+"\n")    
        elif "hypothetical" in function:
            outhandle.write('|'.join(
                            map(str,
                                [acc,clrname,full_evalue,
                                 hmm_st,hmm_end,env_st,env_end,
                                 description,"hypothetical"])
                            )+"\n")
        else:
            outhandle.write('|'.join(
                            map(str,
                                [acc,clrname,full_evalue,
                                 hmm_st,hmm_end,env_st,env_end,
                                 description,"annotated"])
                            )+"\n")
   


    
    
    
    
    
    