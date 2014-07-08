"""
Plot the phylogenetic trees for toxins and immunity genes
"""

import os
import sys
import site
import re
import numpy as np
import numpy.random

base_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
base_path="%s/src"%base_path
for directory_name in os.listdir(base_path):
    site.addsitedir(os.path.join(base_path, directory_name))
    
import fasttree
import muscle
import fasta
from fasttree import FastTree
from fasttree import UnAlignedFastTree
from  acc2species import AccessionToSpecies
from accessionMap import GGAccession
"""
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.cbook as cbook
from matplotlib._png import read_png
from matplotlib.offsetbox import OffsetImage 
import matplotlib.gridspec as gridspec
import matplotlib.image as mpimg
"""
import numpy
import pylab
import argparse
import cPickle
from Bio import SeqIO
from Bio import Phylo

base_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
for directory_name in os.listdir(base_path):
    site.addsitedir(os.path.join(base_path, directory_name))
import quorum

from collections import defaultdict
from collections import Counter
from pandas import *
from itol import *
""" Transforms operon txt file into fasta file """
def getFasta(txt,fastadb,fastaindex,fastaout):
    outhandle = open(fastaout,'w')
    indexer = fasta.Indexer(fastadb,fastaindex)
    indexer.load()
    i = 0
    with open(txt,'r') as handle:
        for ln in handle:
            ln = ln.rstrip()
            toks = ln.split('|')
            acc,clrname,full_evalue,hmm_st,hmm_end,env_st,env_end,description=toks
            full_evalue = float(full_evalue)
            hmm_st,hmm_end,env_st,env_end = map(int,[hmm_st,hmm_end,env_st,env_end])
            seq = indexer.fetch(acc,env_st,env_end)
            seq = fasta.format(seq)
            outhandle.write(">%s:%d %s\n%s\n"%(acc,i,description,seq))
            i+=1
if __name__=="__main__":
    db     = "/home/mortonjt/Projects/Bacfinder/db"
    quorum = "/home/mortonjt/Projects/Bacfinder/workspace/quorum"
    operons= "%s/intermediate/operons.txt"%quorum
    rrnaFile = "%s/rrna.fa"%db
    itol = iTOL(operons,rrnaFile)
    itol.getRRNAs()  #### Note: This will only get the RRNAs for chromosomal bacteriocins
    itol.buildTree()

    #toxins = quorum%"/intermediate/toxin_genes.txt"
    #immunity = quorum%"/intermediate/immunity_genes.txt"
    #allfna = quorum%"/data/all_trans.fna"
    #allfai = quorum%"/data/all_trans.fai"
    #toxinfa = "toxin.fa"
    #immunityfa = "immunity.fa"
    #toxin_treefile = "toxin.tree"
    #immunity_treefile = "immunity.tree"
    #getFasta(toxins,allfna,allfai,toxinfa)
    #getFasta(immunity,allfna,allfai,immunityfa)
    #if not os.path.exists(toxin_treefile):
    #    toxinfasttree = UnAlignedFastTree(toxinfa,toxin_treefile)
    #    toxinfasttree.align()
    #    toxinfasttree.run()
    #if not os.path.exists(immunity_treefile):
    #    immunityfasttree = UnAlignedFastTree(immunityfa,immunity_treefile)
    #    immunityfasttree.align()
    #    immunityfasttree.run()
    #    
    #fig = plt.figure(1)
    #plt.suptitle("Toxin phylogeny",fontsize=22)
    #tree=Phylo.read(toxin_treefile,"newick")
    #Phylo.draw(tree,show_confidence=False,do_show=False)    
    
    #fig = plt.figure(2)
    #plt.suptitle("Immunity phylogeny",fontsize=22)
    #tree=Phylo.read(immunity_treefile,"newick")
    #Phylo.draw(tree,show_confidence=False,do_show=False)
    
    pass



