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
    
            
import cdhit
import fasttree
from fasttree import FastTree
from fasttree import UnAlignedFastTree
from  acc2species import AccessionToSpecies
from accessionMap import GGAccession
import hmmer

import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.cbook as cbook
from matplotlib._png import read_png
from matplotlib.offsetbox import OffsetImage 
import matplotlib.gridspec as gridspec
import matplotlib.image as mpimg

import numpy
import pylab
import argparse
import cPickle
from Bio import SeqIO
from Bio import Phylo

from collections import defaultdict
from collections import Counter
from pandas import *
from clique_filter import *

if __name__=="__main__":
    faidx = "/home/mortonjt/Projects/Bacfinder/workspace/quorum/data/all_trans.fai"
    folder = "/home/mortonjt/Projects/Bacfinder/workspace/quorum/intermediate"
    if os.path.exists("test.pickle"):
        all_hits = cPickle.load(open("test.pickle",'rb'))
    else:
        
        toxin_hits     = hmmer.parse("%s/toxin.out"%folder)
        modifier_hits  = hmmer.parse("%s/modifier.out"%folder)
        immunity_hits  = hmmer.parse("%s/immunity.out"%folder)
        regulator_hits = hmmer.parse("%s/regulator.out"%folder)
        transport_hits = hmmer.parse("%s/transport.out"%folder)
        all_hits = toxin_hits+modifier_hits+immunity_hits+regulator_hits+transport_hits 
        
        #Sort by start/end position
        all_hits=sorted(all_hits,key=lambda x: x[6])   
        all_hits=sorted(all_hits,key=lambda x: x[5])
        
        #Sort by genome name
        all_hits=sorted(all_hits,key=lambda x: x[-1])  
        cPickle.dump(all_hits,open("test.pickle",'wb'))
        print "Sorted"   
    print "All hits",len(all_hits)
    all_hits = collapseOverlaps(all_hits,faidx)
    clusters = findContextGeneClusters(all_hits,faidx)
    outhandle = open('%s/operons.txt'%(folder),'w')
    for cluster in clusters:
        for gene in cluster:
            outhandle.write("%s\n"%gene)
        outhandle.write('----------\n')
    
    
    
    
    
    
    
    