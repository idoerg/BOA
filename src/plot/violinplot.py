import numpy as np
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.cbook as cbook
from matplotlib._png import read_png
from matplotlib.offsetbox import OffsetImage 
import matplotlib.gridspec as gridspec
import matplotlib.image as mpimg
from Bio import Phylo

import os
import sys
import site
import re
import numpy as np
import numpy.random


base_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
for directory_name in os.listdir(base_path):
    site.addsitedir(os.path.join(base_path, directory_name))


import cdhit
import fasttree
import hypertree
from fasttree import FastTree
from fasttree import UnAlignedFastTree
from  acc2species import AccessionToSpecies
from accessionMap import GGAccession

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
import pandas as pd

import rpy2.robjects as robj
#import rpy2.robjects.pandas2ri # for dataframe conversion
import pandas.rpy.common as com
from rpy2.robjects.packages import importr
import rpy2.robjects.lib.ggplot2 as ggplot2
import rpy2.robjects as robjects



#robjects.r("library(vioplot)")

def device(filename):
    robjects.r.pdf(file=filename,width=5,height=5)

def close():
    robjects.r("dev.off()")

cluster_id_reg = re.compile(">(\S+)")
"""Checks the interval overlaps"""
def overlap(bst,bend,ast,aend):
    if bst > ast and bst<aend:
        return True
    elif ast > bst and ast<bend:
        return True
    elif ast > bst and aend < bend:
        return True
    elif bst > ast and bend < aend:
        return True
    else:
        return False

"""Distance distribution of each cluster """
def contextGeneDistances(cdhitProc):
    clusterIDs = []
    dists = []
    #heats = defaultdict( Counter )
    for i,cluster in enumerate(cdhitProc.clusters):
        members = cluster.seqs
        for mem in members:
            tag = cluster_id_reg.findall(mem)[0][:-3]
            bst,bend,ast,aend = tag.split('|')[-4:]
            bst,bend,ast,aend = int(bst),int(bend),int(ast),int(aend)
            if overlap(bst,bend,ast,aend): continue                
            bmid = (bst+bend)/2
            int_st,int_end = ast-bmid,aend-bmid
            interval = xrange(int_st,int_end)
            #heats[i][(int_st+int_end)/2]+=1
            dists+=interval
            clusterIDs+=[i]*len(interval)
            #dists+=[(int_st+int_end)/2]
            #clusterIDs+=[i]
            #dists+= [(int_st+int_end)/2]
            #clusterIDs+=[i]
            
    data = zip(clusterIDs,dists)
    sorted(data,key=lambda x:x[0])
    clusterIDs,dists = zip(*data)
    #heats = pd.DataFrame(heats)
    heats = pd.DataFrame({'clusters':clusterIDs,'distances':dists})
    
    print heats 
    #heats = heats.fillna(0)
    heats_R=com.convert_to_r_dataframe(heats)
    #print heats_R
    importr("ggplot2")
    
    plotViolinFunc = robj.r("""
                            library(ggplot2)
                            function(df){
                            png(filename="violin.png")
                            
                            p <- ggplot(df,aes(x=as.character(clusters),
                                               y=distances)) + 
                                    geom_violin(aes(x=as.character(clusters),
                                               y=distances),
                                               stat="ydensity",
                                               adjust=40,
                                               trim=TRUE,
                                               fill="red") + 
                                    coord_flip()
                            print(p)
                            dev.off()
                            print(p)
                            }
                            """)
    
    plotViolinFunc(heats_R)
    raw_input()
    
  
