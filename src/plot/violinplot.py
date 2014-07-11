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

#import interval_filter
import interval_filter
import fasta

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
#TODO: Better way to manage memory
def countInterval(operon,histogram,fastaindex):
    zero = len(histogram)/2
    indexer = fasta.Indexer("",fastaindex)
    indexer.load()
    intervals = []
    
    operon = sorted(operon,key=lambda x:x[5])
    starts,ends = zip(*operon)[5:7]
    
    gene = operon[0]
    st,end = map(int,gene[5:7])
    name = gene[0]
    st,_ = indexer.sixframe_to_nucleotide( name, st )
    end,_ = indexer.sixframe_to_nucleotide( name, end )
    front = min(st,end)
    for gene in operon:
        name = gene[0]
        cluster = gene[1]
        function = cluster.split('.')[0]
        st,end = map(int,gene[5:7])
        st,_ = indexer.sixframe_to_nucleotide( name, st )
        end,_ = indexer.sixframe_to_nucleotide( name, end )
        st = min([st,0])
        end = max([len(histogram)-zero,end])
        for i in xrange(st,end):
            histogram[i+zero]+=1
    
""" Counts all of the positions in the range of
each functional gene in an operon """
def listIntervals(operon,fastaindex):
    intervals = []
    
    operon = interval_filter.overlaps(operon,fastaindex,False)
    operon = sorted(operon,key=lambda x:x[5])
    toxins = [t for t in operon if t[1].split('.')[0]=='toxin']
    starts,ends = zip(*toxins)[5:7]
    gene = operon[0]
    st,end = map(int,gene[5:7])
    front = st
    name = gene[0]
    print "Functions",[gene[1].split('.')[0] for gene in operon ]
    strand = fasta.strand(fasta.getFrame(gene[0]))
    #if strand== "+":
    #    st,end = min(map(int,starts)),min(map(int,ends)) 
    #    front = min(st,end)
    #else:
    #    st,end = max(map(int,starts)),max(map(int,ends)) 
    #    front = max(st,end)
        
    for i in xrange(0,len(operon)):
        gene = operon[i]
        name = gene[0]
        cluster = gene[1]
        function = cluster.split('.')[0]
        st,end = map(int,gene[5:7])
        if st>end: st,end = end,st
        intervals.append( (function,st-front,end-front) )
    ints = []
    for intv in intervals:
        func,st,end = intv
        interval = xrange(st,end)
        funcs = [func]*len(interval)
        ints+=zip(funcs,interval)
    return ints
""" Count the midpoints of genes of all of the functional genes """
def countMidpoints(operon,fastaindex):
    intervals = []
    operon = interval_filter.overlaps(operon,fastaindex,backtran=False)
    operon = sorted(operon,key=lambda x:x[5])
    #toxins = [t for t in operon if t[1].split('.')[0]=='toxin']
    #starts,ends = zip(*operon)[5:7]
    gene = operon[0]
    strand = fasta.strand(fasta.getFrame(gene[0]))
    name = gene[0]
    st,end = map(int,gene[5:7])
    front = st
    #if strand== "+":
    #    st,end = min(map(int,starts)),min(map(int,ends)) 
    #    front = min(st,end)
    #else:
    #    st,end = max(map(int,starts)),max(map(int,ends)) 
    #    front = max(st,end)
    
    #print "Operon",name,"Strand",strand
    for i in xrange(1,len(operon)):
        gene = operon[i]
        name = gene[0]
        cluster = gene[1]
        function = cluster.split('.')[0]
        st,end = map(int,gene[5:7])
        mid = (st+end)/2 - front
        #print function,'st',st,'end',end,'front',front,'mid',mid
        intervals.append( (function,mid) )
    return intervals

""" Distance distribution of functions within operons"""
def operonDistribution(operonFile,fastaindex,clade="all"):
    if not os.path.exists("intervals.pickle"):
        intervals = []
        with open(operonFile,'r') as handle:
            buf = []
            for ln in handle:
                if ln[0]=="-":
                    intervals+=listIntervals(buf,fastaindex)
                    buf = []
                else:
                    ln = ln.rstrip()
                    toks = ln.split('|')
                    acc,clrname,full_evalue,hmm_st,hmm_end,env_st,env_end,description=toks
                    if clade=="all":
                        buf.append((acc,clrname,full_evalue,hmm_st,hmm_end,env_st,env_end, description))
                    else:
                        currentClade = description.split(' ')[0]
                        if clade==currentClad:
                            buf.append((acc,clrname,full_evalue,hmm_st,hmm_end,env_st,env_end, description))
        print len(intervals)
        cPickle.dump(intervals,open("intervals.pickle",'wb'))
    else:
        intervals = cPickle.load(open('intervals.pickle','rb'))
    print intervals[:10]
    toxin = robjects.IntVector([t[1] for t in intervals if t[0]=="toxin"])
    modifier = robjects.IntVector([t[1] for t in intervals if t[0]=="modifier"])
    immunity = robjects.IntVector([t[1] for t in intervals if t[0]=="immunity"])
    transport = robjects.IntVector([t[1] for t in intervals if t[0]=="transport"])
    regulator = robjects.IntVector([t[1] for t in intervals if t[0]=="regulator"])
    importr("ggplot2")
    
    plotViolinFunc = robj.r("""
                            library(ggplot2)
                            function(toxin,modifier,immunity,transport,regulator){
                            #par(font.axis=2, font.lab=2)
                            #png(filename="violin.png",units="in",width=80,height=80,res=2400)
                            #png(filename="violin.png")
                            pdf("violin_zoom.pdf",width=15,heigh=15)                            
                            df <- rbind(data.frame(x=1:length(toxin),y=sort(toxin),group="toxin"),
                                        data.frame(x=1:length(modifier),y=sort(modifier),group="modifier"),
                                        data.frame(x=1:length(immunity),y=sort(immunity),group="immunity"),
                                        data.frame(x=1:length(transport),y=sort(transport),group="transport"),
                                        data.frame(x=1:length(regulator),y=sort(regulator),group="regulator"))
                            p <- ggplot() + 
                                    geom_violin(data=df,
                                                aes(x=group,y=y,group=group,fill=group),
                                                stat="ydensity",
                                                adjust=0.3,
                                                trim=FALSE) +
                                    scale_fill_discrete("Distribution") + 
                                    coord_flip(ylim=c(-2000,50000)) +
                                    ylab("Distance (base pairs)")+
                                    xlab("Function") +
                                    ggtitle("Distribution of Genes By Function in Operons")+
                                    theme(plot.title=element_text(face="bold", size=30)) +
                                    theme(legend.title = element_text(colour="black", size=24, face="bold")) +
                                    theme(legend.text = element_text(colour="black", size = 18)) +
                                    theme(axis.text.x = element_text(size = 18, colour = 'black',face="plain"),
                                          axis.text.y = element_text(size = 18, colour = 'black',face="plain"),
                                          axis.title.x = element_text(size = 28, colour = 'black',face="plain"),
                                          axis.title.y = element_text(size = 28, colour = 'black', angle = 90,face="plain"))
                                   
                            print(p)
                            dev.off()
                            print(p)
                            }
                            """)
    
    plotViolinFunc(toxin,modifier,immunity,transport,regulator)
    raw_input()

                        
    
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
    
  
