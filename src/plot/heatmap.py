###./imp.py deletedOrg /home/asmariyaz/Desktop/phylo_order /home/asmariyaz/Desktop/txtnames
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
from pandas import *


def produce_tree():
    tree=Phylo.read("/home/asmariyaz/Desktop/reorder.nwk","newick")
    Phylo.draw(tree, axes=phyl_ax)
    tree_f=plt.gcf()
    print type(tree_f)
    return tree_f

""" Lookup species names in tree """
def lookup_by_names(tree):
    names = {}
    for clade in tree.get_terminals():
        if clade.name:
            if clade.name in names:
                raise ValueError("Duplicate key: %s" % clade.name)
            names[clade.name] = clade
    return names
        
""" Remove all leaves that aren't contained in names """
def prune_tree(tree,target_species):
    leafs = tree.get_terminals()
    all_species = set([l.code for l in leafs])
    target_species = set(target_species)
    outlier_species = all_species.difference(target_species)
    name_dict = lookup_by_names(tree)
    for outlier in outlier_species:
        tree.prune(name_dict[outlier])
    return tree

""" Makes a heatmap with a phylotree for axes """
def heatMapTree(heats,treeFile,title,xlabel,ylabel,showX=False,showY=False):
    scoremap = DataFrame(heats).T.fillna(0)
    xlabels = list(scoremap.columns)
    ylabels = list(scoremap.index)
    #if showX: pylab.xticks(xpos, xlabels,rotation=90)
    #if showY: pylab.yticks(ypos, ylabels)
    fig = plt.figure()
    plt.suptitle(title,fontsize=22)
    normalizeTree(treeFile)
    tree=Phylo.read(treeFile,"newick")
    tree = prune_tree(tree,ylabels)
    tree.ladderize()
    scoremap = reorderHeatMap(scoremap,tree)
    
    phyloLabels = zip(map(str,tree.get_terminals()),scoremap.index)
    gs=gridspec.GridSpec(1, 2,width_ratios=[1,1,-2,2],
                         wspace=0)
    
    phyl_ax=plt.subplot(gs[0])
    phyl_ax.set_ylabel(ylabel)
    Phylo.draw(tree,axes=phyl_ax,show_confidence=False,do_show=False)
    plt.rcParams['font.size']=8
    xlabels = list(scoremap.columns)
    ylabels = list(scoremap.index)
    xpos = np.arange(len(xlabels))+1
    xminorpos = np.arange(len(xlabels))+0.5
    ypos = np.arange(len(ylabels)+1)
    #plt.grid()
    
    ht_ax=plt.subplot(gs[1])
    ht_ax.set_xlabel(xlabel)
    ht_ax.set_xlim(0,len(xlabels))
    ht_ax.set_ylim(0,len(ylabels))
    plt.setp(phyl_ax.get_xticklabels(),visible=False)
    plt.setp(phyl_ax.get_yticklabels(),visible=False)
    plt.setp(ht_ax.get_xticklabels(),visible=True)
    plt.setp(ht_ax.get_yticklabels(),visible=False)
    plt.setp(phyl_ax.get_xticklines(),visible=True)
    plt.setp(phyl_ax.get_yticklines(),visible=True)
    plt.setp(ht_ax.get_xticklines(),visible=True)
    plt.setp(ht_ax.get_yticklines(),visible=True)
    ht_ax.xaxis.set_ticks(xminorpos)
    ht_ax.yaxis.set_ticks(ypos)
    ht_ax.set_xticklabels(xlabels,rotation=45,fontsize=10)
    ht_ax.set_yticklabels(ylabels,fontsize=10,alpha=1.0)
    ht_ax.xaxis.set_tick_params(pad=4)
    ht_ax.yaxis.set_tick_params(pad=4)
    
    pylab.xticks(xpos, xlabels)
    pylab.yticks(ypos, ylabels)
    for tick in ht_ax.xaxis.get_major_ticks():
        tick.label1.set_horizontalalignment('right')
    
        
    heatmap = plt.pcolor(scoremap,norm=LogNorm())
    plt.grid(True)
    plt.colorbar()

"""Heatmap of bacteriocins versus species"""
def bacteriocinSpeciesHeatmap(cdhitProc,tree,accTable):
    genomeHeats = defaultdict( Counter )
    plasmidHeats = defaultdict( Counter )
    for cluster in cdhitProc.clusters:
        members = cluster.seqs
        for mem in members:
            bacID=bac_reg.findall(mem)[0]
            accID = acc_reg.findall(mem)[0]
            seqType,speciesName = accTable.lookUp(accID)
            speciesName = truncate(speciesName)
            if seqType=="plasmid":
                plasmidHeats[speciesName][bacID] += 1
            else:
                genomeHeats[speciesName][bacID] += 1
    heatMapTree(genomeHeats,tree,"Bacteriocins vs Species Genome Heatmap",
            xlabel='Bacteriocin ID',ylabel='Species',
            showX=True,showY=True)
    heatMapTree(plasmidHeats,tree,"Bacteriocins vs Species Plasmid Heatmap",
            xlabel='Bacteriocin ID',ylabel='Species',
            showX=True,showY=True)
"""
Bacteriocins versus clusters
"""
def bacteriocinHeatmap(cdhitProc,accTable):
    """ Fill in 2D array """
    genomeHeats = defaultdict( list )
    plasmidHeats = defaultdict( list )
    numClusters = len(cdhitProc.clusters)
    species = []
    for i,cluster in enumerate(cdhitProc.clusters):
        members = cluster.seqs
        head =""
        for mem in members:
            bacID=bac_reg.findall(mem)[0]
            accID = acc_reg.findall(mem)[0]
            seqType,speciesName = accTable.lookUp(accID)
            if mem[-1]=="*": head=mem
            if seqType=="plasmid":
                if len(plasmidHeats[bacID])==0:
                    plasmidHeats[bacID] = [0]*numClusters
                plasmidHeats[bacID][i] += 1
            else:
                if len(genomeHeats[bacID])==0:
                    genomeHeats[bacID] = [0]*numClusters
                genomeHeats[bacID][i] += 1
            
    heatMap(genomeHeats,"Bacteriocins vs Clusters Genome Heatmap",
            xlabel='Cluster number',ylabel='Bacteriocin ID',
            showX=False,showY=True)
    heatMap(plasmidHeats,"Bacteriocins vs Clusters Plasmid Heatmap",
            xlabel='Cluster number',ylabel='Bacteriocin ID',
            showX=False,showY=True)
    
"""
Histogram of bacteriocin based on bacteriocin ID
"""
def bacteriocinIDHistogram(cdhitProc):
    ids = []
    for cluster in cdhitProc.clusters:
        members = cluster.seqs
        for mem in members:
            bacID=bac_reg.findall(mem)[0]
            ids.append(bacID)
    numClusters = len(cdhitProc.clusters)
    hist = Counter()
    for i in ids:
        hist[i]+=1
    counts = hist.values()
    labels = hist.keys()
    ind = numpy.arange(len(labels))
    plt.figure()
    plt.bar(ind,counts)
    plt.xticks(ind+0.5,labels,fontsize=14)
    plt.xlabel("Bagel Bacteriocin IDs",fontsize=18)
    plt.ylabel("Number of homologs",fontsize=18)
    plt.title('Bacteriocin homologs',fontsize=22)
    

"""Histogram of bacteriocin counts in each cluster"""
def bacteriocinClusterHistogram(cdhitProc):
    clrCnts = []
    labels = []
    for cluster in cdhitProc.clusters:
        counts = Counter()
        members = cluster.seqs
        head = ''
        for mem in members:
            if mem[-1]=='*': head = mem
            bacID=bac_reg.findall(mem)[0]
            counts[bacID]+=1
        size = len(counts.keys())
        species = acc_reg.findall(head)[0]
        clrCnts.append(size)
        labels.append(species)
    
    ind = numpy.arange(len(labels))        
    plt.figure()
    plt.bar(ind,clrCnts)
    plt.xticks(ind,labels,rotation=45)
    plt.ylabel("Number of different bacteriocins",fontsize=18)
    plt.title('Bacteriocin Count',fontsize=22)
    

"""Histogram of species diversity in each cluster"""
def speciesClusterHistogram(cdhitProc):
    clrCnts = []
    labels = []
    for cluster in cdhitProc.clusters:
        counts = Counter()
        members = cluster.seqs
        head = ''
        for mem in members:
            if mem[-1]=='*': head = mem
            speciesID=acc_reg.findall(mem)[0]
            counts[speciesID]+=1
        size = len(counts.keys())
        tag = cluster_id_reg.findall(head)[0][:-3]
        anchorgene = tag.split('|')[-1]
        clrCnts.append(size)
        labels.append(anchorgene)
    
    ind = numpy.arange(len(labels))        
    plt.figure()
    plt.barh(ind,clrCnts)
    plt.yticks(ind+0.5,labels)
    plt.ylabel('Anchor genes',fontsize=18)
    plt.xlabel("Number of different species",fontsize=18)
    plt.title('Species Count',fontsize=22)
    
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
def contextGeneDistanceHeatmap(cdhitProc):
    clusterIDs = []
    dists = []
    for i,cluster in enumerate(cdhitProc.clusters):
        members = cluster.seqs
        for mem in members:
            tag = cluster_id_reg.findall(mem)[0][:-3]
            bst,bend,ast,aend = tag.split('|')[-4:]
            bst,bend,ast,aend = int(bst),int(bend),int(ast),int(aend)
            if overlap(bst,bend,ast,aend): continue                
            bmid = (bst+bend)/2
            int_st,int_end = ast-bmid,aend-bmid
            interval = range(int_st,int_end)
            dists+= interval
            clusterIDs+=[i]*len(interval)
    
    plt.figure()
    plt.hist2d(dists,clusterIDs,bins=30,norm=LogNorm())
    plt.colorbar()
    plt.xlabel('Distance from bacteriocin',fontsize=18)
    plt.ylabel("Cluster ID",fontsize=18)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.title('Distance distribution of genes per cluster',fontsize=22)
    
"""Distance distribution of all clusters"""
def anchorGenePosition(cdhitProc):
    dists = []
    i = 0
    test_range = 0
    for cluster in cdhitProc.clusters:
        members = cluster.seqs
        for mem in members:
            tag = cluster_id_reg.findall(mem)[0][:-3]
            bst,bend,ast,aend = tag.split('|')[-4:]
            bst,bend,ast,aend = int(bst),int(bend),int(ast),int(aend)
            if overlap(bst,bend,ast,aend): continue                
            bmid = (bst+bend)/2
            int_st,int_end = ast-bmid,aend-bmid
            interval = range(int_st,int_end)
            #dists+= interval
            dists+=[(int_st+int_end)/2]
            i+=1
            test_range+=(int_st>=0 and int_st<10000)
            
    print "Number of sequences",i
    print "Number of bases",len(dists)
    print "Number of sequences between 0 and 1000",test_range
    
    plt.figure()
    plt.hist(dists,bins=80)
    plt.xlabel('Distance from bacteriocin',fontsize=18)
    plt.ylabel("Number of anchor genes",fontsize=18)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.title('Distance distribution of anchor genes',fontsize=22)

"""Histogram of bacteriocin counts in each cluster"""
def anchorGeneClusterHistogram(cdhitProc):
    clrCnts = []
    labels = []
    for cluster in cdhitProc.clusters:
        counts = Counter()
        members = cluster.seqs
        head = ''
        size = len(members)
        clrCnts.append(size)
        
    plt.figure()
    plt.hist(clrCnts,bins=80)
    plt.xlabel("Number of members in a cluster",fontsize=18)
    plt.ylabel("Number of clusters",fontsize=18)
    plt.title('Cluster Size Count',fontsize=22)
    

"""
Histogram of bacteriocin length
"""
def bacteriocinLength(bacteriocinFile):
    lens = []
    with open(bacteriocinFile,'r') as handle:
        dists = []
        for ln in handle:
            ln = ln.rstrip()
            toks = ln.split('\t')
            seq = toks[6]
            #seq = seq.replace("-",'')
            lens.append(len(seq))
            
    plt.figure()
    plt.hist(lens,bins=100)
    plt.xlabel('Distance from bacteriocin',fontsize=18)
    plt.ylabel("Number of anchor genes",fontsize=18)
    plt.title('Distance distribution of anchor genes',fontsize=22)

"""
Histogram of number of bacteriocins versus species
"""
def bacteriocinSpeciesHistogram(bacteriocinFile,accTable):
    cnts = Counter()
    with open(bacteriocinFile,'r') as handle:
        for ln in handle:
            ln = ln.rstrip()
            toks = ln.split('\t')
            species = acc_reg.findall(toks[1])[0]
            cnts[species]+=1
    bacCnts = cnts.values()
    plt.figure()
    plt.hist(bacCnts,bins=40)
    plt.xlabel("Number of bacteriocins")
    plt.ylabel("Number of species")
    plt.title("Number of bacteriocins per species")
    
def findhead(clrs):#finds representative sequence
    for clr in clrs:
        if "*" in clr:
            return clr
"""
Write all sequence in top scoring clusters to fasta file
"""
def getFASTA(infasta,cdhit,outfasta):
    handle = open("lengths.txt",'w')
    record_dict = SeqIO.to_dict(SeqIO.parse(infasta,'fasta'))
    with open(outfasta,'w') as out:
        for cluster in cdhit.clusters:
            header = findhead(cluster.seqs)
            tag = cluster_id_reg.findall(header)[0][:-3]
            record = record_dict[tag]
            out.write(">%s\n"%record.id)
            out.write("%s\n"%str(record.seq))
            handle.write("%d\n"%len(record.seq))

def removeDuplicates(items):
    uniqueDict = {tuple(x[-5:-1]):x for x in items}
    return uniqueDict.values()

"""
Preprocess fasta file
"""
def preprocessFasta(blastTab,fastaout):
    items = []
    with open(blastTab,'r') as handle:
        for ln in handle:
            ln = ln.rstrip()
            toks = ln.split('\t')
            assert len(toks)>=10
            toks = [tok.replace(' ','') for tok in toks] #remove white space
            items.append(tuple(toks))
    items = removeDuplicates(items)
    with open(fastaout,'w') as handle:
        for item in items:
            bacID,gi,bst,bend,bstrand,species,ast,aend,astrand,seq = item
            seqstr = ">%s|%s|%s|%s|%s|%s\n%s\n"%(bacID,gi,bst,bend,ast,aend,seq)
            handle.write(seqstr) 
def getTree(cdhitProc,rrnaFile):
    
    record_dict = SeqIO.to_dict(SeqIO.parse(open(rrnaFile,'r'), "fasta"))
    rrnas = []
    seen = set()
    for cluster in cdhitProc.clusters:
        members = cluster.seqs
        for mem in members:
            #Obtain accession IDs from cdhitProc 
            acc = acc_reg.findall(mem)[0]
            try:
                if acc not in seen:
                    record = record_dict[acc]
                    rrnas.append(record)
                    seen.add(acc)
            except Exception as k:
                print 'Accession missing',k
                
            #Obtain corresponding 16SRNAs
    #print "Number of rRNAs",len(rrnas)
    basename,_ = os.path.splitext(rrnaFile)
    tmp_rrna = "%s_filtered.fasta"%basename
    tree = "%s_filtered.tree"%basename
    SeqIO.write(rrnas, open(tmp_rrna,'w'), "fasta")
    
    #Run FastTree
    ft = UnAlignedFastTree(tmp_rrna,tree)
    ft.align() #Run multiple sequence alignment and spit out aligned fasta file
    ft.run() #Run fasttree on multiple alignment and spit out newick tree
    ft.cleanUp() #Clean up!
    return tree

