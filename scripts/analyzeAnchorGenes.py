import os
import sys
import site
import re
import numpy as np
import numpy.random
import matplotlib.pyplot as plt
import numpy
import pylab
import argparse
import cPickle
from collections import defaultdict
from matplotlib.colors import LogNorm
from Bio import SeqIO
base_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
site.addsitedir(os.path.join(base_path, 'src'))
import cdhit
from  acc2species import AccessionToSpecies
from collections import Counter

bac_reg = re.compile(">(\d+.\d)")
species_reg = re.compile("([A-Z]+_\d+)")
cluster_id_reg = re.compile(">(\S+)")

def truncate(speciesName):
    toks = speciesName.split(' ')
    
    return '_'.join(toks[:2])
    
def heatMap(heats,speciesLabels,title):
    r,c = len(heats.keys()),len(heats[heats.keys()[0]])
    scoremap = np.zeros((r,c))
    i = 0
    for k,v in heats.iteritems():
        s = '\t'.join(map(str,v))
        scoremap[i,:] = v
        i+=1
    
    fig, ax = plt.subplots()
    ypos = np.arange(r)+0.5
    pylab.yticks(ypos, heats.keys())
    xpos = np.arange(c)+0.5
    #pylab.xticks(xpos, speciesLabels,rotation=90)
    
    #heathandle.close()
    heatmap = plt.pcolor(scoremap, cmap='PuBu_r')
    plt.colorbar()
    plt.title(title,fontsize=22)
    plt.xlabel('Cluster number')
    plt.ylabel('Bacteriocin ID')
    plt.xlim(0,c)
    fig.subplots_adjust(bottom=0.2)
"""Heatmap of bacteriocins versus species"""
def bacteriocinSpeciesHeatmap(cdhitProc,accTable):
    genomeHeats = defaultdict( dict )
    plasmidHeats = defaultdict( dict )
    for i,cluster in enumerate(cdhitProc.clusters):
        members = cluster.seqs
        head =""
        for mem in members:
            bacID=bac_reg.findall(mem)[0]
            accID = species_reg.findall(mem)[0]
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
    pass

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
            accID = species_reg.findall(mem)[0]
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
            
        #accID = species_reg.findall(head)[0]
        #seqType,speciesName = accTable.lookUp(accID)
        #speciesName = truncate(speciesName)
        species.append(str(i))
        
    heatMap(genomeHeats,species,"Whole genome heatmap")
    heatMap(plasmidHeats,species,"Plasmid heatmap")
    
    

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
    

"""Histogram of bacteriocin diversity in each cluster"""
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
        species = species_reg.findall(head)[0]
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
            speciesID=species_reg.findall(mem)[0]
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
    
"""
Checks the interval overlaps
"""
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

def anchorGeneDistanceHeatmap(cdhitProc):
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
    plt.hist2d(dists,clusterIDs,bins=200,norm=LogNorm())
    plt.colorbar()
    plt.xlabel('Distance from bacteriocin',fontsize=18)
    plt.ylabel("Cluster ID",fontsize=18)
    plt.title('Distance distribution of anchor genes per cluster',fontsize=22)
        
    
def anchorGenePosition(cdhitProc):
    dists = []
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
            dists+= interval
    plt.figure()
    plt.hist(dists,bins=300)
    plt.xlabel('Distance from bacteriocin',fontsize=18)
    plt.ylabel("Number of anchor genes",fontsize=18)
    plt.title('Distance distribution of anchor genes',fontsize=22)

"""Histogram of bacteriocin diversity in each cluster"""
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
    plt.title('Bacteriocin Count',fontsize=22)
    

    
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


def findhead(clrs):#finds representative sequence
    for clr in clrs:
        if "*" in clr:
            return clr

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

if __name__=="__main__":
    parser = argparse.ArgumentParser(description=\
        'Format Blast output for iTOL visualization ')
    parser.add_argument(\
        '--accession-table', type=str, required=False,
        help='A table that maps accession ids to species')
    parser.add_argument(\
        '--bacteriocins', type=str, required=False,
        help='Bacteriocins in blast format')
    parser.add_argument(\
        '--anchor-genes', type=str, required=False,
        help='Anchor genes from genbank files in blast format')
    parser.add_argument(\
        '--anchor-genes-fasta', type=str, required=False,
        help='Anchor genes from genbank files in fasta format')
    
    args = parser.parse_args()
    directory = os.path.dirname(args.anchor_genes_fasta)    
    input_file = args.anchor_genes_fasta
    cluster_file = "%s/anchor_gene_cluster"%directory
    outfasta = "%s/anchor_gene_cluster.fa"%directory
    threshold = 0.7
    accTable = AccessionToSpecies(args.accession_table)
    clusterFile = "cluster.pickle"
    if os.path.exists(clusterFile):
        cdhitProc = cPickle.load(open(clusterFile,'rb'))
    else:
        cdhitProc = cdhit.CDHit(input_file,cluster_file,threshold)
        cdhitProc.run()
        cdhitProc.parseClusters()
        cdhitProc.countOut()
        cPickle.dump(cdhitProc,open(clusterFile,'wb'))
    getFASTA(cluster_file,cdhitProc,outfasta)
    cdhitProc.filterSize(50)  #filter out everything less than n 
    cdhitProc.filterSize(30)  #filter out everything less than n
    anchorGeneDistanceHeatmap(cdhitProc)
    anchorGeneClusterHistogram(cdhitProc)
    
    #bacteriocinIDHistogram(cdhitProc)
    #bacteriocinClusterHistogram(cdhitProc)
    #speciesClusterHistogram(cdhitProc)
    
    bacteriocinHeatmap(cdhitProc,accTable)    
    anchorGenePosition(cdhitProc)
    plt.show()    
    
    
