import os
import sys
import site
import re
import numpy as np
import numpy.random
import matplotlib.pyplot as plt
import numpy
import pylab
from collections import defaultdict
from matplotlib.colors import LogNorm
from Bio import SeqIO
base_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
site.addsitedir(os.path.join(base_path, 'src'))
import cdhit
from collections import Counter

bac_reg = re.compile(">(\d+.\d)")
species_reg = re.compile("([A-Z]+_\d+)")
cluster_id_reg = re.compile(">(\S+)")
def bacteriocinHeatmap(cdhitProc):
    """
    Fill in 2D array
    """
    bacHeats = defaultdict( list )
    numClusters = len(cdhitProc.clusters)
    for i,cluster in enumerate(cdhitProc.clusters):
        members = cluster.seqs
        for mem in members:
            bacID=bac_reg.findall(mem)[0]
            if len(bacHeats[bacID])==0:
                bacHeats[bacID] = [0]*numClusters
            bacHeats[bacID][i] += 1
    return bacHeats

"""
Histogram of bacteriocins
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
    for k in hist.keys():
        hist[k] /= float(numClusters)
    return hist

"""Histogram of bacteriocin diversity in each cluster"""
def bacteriocinClusterHistogram(cdhitProc):
    clrCnts = []
    for cluster in cdhitProc.clusters:
        counts = Counter()
        members = cluster.seqs
        for mem in members:
            bacID=bac_reg.findall(mem)[0]
            counts[bacID]+=1
        size = len(counts.keys())
        clrCnts.append(size)
    return clrCnts

"""Histogram of bacteriocin diversity in each cluster"""
def speciesClusterHistogram(cdhitProc):
    clrCnts = []
    for cluster in cdhitProc.clusters:
        counts = Counter()
        members = cluster.seqs
        for mem in members:
            speciesID=species_reg.findall(mem)[0]
            counts[speciesID]+=1
        size = len(counts.keys())
        clrCnts.append(size)
    return clrCnts

def findhead(clrs):
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
    input_file = "all_bacteria_annotations.txt.fa"
    output_file = "anchor_gene_cluster"
    outfasta = "anchor_gene_cluster_top_hits.fa"
    threshold = 0.7
    cdhitProc = cdhit.CDHit(input_file,output_file,threshold)
    cdhitProc.run()
    cdhitProc.parseClusters()
    cdhitProc.countOut()
    cdhitProc.filterSize(10)#filter out everything less than n
    getFASTA(output_file,cdhitProc,outfasta)
    handle = open('bacteriocin_table.txt','w')
    for clr in cdhitProc.clusters:
        locus = clr.seqs[0].split('|')[-1]
        handle.write("%s\t\n"%locus)
    hist = bacteriocinIDHistogram(cdhitProc)
    keys = hist.keys()
    values = hist.values()

    bachandle = open("bacteriocinFrequency.txt",'w')
    clusterhandle = open("bacteriocinDiversity.txt",'w')
    specieshandle = open("speciesClusters.txt",'w')
    heathandle = open("bacteriocinHeatmap.txt",'w')

    bachandle.write('\n'.join(["%s\t%s"%(z) for z in zip(keys,values)]))
    clrcnts = bacteriocinClusterHistogram(cdhitProc)
    clusterhandle.write('\n'.join(map(str,clrcnts)))
    speciescnts = speciesClusterHistogram(cdhitProc)
    specieshandle.write('\n'.join(map(str,speciescnts)))

    bacHeat = bacteriocinHeatmap(cdhitProc)
    r,c = len(bacHeat.keys()),len(bacHeat[bacHeat.keys()[0]])
    scoremap = np.zeros((r,c))
    i = 0
    
    heathandle.write('\t'.join(map(str,range(c+1)))+'\n')
    heathandle.write('-\n')
    for k,v in bacHeat.iteritems():
        s = '\t'.join(map(str,v))
        scoremap[i,:] = v
        i+=1
        heathandle.write("%s|\t%s\n"%(k,s))

    fig, ax = plt.subplots()
    pos = np.arange(r)+0.5
    pylab.yticks(pos, bacHeat.keys())

    heathandle.close()
    heatmap = plt.pcolor(scoremap, cmap='PuBu_r')
    plt.colorbar()
    plt.title('Bacteriocins mapped to anchor gene clusters')
    plt.xlabel('Cluster number')
    plt.ylabel('Bacteriocin ID')
    plt.xlim(0,c)
    plt.show()


