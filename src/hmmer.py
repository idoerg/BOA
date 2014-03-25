"""
A wrapper script for hmmsearch, hmmbuild and the 
multiple alignment step using ClustalW
"""

import cdhit
import numpy
import re
import os
from Bio import SeqIO

cluster_id_reg = re.compile(">(\S+)")

"""Obtain all file handles for each cluster"""
def obtainFileHandles(clusterProc):
    infasta = clusterProc.input    
    directory = os.path.dirname(infasta)
    return ["%s/cluster%d.fa"%(directory,i) for i in range(0,len(clusterProc.clusters))]

"""Remove all useless files"""
def cleanUp(clusterProc):
    fnames = obtainFileHandles(clusterProc)
    for handle in fnames:
        os.remove(handle)
        
"""Obtains ids for all of the sequences and write them into separate cluster files"""
def writeClusters(clusterProc):
    i = 0
    infasta = clusterProc.input    
    directory = os.path.dirname(infasta)
    record_dict = SeqIO.to_dict(SeqIO.parse(infasta,'fasta'))    
    for cluster in clusterProc.clusters:
        outfile = '%s/cluster%d.fa'%(directory,i)
        handle = open(outfile,'w')
        for subc in cluster.seqs:
            subtitle = cluster_id_reg.findall(subc)[0][:-3]
            record = record_dict[subtitle]
            handle.write(">%s\n"%record.id)
            handle.write("%s\n"%str(record.seq))
        handle.close()
        i+=1
    
"""First separate all of the sequences into clusters of 70%"""
def go(infilepath,outfilepath):
    clusterProc = cdhit.CDHit(infilepath,outfilepath,0.7)
    clusterProc.run()
    clusterProc.parseClusters()
    writeClusters(clusterProc)
    #cleanUp(clusterProc)

if __name__=="__main__":
    directory = '/media/HD/Documents/Jamie/MiamiBio/Bacfinder/workspace'
    inputfile = 'bacteriocin.fa'
    outputfile = 'bacteriocinCluster'
    infilepath = '%s/%s'%(directory,inputfile)
    outfilepath = '%s/%s'%(directory,outputfile)
    go(infilepath,outfilepath)
    