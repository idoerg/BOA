"""
Plot histogram of all context genes
"""


import os
import sys
import site
import re
import numpy as np
import numpy.random
import matplotlib.pyplot as plt

base_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
base_path="%s/src"%base_path
for directory_name in os.listdir(base_path):
    site.addsitedir(os.path.join(base_path, directory_name))

import clique_filter
import hmmer
from collections import Counter
import itertools
import cPickle
import gff
import interval_filter
from Bio import SeqIO

"""
Counting clusters with a set number of function

example of a sequence id:
>accession=AE001437.1|function=transport|start=435542|end=436306|strand=+|score=90.1|protein_id=AAK78353.1|cluster_0|-
"""
def operonCounts(fname,numFuncs):
    prevCluster = None
    cluster_counts = 0 #number of clusters with specified number of functions
    operon_cnts = Counter() #['toxin','modifier','immunity','transport','regulator']
    function_cnts = Counter()
        
    for record in SeqIO.parse(open(fname,'r'),"fasta"):
        seqid = record.id
        accession,function,st,end,strand,score,protid,cluster,description = seqid.split("|")
        function = function.split("=")[1]
        if prevCluster==None:
            prevCluster = cluster
            function_cnts[function]+=1
        elif prevCluster==cluster:
            function_cnts[function]+=1
        else:
            prevCluster = cluster
            if len(function_cnts.keys())==numFuncs:
                operon_cnts=operon_cnts+function_cnts
                cluster_counts+=1
            function_cnts = Counter()
            function_cnts[function]+=1
    
    return operon_cnts,cluster_counts 

def functionCounts(fname,funcs):
    prevCluster = None
    cluster_counts = 0 #number of clusters with specified number of functions
    function_cnts = Counter() ###['toxin','modifier','immunity','transport','regulator']
    operon_cnts = Counter()
    for record in SeqIO.parse(open(fname,'r'),"fasta"):
        seqid = record.id
        accession,function,st,end,strand,score,protid,cluster,description = seqid.split("|")
        function = function.split("=")[1]
        if prevCluster==None:
            prevCluster = cluster
            function_cnts[function]+=1
        elif prevCluster==cluster:
            function_cnts[function]+=1
        else:
            prevCluster = cluster
            print function_cnts

            if set(function_cnts.keys())==set(funcs):
                cluster_counts+=1
                operon_cnts = operon_cnts+function_cnts
            function_cnts = Counter()
            function_cnts[function]+=1

    return operon_cnts,cluster_counts 

            
def bargraph(clusters,numFuncs):
	operon_cnts = Counter(['toxin','modifier','immunity','transport','regulator'])
	funcs = set()
	for cluster in clusters:
		cnts = Counter()
		for node in cluster:
			toks = node.split('|')
			acc,clrname,full_evalue,hmm_st,hmm_end,env_st,env_end,description=toks
			function = clrname.split('.')[0]
			funcs.add(function)
			cnts[function]+=1

		if numFuncs==len(funcs): operon_cnts = operon_cnts+cnts
	return operon_cnts
    
    

if __name__=="__main__":
	#faidx = "/Users/mortonyt/Documents/MiamiBio/workspace/all_trans.fai"
	#gffFile = "/Users/mortonyt/Documents/MiamiBio/workspace/all.gff"
	#folder = "/Users/mortonyt/Documents/MiamiBio/workspace"
#     faidx = "/home/mortonjt/Projects/Bacfinder/workspace/quorum/data/all_trans.fai"
#     gffFile = "/home/mortonjt/Projects/Bacfinder/workspace/quorum/data/all.gff"
#     folder = "/home/mortonjt/Projects/Bacfinder/workspace/quorum/intermediate"
#     if os.path.exists("clusters.pickle"):
# 		all_hits = cPickle.load(open("all_hits.pickle",'rb'))
# 		clusters = cPickle.load(open("clusters.pickle",'rb'))
#     else:
# 		if not os.path.exists("all_hits.pickle"):
# 			toxin_hits     = hmmer.parse("%s/toxin.out"%folder)
# 			modifier_hits  = hmmer.parse("%s/modifier.out"%folder)
# 			immunity_hits  = hmmer.parse("%s/immunity.out"%folder)
# 			regulator_hits = hmmer.parse("%s/regulator.out"%folder)
# 			transport_hits = hmmer.parse("%s/transport.out"%folder)
# 
# 			gff = gff.GFF(gff_file=gffFile,fasta_index=faidx)
# 			toxin_hits     = gff.call_orfs(toxin_hits    )
# 			modifier_hits  = gff.call_orfs(modifier_hits )
# 			immunity_hits  = gff.call_orfs(immunity_hits )
# 			regulator_hits = gff.call_orfs(regulator_hits)
# 			transport_hits = gff.call_orfs(transport_hits)
# 			all_hits = toxin_hits+modifier_hits+immunity_hits+regulator_hits+transport_hits
# 			all_hits=sorted(all_hits,key=lambda x: x[6])
# 			all_hits=sorted(all_hits,key=lambda x: x[5])
# 			cPickle.dump(all_hits,open("all_hits.pickle",'wb'))
# 		else:
# 			all_hits = cPickle.load(open("all_hits.pickle",'rb'))
# 
# 		all_hits = sorted(all_hits,key=lambda x: x[0] )
# 		all_hits = interval_filter.overlaps(all_hits,faidx)
# 		all_hits = interval_filter.unique(all_hits)
# 		clusters = clique_filter.findContextGeneClusters(all_hits,faidx,radius=50000,functions = ['toxin','transport'],backtrans=False)
#         
#         #cPickle.dump((all_hits,clusters),open("clusters.pickle",'wb'))

    #func5 = bargraph(clusters,5)
    #func4 = bargraph(clusters,5)
    #func3 = bargraph(clusters,5)
    #func2 = bargraph(clusters,5)
    workspace = "/Users/mortonyt/Documents/MiamiBio/workspace"
    operons = "%s/operons.txt"%workspace
    predicted = "%s/predicted_operons.txt"%workspace
    
    func5,num5 = operonCounts(operons, 5)
    func4,num4 = operonCounts(operons, 4)
    func3,num3 = operonCounts(operons, 3)
    func2,num2 = operonCounts(operons, 2)
    pred,numPred  = functionCounts(predicted,funcs=['modifier','immunity','transport','regulator'])
    
    toxinCnts = np.array([func5['toxin'],func4['toxin'],func3['toxin'],func2['toxin'],pred['toxin']],dtype='float')
    modifierCnts = np.array([func5['modifier'],func4['modifier'],func3['modifier'],func2['modifier'],pred['modifier']],dtype='float')
    immunityCnts = np.array([func5['immunity'],func4['immunity'],func3['immunity'],func2['immunity'],pred['immunity']],dtype='float')
    transportCnts = np.array([func5['transport'],func4['transport'],func3['transport'],func2['transport'],pred['transport']],dtype='float')
    regulatorCnts = np.array([func5['regulator'],func4['regulator'],func3['regulator'],func2['regulator'],pred['regulator']],dtype='float')
    operonCnts = np.array([num5, num4, num3, num2, numPred],dtype='float')
    print toxinCnts
    print modifierCnts
    print immunityCnts
    print transportCnts
    print regulatorCnts
    print operonCnts
    plt.figure(1)
    N = 5
    width = 0.1
    ind = np.arange(N)
    print "Read data"
    plt.bar(ind,operonCnts,width,color='k',label='operons')
    plt.bar(ind+width,toxinCnts,width,color='r',label='toxins')
    plt.bar(ind+2*width,transportCnts,width,color='y',label='transport')
    plt.bar(ind+3*width,immunityCnts,width,color='b',label='immunity')
    plt.bar(ind+4*width,modifierCnts,width,color='g',label='modifier')
    plt.bar(ind+5*width,regulatorCnts,width,color='m',label='regulator')
    plt.legend(loc=1)
    plt.ylabel("Number of genes")
    plt.title("Number of gene function counts")
    plt.xlabel("Type of operons")
    plt.xticks(ind+3*width, ('Five\nfunctions', 'Four\nfunctions', 'Three\nfunctions', 'Two\nfunctions', 'Toxins\nUndetected'))
    
    plt.figure(2)
    toxinCnts = np.divide(toxinCnts,operonCnts)
    modifierCnts = np.divide(modifierCnts,operonCnts)
    immunityCnts = np.divide(immunityCnts,operonCnts)
    regulatorCnts = np.divide(regulatorCnts,operonCnts)
    transportCnts = np.divide(transportCnts,operonCnts)
    print toxinCnts
    print modifierCnts
    print immunityCnts
    print transportCnts
    print regulatorCnts
    print operonCnts
    N = 5
    width = 0.1
    ind = np.arange(N)
    print "Read data"
    plt.bar(ind+width,toxinCnts,width,color='r',label='toxins')
    plt.bar(ind+2*width,transportCnts,width,color='y',label='transport')
    plt.bar(ind+3*width,immunityCnts,width,color='b',label='immunity')
    plt.bar(ind+4*width,modifierCnts,width,color='g',label='modifier')
    plt.bar(ind+5*width,regulatorCnts,width,color='m',label='regulator')
    plt.legend(loc=1)
    plt.ylabel("Genes per operon")
    plt.title("Ratio of gene function counts")
    plt.xlabel("Type of operons")
    plt.xticks(ind+3*width, ('Five\nfunctions', 'Four\nfunctions', 'Three\nfunctions', 'Two\nfunctions', 'Toxins\nUndetected'))
    
    plt.show()
        
