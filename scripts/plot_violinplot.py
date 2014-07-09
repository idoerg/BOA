
import os,sys,site
import argparse
import cPickle

base_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
base_path="%s/src"%base_path
for directory_name in os.listdir(base_path):
    site.addsitedir(os.path.join(base_path, directory_name))
import cdhit
import violinplot

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
if __name__=="__main__":
    quorum_data = "/media/HD/Documents/Jamie/MiamiBio/Bacfinder/workspace/quorum_data"
    operons = "%s/big_operons.txt"%quorum_data
    fastaindex = "%s/all_trans.fai"%quorum_data
    violinplot.operonDistribution(operons,fastaindex)
# 
#     
#     parser = argparse.ArgumentParser(description=\
#         'Format Blast output for iTOL visualization ')
#     parser.add_argument(\
#         '--accession-table', type=str, required=False,
#         help='A table that maps accession ids to species')
#     parser.add_argument(\
#         '--bacteriocins', type=str, required=False,
#         help='Bacteriocins in blast format')
#     parser.add_argument(\
#         '--anchor-genes', type=str, required=False,
#         help='Anchor genes from genbank files in blast format')
#     parser.add_argument(\
#         '--tree', type=str, required=False,default=None,
#         help='Newick tree')
#     parser.add_argument(\
#         '--rRNA', type=str, required=False,
#         help='16SRNA sequences')
#     args = parser.parse_args()
#     #root = os.environ['BACFINDER_HOME']
#     directory = os.path.dirname(args.anchor_genes)  
#     cluster_file = "%s/context_gene_cluster"%directory
#     outfasta = "%s/context_gene_cluster.fa"%directory
#     threshold = 0.7
#     clusterFile = "cluster.pickle"
#     
#     if os.path.exists(clusterFile):
#         cdhitProc = cPickle.load(open(clusterFile,'rb'))
#     else:
#         clrfasta = "anchorgenes.fa"
#         preprocessFasta(args.anchor_genes,clrfasta)
#         cdhitProc = cdhit.CDHit(clrfasta,cluster_file,threshold)
#         cdhitProc.run()
#         cdhitProc.parseClusters()
#         cdhitProc.countOut()
#         cPickle.dump(cdhitProc,open(clusterFile,'wb'))
#         os.remove(clrfasta)
#     cdhitProc.filterSize(55) 
#     violinplot.contextGeneDistances(cdhitProc)
#     
# 
# 
# 

