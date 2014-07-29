"""
Runs fasttree and generates a phylogenetic tree
"""

import Bio
from Bio import SeqIO, SeqFeature
from Bio.SeqRecord import SeqRecord

from Bio.Blast import NCBIXML
from Bio.Blast import NCBIStandalone
import sys
import os
import site
import argparse
import string
import numpy
import re
from collections import defaultdict
base_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
for directory_name in os.listdir(base_path):
    site.addsitedir(os.path.join(base_path, directory_name))
import quorum

import subprocess
import cdhit
from accessionMap import GGAccession
from clustalw import ClustalW
from muscle import Muscle
from  acc2species import AccessionToSpecies
from muscle import Muscle
from mafft import MAFFT



class FastTree(object):
    def __init__(self,algnFile,treeFile):
        self.algnFile = algnFile
        self.treeFile = treeFile
    def run(self,module=subprocess):
        cmd = "fasttree -nt %s"%self.algnFile
        child = module.Popen(cmd,stdout=open(self.treeFile,'w+'),shell=True)
        if module==quorum: child.submit()
        return child
    
class UnAlignedFastTree(FastTree):
    def __init__(self,rawSeqs,treeFile):
        self.rawSeqs=rawSeqs
        basename,_ = os.path.splitext(rawSeqs)
        algnFile = "%s.align"%basename
        #treeFile = "%s.tree"%basename
        #self.accSeqs = "%s.acc"%basename
        super(UnAlignedFastTree,self).__init__(algnFile,treeFile)
    """ Run fasta tree """
    def run(self,module=subprocess):
        proc = super(UnAlignedFastTree,self).run(module)
        return proc
    """ Run multiple alignment """
    def align(self,module=subprocess,MSA=Muscle,iters=4,threads=8,hours=12):
        #assert os.path.exists(self.accSeqs)
        #assert os.stat(self.accSeqs).st_size!=0 #assert not empty
        #cw = ClustalW(self.rawSeqs,self.algnFile)
        self.cw = MSA(self.rawSeqs,self.algnFile,module=module)
        if MSA==Muscle:
            self.cw.run(fasta=True,maxiters=iters,maxhours=hours)
        elif MSA==MAFFT:
            self.cw.run(fasta=True,maxiters=iters,threads=threads)
        self.cw.wait()
        #print self.accSeqs
        #cw.cleanUp()
    def cleanUp(self):
        #os.remove(self.accSeqs)
        #os.remove(self.algnFile)
        self.cw.cleanUp()
        pass

acc_reg = re.compile("([A-Z]+_\d+)")
def truncate(speciesName):
    toks = speciesName.split(' ')    
    return '_'.join(toks[:2])

def removeDuplicates(items):
    uniqueDict = {tuple(x[-5:-1]):x for x in items}
    return uniqueDict.values()

"""
Turns blast tabbed output into a fasta file
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
"""
Replaces all accession ids with species names
and matches 16SRNA sequences with species
"""
def accTospeciesSwap(cdhitProc, #Clustering object
                     rrnaFasta, #Fasta from 16SRNA
                     accTable,
                     fastaout):  
    rrna_dict = SeqIO.to_dict(SeqIO.parse(rrnaFasta,'fasta'))
    outHandle = open(fastaout,'w')
    seen = set()
    for cluster in cdhitProc.clusters:
        members = cluster.seqs
        for mem in members:
            accID = acc_reg.findall(mem)[0]
            seqType,speciesName = accTable.lookUp(accID)
            speciesName = truncate(speciesName)
            if speciesName!=None and speciesName not in seen:
                try:
                    rrnaRecord = rrna_dict[accID]
                    rrnaRecord.id = speciesName
                    outHandle.write(rrnaRecord.format('fasta'))
                    seen.add(speciesName)
                    print "Printed record",speciesName
                except KeyError as k:
                    print "16SRNA not found",k
            else:
                print "Already in set",speciesName
    outHandle.close()
    
def go(refTabs,
       in16SRNA,
       treeFile,
       accTable,
       threshold):
    basename,_ = os.path.splitext(refTabs)
    refFasta = "%s_ref.fasta"%basename
    clrFasta = "%s_clr.fasta"%basename
    filtered16SRNA = "%s_16SRNA.fasta"%basename
    
    preprocessFasta(refTabs,refFasta)
    cdhitProc = cdhit.CDHit(refFasta,clrFasta,0.7)
    cdhitProc.run()
    cdhitProc.parseClusters()
    cdhitProc.filterSize(threshold)     
    
    accTospeciesSwap(cdhitProc,
                     in16SRNA,
                     accTable,
                     filtered16SRNA)
    ft = UnAlignedFastTree(filtered16SRNA,treeFile)
    ft.align() #Run multiple sequence alignment and spit out aligned fasta file
    ft.run() #Run fasttree on multiple alignment and spit out newick tree
    ft.cleanUp() #Clean up!
    #os.remove(clrFasta)
    #os.remove(refFile)
    #os.remove(filtered16SRNA)
    
if __name__=="__main__":
    parser = argparse.ArgumentParser(description=\
        'Runs fasttree on unaligned 16SRNA sequences from greengenes')
    parser.add_argument(\
        '--blast-tab', type=str, required=False,
        help='Blast output')
    parser.add_argument(\
        '--accession-table', type=str, required=False,
        help='A table that maps accession ids to species')
    parser.add_argument(\
        '--rRNA', type=str, required=False,
        help='FASTA files containing 16SRNA sequences')
    parser.add_argument(\
        '--tree', type=str, required=False,
        help='Output newick tree')
    parser.add_argument(\
        '--cluster-threshold', type=int, required=False, default=60,
        help='Threshold for filtering clusters by size')
    parser.add_argument(\
        '--test', action='store_const', const=True, default=False,
        help='Run unittests')
    args = parser.parse_args()

    if not args.test:
        accTable = AccessionToSpecies(args.accession_table)

        go(args.blast_tab,
           args.rRNA,
           args.tree,
           accTable,
           args.cluster_threshold)
        
    else:
        del sys.argv[1:]
        import unittest
        class TestRun(unittest.TestCase):
            def setUp(self):
                self.infile="test.fa"
                self.ggtable = "../data/gg_13_5_accessions.txt"
                self.treefile = "test.tree"
                #self.accSeqs= "acc.fa"
                entries = [">1111886",
                           "AACGAACGCTGGCGGCATGCCTAACACATGCAAGTCGAACGAGACCTTCGGGTCTAGTGGCGCACGGGTGCGTAACGCGTGGGAATCTGCCCTTGGGTACGGAATAACAGTTAGAAATGACTGCTAATACCGTATAATGACTTCGGTCCAAAGATTTATCGCCCAGGGATGAGCCCGCGTAGGATTAGCTTGTTGGTGAGGTAAAGGCTCACCAAGGCGACGATCCTTAGCTGGTCTGAGAGGATGATCAGCCACACTGGGACTGAGACATGGCCCAGACTCCTACGGGAGGCAGCAGTGGGGAATATTGGACAATGGGCGAAAGCCTGATCCAGCAATGCCGCGTGAGTGATGAAGGCCTTAGGGTTGTAAAGCTCTTTTACCCGGGATGATAATGACAGTACCGGGAGAATAAGCCCCGGCTAACTCCGTGCCAGCAGCCGCGGTAATACGGAGGGGGCTAGCGTTGTTCGGAATTACTGGGCGTAAAGCGCACGTAGGCGGCTTTGTAAGTTAGAGGTGAAAGCCCGGGGCTCAACTCCGGAATTGCCTTTAAGACTGCATCGCTAGAATTGTGGAGAGGTGAGTGGAATTCCGAGTGTAGAGGTGAAATTCGTAGATATTCGGAAGAACACCAGTGGCGAAGGCGACTCACTGGACACATATTGACGCTGAGGTGCGAAAGCGTGGGGAGCAAACAGGATTAGATACCCTGGTAGTCCACGCCGTAAACGATGATGACTAGCTGTCGGGGCGCTTAGCGTTTCGGTGGCGCAGCTAACGCGTTAAGTCATCCGCCTGGGGAGTACGGCCGCAAGGTTAAAACTCAAAGAAATTGACGGGGGCCTGCACAAGCGGTGGAGCATGTGGTTTAATTCGAAGCAACGCGCAGAACCTTACCAGCGTTTGACATGGTAGGACGGTTTCCAGAGATGGATTCCTACCCTTACGGGACCTACACACAGGTGCTGCATGGCTGTCGTCAGCTCGTGTCGTGAGATGTTGGGTTAAGTCCCGCAACGAGCGCAACCCTCGTCTTTGGTTGCTACCATTTAGTTGAGCACTCTAAAAAAACTGCCGGTGATAAGCCGGAGGAAGGTGGGGATGACGTCAAGTCCTCATAGCCCTTACGCGCTGGGCTACACACGTGCTACAATGGCGGTGACAGAGGGCAGCAAACCCGCGAGGGTGAGCTAATCTCCAAAAGCCGTCTCAGTTCGGATTGTTCTCTGCAACTCGAGAGCATGAAGGCGGAATCGCTAGTAATCGCGGATCAGCACGCCGCGGTGAATACGTTCCCAGGCCTTGTACACACCGCCCGTCACATCACGAAAGTCGGTTGCACTAGAAGTCGGTGGGCTAACCCGCAAGGGAGGCAGCCGCCTAAAGTGTGATCGGTAATTGGGGTG",
                           ">1111885",
                           "AGAGTTTGATCCTGGCTCAGAATGAACGCTGGCGGCGTGCCTAACACATGCAAGTCGTACGAGAAATCCCGAGCTTGCTTGGGAAAGTAAAGTGGCGCACGGGTGAGTAACGCGTGGGTAACCCACCCCCGAATTCGGGATAACTCCGCGAAAGCGGTGCTAATACCGGATAAGACCCCTACCGCTTCGGCGGCAGAGGTAAAAGCTGACCTCTCCATGGAAGTTAGCGTTTGGGGACGGGCCCGCGTCCTATCAGCTTGTTGGTGGGGTAACAGCCCACCAAGGCAACGACGGGTAACTGGTCTGAGAGGATGATCAGTCACACTGGAACTGGAACACGGTCCAGACTCCTACGGGAGGCAGCAGTGAGGAATTTTGCGCAATGGGCGAAAGCCTGACGCAGCAACGCCGCGTGGGTGAAGAAGGCTTTCGGGTCGTAAAGCCCTGTCAGGTGGGAAGAAACCTTTCCGGTACTAATAATGCCGGAAATTGACGGTACCACCAAAGGAAGCACCGGCCAACTCCGTGCCAGCAGCCGCGGTAATACGGAGGGTGCAAGCGTTGTTCGGAATTATGGGGCGTAAAGAGCGTGTGGGCGGTTAGGAAAGTCAGATGTGAAAGCCCTGGGCTCAACCCAGGAAGTGCATTTGAAACTGCCTAACTTGAGTACGGGAGAGGAAGGGGGAATTCCCGGTGTAGAGGTGAAATTCGTAGATATCGGGAGGAATACCGGTGGCGAAGGCGCCCTTCTGGACCGATACTGACGCTGAGACGCGAAAGCGTGGGGAGCAAACAGGATTAGATACCCTGGTAGTCCACGCCGTAAACGATGAGCACTAGGTGTAGCGGGTATTGACCCCTGCTGTGCCGTAGCTAACGCATTAAGTGCTCCGCCTGGGGATTACGGTCGCAAGACTAAAACTCAAAGGAATTGACGGGGGCCCGCACAAGCGGTGGAGCATGTGGTTTAATTCGACGCAACGCGAAGAACCTTACCTGGGCTTGACATCCCCGGACAGCCCTGGAAACAGGGTCTCCCACTTCGGTGGGCTGGGTGACAGGTGCTGCATGGCTGTCGTCAGCTCGTGTCGTGAGATGTTGGGTTAAGTCCCGCAACGAGCGCAACCCCTGCCTTTAGTTGCCATCATTTAGCTGGGCACTCTAAAGGGACTGCCGGTGTTAAACCGGAGGAAGGTGGGGACGACGTCAAGTCCTCATGGCCTTTATGCCCAGGGCTACACACGTGCTACAATGGGCGGTACAAAGGGCAGCGACATCGTGAGGTGAAGCAAATCCCAAAAAACCGCTCTCAGTTCGGATCGGAGTCTGCAACTCGACTTCGTGAAGGTGGAATCACTAGTAATCGTGGATCAGCATGCCACGGTGAATACGTTCCCGGGCCTTGTACACACCGCCCGTCACACCACGAAAGTCTGCTGTACCAGAAGTCGCTGGGCTAACCCGCCCTAGGCGGGAGGTAGGCGCCTAAGGTACGGCCGGTAATTGGGGTGAAGTCGTAACAAGGTAACC",
                           ">1111883",
                           "GCTGGCGGCGTGCCTAACACATGTAAGTCGAACGGGACTGGGGGCAACTCCAGTTCAGTGGCAGACGGGTGCGTAACACGTGAGCAACTTGTCCGACGGCGGGGGATAGCCGGCCCAACGGCCGGGTAATACCGCGTACGCTCGTTTAGGGACATCCCTGAATGAGGAAAGCCGTAAGGCACCGACGGAGAGGCTCGCGGCCTATCAGCTAGTTGGCGGGGTAACGGCCCACCAAGGCGACGACGGGTAGCTGGTCTGAGAGGATGGCCAGCCACATTGGGACTGAGACACGGCCCAGACTCCTACGGGAGGCAGCAGTGGGGAATCTTGCGCAATGGCCGCAAGGCTGACGCAGCGACGCCGCGTGTGGGATGACGGCCTTCGGGTTGTAAACCACTGTCGGGAGGAACGAATACTCGGCTAGTCCGAGGGTGACGGTACCTCCAAAGGAAGCACCGGCTAACTCCGTGCCAGCAGCCGCGGTAATACGGAGGGTGCGAGCGTTGTCCGGAATCACTGGGCGTAAAGGGCGCGTAGGTGGCCCGTTAAGTGGCTGGTGAAATCCCGGGGCTCAACTCCGGGGCTGCCGGTCAGACTGGCGAGCTAGAGCACGGTAGGGGCAGATGGAATTCCCGGTGTAGCGGTGGAATGCGTAGATATCGGGAAGAATACCAGTGGCGAAGGCGTTCTGCTGGGCCGTTGCTGACACTGAGGCGCGACAGCGTGGGGAGCAAACAGGATTAGATACCCTGGTAGTCCACGCCGTAAACGATGGACACTAGACGTCGGGGGGAGCGACCCTCCCGGTGTCGTCGCTAACGCAGTAAGTGTCCCGCCTGGGGAGTACGGCCGCAAGGCTGAAACTCAAAGGAATTGACGGGGGCCCGCACAAGCGGTGGAGCATGTGGTTTAATTCGAAGCAACGCGAAGAACCTTACCTGGGCTTGACATGCTGGTGCAAGCCGGTGGAAACATCGGCCCCTCTTCGGAGCGCCAGCACAGGTGCTGCATGGCTGTCGTCAGCTCGTGTCGTGAGATGTTGGGTTAAGTCCCGCAACGAGCGCAACTCTCGCTCCCAGTTGCCAGCGGTTCGGCCGGGGACTCTGGGGGGACTGCCGGCGTTAAGCCGGAGGAAGGTGGGGACGACGTCAAGTCATCATGGCCCTTACGTCCAGGGCGACACACGTGCTACAATGCCTGGTACAGCGCGTCGCGAACTCGCAAGAGGGAGCCAATCGCCAAAAGCCGGGCTAAGTTCGGATTGTCGTCTGCAACTCGACGGCATGAAGCCGGAATCGCTAGTAATCGCGGATCAGCCACGCCGCGGTGAATACGTTCCCGGGCCTTGTACACACCGCCCGTCACGCCATGGAAGCCGGAGGGACCCGAAACCGGTGGGCCAACCGCAAGGGGGCAGCCGTCTAAGGT",
                           ">1111882",
                           "AGAGTTTGATCATGGCTCAGGATGAACGCTAGCGGCAGGCCTAACACATGCAAGTCGAGGGGTAGAGGCTTTCGGGCCTTGAGACCGGCGCACGGGTGCGTAACGCGTATGCAATCTGCCTTGTACTAAGGGATAGCCCAGAGAAATTTGGATTAATACCTTATAGTATATAGATGTGGCATCACATTTCTATTAAAGATTTATCGGTACAAGATGAGCATGCGTCCCATTAGCTAGTTGGTATGGTAACGGCATACCAAGGCAATGATGGGTAGGGGTCCTGAGAGGGAGATCCCCCACACTGGTACTGAGACACGGACCAGACTCCTACGGGAGGCAGCAGTGAGGAATATTGGTCAATGGGCGCAAGCCTGAACCAGCCATGCCGCGTGCAGGATGACGGTCCTATGGATTGTAAACTGCTTTTGTACGGGAAGAAACACTCCTACGTGTAGGGGCTTGACGGTACCGTAAGAATAAGGATCGGCTAACTCCGTGCCAGCAGCCGCGGTAATACGGAGGATCCAAGCGTTATCCGGAATCATTGGGTTTAAAGGGTCCGTAGGCGGTTTTATAAGTCAGTGGTGAAATCCGGCAGCTCAACTGTCGAACTGCCATTGATACTGTAGAACTTGAATTACTGTGAAGTAACTAGAATATGTAGTGTAGCGGTGAAATGCTTAGATATTACATGGAATACCAATTGCGAAGGCAGGTTACTAACAGTATATTGACGCTGATGGACGAAAGCGTGGGGAGCGAACAGGATTAGATACCCTGGTAGTCCACGCCGTAAACGATGGATACTAGCTGTTTGGCAGCAATGCTGAGTGGCTAAGCGAAAGTGTTAAGTATCCCACCTGGGGAGTACGAACGCAAGTTTGAAACTCAAAGGAATTGACGGGGGCCCGCACAAGCGGTGGAGCATGTGGTTTAATTCGATGATACGCGAGGAACCTTACCAGGGCTTAAATGTAGAGTGACAGGACTGGAAACAGTTTTTTCTTCGGACACTTTACAAGGTGCTGCATGGTTGTCGTCAGCTCGTGCCGTGAGGTGTCAGGTTAAGTCCTATAACGAGCGCAACCCCTGTTGTTAGTTGCCAGCGAGTAATGTCGGGAACTCTAACAAGACTGCCGGTGCAAACCGTGAGGAAGGTGGGGATGACGTCAAATCATCACGGCCCTTACGTCCTGGGCTACACACGTGCTACAATGGCCGGTACAGAGAGCAGCCACCTCGCGAGGGGGAGCGAATCTATAAAGCCGGTCACAGTTCGGATTGGAGTCTGCAACCCGACTCCATGAAGCTGGAATCGCTAGTAATCGGATATCAGCCATGATCCGGTGAATACGTTCCCGGGCCTTGTACACACCGCCCGTCAAGCCATGGAAGCTGGGGGTACCTGAAGTCGGTGACCGCAAGGAGCTGCCTAGGGTAAAACTGGTAACTGGGGCTAAGTCGTACAAGGTAGCCGTA",
                           ">1111879",
                           "CCTAATGCATGCAAGTCGAACGCAGCAGGCGTGCCTGGCTGCGTGGCGAACGGCTGACGAACACGTGGGTGACCTGCCCCGGAGTGGGGGATACCCCGTCGAAAGACGGGACAATCACGCATACGCTCTTTGGAGGAAAGCCATCCGGCGCTCTGGGAGGGGCCTGCGGCCCATCAGGTAGTTGGTGTGGTAACGGCGCACCAAGCCAATGACGGGTACCCGGTCTGAGAGGACGACCGGCCAGACTGGAACTGCGACACGGCCCAGACTCCTACGGGAGGCAGCAGCAAGGAATTTTCCCCAATGGGCGCAAGCCTGAGGCAGCAACGCCGCGTGCGGGATGACGGACTTCGGGTTGTAAACCGCTTTTCGGGGGGACAACCCTGACGGTACCCCCGGAACAAGCCCCGGCTAACTCTGTGCCAGCAGCCGCGGTAAGACAGAGGGGGCAAGCGTTGTCCGGAGTCACTGGGCGTAAAGCGCGCGCAGGCGGCTGCCTAAGTGTCGTGTGAAAGCCCCCGGCTCAACCGGGGGAGGCCATGGCAAACTGGGTGGCTCGAGCGGCGGAGAGGTCCCTCGAATTGCCGGTGTAGCGGTGAAATGCGTAGAGATCGGCAGGAAGACCAAGGGGGAAGCCAGGGGGCTGGCCGCCGGCTGACGCTGAGGCGCGACAGCGTGGGGAGCAAACCGGATTAGATACCCGGGTAGTCCACGCCGTAAACGATGACCACTCGGCGTGTGGCGACTATTAACGTCGCGGCGCGCCCTAGCTCACGCGATAAGTGGTCCGCCTGGGAACTACGAGCGCAAGCTTAAAACTCAAAGGAATTGACGGGGGCCCGCACAAGCAGCGGAGCGTGTGGTTTAATTCGACGCAACCCGCAGAACCTTACCCAGACTGGACATGACGGTGCAGACGGCGGAAACGTCGTCGCCTGCGAGGGTCCGTCACAGGTGCTGCATGGCTGTCGTCAGCTCGTGTCGTGAGATGTTGGGTTAAGTCCCGCAACGAGCGCAACCCCTGCGGTTAGTTACCCGTGTCTAACCGGACTGCCCTTCGGGGAGGAAGGCGGGGATGACGTCAAGTCCGCATGGCCCTTACGTCTGGGGCGACACACACGCTACAATGGCGCCGACAATGCGTCGCTCCCGCGCAAGCGGATGCTAATCGCCAAACGGCGCCCCAGTGCAGATCGGGGGCTGCAACTCGCCCCCGTGAAGGCGGAGTTGCTAGTAACCGCGTATCAGCCATGGCGCGGTGAATACGTACCCGGGCCTTGTACACACCGCCCGTCACGTCATGGAGTTGTCAATGCCTGAAGTCCGCCAGCTAACC"
                           ]
                open(self.infile,'w').write('\n'.join(entries))
            def tearDown(self):
                os.remove(self.infile)
                import glob
                for file in glob.glob("*.tree"):
                    os.remove(file)
                
            def testAlign(self):
                ft = UnAlignedFastTree(self.infile,self.treefile)
                ft.align()
                self.assertTrue( os.path.exists(ft.algnFile) )
                self.assertTrue( os.stat(ft.algnFile).st_size!=0 )
                ft.cleanUp()
                os.remove(ft.treeFile)
            def testRun(self):
                ft = UnAlignedFastTree(self.infile,self.treefile)
                ft.align()
                ft.run()
                self.assertTrue( os.path.exists(ft.treeFile) )
                self.assertTrue( os.stat(ft.treeFile).st_size!=0 )
                ft.cleanUp() 
                os.remove(ft.treeFile)
        unittest.main()
        
        
        
        
        
        
        