

"""
1. Gather blasted bacteriocins and context genes
2. Cluster bacteriocins and context genes
3. Run naive bayes to classifiy individual bacteriocins/context genes
4. Run majority vote to classify entire clustesr
"""
import Bio
from Bio import SeqIO, SeqFeature
from Bio.SeqRecord import SeqRecord

from collections import defaultdict

import sys
import os,shutil
import site
import argparse
import string
import numpy
import re
import subprocess
import pickle

base_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
for directory_name in os.listdir(base_path):
    site.addsitedir(os.path.join(base_path, directory_name))
print base_path
import genbank
import blast
import intergene
import genome
import annotated_genes
import bacteriocin
import nbayes
import cdhit

""" Remove duplicate entries"""
def removeDuplicates(items):
    uniqueDict = {tuple(x[-5:-1]):x for x in items}
    return uniqueDict.values()
""" Preprocess fasta file """
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

class PipelineHandler(object):
    def __init__(self,
                 rootdir,
                 genome_files,
                 intergenes,
                 annotated_genes,
                 bacteriocins,
                 bacteriocin_radius,
                 similarity,
                 bac_evalue,
                 training_labels,
                 intermediate,
                 output,
                 numThreads,
                 formatdb,
                 verbose                
                 ):
        
        #Declare global vars
        self.rootdir            =     rootdir   			
        self.genome_files       =     genome_files
        self.bacteriocins   	=	  bacteriocins  		
        self.bacteriocin_radius =	  bacteriocin_radius	
        self.similarity	        =	  similarity			
        self.bac_evalue	        =	  bac_evalue
        self.training_labels    =     training_labels			
        self.intermediate   	=	  intermediate  		
        self.output	         	=	  output
        self.numThreads         =     numThreads
        self.formatdb           =     formatdb				
        self.verbose			=	  verbose
        if not os.path.exists(self.intermediate):
            #shutil.rmtree(self.intermediate)
            os.mkdir(self.intermediate)
            
        if intergenes!=None and annotated_genes!=None: 
            self.intergenes	        =	  intergenes			
            self.annotated_genes    =	  annotated_genes   	
        else:
            self.intergenes = "%s/intergeneDB.fa"%self.intermediate
            self.annotated_genes = "%s/annotated_genesDB.fa"%self.intermediate
        
        self.blasted_tab_bacteriocins = "%s/blasted_bacteriocins.txt"%self.intermediate
        self.blasted_fasta_bacteriocins = "%s/blasted_bacteriocins.fa"%self.intermediate
        self.cand_context_genes_tab = "%s/cand_context_genes.txt"%self.intermediate
        self.cand_context_genes_fasta = "%s/cand_context_genes.fa"%self.intermediate
        self.cand_context_cluster = "%s/cand_context_cluster"%self.intermediate
        #Declare object handlers
        self.clusterer = None                       
        self.textClassifier = None           
        self.nbpickle = "nb.zip"
        self.clusterpickle = "cluster.zip"
        self.textout = "text_out.txt"
    def load(self):
        self.textClassifier = nbayes.NBayes(self.training_labels)
        self.textClassifier.load(self.nbpickle)
        self.clusterer = cdhit.CDHit(self.cand_context_genes_fasta,
                                     self.cand_context_cluster,
                                     self.similarity)
        self.clusterer.load(self.clusterpickle)
    """ Builds database such as the 
        intergenic database
        annotated genes database 
        naive bayes model"""
    def preprocess(self):
        print "Preprocessing"
        annotated_genes.go(self.rootdir,self.annotated_genes) 
        intergene.go(self.rootdir,self.intergenes)
        self.textClassifier = nbayes.NBayes(self.training_labels)
        self.textClassifier.train()
        self.textClassifier.dump(self.nbpickle)
        print "Dumped pickle file"
    """  Runs blast to identify bacteriocins and context genes"""
    def blast(self):
        print "Blasting"
        bacteriocin.main(self.genome_files,
                         self.bacteriocins,
                         self.intergenes,
                         self.annotated_genes,
                         open(self.blasted_tab_bacteriocins,'w'),
                         open(self.cand_context_genes_tab,'w'),
                         self.intermediate,
                         self.bac_evalue,
                         self.numThreads,
                         self.formatdb,
                         self.bacteriocin_radius,
                         self.verbose,
                         False)
    """ Clusters bacteriocins and context genes together"""
    def cluster(self):
        print "Clustering"
        preprocessFasta(self.cand_context_genes_tab,
                        self.cand_context_genes_fasta)
        
        
        self.clusterer = cdhit.CDHit(self.cand_context_genes_fasta,
                                     self.cand_context_cluster,
                                     self.similarity)
        self.clusterer.run()
        self.clusterer.parseClusters()
        self.clusterer.dump(self.clusterpickle)
    """ Classifies individual bacteriocins and context genes based on their text"""
    def naiveBayes(self):
        sets = self.textClassifier.classify(self.cand_context_cluster)
        titles,labels = zip(*sets)
        print sets
        open(self.textout,'w').write('\n'.join(map(str,sets)))
        pass

    """ Classifies entire clusters based on a majority vote"""
    def majorityVote(self):
        pass
    def main(self):
        pass

if __name__=="__main__":
    parser = argparse.ArgumentParser(description=\
        'Finds bacteriocins and context genes')
    parser.add_argument(\
        '--pipeline-section', type=str, required=False, default="all",
        help='Section of the pipeline to run (all, preprocess, blast, cluster, nbayes, mvote)')
    parser.add_argument(\
        '--root-dir',type=str, required=False,
        help='Root directory')
    parser.add_argument(\
        '--genome-files', type=str, nargs="+", required=False,
        help='FASTA files containing bacterial genomes')
    parser.add_argument(\
        '--intergenes', type=str, required=False,default=None,
        help='FASTA files containing intergenic regions')
    parser.add_argument(\
        '--annotated-genes', type=str, required=False,default=None,
        help='FASTA files containing annotated genetic regions')
    parser.add_argument(\
        '--bacteriocins', type=str, required=False,default=None,
        help='The bacteriocin proteins that are to be blasted')
    parser.add_argument(\
        '--bacteriocin-radius', type=int, required=False, default=5000,
        help='The search radius around every specified bacteriocin')
    parser.add_argument(\
        '--similarity', type=int, required=False, default=0.7,
        help='Clustering similarity')    
    parser.add_argument(\
        '--bac-evalue', type=float, required=False, default=0.00001,
        help='The evalue for bacteriocin hits')
    parser.add_argument(\
        '--training-labels',type=str,required=False,default=None,
        help='Training labels for Naive Bayes')
    parser.add_argument(\
        '--intermediate', type=str, required=False,default='.',
        help='Directory for storing intermediate files')
    parser.add_argument(\
        '--output', type=str, required=False,
        help='The output file basename for filtered annotationed regions and bacteriocins')
    parser.add_argument(\
        '--formatdb', action='store_const', const=True, default=False,
        help='Indicates if formatdb should be run')
    parser.add_argument(\
        '--num-threads', type=int, required=False, default=1,
        help='The number of threads to be run by BLAST')
    parser.add_argument(\
        '--verbose', action='store_const', const=True, default=False,
        help='Messages for debugging')
    parser.add_argument(\
        '--test', action='store_const', const=True, default=False,
        help='Run unittests')
    args = parser.parse_args()
    
    if not args.test:
        proc = PipelineHandler(  args.root_dir,
                                 args.genome_files,
                                 args.intergenes,
                                 args.annotated_genes,
                                 args.bacteriocins,
                                 args.bacteriocin_radius,
                                 args.similarity,
                                 args.bac_evalue,
                                 args.training_labels,
                                 args.intermediate,
                                 args.output,
                                 args.num_threads,
                                 args.formatdb,
                                 args.verbose                
                                ) 
    
        proc.preprocess()
        proc.blast()
        proc.cluster()
        proc.naiveBayes()
        
    else:
        del sys.argv[1:]
        import unittest
        import test_modules
        class TestRun(unittest.TestCase):
            def setUp(self):
                self.root = os.environ['BACFINDER_HOME']
                self.exampledir = "%s/example/Streptococcus_pyogenes"%self.root
                self.bacdir = "%s/bacteriocins"%self.root
                self.genome_files = test_modules.getFNA(self.exampledir)
                self.bacteriocins = "%s/bagel.fa"%self.bacdir
                self.intergenes = "test_intergenes.fa"
                self.annotated_genes = "test_genes.fa"
                self.intermediate = "intermediate"
                self.training_labels = "%s/data/training/training.txt"%self.root
                if not os.path.exists(self.intermediate):
                    os.mkdir(self.intermediate)
                self.bac_evalue = 0.000001
                self.formatdb = True
                self.bacteriocin_radius = 50000
                self.verbose = True
                self.keep_tmp = False
                self.similarity = 0.65
                self.numThreads = 1
                self.output = "out"
                self.keep_tmp = True
            def testrun(self):
                print "Test Run"
                proc = PipelineHandler(  self.root,
                                         self.genome_files,
                                         self.intergenes,
                                         self.annotated_genes,
                                         self.bacteriocins,
                                         self.bacteriocin_radius,
                                         self.similarity,
                                         self.bac_evalue,
                                         self.training_labels,
                                         self.intermediate,
                                         self.output,
                                         self.numThreads,
                                         self.formatdb,
                                         self.verbose                
                                        ) 
                
                proc.preprocess()
                self.assertTrue(os.path.getsize(self.annotated_genes) > 0)
                self.assertTrue(os.path.getsize(self.intergenes) > 0)
                
                #proc.blast()
                
                self.assertTrue(os.path.getsize(proc.blasted_tab_bacteriocins) > 0)
                self.assertTrue(os.path.getsize(proc.cand_context_genes_tab) > 0)       
                
                #proc.cluster()
                
                self.assertTrue(os.path.getsize("%s"%(proc.cand_context_genes_fasta)) > 0)                
                self.assertTrue(os.path.getsize(proc.cand_context_cluster) > 0)
                
                proc.load()
                proc.naiveBayes()
                self.assertTrue(os.path.getsize(proc.textout) > 0)
                
        unittest.main()


