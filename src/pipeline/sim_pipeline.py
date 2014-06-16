

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
import nbayes_sim
import pipeline

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
        '--training-directory',type=str,required=False,default=None,
        help='''Training directory containing all 
                genbank files required for training the Naive Bayes model''')
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
                                 args.training_directory,
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
                #self.training_labels = "%s/data/training/training.txt"%self.root
                
                self.trainDir = "%s/data/training/protein"%self.root
                self.training_labels = "%s/data/training/training_proteins.txt"%self.root
                
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
                
                self.db = "testdb"
                if os.path.exists(self.db):
                    os.remove(self.db)
                
                self.fasta = "test.fa"
                self.function_check = "correct_functions.txt"
            def testNBayes(self):
                nsim = nbayes_sim.NBayesSim(self.trainDir,self.training_labels)
                nsim.simulationOutput(10,self.function_check,self.db,self.fasta)
                
                print "Test Run"
                proc = pipeline.PipelineHandler( self.root,
                                                 self.genome_files,
                                                 self.intergenes,
                                                 self.annotated_genes,
                                                 self.bacteriocins,
                                                 self.bacteriocin_radius,
                                                 self.similarity,
                                                 self.bac_evalue,
                                                 self.training_labels,
                                                 self.trainDir,
                                                 self.intermediate,
                                                 self.output,
                                                 self.numThreads,
                                                 self.formatdb,
                                                 self.verbose                
                                                ) 

                proc.cand_context_genes_fasta = self.fasta
                proc.preprocess(buildAnnotations=False)
                #self.assertTrue(os.path.getsize(self.annotated_genes) > 0)
                #self.assertTrue(os.path.getsize(self.intergenes) > 0)
                
                #proc.blast()
                
                #self.assertTrue(os.path.getsize(proc.blasted_tab_bacteriocins) > 0)
                #self.assertTrue(os.path.getsize(proc.cand_context_genes_tab) > 0)       
                
                #proc.cluster(preprocess=False)
                
                #self.assertTrue(os.path.getsize("%s"%(proc.cand_context_genes_fasta)) > 0)                
                #self.assertTrue(os.path.getsize(proc.cand_context_cluster) > 0)
                
                proc.load()
                p = proc.textClassifier.crossvalidation()
                print "Accuracy:",p
                proc.naiveBayes(self.db)
                self.assertTrue(os.path.getsize(proc.textout) > 0)
    
                rows = genbank.proteinQueryAll(self.db)
                open("testdb.txt",'w').write('\n'.join(map(str,rows)))
        unittest.main()


