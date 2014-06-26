""" 
A pipeline to launch parallelize computation on the 
quorum computing cluster
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
import tempfile

base_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
for directory_name in os.listdir(base_path):
    site.addsitedir(os.path.join(base_path, directory_name))

import genbank
import blast
import intergene
import genome
import annotated_genes
import bacteriocin
import nbayes
import rforests
import cdhit
import fasta
from quorum import *

class QuorumPipelineHandler(object):
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
                 training_directory,
                 text_database,
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
        self.training_directory =     training_directory
        self.text_database      =     text_database			
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
        self.output = "%s/blasted"%self.intermediate
        self.blasted_tab_bacteriocins = "%s/blasted_bacteriocins.txt"%self.intermediate
        self.blasted_fasta_bacteriocins = "%s/blasted_bacteriocins.fa"%self.intermediate
        self.cand_context_genes_tab = "%s/cand_context_genes.txt"%self.intermediate
        self.cand_context_genes_fasta = "%s/cand_context_genes.fa"%self.intermediate
        self.cand_context_cluster = "%s/cand_context_cluster"%self.intermediate
        self.classifier_out = "%s/classify"%self.intermediate
        #Declare object handlers
        self.clusterer = None                       
        self.textClassifier = None           
        self.nbpickle = "nb.zip"
        self.clusterpickle = "cluster.zip"
        self.textout = "text_out.txt"
        self.jobs  = []
        self.split_files = [] 
        
    def cleanup(self):
        print "Cleaning up"
        for job in self.jobs:
            job.erase_files()
        for fname in self.split_files:
            os.remove(fname) 
    """ Builds database such as the 
        intergenic database
        annotated genes database 
        naive bayes model"""
    def preprocess(self,root,buildAnnotations=True):
        print "Preprocessing"
        if buildAnnotations:
            annotated_genes.go(root,self.annotated_genes) 
            intergene.go(root,self.intergenes)
        
        
    """  Runs blast to identify bacteriocins and context genes"""
    def blast(self,njobs=1):
        print "Blasting"
        
        """ First split up the main bacteriocin file into a bunch of smaller files"""
        split_bacfiles = ["%s/bacteriocin.%d"%(self.intermediate,i)
                          for i in xrange(njobs)]
        self.split_files = split_bacfiles
        split_bachandles = [open(f,'w') for f in split_bacfiles]
        out_fnames = ["%s/blasted.%d"%(self.intermediate,i) for i in xrange(njobs)]
        out_bac = ["%s.bacteriocins.txt"%(out) for out in out_fnames]
        out_genes = ["%s.annotated.txt"%(out) for out in out_fnames]
        index=0
        for record in SeqIO.parse(self.bacteriocins,"fasta"):
            split_bachandles[index].write(">%s\n%s\n"%(str(record.id),
                                                       str(record.seq)))
            index=(index+1)%njobs
        #Close files
        for handle in split_bachandles: handle.close()
            
        blast_cmd = ' '.join([
                    """module load anaconda; module load blast;"""
                    """python %s/src/genome/bacteriocin.py  """,
                    """ --genome-files %s        """,
                    """ --annotated-genes=%s     """,
                    """ --intergenes=%s          """,
                    """ --bacteriocins=%s   	 """,
                    """ --bacteriocin-radius=%d  """,
                    """ --bac-evalue=%f 		 """,
                    """ --num-threads=%d    	 """,
                    """ --intermediate=%s   	 """,
                    """ --output=%s     		 """,
                    """ --formatdb               """,
                    """ --verbose                """
                    ])
        
        """ Release jobs """
        jobs = []
        for i in xrange(njobs):
            cmd = blast_cmd%(self.rootdir,
                            " ".join(map(str,self.genome_files)),
                            self.annotated_genes,
                            self.intergenes,
                            split_bacfiles[i],
                            self.bacteriocin_radius,
                            self.bac_evalue,
                            self.numThreads,
                            self.intermediate,
                            out_fnames[i])
            
            batch_file = "%s/blast%i.%d.job"%(os.getcwd(),i,os.getpid())
            proc = Popen(cmd,shell=True,batch_file=batch_file,
                         stdin=PIPE,stdout=PIPE )
            proc.submit()
            #proc.output = out_fnames[i]
            jobs.append(proc)
            self.jobs.append(proc)
        
        for job in jobs: job.wait() 
        
        """ Collect all of the results from the jobs"""
        bacteriocins_out = open(self.blasted_tab_bacteriocins,'w')
        context_genes_out = open(self.cand_context_genes_tab,'w')
        out_bac = ["%s.bacteriocins.txt"%(out) for out in out_fnames]
        out_genes = ["%s.annotated.txt"%(out) for out in out_fnames]
        for i in xrange(njobs):
            shutil.copyfileobj(open(out_bac[i]),bacteriocins_out)
            shutil.copyfileobj(open(out_genes[i]),context_genes_out)
        bacteriocins_out.close()
        context_genes_out.close()
        
        
    """ Clusters bacteriocins and context genes together"""
    def cluster(self,preprocess=True,numThreads=8,mem=6000):
        print "Clustering"
        if preprocess:
            fasta.preprocess(self.cand_context_genes_tab,
                            self.cand_context_genes_fasta)
            
        
        self.clusterer = cdhit.CDHit(self.cand_context_genes_fasta,
                                     self.cand_context_cluster,
                                     self.similarity,
                                     numThreads,mem)
        self.clusterer.run()
        self.clusterer.parseClusters()
        self.clusterer.dump(self.clusterpickle)
        
    """ Classifies individual bacteriocins and context genes based on their text"""
    def textmine(self,njobs=1):
        
        print "Classifying"
        """ First split up the main bacteriocin file into a bunch of smaller files"""
        split_fastafiles = ["%s/cluster.%d"%(self.intermediate,i)
                          for i in xrange(njobs)]
        self.split_files += split_fastafiles
        split_fastahandles = [open(f,'w') for f in split_fastafiles]
        out_classes = ["%s/classify.%d"%(self.intermediate,i) for i in xrange(njobs)]
        
        index=0
        for record in SeqIO.parse(self.cand_context_cluster,"fasta"):
            split_fastahandles[index].write(">%s\n%s\n"%(str(record.id),
                                                       str(record.seq)))
            index=(index+1)%njobs
        #Close files
        for handle in split_fastahandles: handle.close()
        classify_cmd = ' '.join([
                                 """module load anaconda; module load blast;""",
                                 """python %s/src/classify/rforests.py""",
                                 """--training-directory=%s""",
                                 """--training-labels=%s""",
                                 """--text-database=%s""",
                                 """--num-trees=1000""",
                                 """--fasta=%s""",
                                 """--output=%s"""         
                                 ])    
        print out_classes
        
        """ Release jobs """
        jobs = []
        for i in xrange(njobs):
            cmd = classify_cmd%(self.rootdir,
                                self.training_directory,
                                self.training_labels,
                                self.text_database,
                                split_fastafiles[i],
                                out_classes[i]
                                )
            
            batch_file = "%s/classify%i.%d.job"%(os.getcwd(),i,os.getpid())
            
            proc = Popen( cmd,shell=True,batch_file=batch_file,stdin=PIPE,stdout=PIPE )
            proc.submit()
            #proc.output = out_classes[i]
            jobs.append(proc)
            self.jobs.append(proc)
        
        for job in jobs: job.wait()
        
        """ Collect all of the results from the jobs"""
        classify_out = open(self.classifier_out,'w')
        for i in xrange(njobs):
            shutil.copyfileobj(open(out_classes[i]),classify_out)
        classify_out.close()
        
        
        
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
        '--text-database', type=str, required=False,
        help='SQL database containing text annotations')
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
        proc = QuorumPipelineHandler(args.root_dir,
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
                self.textdb = "%s/db/bacteria_database"%(self.root)
                #self.intergenes = "%s/db/intergenes.txt"%(self.root)
                #self.annotated_genes = "%s/db/annotated_genes.txt"%(self.root)
                self.intermediate = "intermediate"
                #self.training_labels = "%s/data/training/training.txt"%self.root
                self.training_directory = "%s/data/training/protein"%self.root
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
                self.proc = None
            def tearDown(self):
                #self.proc.cleanup()
                pass
            def testrun(self):
                print "Test Run"
                self.proc = QuorumPipelineHandler(  
                                         self.root,
                                         self.genome_files,
                                         self.intergenes,
                                         self.annotated_genes,
                                         self.bacteriocins,
                                         self.bacteriocin_radius,
                                         self.similarity,
                                         self.bac_evalue,
                                         self.training_labels,
                                         self.training_directory,
                                         self.textdb,
                                         self.intermediate,
                                         self.output,
                                         self.numThreads,
                                         self.formatdb,
                                         self.verbose                
                                        ) 
                
                #self.proc.preprocess(root=self.exampledir,buildAnnotations=True)
                self.assertTrue(os.path.getsize(self.annotated_genes) > 0)
                self.assertTrue(os.path.getsize(self.intergenes) > 0)
                
                self.proc.blast(njobs=2)
                
                self.assertTrue(os.path.getsize(self.proc.blasted_tab_bacteriocins) > 0)
                self.assertTrue(os.path.getsize(self.proc.cand_context_genes_tab) > 0)
                
                self.proc.cluster()
                
                self.assertTrue(os.path.getsize("%s"%(self.proc.cand_context_genes_fasta)) > 0)                
                self.assertTrue(os.path.getsize(self.proc.cand_context_cluster) > 0)
                
                self.proc.textmine(njobs=10)
                
                self.assertTrue(os.path.getsize(self.proc.classifier_out ) > 0)
                self.proc.cleanup()
                
        unittest.main()       
        
        
        
        
        
