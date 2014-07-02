"""
A module for finding context genes associated with bacteriocins
using a training dataset and blast

Output
1.  bacteriocin ID 
2.  organism 
3.  bacteriocin start 
4.  bacteriocin end 
5.  bacteriocin strand  
6.  annotated gene organism 
7.  annotated gene locus 
8.  annotated gene protein id 
9.  annotated gene start 
10. annotated gene end 
11. annotated gene strand 
12. context gene organism
13. context gene protein id
14. context gene class label
"""
import os,sys,site,shutil
base_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
for directory_name in os.listdir(base_path):
    site.addsitedir(os.path.join(base_path, directory_name))
import blast
import hmmer
import training
import argparse 
import fasta
from Bio import Seq
from Bio import SeqIO
from Bio import SeqRecord
from Bio.Alphabet import IUPAC
from Bio.SeqRecord import SeqRecord


class ContextGeneBLAST(object):
    def __init__(self,trainingDir,
                 training_labels,
                 intermediate,
                 evalue=0.000001,
                 threads):                 
        self.labels = training.Labels(trainingDir,training_labels)
        self.blast_proc = None
        self.evalue = evalue
        self.intermediate = intermediate
        self.blastdb = "training_%d.fa"%os.getpid()
        self.threads = threads
    def getTrainingFasta(self):
        return self.blastdb
    def cleanup(self):
        os.remove(self.blastdb)
        self.blast_proc.cleanup()
    """ Build blast database and find context genes"""
    def find(self,queryFile,formatdb=True):
        trainingEntries = self.labels.getTrainingSequences()
        handle = open(self.blastdb,'w')
        #with open(self.blastdb,'w') as handle:
        for entry in trainingEntries:
            org,protID,label,seq = entry
            handle.write(">%s|%s|%s\n%s\n"%
                         (org,protID,label,seq))
        handle.close()
        self.blast_proc = blast.BLAST(self.blastdb,
                                      queryFile,
                                      self.intermediate,
                                      self.evalue)
        if formatdb:
            self.blast_proc.buildDatabase(base="protein")
        self.blast_proc.run(blast_cmd='blastp',mode='xml',
                            num_threads=self.threads)
        hits = self.blast_proc.parseBLAST("xml")
        return hits
    def write(self,hits,output):
        #with open(output,'w') as handle:
        handle = open(output,'w')
        for hit in hits:
            handle.write(">%s|%s\n%s\n"%(hit.query_id,
                                         hit.sbjct_id,
                                         hit.sbjct))
            #print hit.sbjct
            #record = SeqRecord(Seq(hit.sbjct,IUPAC.protein),id="%s|%s"%(hit.query_id,hit.sbjct_id))
            #SeqIO.write(record,handle,"fasta")
        handle.close()
            
if __name__=="__main__":
     parser = argparse.ArgumentParser(description=\
        'Finds bacteriocins and context genes')
     parser.add_argument(\
        '--evalue', type=float, required=False, default=0.00001,
        help='The evalue for context gene hits')
     parser.add_argument(\
        '--training-labels',type=str,required=False,default=None,
        help='Training labels for Naive Bayes')
     parser.add_argument(\
        '--training-directory',type=str,required=False,default=None,
        help='''Training directory containing all 
                genbank files required for training''')
     parser.add_argument(\
        '--query', type=str, required=False,default='.',
        help='Query files used to find context genes')
     parser.add_argument(\
        '--output', type=str, required=False,default='.',
        help='Output context genes')
     parser.add_argument(\
        '--intermediate', type=str, required=False,default='.',
        help='Directory for storing intermediate files')
     parser.add_argument(\
        '--test', action='store_const', const=True, default=False,
        help='Run unittests')
     args = parser.parse_args()
     if not args.test:
         cg = ContextGeneBLAST(args.training_directory,
                               args.training_labels,
                               args.intermediate,
                               args.evalue)
         hits = cg.find(args.query)
         cg.write(hits,args.output)
         cg.cleanup()
     else:
         del sys.argv[1:]
         import unittest
         import test_fasta
         class TestContextGeneBLAST(unittest.TestCase):
             def setUp(self):
                 self.cluster = test_fasta.cluster
                 self.root = os.environ['BACFINDER_HOME']
                 self.trainDir = "%s/data/training/protein"%self.root
                 self.labelFile = "%s/data/training/training_proteins.txt"%self.root
                 self.intermediate="./intermediate"
                 if not os.path.exists(self.intermediate):
                     os.mkdir(self.intermediate)
                 self.evalue=0.000001
             def tearDown(self):
                 shutil.rmtree(self.intermediate)
             def test1(self):
                 cg = ContextGeneBLAST(self.trainDir,
                                       self.labelFile,
                                       self.intermediate,
                                       self.evalue)
                 hits = cg.find(cg.getTrainingFasta())
                 self.assertTrue(len(hits)>0)
                 cg.cleanup()
                 
         unittest.main()
         
         
         
         
         
         
         
         
         
         
         
         
         
         
         
         
    