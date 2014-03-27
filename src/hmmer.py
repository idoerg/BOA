"""
A wrapper script for hmmsearch, hmmbuild and the 
multiple alignment step using ClustalW
"""

import numpy
import re
import os,sys
from Bio import SeqIO
from cdhit import CDHit
from clustalw import ClustalW
import argparse
import subprocess
cluster_id_reg = re.compile(">(\S+)")

"""Runs HMMER on a single cluster"""
class HMM(object):
    
    def __init__(self,clusterfa):
        directory = os.path.dirname(clusterfa)
        basename = os.path.splitext(os.path.basename(clusterfa))[0]
        self.clustal = None        
        self.clusterfa = clusterfa
        self.sto = "%s/%s.sto"%(directory,basename)
        self.hmm = "%s/%s.hmm"%(directory,basename)
        self.table = "%s/%s.table"%(directory,basename) #results from hmmsearch
        self.logs = ["%s/%s_hmmbuild.log"%(directory,basename),
                     "%s/%s_hmmsearch.log"%(directory,basename)]
        self.multipleAlignment()
        self.hmmbuild()
    """Clean up useless files"""        
    def cleanUp(self):
        self.clustal.cleanUp()
        os.remove(self.clusterfa)
        os.remove(self.sto)
        os.remove(self.hmm)
        os.remove(self.table)
        for log in self.logs:
            os.remove(log)
    """Runs Viterbi algorithm on an input fasta"""
    def hmmsearch(self,infasta):
        cmd = "hmmsearch --noali --domtblout %s %s %s"%(self.table,self.hmm,infasta)
        proc = subprocess.Popen(cmd,stderr=open(self.logs[0],'w+'),shell=True)
        proc.wait()
        
    """Builds an HMM for a cluster"""
    def hmmbuild(self):
        cmd = "hmmbuild %s %s "%(self.hmm,self.sto)
        proc = subprocess.Popen(cmd,stderr=open(self.logs[1],'w+'),shell=True)
        proc.wait()
        
    """Perform multiple alignment with ClustalW on a single cluster"""
    def multipleAlignment(self):
        cw = ClustalW(self.clusterfa,self.sto)
        cw.run()
        cw.outputSTO()
        self.clustal = cw
"""Performs clustering and runs HMMER on every cluster"""
class HMMER(object):
    def __init__(self,fasta,minClusters=2):
        directory = os.path.dirname(fasta)
        basename = os.path.splitext(os.path.basename(fasta))[0]
        self.clrreps = "%s/%s_cluster"%(directory,basename)
        self.cluster = "%s.clstr"%(self.clrreps)
        self.fasta = fasta
        self.clusterfas = []
        self.hmms = []
        self.tables = []
        self.minClusters = minClusters #The mininum number of members in a cluster
        
    """Remove all useless files"""
    def cleanUp(self):
        for hmm in self.hmms:
            hmm.cleanUp()
        os.remove(self.clrreps)
        os.remove(self.cluster)
        
    """Get file names for all of the cluster files"""
    def getClusters(self):
        return self.clusterfas
    """Spawn hmms from each cluster"""
    def HMMspawn(self):
        for clrfa in self.clusterfas:
            hmm = HMM(clrfa)
            self.hmms.append(hmm)
    """Performs HMMER using all clusters on infasta"""
    def search(self,infasta,out):
        for hmm in self.hmms:
            hmm.hmmsearch(infasta)
            self.tables.append(hmm.table)
        with open(out, 'w') as outfile:
            for fname in self.tables:
                with open(fname) as infile:
                    for line in infile:
                        outfile.write(line)
                    
    """Obtains ids for all of the sequences and write them into separate cluster files"""
    def writeClusters(self):
        clusterProc = CDHit(self.fasta,self.clrreps,0.7)
        clusterProc.run()
        clusterProc.parseClusters()    
        i = 0
        infasta = clusterProc.input    
        directory = os.path.dirname(infasta)
        record_dict = SeqIO.to_dict(SeqIO.parse(infasta,'fasta'))    
        for cluster in clusterProc.clusters:
            if len(cluster.seqs)<self.minClusters:#filter out small clusters
                continue
            outfile = '%s/cluster%d.fa'%(directory,i)
            self.clusterfas.append(outfile)
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
    parser = argparse.ArgumentParser(description=\
        'Runs HMMER on CDHit clusters')
    parser.add_argument(\
        '--fasta', type=str, required=False,
        help='Input sequences to cluster')
    parser.add_argument(\
        '--test', action='store_const', const=True, default=False,
        help='Run unittests')
    args = parser.parse_args()
    if not args.test:
        directory = os.path.dirname(args.fasta)
        basename = os.path.basename(args.fasta)
        outputfile = "%s/%s.fa"%(directory,basename)
        go(args.fasta,outputfile)
    else:
        del sys.argv[1:]
        import unittest
                
        class TestHMM(unittest.TestCase):
            def setUp(self):
                seqs = ['>20.1',
                        'CKQSCSFGPFTFVCDGNTK',
                        '>21.1',
                        'CRQSCSFGPLTFVCDGNTK',
                        '>22.1', 
                        'CANSCSYGPLTWSCDGNTK']
                self.fasta = "./test.fa"
                self.genome = "./testGenome.fa"
                self.results = './results.txt'
                handle = open(self.fasta,'w')
                handle.write('\n'.join(seqs))
                handle.close()
                genome = ['>genome',
                          'CKQSCSFGPFTFVCDGNTK',
                          'CRQSCSFGPLTFVCDGNTK',
                          'CANSCSYGPLTWSCDGNTK']
                handle = open(self.genome,'w')
                handle.write('\n'.join(genome))
                handle.close()
                self.testhmmer = HMMER(self.fasta)
             
            def tearDown(self):
                os.remove(self.fasta)
                os.remove(self.genome)
                self.testhmmer.cleanUp()
                self.testhmm.cleanUp()
            
            def testHMMSearch(self):
                self.testhmmer.writeClusters()
                clrfnames = self.testhmmer.getClusters()
                testclr = clrfnames[0]
                self.testhmm = HMM(testclr)
                self.testhmm.multipleAlignment()
                self.assertTrue(os.path.exists("cluster0.aln"))
                self.assertTrue(os.path.exists("cluster0.sto"))
                self.testhmm.hmmbuild()
                self.assertTrue(os.path.exists("cluster0.hmm"))
                self.testhmm.hmmsearch(self.genome)
                self.assertTrue(os.path.exists("cluster0.table"))
        
        
        class TestHMMER(unittest.TestCase):

            def setUp(self):
                seqs = ['>20.1',
                        'CKQSCSFGPFTFVCDGNTK',
                        '>21.1',
                        'CRQSCSFGPLTFVCDGNTK',
                        '>22.1', 
                        'CANSCSYGPLTWSCDGNTK']
                self.fasta = "./test.fa"
                self.genome = "./testGenome.fa"
                handle = open(self.fasta,'w')
                handle.write('\n'.join(seqs))
                handle.close()
                genome = ['>genome',
                          'CKQSCSFGPFTFVCDGNTK',
                          'CRQSCSFGPLTFVCDGNTK',
                          'CANSCSYGPLTWSCDGNTK']
                handle = open(self.genome,'w')
                handle.write('\n'.join(genome))
                handle.close()
                self.testhmm = HMMER(self.fasta)
             
            def tearDown(self):
                os.remove(self.fasta)
                os.remove(self.genome)
                self.testhmm.cleanUp()
            
            def testSearch(self):
                self.testhmm.writeClusters()
                self.assertTrue(os.path.exists("test_cluster"))
                self.assertTrue(os.path.exists("test_cluster.clstr"))
                self.assertTrue(os.path.exists("cluster0.fa"))
                self.testhmm.HMMspawn()
                self.assertTrue(os.path.exists("cluster0.hmm"))
                self.testhmm.search(self.genome,"results.txt")
                self.assertTrue(os.path.exists("cluster0.table"))
                self.assertTrue(os.path.exists("results.txt"))
                
        unittest.main()
    
    
    
    
    
    