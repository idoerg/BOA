"""
A wrapper script for hmmsearch, hmmbuild and the 
multiple alignment step using ClustalW

TODO: fix the unittest tests here
"""

import numpy
import re
import os,sys, site


from Bio import SeqIO
from cdhit import CDHit
#from clustalw import ClustalW
from muscle import Muscle
from mafft import MAFFT
import argparse
import subprocess
cluster_id_reg = re.compile(">(\S+)")

base_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
for directory_name in os.listdir(base_path):
    site.addsitedir(os.path.join(base_path, directory_name))
import quorum
import fasta

"""Runs HMMER on a single cluster"""
class HMM(object):
    """ Make import variable. Could either import subprocess or quorum """
    def __init__(self,
                 clusterfa,
                 module=subprocess,
                 threads=8,
                 reverseComplement=False):
        self.module=module
        directory = os.path.dirname(clusterfa)
        #basename = os.path.splitext(os.path.basename(clusterfa))[0]
        basename = os.path.basename(clusterfa)
        self.clustal = None        
        self.clusterfa = clusterfa
        self.sto = "%s/%s.sto"%(directory,basename)
        self.hmm = "%s/%s.hmm"%(directory,basename)
        self.table = "%s/%s.table"%(directory,basename) #results from hmmsearch
        self.logs = ["%s/%s_hmmbuild.log"%(directory,basename),
                     "%s/%s_hmmsearch.log"%(directory,basename)]
        self.threads = threads
    """Clean up useless files"""        
    def cleanUp(self):
        self.clustal.erase_files()
        #if os.path.exists(self.clusterfa): os.remove(self.clusterfa)
        if os.path.exists(self.sto): os.remove(self.sto)
        if os.path.exists(self.hmm): os.remove(self.hmm)
        if os.path.exists(self.table): os.remove(self.table)
        for log in self.logs:
            if os.path.exists(log): os.remove(log)
    """Runs Viterbi algorithm on an input protein fasta"""
    def hmmsearch(self,infasta,maxpower=False):
        if maxpower:
            cmd = "hmmsearch --noali --notextw --max --domtblout %s %s %s"%(self.table,self.hmm,infasta)
        else:
            cmd = "hmmsearch --noali --notextw --domtblout %s %s %s"%(self.table,self.hmm,infasta)
        print cmd
        proc = self.module.Popen(cmd,stderr=open(self.logs[0],'w+'),shell=True,threads=self.threads)
        if self.module==quorum: proc.submit()
        return proc
    """Runs Viterbi algorithm on an input nucleotide fasta"""
    def nhmmsearch(self,infasta):
        cmd = "nhmmer --noali --domtblout --max %s %s %s"%(self.table,self.hmm,infasta)
        proc = self.module.Popen(cmd,stderr=open(self.logs[0],'w+'),shell=True,threads=self.threads)
        if self.module==quorum: proc.submit()
        return proc 
    """Builds an HMM for a cluster"""
    def hmmbuild(self,module=subprocess):
        cmd = "hmmbuild %s %s "%(self.hmm,self.sto)
        proc = self.module.Popen(cmd,stderr=open(self.logs[1],'w+'),shell=True,threads=self.threads)
        if self.module==quorum: proc.submit()
        #proc.wait()
        return proc
    """Perform multiple alignment with Muscle on a single cluster"""
    def multipleAlignment(self,module=subprocess,msa=MAFFT,maxiters=20):
        cw = msa(self.clusterfa,self.sto,module)
        self.clustal = cw
        if msa==Muscle:
            proc = cw.run(fasta=True,maxiters=maxiters)
        else:
            proc = cw.run(fasta=True,maxiters=maxiters,threads=self.threads)
        #cw.outputSTO()
        return cw
        
"""Performs clustering and runs HMMER on every cluster"""
class HMMER(object):
    def __init__(self,fasta="",
                 module=subprocess,
                 threads=8,
                 minClusters=2):
        self.module = module
        directory = os.path.dirname(fasta)
        #basename = os.path.splitext(os.path.basename(fasta))[0]
        self.basename = os.path.basename(fasta)
        self.clrreps = "%s/%s_cluster"%(directory,self.basename)
        self.cluster = "%s.clstr"%(self.clrreps)
        self.fasta = fasta
        self.clusterfas = []
        self.hmms = []
        self.tables = []
        self.minClusters = minClusters #The mininum number of members in a cluster
        self.threads = threads
    """Remove all useless files"""
    def cleanUp(self):
        for hmm in self.hmms:
            hmm.cleanUp()
        if os.path.exists(self.clrreps): os.remove(self.clrreps)
        if os.path.exists(self.cluster): os.remove(self.cluster)
        
    """Get file names for all of the cluster files"""
    def getClusters(self):
        return self.clusterfas
    """Wait for all jobs to complete"""
    def wait(self,procs,fasta=False):
        for p in procs: 
            p.wait()
            if fasta: 
                #p.outputFASTA()
                assert os.path.exists(p.output)  
            if self.module==quorum: p.erase_files()
        
    """Spawn hmms from each cluster"""
    def HMMspawn(self,msa=MAFFT,njobs=4,maxiters=20):
        procs = []
        i = 0
        for clrfa in self.clusterfas:
            hmm = HMM(clrfa,module=self.module,threads=self.threads)
            self.hmms.append(hmm)
        for hmm in self.hmms:
            proc = hmm.multipleAlignment(msa=msa,maxiters=maxiters)
            procs.append(proc)
            i+=1
            if i==njobs: #make sure jobs don't overload
                self.wait(procs,fasta=True)
                procs,i = [],0
        self.wait(procs,fasta=True)
        
        procs,i = [],0
        for hmm in self.hmms:
            proc = hmm.hmmbuild()
            procs.append(proc)
            i+=1
            if i==njobs: #make sure jobs don't overload
                self.wait(procs)
                procs,i = [],0
        self.wait(procs)
        
            
            
    """Performs HMMER using all clusters on infasta"""
    def search(self,infasta,out,njobs=4):
        procs = []
        i = 0
        for hmm in self.hmms:
            proc = hmm.hmmsearch(infasta)
            self.tables.append(hmm.table)
            procs.append(proc)
            i+=1
            if i==njobs: #make sure jobs don't overload
                self.wait(procs)
                procs,i = [],0
        self.wait(procs)
        procs,i = [],0     
            
        with open(out, 'w') as outfile:
            for fname in self.tables:
                if os.path.exists(fname):
                    with open(fname) as infile:
                        for line in infile:
                            outfile.write(line)
                        
    """Obtains ids for all of the sequences and write them into separate cluster files"""
    def writeClusters(self,similarity=0.7,memory=800):
        clusterProc = CDHit(self.fasta,self.clrreps,similarity,self.threads,memory)
        clusterProc.run()
        clusterProc.parseClusters()    
        i = 0
        infasta = clusterProc.input    
        directory = os.path.dirname(infasta)
        record_dict = SeqIO.to_dict(SeqIO.parse(infasta,'fasta'))    
        
        for cluster in clusterProc.clusters:
            print cluster
            if len(cluster.seqs)<self.minClusters:#filter out small clusters
                continue
            outfile = '%s/%s.cluster%d.fa'%(directory,self.basename,i)
            self.clusterfas.append(outfile)
            handle = open(outfile,'w')
            for subc in cluster.seqs:
                subtitle = cluster_id_reg.findall(subc)[0][:-3]
                print subtitle
                record = record_dict[subtitle]
                handle.write(">%s\n"%record.id)
                handle.write("%s\n"%str(record.seq))
            handle.close()
            i+=1
""" Parse the output of hmmsearch output"""            
def parse(hmmerout):
    entries = []
    with open(hmmerout,'r') as handle:
        for ln in handle:
            if ln[0]=="#": continue
            ln = ln.rstrip()
            toks = re.split("\s+",ln)
            target,query = toks[0],toks[3]
            try:
                #full_evalue = float(toks[6])
                score = float(toks[13]) 
                hmm_st,hmm_end,env_st,env_end = map(int,[ toks[15],toks[16],toks[19],toks[20] ]) 
                description = ' '.join(toks[22:])
                entries.append( (
                       target,query,
                       score,
                       hmm_st,hmm_end,env_st,env_end,
                       description 
                       ) )
            except Exception as e:
                #Ignore HMMER dumbass parsing exceptions
                continue
    return entries

def parse_scores(hmmerout):
    scores = []
    with open(hmmerout,'r') as handle:
        for ln in handle:
            if ln[0]=="#": continue
            ln = ln.rstrip()
            toks = re.split("\s+",ln)
            try:
                score = float(toks[6])
                scores.append(score)
            except Exception as e:
                    continue
    return scores
    
def hmmerstr(hits):
    all_hits = []
    for hit in hits:
        s = '|'.join(map(str,hit))
        all_hits.append(s)
    return '\n'.join(all_hits)
"""First separate all of the sequences into clusters of 70%"""
def go(trainingSet,genome,results):
    hmmer = HMMER(trainingSet)
    hmmer.HMMspawn()
    hmmer.search(genome,results)
    hmmer.cleanUp()
    
if __name__=="__main__":
    parser = argparse.ArgumentParser(description=\
        'Runs HMMER on CDHit clusters')
    parser.add_argument(\
        '--training', type=str, required=False,
        help='Input sequences to cluster')
    parser.add_argument(\
        '--genome', type=str, required=False,
        help='Genome to align against')
    parser.add_argument(\
        '--output', type=str, required=False,
        help='All HMMER alignments')
    parser.add_argument(\
        '--test', action='store_const', const=True, default=False,
        help='Run unittests')
    args = parser.parse_args()
    if not args.test:
        directory = os.path.dirname(args.output)
        basename = os.path.basename(args.output)
        go(args.traing,args.genome,args.output)
    else:
        del sys.argv[1:]
        import unittest
        class TestStr(unittest.TestCase):
            def test(self):
                hits = [('CP002279.1_3','toxin.fa.cluster2.fa',0,0,100,25000,25100,
                        'Mesorhizobium opportunistum WSM2075, complete genome')]
                result = hmmerstr(hits)
                self.assertEquals(result,
                                  'CP002279.1_3|toxin.fa.cluster2.fa|0|0|100|25000|25100|Mesorhizobium opportunistum WSM2075, complete genome')
        
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
                proc = self.testhmm.multipleAlignment(msa=Muscle,maxiters=4) 
                proc.wait()
                self.assertTrue(os.path.exists(self.testhmm.clusterfa))
                self.assertTrue(os.path.getsize(self.testhmm.sto)>0)
                proc = self.testhmm.hmmbuild()
                proc.wait()
                self.assertTrue(os.path.exists(self.testhmm.hmm))
                proc = self.testhmm.hmmsearch(self.genome)
                proc.wait()
                self.assertTrue(os.path.exists(self.testhmm.table))
            
            def testFFTHMMSearch(self):
                self.testhmmer.writeClusters()
                clrfnames = self.testhmmer.getClusters()
                testclr = clrfnames[0]
                self.testhmm = HMM(testclr)
                proc = self.testhmm.multipleAlignment(msa=MAFFT,maxiters=4)
                proc.wait()
                self.assertTrue(os.path.exists(self.testhmm.clusterfa))
                self.assertTrue(os.path.getsize(self.testhmm.sto)>0)
                proc = self.testhmm.hmmbuild()
                proc.wait()
                self.assertTrue(os.path.exists(self.testhmm.hmm))
                proc = self.testhmm.hmmsearch(self.genome)
                proc.wait()
                self.assertTrue(os.path.exists(self.testhmm.table))
             
        
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
                pass
            def testSearch(self):
                self.testhmm.writeClusters()
                print "Cluster fastas",self.testhmm.clusterfas
                
                self.assertTrue(os.path.exists("test.fa_cluster"))
                self.assertTrue(os.path.exists("test.fa.cluster0.fa"))
               
                self.testhmm.HMMspawn()
                self.assertTrue(os.path.exists("test.fa.cluster0.fa.hmm"))
                self.testhmm.search(self.genome,"results.txt")
                self.assertTrue(os.path.exists("test.fa.cluster0.fa.table"))
                self.assertTrue(os.path.exists("results.txt"))
        
        class TestHMMERQuorum(unittest.TestCase):

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
                self.testhmm = HMMER(self.fasta,quorum)
            
            def tearDown(self):
                os.remove(self.fasta)
                os.remove(self.genome)
                self.testhmm.cleanUp()
                os.remove("quorum_epilogue.py")
                pass
            def testSearch(self):
                self.testhmm.writeClusters()
                self.assertTrue(os.path.exists("test.fa_cluster"))
                self.assertTrue(os.path.exists("test.fa.cluster0.fa"))
               
                self.testhmm.HMMspawn(msa=MAFFT,maxiters=4) 
                
                self.assertTrue(os.path.exists("test.fa.cluster0.fa.hmm"))
                self.testhmm.search(self.genome,"results.txt")
                self.assertTrue(os.path.exists("test.fa.cluster0.fa.table"))
                self.assertTrue(os.path.exists("results.txt"))
            
        
        unittest.main()
    
    
    
    
    
    
