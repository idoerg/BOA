"""
Perform mafft multiple alignment on a single file
"""
import sys,os, subprocess
from Bio import AlignIO
import os,site,sys
base_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
for directory_name in os.listdir(base_path):
    site.addsitedir(os.path.join(base_path, directory_name))
import quorum
import glob

class MAFFT(object):
    def __init__(self,input_file,output_file):
        self.input = input_file #input fasta file
        self.output = output_file
        basename,_ = os.path.splitext(input_file)
        #output from clustalw
        self.aln = "%s.aln"%basename
        
    def run(self,fasta=False,module=subprocess,maxiters=4,threads=1):
        
        if fasta:
            cline = "mafft --maxiterate %d --thread %d %s"%(maxiters,threads,self.input)
            out = self.output
        else:
            cline = "mafft --maxiterate %d --thread %d --clustalout %s"%(maxiters,threads,self.input)
            out = self.aln
        print cline
        child = module.Popen(str(cline),
                             stdout=open(out,'w'),
                             shell=True)
        if module==quorum: child.submit()
        return child
        #child.wait()
    def outputSTO(self):
        handle = open(self.output, 'w+')
        align = AlignIO.read(self.aln, "clustal")
        AlignIO.write([align], handle, 'stockholm')
        handle.close()
    def cleanUp(self):
        if os.path.exists(self.aln):
            os.remove(self.aln)
        
        
if __name__=='__main__':
    import unittest
    class TestRun(unittest.TestCase):
        def setUp(self):
            self.infile = "test.fa"
            entries = ['>testseq1',
                       'AGCTACT',
                       '>testseq2',
                       'AGCTAGCT']
            open(self.infile,'w').write('\n'.join(entries))
            self.outsto = "out.sto"
            self.outfasta = "out.fasta"
            basename,_ = os.path.splitext(self.infile)
            #output from clustalw
            self.aln = "%s.aln"%basename
            
        def tearDown(self):
            #if os.path.exists(self.infile):
            #    os.remove(self.infile)
            #if os.path.exists(self.infile):
            #    os.remove(self.aln)
            pass
        def testRunSto(self):
           cw = MAFFT(self.infile,self.outsto)
           proc = cw.run()
           proc.wait()
           cw.outputSTO()
           handle = open(self.outsto,'r')
           lines = handle.readlines()
           line1 = lines[0].rstrip()
           self.assertEquals(line1,"# STOCKHOLM 1.0")
           os.remove(self.outsto)
           
    class TestQuorum(unittest.TestCase):
        def setUp(self):
            self.infile = "test.fa"
            entries = ['>testseq1',
                       'AGCTACT',
                       '>testseq2',
                       'AGCTAGCT']
            open(self.infile,'w').write('\n'.join(entries))
            self.outsto = "out.sto"
            self.outfasta = "out.fasta"
            basename,_ = os.path.splitext(self.infile)
            #output from clustalw
            self.aln = "%s.aln"%basename
            
        def tearDown(self):
            if os.path.exists(self.infile):
                os.remove(self.infile)
            if os.path.exists(self.infile):
                os.remove(self.aln)
            for file in glob.glob("tmp*"):
                os.remove(file)
        def testQuorum(self):
           cw = MAFFT(self.infile,self.outfasta)
           proc = cw.run(fasta=True,module=quorum)
           proc.wait()
           result = proc.ofile_string()
           print result
           self.assertEquals(result[0],">")
           
    unittest.main()
    
    
    
