"""
Perform clustalw multiple alignment on a single file
"""
import sys,os, subprocess
from Bio import AlignIO
import os,site,sys
base_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
for directory_name in os.listdir(base_path):
    site.addsitedir(os.path.join(base_path, directory_name))
import quorum
class Muscle(object):
    def __init__(self,input_file,output_file):
        self.input = input_file #input fasta file
        self.output = output_file
        basename,_ = os.path.splitext(input_file)
        #output from clustalw
        self.aln = "%s.aln"%basename
        
        
    def run(self,fasta=False,module=subprocess,maxiters=4,maxhours=12):
        if fasta:
            cline = "muscle -in %s -out %s -maxiters %d -maxhours %d"%(self.input,self.output,maxiters,maxhours)
        else:
            cline = "muscle -in %s -out %s -maxiters %d -maxhours %d -clw"%(self.input,self.aln,maxiters,maxhours)
        print cline
        child = module.Popen(str(cline),
                             #stdout=subprocess.PIPE,
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
            os.remove(self.infile)
            os.remove(self.aln)
            
        def testRunSto(self):
           cw = Muscle(self.infile,self.outsto)
           cw.run()
           cw.outputSTO()
           handle = open(self.outsto,'r')
           lines = handle.readlines()
           line1 = lines[0].rstrip()
           self.assertEquals(line1,"# STOCKHOLM 1.0")
           os.remove(self.outsto)
        def testRunFasta(self):
           cw = Muscle(self.infile,self.outfasta)
           cw.run()
           cw.outputFASTA()
           handle = open(self.outfasta,'r')
           lines = handle.readlines()
           line1 = lines[0].rstrip()
           self.assertEquals(line1[0],">")
           os.remove(self.outfasta)
           
        def testQuorum(self):
           cw = Muscle(self.infile,self.outfasta)
           cw.run(quorum)
           cw.outputFASTA()
           handle = open(self.outfasta,'r')
           lines = handle.readlines()
           line1 = lines[0].rstrip()
           self.assertEquals(line1[0],">")
           os.remove(self.outfasta)
         
    unittest.main()
    
    
    
