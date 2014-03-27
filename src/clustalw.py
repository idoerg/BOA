
"""
Perform clustalw multiple alignment on a single file
"""
import sys,os, subprocess
from Bio import AlignIO
from Bio.Align.Applications import ClustalwCommandline

class ClustalW(object):
    def __init__(self,input_file,output_file):
        self.input = input_file #input fasta file
        self.output = output_file
        basename,_ = os.path.splitext(input_file)
        #output from clustalw
        self.aln = "%s.aln"%basename
        self.dnd = "%s.dnd"%basename
    def run(self):
        cline = ClustalwCommandline(infile=self.input)
        child = subprocess.Popen(str(cline),
                                 stdout=subprocess.PIPE,
                                 shell=True)
        child.wait()
    def outputSTO(self):
        handle = open(self.output, 'w+')
        align = AlignIO.read(self.aln, "clustal")
        AlignIO.write([align], handle, 'stockholm')
        handle.close()
    def cleanUp(self):
        os.remove(self.aln)
        os.remove(self.dnd)
        
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
            self.outfile = "out.sto"
            basename,_ = os.path.splitext(self.infile)
            #output from clustalw
            self.aln = "%s.aln"%basename
            self.dnd = "%s.dnd"%basename 
            
        def tearDown(self):
            os.remove(self.infile)
            os.remove(self.outfile)
            os.remove(self.aln)
            os.remove(self.dnd)
        
        def testrun(self):
           cw = ClustalW(self.infile,self.outfile)
           cw.run()
           cw.outputSTO()
           handle = open(self.outfile,'r')
           lines = handle.readlines()
           line1 = lines[0].rstrip()
           self.assertEquals(line1,"# STOCKHOLM 1.0") 
    unittest.main()
    
    
    