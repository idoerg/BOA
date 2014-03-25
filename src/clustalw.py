
"""
Perform clustalw multiple alignment on a single file
"""
import sys,os, subprocess
from Bio import AlignIO
from Bio.Align.Applications import ClustalwCommandline

class ClustalW(object):
    def __init__(self,input_file,output_file):
        self.input = input_file
        self.output = output_file
        basename,_ = os.path.splitext(input_file)
        #output from clustalw
        self.aln = "%s.aln"%basename
        self.dnd = "%s.dnd"%basename
    def run(self):
        cline = ClustalwCommandline(infile=self.input)
        print cline
        child = subprocess.Popen(str(cline),
                                 stdout=subprocess.PIPE,
                                 shell=True)
    def parse(self):
        align = AlignIO.read(self.aln, "clustal")
        AlignIO.write([align], open(self.output, 'w+'), 'stockholm')
