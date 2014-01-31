"""
Given a directory, walks into all of the subdirectorys and looks for intergenic regions
"""
import sys
import os
import site
import argparse
import string
import numpy
import re
import subprocess
import intergene
parser = argparse.ArgumentParser(description=\
    'Finds intergenic regions from genback files')
parser.add_argument(\
    '--root-dir', type=str,required=True,default="",
    help='Root directory of all of the files of interest')
parser.add_argument(\
    '--output-file', type=str, required=True,
    help='The output file containing the BLAST output')
args = parser.parse_args()


def go():
    outHandle = open(args.output_file,'w')
    for root, subFolders, files in os.walk(args.root_dir):
        for fname in files:
            genome_files = []
            ext = os.path.splitext(os.path.basename(fname))[1]
            if ext==".gbk":
                absfile=os.path.join(root,fname)
                intergene.get_interregions(absfile,args.output_file,
                                           intergene_length=1,
                                           write_flag='a')
                
if __name__=="__main__":
    go()
