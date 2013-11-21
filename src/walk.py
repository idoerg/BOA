"""
Given a directory, walks into all of the subdirectorys and looks for bacteriocins
"""
import sys
import os
import site
import argparse
import string
import numpy
import re
import subprocess
import bacteriocin

parser = argparse.ArgumentParser(description=\
    'Finds intergenic regions from genback file')
parser.add_argument(\
    '--root-dir', type=str,required=True,default="",
    help='Root directory of all of the files of interest')
parser.add_argument(\
    '--genes', type=str,required=True,default="",
    help='A FASTA file containing all of the target genes of interest')
parser.add_argument(\
    '--bacteriocins', type=str, required=True,
    help='The bacteriocin proteins that are to be blasted')
parser.add_argument(\
    '--radius', type=int, required=False, default=10000,
    help='The search radius around every specified gene')
parser.add_argument(\
    '--bac-evalue', type=float, required=False, default=0.0000000001,
    help='The search radius around every specified gene')
parser.add_argument(\
    '--gene-evalue', type=float, required=False, default=0.00001,
    help='The search radius around every specified gene')
parser.add_argument(\
    '--intermediate', type=str, required=True,
    help='Directory for storing intermediate files')
parser.add_argument(\
    '--output-file', type=str, required=True,
    help='The output file containing the BLAST output')
parser.add_argument(\
    '--num-threads', type=int, required=False, default=1,
    help='The number of threads to be run on blast')
parser.add_argument(\
    '--keep-tmp', action='store_const', const=True, default=False,
    help='Keeps temporary files such as blast database and blast output xml')
parser.add_argument(\
    '--verbose', action='store_const', const=True, default=False,
    help='Messages for debugging')

args = parser.parse_args()
def go():
    outHandle = open(args.output_file,'w')
    outHandle.write("bacteriocin\torganism\torganismSt\torganismEnd\tgene\tgeneSt\tbacterciocin_sequence\n")
    for root, subFolders, files in os.walk(args.root_dir):
        for fname in files:
            genbank_files = []
            ext = os.path.splitext(os.path.basename(fname))[1]
            if ext==".gbk":
                absfile=os.path.join(root,fname)
                genbank_files.append(absfile)
                
            if len(genbank_files)>0:
                print "genbank_files",genbank_files
                print "bacteriocins",args.bacteriocins
                print "genes",args.genes
                print "intermediate",args.intermediate

                bacteriocin.main(genbank_files,
                                 args.bacteriocins,
                                 args.genes,
                                 outHandle,
                                 args.intermediate,
                                 args.gene_evalue,
                                 args.bac_evalue,
                                 args.num_threads,
                                 args.radius,
                                 args.verbose,
                                 args.keep_tmp)


if __name__=="__main__":
    go()
