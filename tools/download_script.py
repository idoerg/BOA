#!/usr/bin/python

# Copyright(C) 2013 David Ream
# Released under Biopython license. http://www.biopython.org/DIST/LICENSE
# Do not remove this comment

import os
import argparse
import time
import sys
from Bio import Entrez

# This exists to  make the main function easier to read. It contains code to run the argument parser, and does nothing else.
def parser_code():

    parser = argparse.ArgumentParser(description='Determine information about the genbank files under study, and report this information for use by other modules.')

    parser.add_argument("-i", "--infile", dest="infile", metavar="FILE", default='./viral_info.txt',
                help="Folder containing all genbank files for use by the program.")
                 
    parser.add_argument("-o", "--outfolder", dest="outfolder", metavar="FOLDER", default='./sequences/',
                help="Folder where results will be stored.")
                
    return parser.parse_args()
    
# A function that will check the command line input for errors. If serious errors exist, it will exit.
def check_options(parsed_args):

    if os.path.exists(parsed_args.infile):
        infile = parsed_args.infile
    else:
        print "The file %s does not exist." % parsed_args.infile
        sys.exit()
    
    # if the directory that the user specifies does not exist, then the program makes it for them. 
    if not os.path.isdir(parsed_args.outfolder):
        os.makedirs(parsed_args.outfolder)
    outfolder = parsed_args.outfolder

    return infile, outfolder

def download_genbank_files(inlist, outfolder):
    for name, accession in inlist:
        net_handle = Entrez.efetch(db="nucleotide", id=accession, rettype="gb")
        out_handle = open(outfolder + name + '.gbk', "w")
        out_handle.write(net_handle.read())
        out_handle.close()
        net_handle.close()
    
    
def main():

    # Timer, used during debug to determine the fastest implementation for a code block
    start = time.time()
    
    Entrez.email = 'jamietmorton@gmail.com'

    parsed_args = parser_code()
    
    infile, outfolder = check_options(parsed_args)
    
    in_list = []
    for item1,item2 in [(i.strip().split('\t')) for i in open(infile).readlines()]:
        in_list.append((item1, item2))
    
    download_genbank_files(in_list, outfolder)
    
    
    print time.time() - start

if __name__ == '__main__':
    main()
