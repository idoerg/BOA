import Bio
from Bio import SeqIO, SeqFeature
from Bio.SeqRecord import SeqRecord
import sys
import os
import site
import argparse
import string
import numpy
import re
import subprocess

accession_reg = re.compile("(NC_\d\d\d\d\d\d.1)")

def go(bacteriocin_file,RNA_file):

    record_dict = {}
    for record in SeqIO.parse( open(bacteriocin_file), "fasta") :
        descriptions = accession_reg.findall(record.description)
        if len(descriptions)>0:
            description = descriptions[0]
            record_dict[description] = record
    for record in SeqIO.parse( open(RNA_file), "fasta") :
        descriptions = accession_reg.findall(record.description)
        if len(descriptions)>0:
            if descriptions[0] in record_dict:
                print record.format("fasta")

if __name__=="__main__":
    parser = argparse.ArgumentParser(description=\
        'Filters out sequences without bacteriocins')
    parser.add_argument(\
        '--bacteriocin-file', type=str, required=False,
        help='FASTA files containing bacteriocin genes')
    parser.add_argument(\
        '--bacteria-16SRNA', type=str, required=False,
        help='FASTA files containing 16SRNA genes')
    args = parser.parse_args()
    go(args.bacteriocin_file,args.bacteria_16SRNA)
    
