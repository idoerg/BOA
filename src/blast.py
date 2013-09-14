import sys
import Bio
from Bio import SeqIO, SeqFeature
from Bio.SeqRecord import SeqRecord
import os

import site
import argparse
import string
import numpy

parser = argparse.ArgumentParser(description=\
    'Finds intergenic regions from genback file')
parser.add_argument(\
    '--genbank-path', type=str, required=True,
    help='The path of the genbank file')
parser.add_argument(\
    '--intergene-length', type=int, default=1,required=False,
    help='The path of the genbank file')
parser.add_argument(\
    '--output-dir', type=str, default='.',required=False,
    help='The output directory of the fasta files')
args = parser.parse_args()

def findGene(genbank_path,gname):
    seq_record = SeqIO.parse(open(genbank_path), "genbank").next()
    for feature in seq_record.features:
        if feature.type == 'CDS' and feature.gene==gname:
            return feature
    return None

if __name__=="__main__":
    sagB_gene = find(args.genbank_path,"sagB")

