import Bio
from Bio import SeqIO, SeqFeature
from Bio.SeqRecord import SeqRecord

import sys
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
    '--intergenes', type=str, required=True,
    help='The path of the fasta file containing the intergenes')
parser.add_argument(\
    '--genes', type=str, required=False, nargs='+',
    help='The genes whose neighbors will be considered by BLAST')
parser.add_argument(\
    '--radius', type=str, required=False, default=10000,
    help='The search radius around every specified gene')

args = parser.parse_args()

import genbank

#May want to consider multiple genes eventually
def getIntergenes(inGeneFile,gene,radius):
    records = []
    for interGene in SeqIO.parse(inGeneFile,"fasta"):
        toks = interGene.id.split(" ")
        start,end = toks[3].split('-')
    pass

def buildDatabase(genes,radius):
    pass


if __name__=="__main__":
    geneDict = genbank.GenBank(args.genbank_path)
    sagB_gene = geneDict.findGene("sagB")
    print sagB_gene

