import Bio
from Bio import SeqIO, SeqFeature
from Bio.SeqRecord import SeqRecord

from Bio.Blast import NCBIXML
from Bio.Blast import NCBIStandalone

import sys
import os
import site
import argparse
import string
import numpy
import re
import subprocess

parser = argparse.ArgumentParser(description=\
    'Finds intergenic regions from genback file')
parser.add_argument(\
    '--genes', type=str,required=True,default="",
    help='A FASTA file containing all of the genes of interest')
parser.add_argument(\
    '--genbank-files', type=str, nargs="+", required=False,
    help='The genbank files containing annotated genes')
parser.add_argument(\
    '--genomes', type=str,nargs="+",required=False,
    help='The bacterial genomes')
parser.add_argument(\
    '--intergenes', type=str, required=True,
    help='The path of the fasta file containing the intergenes')
parser.add_argument(\
    '--bacteriocins', type=str, required=True,
    help='The bacteriocin proteins that are to be blasted')
parser.add_argument(\
    '--radius', type=int, required=False, default=10000,
    help='The search radius around every specified gene')
parser.add_argument(\
    '--test', action='store_const', const=True, default=False,
    help='Run unittests')

import genbank
import blast

blast.addArgs(parser)
args = parser.parse_args()

loc_reg = re.compile("(\d+):(\d+)")
#May want to consider multiple genes eventually
def getIntergenes(inGeneFile,gene,radius):
    records = []
    geneSt, geneEnd = loc_reg.findall(str(gene.location))[0]
    geneSt, geneEnd = int(geneSt), int(geneEnd)
    for interGene in SeqIO.parse(inGeneFile,"fasta"):
        toks = interGene.description.split(" ")
        start,end = toks[2].split('-')
        start,end = int(start),int(end)
        if (start>geneSt-radius and start<geneSt+radius and
            end>geneSt-radius and end<geneSt+radius):
            records.append(interGene)
    return records

"""
Use the Genbank files to look up genes
"""
def handleAnnotatedGenes(genbank_files,intergene_file,input_genes):
    intergenes = []
    for gbk_file in genbank_files:
        geneDict = genbank.GenBank(gbk_file)
        for gene in input_genes:
            gene_name = gene.id
            gene_record = geneDict.findGene(gene_name)
            intergenes += getIntergenes(intergene_file,gene_record,args.radius)
    return intergenes

"""
Use BLAST to look up genes in genomes

Right now, only consider the first hit
"""
def handleUnannotatedGenes():
    pass

def go():
    intergenes = []
    input_genes = [x for x in SeqIO.parse(args.genes,"fasta")]
    intergenes = handleAnnotatedGenes(args.genbank_files,args.intergenes,input_genes)

    blast_obj = blast.BLAST(args.bacteriocins,intergenes)
    blast_obj.buildDatabase()
    blast_obj.run(args.num_threads)
    hits = blast_obj.parseBLAST()
    outHandle = open(args.output_file,'w')

    outHandle.write("\n".join( map( str, hits)))
    if not args.keep_tmp: blast_obj.cleanup()

if __name__=="__main__":
    if len(args.genbank_files)==0 and len(args.genomes)==0:
        print sys.stderr,"Neither a genbank files or bacterial genomes have been specified"
        sys.close(0)
    go()
