import Bio
from Bio import SeqIO, SeqFeature
from Bio import SearchIO
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
    '--bacteriocins', type=str, required=True,
    help='The bacteriocin proteins that are to be blasted')
parser.add_argument(\
    '--radius', type=int, required=False, default=10000,
    help='The search radius around every specified gene')
parser.add_argument(\
    '--intermediate', type=str, required=True,
    help='Directory for storing intermediate files')
parser.add_argument(\
    '--test', action='store_const', const=True, default=False,
    help='Run unittests')

import genbank
import blast
import intergene

blast.addArgs(parser)
args = parser.parse_args()

loc_reg = re.compile("(\d+):(\d+)")

"""
Narrows down search area for intergenes
"""
def filterIntergenes(intergenic_file,geneSt,geneEnd,radius):
    records = []
    for interGene in SeqIO.parse(intergenic_file,"fasta"):
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
def handleAnnotatedIntergenes(genbank_file,intergene_file,input_genes):
    intergenes = []
    geneDict = genbank.GenBank(genbank_file)
    for gene in input_genes:
        gene_name = gene.id
        gene_record = geneDict.findGene(gene_name)
        geneSt, geneEnd = loc_reg.findall(str(gene_record.location))[0]
        geneSt, geneEnd = int(geneSt), int(geneEnd)
        intergenes += filterIntergenes(intergene_file,geneSt,geneEnd,args.radius)
    return intergenes

"""
1) BLAST sagB against all proteins
2) Use protein ID to find the location of the protein in the genbank file
3) Return intergenic regions

Right it only considers the first hit
"""

def handleUnannotatedIntergenes(genbank_file,protein_file,intergene_file,input_genes):
    proteinDict =  genbank.GenBank(genbank_file,"protein")
    blast_obj = blast.BLAST(protein_file,input_genes)
    blast_obj.buildDatabase("protein")
    blast_obj.run(blast_cmd="blastp",num_threads=args.num_threads)
    blast_file = blast_obj.getFile()
    blast_qresults = SearchIO.read(blast_file,'blast-xml')

    for qresult in blast_qresults:
        query_id = qresult.id
        toks = query_id.split('|')
        protein_id = toks[3]
        protein_record = proteinDict.findProtein(protein_id)
        geneSt, geneEnd = loc_reg.findall(str(protein_record.location))[0]
        geneSt,geneEnd = int(geneSt), int(geneEnd)
        return filterIntergenes(intergene_file,geneSt,geneEnd,args.radius)


"""
Finds all intergenic regions
"""
def getAllIntergenicRegions(genbank_files,input_genes,intermediate):
    intergenic_regions = []
    for gbk in genbank_files:
        base = os.path.splitext(os.path.basename(gbk))[0]
        path = os.path.dirname(os.path.abspath(gbk))
        intergene_file = "%s/%s.fasta"%(intermediate,base)
        protein_file = "%s/%s.faa"%(path,base)
        sequences = SeqIO.parse(gbk, "genbank")

        intergene.get_interregions(gbk,intergene_file)
        intergenic_regions+=handleUnannotatedIntergenes(gbk,protein_file,intergene_file,input_genes)

        if not args.keep_tmp:
            os.remove(intergene_file)

    return intergenic_regions

def go():
    intergenes = []
    input_genes = [x for x in SeqIO.parse(args.genes,"fasta")]
    intergenes = getAllIntergenicRegions(args.genbank_files,input_genes,args.intermediate)

    blast_obj = blast.BLAST(args.bacteriocins,intergenes)
    blast_obj.buildDatabase("protein")
    blast_obj.run(blast_cmd="blastx",num_threads=args.num_threads)
    hits = blast_obj.parseBLAST()
    outHandle = open(args.output_file,'w')

    outHandle.write("\n".join( map( str, hits)))
    if not args.keep_tmp: blast_obj.cleanup()

if __name__=="__main__":
    go()
