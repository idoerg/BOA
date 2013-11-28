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

import genbank
import blast
import intergene

loc_reg = re.compile("(\d+):(\d+)\S\((\S)\)")
class GenomeHandler:
    def __init__(self,genome,intermediate,evalue,num_threads,radius,verbose,keep_tmp):
        self.pid = os.getpid() #Use current pid to name temporary files
        #self.genbank = genbank
        #self.genome_file = "%s.fna"%(os.path.splitext(genbank)[0])
        self.genome_file = genome
        # self.genomic_query = "%s/genomicQuery.%d.fa"%(intermediate,self.pid)
        #self.genome_file = "%s/genome.%d.fasta"%(intermediate,self.pid)
        # self.target_genes = "%s/target_genes.%d.fasta"%(intermediate,self.pid)
        self.evalue = evalue
        self.num_threads = num_threads
        self.radius = radius
        self.verbose = verbose
        self.keep_tmp = keep_tmp

        # SeqIO.write(input_genes,self.target_genes,"fasta")
        self.noHits = True
        self.intermediate = intermediate

    def cleanup(self):
        #os.remove(self.genomic_query)
        #os.remove(self.genome_file)
        pass

    """Gets the filtered intergenic regions"""
    def getGenomicQuery(self):
        return self.genomic_query

    def getGenomeFile(self):
        return self.genome_file

    def getAlignedGenes(self,genes,gene_evalue,num_threads,formatdb):
        # try:
        geneBlast = blast.BLAST(self.genome_file,genes,self.intermediate,gene_evalue)
        if formatdb: geneBlast.buildDatabase("nucleotide")
        geneBlast.run(blast_cmd="tblastn",mode="xml",num_threads=num_threads)
        hits = geneBlast.parseBLAST("xml")
        return hits
        # except Exception as ew:
        #     print ew
        #     return None
def main(genome_files,bacteriocins,genes,outHandle,intermediate,gene_evalue,bac_evalue,num_threads,radius,verbose,keep_tmp):
    for gbk in genome_files:
        ghr = GenomeHandler(gbk,intermediate,gene_evalue,num_threads,radius,verbose,keep_tmp)
        hits = ghr.getAlignedGenes(genes,gene_evalue,num_threads)
        outHandle.write("\n".join( map( str, hits))+"\n")

if __name__=="__main__":
    parser = argparse.ArgumentParser(description=\
        'Finds target genes in a FASTA file')
    parser.add_argument(\
        '--genome-files', type=str, nargs="+", required=False,
        help='The FASTA containg the bacterial genome')
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
        '--gene-evalue', type=float, required=False, default=0.00001,
        help='The evalue for gene hits')
    parser.add_argument(\
        '--bac-evalue', type=float, required=False, default=0.000000001,
        help='The evalue for bacteriocin hits')
    parser.add_argument(\
        '--intermediate', type=str, required=True,
        help='Directory for storing intermediate files')
    parser.add_argument(\
        '--verbose', action='store_const', const=True, default=False,
        help='Messages for debugging')
    parser.add_argument(\
        '--test', action='store_const', const=True, default=False,
        help='Run unittests')

    blast.addArgs(parser)
    args = parser.parse_args()
    outHandle = open(args.output_file,'w')

    main(args.genome_files,
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
