"""
Finds bacteriocins

Todo
1) Fix the unittests
2) Deprecate the BLAST XML parser
"""

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


import genbank
import blast
import intergene


loc_reg = re.compile("(\d+):(\d+)")

class IntergeneHandler:
    def __init__(self,genbank_files,input_genes,intermediate,evalue,num_threads,radius,verbose,keep_tmp):
        self.pid = os.getpid() #Use current pid to name temporary files
        self.genbank_files = genbank_files
        self.genomic_query = "%s/genomicQuery.%d.fa"%(intermediate,self.pid)
        self.intergene_file = "%s/intergenes.%d.fasta"%(intermediate,self.pid)
        self.target_genes = "%s/target_genes.%d.fasta"%(intermediate,self.pid)
        self.evalue = evalue
        self.num_threads = num_threads
        self.radius = radius
        self.verbose = verbose
        self.keep_tmp = keep_tmp
        SeqIO.write(input_genes,self.target_genes,"fasta")

        self.intermediate = intermediate

    def cleanup(self):
        os.remove(self.genomic_query)
        os.remove(self.intergene_file)

    """Gets the filtered intergenic regions"""
    def getGenomicQuery(self):
        return self.genomic_query

    """Gets all of the intergenic regions"""
    def getIntergenicFile(self):
        return self.intergene_file

    """
    Narrows down search area for intergenes
    """
    def filterIntergenes(self,intergenic_file,geneSt,geneEnd,radius):
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
    1) BLAST sagB against all proteins
    2) Use protein ID to find the location of the protein in the genbank file
    3) Return intergenic regions

    Right it only considers the first hit
    """
    def getUnannotatedIntergenes(self,genbank_file,protein_file,evalue,num_threads,radius,verbose,keep_tmp):
        intergenes = []
        proteinDict =  genbank.GenBank(genbank_file,"protein")
        blast_obj = blast.BLAST(protein_file,self.target_genes,self.intermediate,evalue)
        blast_obj.buildDatabase("protein")
        blast_obj.run(blast_cmd="blastp",mode="xml",num_threads=num_threads)
        if verbose: print >> sys.stderr,blast_obj.formatDBCommand("protein")
        if verbose: print >> sys.stderr,blast_obj.BLASTCommand("blastp",num_threads=num_threads)
        blast_file = blast_obj.getFile()
        blast_qresults = SearchIO.read(blast_file,'blast-xml')
        organism = os.path.dirname(genbank_file)
        organism_id = os.path.splitext(os.path.basename(genbank_file))[0]

        for qresult in blast_qresults:
            query_id = qresult.id
            toks = query_id.split('|')
            protein_id = toks[3]
            protein_record = proteinDict.findProtein(protein_id)
            geneSt, geneEnd = loc_reg.findall(str(protein_record.location))[0]
            geneSt,geneEnd = int(geneSt), int(geneEnd)
            intergenes += self.filterIntergenes(self.intergene_file,geneSt,geneEnd,radius)
        if not keep_tmp:
            blast_obj.cleanup()
        return intergenes

    """
    Finds all intergenic regions
    """
    def buildIntergenicDatabase(self):
        intergenic_regions = []
        intergene_handle = open(self.genomic_query,'a')

        for gbk in self.genbank_files:
            base = os.path.splitext(os.path.basename(gbk))[0]
            path = os.path.dirname(os.path.abspath(gbk))
            protein_file = "%s/%s.faa"%(path,base) #assuming that protein file in is in the same folder
            try:
                sequences = SeqIO.parse(gbk, "genbank")
                intergene.get_interregions(gbk,self.intergene_file)
                intergenic_regions=self.getUnannotatedIntergenes(gbk,protein_file,
                                                                 self.evalue,
                                                                 self.num_threads,
                                                                 self.radius,
                                                                 self.verbose,
                                                                 self.keep_tmp)
                SeqIO.write(intergenic_regions,intergene_handle,"fasta")
            except IOError:
                print>> sys.stderr,".faa folder must be in the same folder as the genbank file"
        intergene_handle.close()

def main(genbank_files,bacteriocins,genes,outHandle,intermediate,evalue,num_threads,radius,verbose,keep_tmp):
    intergenes = []
    input_genes = [x for x in SeqIO.parse(genes,"fasta")]
    intergene_obj = IntergeneHandler(genbank_files,input_genes,
                                     intermediate,evalue,
                                     num_threads,radius,
                                     verbose,keep_tmp)
    intergene_obj.buildIntergenicDatabase()
    intergeneFile = intergene_obj.getGenomicQuery()
    blast_obj = blast.BLAST(intergeneFile,bacteriocins,intermediate,evalue)
    blast_obj.buildDatabase("nucleotide")
    blast_obj.run(blast_cmd="tblastn",mode="tab",num_threads=num_threads)
    if verbose: print >> sys.stderr,blast_obj.formatDBCommand("nucleotide")
    if verbose: print >> sys.stderr,blast_obj.BLASTCommand("tblastn",num_threads=num_threads)
    hits = blast_obj.parseBLAST("tab")
    outHandle.write("\n".join( map( str, hits)))
    if not keep_tmp:
        blast_obj.cleanup()
        intergene_obj.cleanup()

if __name__=="__main__":
    parser = argparse.ArgumentParser(description=\
        'Finds intergenic regions from genback file')
    parser.add_argument(\
        '--genbank-files', type=str, nargs="+", required=False,
        help='The genbank files containing annotated genes')
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
        '--evalue', type=float, required=False, default=0.00001,
        help='The search radius around every specified gene')
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

    main(args.genbank_files,
         args.bacteriocins,
         args.genes,
         outHandle,
         args.intermediate,
         args.evalue,
         args.num_threads,
         args.radius,
         args.verbose,
         args.keep_tmp)
