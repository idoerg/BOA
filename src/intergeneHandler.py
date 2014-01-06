"""
Finds bacteriocins

Todo
1) Fix the unittests
2) Deprecate the BLAST XML parser
"""

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
        self.noHits = True
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
    def filterIntergenes(self,intergenic_file,geneSt,geneEnd,geneStrand,radius):
        records = []
        for interGene in SeqIO.parse(intergenic_file,"fasta"):
            toks = interGene.description.split(" ")
            start,end = toks[2].split('-')
            start,end = int(start),int(end)
            if (start>geneSt-radius and start<geneSt+radius and
                end>geneSt-radius and end<geneSt+radius):
                interGene.id+="|locus:%d-%d%s"%(geneSt,geneEnd,geneStrand)
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
        blast_obj.buildDatabase(base="protein")
        blast_obj.run(blast_cmd="blastp",mode="xml",num_threads=num_threads)
        if verbose: print >> sys.stderr,blast_obj.formatDBCommand()
        if verbose: print >> sys.stderr,blast_obj.BLASTCommand()
        blast_file = blast_obj.getFile()
        hits = blast_obj.parseBLAST("xml")
        if len(hits)>0: self.noHits=False
        for qresult in hits:
            description = qresult.description
            query_id = description.split(' ')[1]
            toks = query_id.split('|')
            protein_id = toks[3]
            protein_record = proteinDict.findProtein(protein_id)
            geneSt, geneEnd, geneStrand = loc_reg.findall(str(protein_record.location))[0]
            geneSt,geneEnd = int(geneSt), int(geneEnd)
            intergenes += self.filterIntergenes(self.intergene_file,geneSt,geneEnd,geneStrand,radius)
        if not keep_tmp:
            blast_obj.cleanup()
        return intergenes

    """
    Finds all intergenic regions
    """
    def buildIntergenicDatabase(self):
        intergenic_regions = []
        intergene_handle = open(self.genomic_query,'w')
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
                if len(intergenic_regions)==0:
                    if self.verbose: print >>sys.stderr,"No target genes detected"
                else:
                    SeqIO.write(intergenic_regions,intergene_handle,"fasta")
            except IOError:
                print>> sys.stderr,".faa folder must be in the same folder as the genbank file"
        intergene_handle.close()


def main(genbank_files,bacteriocins,genes,outHandle,intermediate,gene_evalue,bac_evalue,num_threads,radius,verbose,keep_tmp):
    try:
        intergenes = []
        input_genes = [x for x in SeqIO.parse(genes,"fasta")]
        intergene_obj = IntergeneHandler(genbank_files,input_genes,
                                         intermediate,gene_evalue,
                                         num_threads,radius,
                                         verbose,keep_tmp)
        intergene_obj.buildIntergenicDatabase()
        intergeneFile = intergene_obj.getGenomicQuery()
        if intergene_obj.noHits==False:
            blast_obj = blast.BLAST(intergeneFile,bacteriocins,intermediate,bac_evalue)
            blast_obj.buildDatabase("nucleotide")
            blast_obj.run(blast_cmd="tblastn",mode="xml",num_threads=num_threads)
            if verbose: print >> sys.stderr,blast_obj.formatDBCommand()
            if verbose: print >> sys.stderr,blast_obj.BLASTCommand()
            hits = blast_obj.parseBLAST("coord")
            outHandle.write("organism\tlocus\t\thit\t\tbacteriocin\n")
            outHandle.write("\n".join( map( str, hits))+"\n")
            if not keep_tmp: blast_obj.cleanup()
        if not keep_tmp: intergene_obj.cleanup()
    except Exception as e:
        print "Error",e

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

    main(args.genbank_files,
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