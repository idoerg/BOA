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
import genome
import intervals

loc_reg = re.compile("(\d+):(\d+)\S\((\S)\)")
class BacteriocinHandler:
    def __init__(self,genome,intermediate,evalue,num_threads,radius,verbose,keep_tmp):
        self.pid = os.getpid() #Use current pid to name temporary files
        self.genbank = genbank
        #self.genome_file = "%s.fna"%(os.path.splitext(genbank)[0])
        self.genome_file = genome
        #self.genome_file = "%s/genome.%d.fasta"%(intermediate,self.pid)
        self.evalue = evalue
        self.num_threads = num_threads
        self.radius = radius
        self.verbose = verbose
        self.keep_tmp = keep_tmp
        self.noHits = True
        self.intermediate = intermediate

    def cleanup(self):
        os.remove(self.genomic_query)
        os.remove(self.genome_file)

    """Gets the filtered intergenic regions"""
    def getGenomicQuery(self):
        return self.genomic_query

    def getGenomeFile(self):
        return self.genome_file

    def getAlignedBacteriocins(self,bacteriocins,bac_evalue,num_threads,formatdb):
        try:
            bacBlast = blast.BLAST(self.genome_file,bacteriocins,self.intermediate,bac_evalue)
            if formatdb: bacBlast.buildDatabase("nucleotide")
            bacBlast.run(blast_cmd="tblastn",mode="xml",num_threads=num_threads)
            hits = bacBlast.parseBLAST("xml")
            return hits
        except Exception as ew:
            print ew
            return None

#Filters out bacteriocins not contained in a gene neighborhood
def filterBacteriocins(bacteriocins,genes,radius):
    ints = intervals.Intervals()
    for gene in genes:
        start,end,refid = gene.sbjct_start,gene.sbjct_end,gene.query_id
        ints.append((start-radius,end+radius,refid,gene))
    filtered = []
    geneNeighborhoods = []
    for bact in bacteriocins:
        start,end = bact.sbjct_start,bact.sbjct_end
        nearestGene = ints.search( (start,end) )
        if nearestGene!=None:
            filtered.append( bact )
            geneNeighborhoods.append(nearestGene)
    return filtered,geneNeighborhoods


def main(genome_files,bacteriocins,genes,outHandle,intermediate,gene_evalue,bac_evalue,num_threads,formatdb,radius,verbose,keep_tmp):
    for gnome in genome_files:
        gnomehr = genome.GenomeHandler(gnome,
                                       intermediate,
                                       gene_evalue,
                                       num_threads,
                                       radius,
                                       verbose,
                                       keep_tmp)
        genes = gnomehr.getAlignedGenes(genes,gene_evalue,num_threads,formatdb)
        bacthr = BacteriocinHandler(gnome,
                                    intermediate,
                                    bac_evalue,
                                    num_threads,
                                    radius,
                                    verbose,
                                    keep_tmp)
        bacteriocins = bacthr.getAlignedBacteriocins(bacteriocins,
                                                     bac_evalue,
                                                     num_threads,
                                                     formatdb)
        if verbose and genes == None: print "No genes found"
        if verbose and bacteriocins == None: print "No bacteriocins found"
        if genes == None or bacteriocins == None: continue
        if verbose: print "Genes found\n","\n".join(map(str,genes))
        if verbose: print "Bacteriocins found\n","\n".join(map(str,bacteriocins))
        bacteriocins,geneNeighborhoods = filterBacteriocins(bacteriocins,
                                                            genes,
                                                            radius)
        records = zip(bacteriocins,geneNeighborhoods)
        for record in records:
            bacteriocin,gene = record
            bac_loc = "%s-%s"%(bacteriocin.sbjct_start,
                               bacteriocin.sbjct_end)
            geneStart,geneEnd,geneName,geneRecord = gene[0],gene[1],gene[2],gene[3]
            gene_loc = "%s-%s"%(geneStart,geneEnd)
            bacID = bacteriocin.query_id.split(' ')[0]
            organism = bacteriocin.sbjct_id.split(' ')[0]
            mid = (geneStart+geneEnd)/2
            if verbose:  print "%s\t %s\t %s\t %s\t %s\t %s\t %s\t %s\t %s"%(bacID,
                                                                             organism,
                                                                             bacteriocin.sbjct_start,
                                                                             bacteriocin.sbjct_end,
                                                                             bacteriocin.strand,
                                                                             geneName,
                                                                             mid,
                                                                             geneRecord.strand,
                                                                             bacteriocin.sbjct)
            outHandle.write("%s\t %s\t %s\t %s\t %s\t %s\t %s\t %s\t %s\n"%(bacID,
                                                                            organism,
                                                                            bacteriocin.sbjct_start,
                                                                            bacteriocin.sbjct_end,
                                                                            bacteriocin.strand,
                                                                            geneName,
                                                                            mid,
                                                                            geneRecord.strand,
                                                                            bacteriocin.sbjct))

if __name__=="__main__":
    parser = argparse.ArgumentParser(description=\
        'Finds bacteriocins surrounding target genes')
    parser.add_argument(\
        '--genome-files', type=str, nargs="+", required=False,
        help='FASTA files containing bacterial genomes')
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
        '--bac-evalue', type=float, required=False, default=0.000001,
        help='The evalue for bacteriocin hits')
    parser.add_argument(\
        '--intermediate', type=str, required=True,
        help='Directory for storing intermediate files')
    parser.add_argument(\
        '--verbose', action='store_const', const=True, default=False,
        help='Messages for debugging')
    parser.add_argument(\
        '--formatdb', action='store_const', const=True, default=False,
        help='Indicates if formatdb should be run')
    parser.add_argument(\
        '--test', action='store_const', const=True, default=False,
        help='Run unittests')

    blast.addArgs(parser)
    args = parser.parse_args()
    outHandle = open(args.output_file,'w')
    outHandle.write("bacteriocin\tbacteriocin_location\torganism\tgene\tbacterciocin_sequence\n")

    main(args.genome_files,
         args.bacteriocins,
         args.genes,
         outHandle,
         args.intermediate,
         args.gene_evalue,
         args.bac_evalue,
         args.num_threads,
         args.formatdb,
         args.radius,
         args.verbose,
         args.keep_tmp)
