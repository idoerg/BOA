import Bio
from Bio import SeqIO, SeqFeature
from Bio.SeqRecord import SeqRecord

from Bio.Blast import NCBIXML
from Bio.Blast import NCBIStandalone

from collections import defaultdict

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
import intergeneHandler

loc_reg = re.compile("(\d+):(\d+)\S\((\S)\)")
class BacteriocinHandler:
    def __init__(self,genome,intermediate,evalue,num_threads,radius,verbose,keep_tmp):
        self.pid = os.getpid() #Use current pid to name temporary files
        self.genbank = genbank
        self.genome_file = genome
        self.evalue = evalue
        self.num_threads = num_threads
        self.radius = radius
        self.verbose = verbose
        self.keep_tmp = keep_tmp
        self.noHits = True
        self.intermediate = intermediate

    def cleanup(self):
        os.remove(self.genomic_query)
        #os.remove(self.genome_file)

    """Gets the filtered intergenic regions"""
    def getGenomicQuery(self):
        return self.genomic_query

    def getGenomeFile(self):
        return self.genome_file

    """Get all of the bacteriocins in all of the bacterial genomes"""
    def getAlignedBacteriocins(self,bacteriocins,bac_evalue,num_threads,formatdb):
        bacBlast = blast.BLAST(self.genome_file,bacteriocins,self.intermediate,bac_evalue)
        if formatdb:
            bacBlast.buildDatabase("nucleotide")
            print bacBlast.formatDBCommand()
        bacBlast.run(blast_cmd="tblastn",mode="xml",num_threads=num_threads)
        hits = bacBlast.parseBLAST("xml")
        return hits

"""
Given a set of genomes and bacteriocins, determine which bacteriocins are in intergenic regions
"""
def identifyIntergenic(bacteriocins,intergene_file):
    print "filename",intergene_file
    intergeneObj = intergeneHandler.IntergeneHandler(intergene_file)
    intergeneObj.getIntervals()
    intergeneDict = dict()
    for bact in bacteriocins:
        gene = bact.sbjct_id,bact.sbjct_start,bact.sbjct_end,bact.strand
        intergeneDict[(bact.sbjct_id,bact.sbjct_start,bact.sbjct_end,bact.strand)] = intergeneObj.overlapIntergene(gene)
    return intergeneDict
        
"""
Filters out bacteriocins not contained in a gene neighborhood
bacteriocins: list of bacteriocins blast hits
genes: list of target genes blast hits 
radius: radius around target gene
filtered: list of target genes blast hits that make it
          through the filtering criteria imposed by the radius
geneNeighborhoods: a list of intervals with a radius surrounding the
          target genes
"""
def filterBacteriocins(bacteriocins,genes,radius):
    #ints = intervals.Intervals()
    intervalDict = defaultdict( intervals.Intervals )
    for gene in genes:
        start,end,refid,orgid = gene.sbjct_start,gene.sbjct_end,gene.query_id,gene.sbjct_id
        intervalDict[orgid].append( (start-radius,end+radius,refid,gene) )
    filtered = []
    geneNeighborhoods = []
    for bact in bacteriocins:
        start,end,orgid = bact.sbjct_start,bact.sbjct_end,bact.sbjct_id
        if orgid not in intervalDict:
            continue
        nearestGene = intervalDict[orgid].search( (start,end) )
        if nearestGene!=None:
            filtered.append( bact )
            nearestGene = intervals.reformat(nearestGene,radius)
            gstart,gend,refid,gene = nearestGene
            print "Radius: ",radius
            print "Bacteriocin: (%d,%d,%s)"%(start,end,orgid)
            print "Gene: (%d,%d,%s)"%(gstart,gend,gene.sbjct_id)
            geneNeighborhoods.append(nearestGene)
    return filtered,geneNeighborhoods


def main(genome_files,bacteriocins,
         genes,intergene_file,
         outHandle,intermediate,
         gene_evalue,bac_evalue,
         num_threads,formatdb,
         radius,verbose,keep_tmp):
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
        bact_gene_pairs = zip(bacteriocins,geneNeighborhoods)
        intergeneDict = identifyIntergenic(bacteriocins,intergene_file)        
        for bact_gene in bact_gene_pairs:
            bacteriocin,gene = bact_gene
            inIntergene = intergeneDict[(bacteriocin.sbjct_id,
                                         bacteriocin.sbjct_start,
                                         bacteriocin.sbjct_end,
                                         bacteriocin.strand)]
            regionType = "intergene" if inIntergene else "gene"
            bac_loc  = "%s-%s"%(bacteriocin.sbjct_start,
                                bacteriocin.sbjct_end)
            geneStart,geneEnd,geneName,geneRecord = gene[0],gene[1],gene[2],gene[3]
            gene_loc = "%s-%s"%(geneStart,geneEnd) 
            bacID    = bacteriocin.query_id
            organism = bacteriocin.sbjct_id
            bacID    = bacteriocin.query_id.split(' ')[0]
            organism = bacteriocin.sbjct_id.split(' ')[0]
            argstr   = "%s\t %s\t %s\t %s\t %s\t %s\t %d\t %d\t %s\t %s\t %"
            if verbose:  print argstr%(bacID,
                                       organism,
                                       bacteriocin.sbjct_start,
                                       bacteriocin.sbjct_end,
                                       bacteriocin.strand,
                                       geneName,
                                       geneStart,
                                       geneEnd,
                                       geneRecord.strand,
                                       regionType,
                                       bacteriocin.sbjct)
            outHandle.write(argstr%(bacID,
                                    organism,
                                    bacteriocin.sbjct_start,
                                    bacteriocin.sbjct_end,
                                    bacteriocin.strand,
                                    geneName,
                                    geneStart,
                                    geneEnd,
                                    geneRecord.strand,
                                    regionType,
                                    bacteriocin.sbjct))

if __name__=="__main__":
    parser = argparse.ArgumentParser(description=\
        'Finds bacteriocins surrounding target genes')
    parser.add_argument(\
        '--genome-files', type=str, nargs="+", required=False,
        help='FASTA files containing bacterial genomes')
    parser.add_argument(\
        '--genes', type=str,required=False,default="",
        help='A FASTA file containing all of the target genes of interest')
    parser.add_argument(\
        '--intergenes', type=str, required=False,
        help='FASTA files containing intergenic regions')
    parser.add_argument(\
        '--bacteriocins', type=str, required=False,
        help='The bacteriocin proteins that are to be blasted')
    parser.add_argument(\
        '--radius', type=int, required=False, default=10000,
        help='The search radius around every specified gene')
    parser.add_argument(\
        '--gene-evalue', type=float, required=False, default=0.00001,
        help='The evalue for gene hits')
    parser.add_argument(\
        '--bac-evalue', type=float, required=False, default=0.00001,
        help='The evalue for bacteriocin hits')
    parser.add_argument(\
        '--intermediate', type=str, required=False,
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
    #outHandle.write("bacteriocin\tbacteriocin_location\torganism\tgene\tbacterciocin_sequence\n")
    if not args.test:
        outHandle = open(args.output_file,'w')
        main(args.genome_files,
             args.bacteriocins,
             args.genes,
             args.intergenes,
             outHandle,
             args.intermediate,
             args.gene_evalue,
             args.bac_evalue,
             args.num_threads,
             args.formatdb,
             args.radius,
             args.verbose,
             args.keep_tmp)
    else:
        del sys.argv[1:]
        import unittest        
        class TestFilters(unittest.TestCase):            
            def test_filter_bacteriocins_1(self):
                bacteriocins = [blast.XMLRecord(description="",
                                                expected_value=0.00001,
                                                score = 0,
                                                query_id = "bacteriocin1",
                                                query="ACGTACGTTT",
                                                query_start = 1,
                                                query_end   = 10,
                                                sbjct_id = "org1",
                                                sbjct="ACGTACGTTT",
                                                sbjct_start = 100,
                                                sbjct_end   = 110,
                                                strand = "-"),
                                blast.XMLRecord(description="",
                                                expected_value=0.00001,
                                                score = 0,
                                                query_id = "bacteriocin2",
                                                query="ACGTACGTTT",
                                                query_start = 1,
                                                query_end   = 10,
                                                sbjct_id = "org1",
                                                sbjct="ACGTACGTTT",
                                                sbjct_start = 1000,
                                                sbjct_end   = 1010,
                                                strand = "-"),
                                blast.XMLRecord(description="",
                                                expected_value=0.00001,
                                                score = 0,
                                                query_id = "bacteriocin3",
                                                query="ACGTACGTTT",
                                                query_start = 1,
                                                query_end   = 10,
                                                sbjct_id = "org2",
                                                sbjct="ACGTACGTTT",
                                                sbjct_start = 100,
                                                sbjct_end   = 110,
                                                strand = "-")]
                genes = [blast.XMLRecord(description="",
                                         expected_value=0.00001,
                                         score = 0,
                                         query_id = "gene1",
                                         query="AAAAAAAAAA",
                                         query_start = 1,
                                         query_end   = 10,
                                         sbjct_id = "org1",
                                         sbjct="AAAAAAAAAA",
                                         sbjct_start = 150,
                                         sbjct_end   = 160,
                                         strand = "-"),
                         blast.XMLRecord(description="",
                                         expected_value=0.00001,
                                         score = 0,
                                         query_id = "gene2",
                                         query="AAAAAAAAAA",
                                         query_start = 1,
                                         query_end   = 10,
                                         sbjct_id = "org2",
                                         sbjct="AAAAAAAAAA",
                                         sbjct_start = 1050,
                                         sbjct_end   = 1060,
                                         strand = "-")]
                radius = 100
                filtered,hoods = filterBacteriocins(bacteriocins,genes,radius)
                print map( str, filtered)
                self.assertEquals(1,len(filtered))
                self.assertEquals(1,len(hoods))
                self.assertTrue(bacteriocins[0] in filtered)
                start,end,refid,gene = hoods[0]
                self.assertEquals(start,150)
                self.assertEquals(end,160)
                self.assertEquals(refid,"gene1")
        unittest.main()
