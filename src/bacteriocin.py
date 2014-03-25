
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
#import intervals
import annotations
import intergeneHandler
import pickle

from bx.intervals import *

loc_reg = re.compile("(\d+):(\d+)\S\((\S)\)")
class BacteriocinHandler:
    def __init__(self,genome,intermediate,evalue,num_threads,verbose,keep_tmp):
        self.pid = os.getpid() #Use current pid to name temporary files
        self.genbank = genbank
        self.genome_file = genome
        self.evalue = evalue
        self.num_threads = num_threads
        self.verbose = verbose
        self.keep_tmp = keep_tmp
        self.noHits = True
        self.intermediate = intermediate
 
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
    print "Building intergenic dictionary"
    print "filename",intergene_file
    print "Number of bacteriocins",len(bacteriocins)
    intergeneObj = intergeneHandler.IntergeneHandler(intergene_file)    
    intergeneObj.getIntervals()
    intergeneDict = dict()
    for bact in bacteriocins:
        gene = bact.sbjct_id,bact.sbjct_start,bact.sbjct_end,bact.strand
        overlaps = intergeneObj.overlapIntergene(gene)
        intergeneDict[(bact.sbjct_id,bact.sbjct_start,bact.sbjct_end,bact.strand)] = overlaps
    return intergeneDict

"""
Filters out all query sequences that don't fall in one of the intervals
Note: returns a list of all of the overlapping intervals
intervalDict: A dictionary of intervals indexed by (organismID,strand)
queries: A list of query objects with
        (start,end,orgid,strand,query_obj)
"""
def spatialFilter(queries,intervalDict,radius):
    filtered = []
    geneNeighborhoods = []
    for query in queries:
        header, query_obj = query
        start,end,orgid,strand = header
        if (orgid,strand) not in intervalDict: continue
        nearTargets = intervalDict[(orgid,strand)].find( start,end )
        if len(nearTargets)>0:
            for geneobj in nearTargets:
                filtered.append(query_obj)
                geneNeighborhoods.append(geneobj)
    return filtered,geneNeighborhoods
               

"""
Filter out all annotations that are within the radius of a bacteriocin
"""
annot_reg = re.compile("([A-Z]+_[0-9]+).")
def filterAnnotations(annots,bacteriocins,radius):
    intervalDict = defaultdict( IntervalTree )
    for bact in bacteriocins:
        start,end,orgid,strand = bact.sbjct_start,bact.sbjct_end,bact.sbjct_id,bact.strand
        orgid = annot_reg.findall(orgid)[0]
        stBound,endBound = start-radius,end+radius
        intervalDict[(orgid,strand)].add( stBound,endBound,bact )                                                 
    headers = []
    for a in annots: 
        orgid,start,end,strand = a[0],a[1],a[2],a[3];
        headers.append( (orgid,start,end,strand) )
    queries = zip(headers,annots) 
    filtered,annotNeighborhoods = spatialFilter(queries,intervalDict,radius)    
    return filtered,annotNeighborhoods

"""
Filters out bacteriocins not contained in a gene neighborhood
bacteriocins: list of bacteriocins blast hits
genes: list of target genes blast hits 
radius: radius around target gene
filtered: list of target genes blast hits that make it
          through the filtering criteria imposed by the radius
geneNeighborhoods: a list of intervals with a radius surrounding the
          target genes
header = (orgid,start,end,strand)
"""
def filterBacteriocins(bacteriocins,genes,radius):
    if radius<0: radius = 50000000 #largest bacterial genome
    intervalDict = defaultdict( IntervalTree )
    for gene in genes:    
        start,end,refid,orgid,strand = gene.sbjct_start,gene.sbjct_end,gene.query_id,gene.sbjct_id,gene.strand
        stBound,endBound = start-radius,end+radius
        intervalDict[(orgid,strand)].add( stBound,endBound,gene )
    headers = []
    for b in bacteriocins: 
        headers.append( (b.sbjct_start,b.sbjct_end,b.sbjct_id,b.strand) )
    queries = zip(headers,bacteriocins) 
    filtered,geneNeighborhoods = spatialFilter(queries,intervalDict,radius)
    return filtered,geneNeighborhoods

def writeBacteriocins(bacteriocins,intergeneDict,outHandle,genes=False):
    for bacteriocin in bacteriocins:
        if genes:  bacteriocin,gene = bacteriocin
        inIntergene = intergeneDict[(bacteriocin.sbjct_id,
                                     bacteriocin.sbjct_start,
                                     bacteriocin.sbjct_end,
                                     bacteriocin.strand)]
        regionType = "intergene" if inIntergene else "gene"        
        bacID    = bacteriocin.query_id
        organism = bacteriocin.sbjct_id
        bacID    = bacteriocin.query_id.split(' ')[0]
        organism = bacteriocin.sbjct_id.split(' ')[0]
        if genes:
            geneStart,geneEnd,geneName = gene.sbjct_start,gene.sbjct_end,gene.query_id
            argstr   = "%s\t %s\t %s\t %s\t %s\t %s\t %d\t %d\t %s\t %s\t %s \n"
            result_str = argstr%(bacID,
                                 organism,
                                 bacteriocin.sbjct_start,
                                 bacteriocin.sbjct_end,
                                 bacteriocin.strand,
                                 geneName,
                                 geneStart,
                                 geneEnd,
                                 gene.strand,
                                 regionType,
                                 bacteriocin.sbjct)
        else:
            argstr   = "%s\t %s\t %d\t %d\t %s\t %s\t %s \n"
            result_str = argstr%(bacID,
                                 organism,
                                 bacteriocin.sbjct_start,
                                 bacteriocin.sbjct_end,
                                 bacteriocin.strand,
                                 regionType,
                                 bacteriocin.sbjct)
        outHandle.write(result_str)


def writeAnnotations(annot_bact_pairs,outHandle):
    for annot_bact in annot_bact_pairs:
        annot,bacteriocin = annot_bact
        annot_st,annot_end,annot_org,annot_strand,annot_locus,annot_seq = annot
        bacID    = bacteriocin.query_id
        organism = bacteriocin.sbjct_id
        bacID    = bacteriocin.query_id.split(' ')[0]
        organism = bacteriocin.sbjct_id.split(' ')[0]
        organism+= annot_locus
        argstr   = "%s\t %s\t %s\t %s\t %s\t %s\t %d\t %d\t %s\t %s \n"
        result_str = argstr%(bacID,
                             organism,
                             bacteriocin.sbjct_start,
                             bacteriocin.sbjct_end,
                             bacteriocin.strand,
                             annot_org,
                             annot_st,
                             annot_end,
                             annot_strand,
                             annot_seq)        
        outHandle.write(result_str)

def main(genome_files,bacteriocins,
         genes,intergene_file,
         annotations_file,
         bacteriocinsOut,
         filteredOut,
         annotationsOut,
         intermediate,
         gene_evalue,bac_evalue,
         num_threads,formatdb,
         gene_radius,
         bacteriocin_radius,
         verbose,keep_tmp):
    for gnome in genome_files:
        gnomehr = genome.GenomeHandler(gnome,
                                       intermediate,
                                       gene_evalue,
                                       num_threads,
                                       verbose,
                                       keep_tmp)
        genes = gnomehr.getAlignedGenes(genes,gene_evalue,num_threads,formatdb)
        bacthr = BacteriocinHandler(gnome,
                                    intermediate,
                                    bac_evalue,
                                    num_threads,
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
        if verbose: print "Number of original bacteriocins",len(bacteriocins)
        intergeneDict = identifyIntergenic(bacteriocins,intergene_file)
        writeBacteriocins(bacteriocins,intergeneDict,bacteriocinsOut)
        annots = [annot for annot in annotations.Annotations(annotations_file)]
        annots,bacteriocinNeighborhoods = filterAnnotations(annots,bacteriocins,bacteriocin_radius)
        annot_bact_pairs = zip(annots,bacteriocinNeighborhoods)
        writeAnnotations(annot_bact_pairs, annotationsOut)
        
        bacteriocins,geneNeighborhoods = filterBacteriocins(bacteriocins,genes,gene_radius)
        if verbose: print "Number of filtered bacteriocins",len(bacteriocins)
        bact_gene_pairs = zip(bacteriocins,geneNeighborhoods)
        intergeneDict = identifyIntergenic(bacteriocins,intergene_file)
        writeBacteriocins(bact_gene_pairs,intergeneDict, filteredOut,genes="True")

        #pickle.dump(intergeneDict,open("intergene.dict",'w'))

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
        '--annotations', type=str, required=False,
        help='FASTA files containing annotated regions')
    parser.add_argument(\
        '--bacteriocins', type=str, required=False,
        help='The bacteriocin proteins that are to be blasted')
    parser.add_argument(\
        '--gene_radius', type=int, required=False, default=-1,
        help='The search radius around every specified gene. By default, this filter option is off')
    parser.add_argument(\
        '--bacteriocin-radius', type=int, required=False, default=5000,
        help='The search radius around every specified bacteriocin')
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
        '--output', type=str, required=False,
        help='The output file basename for filtered annotations and bacteriocins')
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
        bacteriocinsOut = open("%s_bacteriocins.txt"%(args.output),'w')
        filteredOut = open("%s_filtered.txt"%(args.output),'w')
        annotationsOut  = open("%s_annotations.txt"%(args.output),'w')
        main(args.genome_files,
             args.bacteriocins,
             args.genes,
             args.intergenes,
             args.annotations,
             bacteriocinsOut,
             filteredOut,
             annotationsOut,
             args.intermediate,
             args.gene_evalue,
             args.bac_evalue,
             args.num_threads,
             args.formatdb,
             args.gene_radius,
             args.bacteriocin_radius,
             args.verbose,
             args.keep_tmp)
    else:
        del sys.argv[1:]
        import unittest
        import test_genbank
        class TestFilters(unittest.TestCase):
            def setUp(self):
                test_input = test_genbank.yeast
                self.test_file = "test.gbk"
                self.out_file = "out.fa"
                handle = open(self.test_file,'w')
                handle.write(test_input)
                handle.close()
                annotations.parseAnnotations("NC_12345",self.test_file,open(self.out_file,'w'))
            def test_filter_bacteriocins_1(self):
                bacteriocins = [blast.XMLRecord(description="",
                                                expected_value=0.00001,
                                                score = 0,
                                                query_id = "bacteriocin1",
                                                query="ACGTACGTTT",
                                                query_start = 1,
                                                query_end   = 10,
                                                sbjct_id = "NC_12345",
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
                                                sbjct_id = "NC_12345",
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
                                            sbjct_id = "NC_12346",
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
                                         sbjct_id = "NC_12345",
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
                self.assertEquals(1,len(filtered))
                self.assertEquals(1,len(hoods))
                self.assertTrue(bacteriocins[0] in filtered)
                record = hoods[0]
                start,end,refid,gene = record.sbjct_start,record.sbjct_end,record.query_id,record.strand
                
                self.assertEquals(start,150)
                self.assertEquals(end,160)
                self.assertEquals(refid,"gene1")
            def test_filter_bacteriocins_2(self):
                bacteriocins = [blast.XMLRecord(description="",
                                                expected_value=0.00001,
                                                score = 0,
                                                query_id = "bacteriocin1",
                                                query="ACGTACGTTT",
                                                query_start = 1,
                                                query_end   = 10,
                                                sbjct_id = "NC_12345",
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
                                                sbjct_id = "NC_12345",
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
                                            sbjct_id = "NC_12346",
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
                                         sbjct_id = "NC_12345",
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
                                         sbjct_id = "NC_12346",
                                         sbjct="AAAAAAAAAA",
                                         sbjct_start = 1050,
                                         sbjct_end   = 1060,
                                         strand = "-")]
                radius = 10000
                filtered,hoods = filterBacteriocins(bacteriocins,genes,radius)
                print '\n'.join(map(str,filtered))
                print '\n'.join(map(str,hoods))
                self.assertEquals(3,len(filtered))
                self.assertEquals(3,len(hoods))
                self.assertTrue(bacteriocins[0] in filtered)

            def test_filter_annotations_1(self):
                bacteriocins = [blast.XMLRecord(description="",
                                                expected_value=0.00001,
                                                score = 0,
                                                query_id = "bacteriocin1",
                                                query="ACGTACGTTT",
                                                query_start = 250,
                                                query_end   = 260,
                                                sbjct_id = "NC_12345",
                                                sbjct="ACGTACGTTT",
                                                sbjct_start = 250,
                                                sbjct_end   = 260,
                                                strand = "-"),
                                blast.XMLRecord(description="",
                                                expected_value=0.00001,
                                                score = 0,
                                                query_id = "bacteriocin2",
                                                query="ACGTACGTTT",
                                                query_start = 450,
                                                query_end   = 460,
                                                sbjct_id = "NC_12345",
                                                sbjct="ACGTACGTTT",
                                                sbjct_start = 450,
                                                sbjct_end   = 460,
                                                strand = "-")]
                radius = 100
                annots = [A for A in annotations.Annotations(self.out_file)]
                filtered,hoods = filterAnnotations(annots,bacteriocins,radius)
                # self.assertEquals(1,len(filtered))
                # self.assertEquals(1,len(hoods))
                # self.assertTrue(bacteriocins[0] in filtered)
                # start,end,refid,gene = hoods[0]
                # self.assertEquals(start,150)
                # self.assertEquals(end,160)
                # self.assertEquals(refid,"gene1")
                
                
        unittest.main()
