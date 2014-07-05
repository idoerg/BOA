""" 
OUTPUT:
Bacteriocins
1. bacterion ID
2. organism
3. bacteriocin start
4. bacteriocin end
5. bacteriocin strand
6. regionType


Annotated Genes
1.  bacteriocin ID 
2.  organism 
3.  bacteriocin start 
4.  bacteriocin end 
5.  bacteriocin strand  
6.  annotated gene organism 
7.  annotated gene locus 
8.  annotated gene protein id 
9.  annotated gene start 
10. annotated gene end 
11. annotated gene strand 
12. annotated gene sequence

TODO: Need to reorganize all regex commansd
"""
import Bio
from Bio import SeqIO, SeqFeature
from Bio.SeqRecord import SeqRecord

from Bio.Blast import NCBIXML

from collections import defaultdict

import sys
import os,shutil
import site
import argparse
import string
import numpy
import re
import subprocess
import pickle

base_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
for directory_name in os.listdir(base_path):
    site.addsitedir(os.path.join(base_path, directory_name))
print base_path
import genbank
import blast
import intergene
import genome
import filter
#import intervals
import annotated_genes


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
    intergeneObj = intergene.IntergeneHandler(intergene_file)    
    intergeneObj.getIntervals()
    intergeneDict = dict()
    for bact in bacteriocins:
        gene = bact.sbjct_id,bact.sbjct_start,bact.sbjct_end,bact.strand
        overlaps = intergeneObj.overlapIntergene(gene)
        intergeneDict[(bact.sbjct_id,bact.sbjct_start,bact.sbjct_end,bact.strand)] = overlaps
    return intergeneDict

""" Write bacteriocins to a tab delimited file """
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
            argstr   = ">%s|%s%s|%s|%s|%s|%d|%d|%s|%s\n%s\n"
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
            argstr   = ">%s|%s%d|%d|%s|%s\n%s\n"
            result_str = argstr%(bacID,
                                 organism,
                                 bacteriocin.sbjct_start,
                                 bacteriocin.sbjct_end,
                                 bacteriocin.strand,
                                 regionType,
                                 bacteriocin.sbjct)
        outHandle.write(result_str)


def writeAnnotatedGenes(annot_bact_pairs,outHandle):
    #print "Writing Annotated Genes"
    #print "Annotations",annot_bact_pairs
    for annot_bact in annot_bact_pairs:
        annot,bacteriocin = annot_bact
        annot_st,annot_end,annot_org,annot_strand,annot_locus,annot_protid,annot_seq = annot
        bacID    = bacteriocin.query_id
        organism = bacteriocin.sbjct_id
        bacID    = bacteriocin.query_id.split(' ')[0]
        organism = bacteriocin.sbjct_id.split(' ')[0]
        #organism+= annot_locus
        argstr   = ">%s|%s%s|%s|%s|%s|%s|%s|%d|%d|%s\n%s\n"
        result_str = argstr%(bacID,
                             organism,
                             bacteriocin.sbjct_start,
                             bacteriocin.sbjct_end,
                             bacteriocin.strand,
                             annot_org,
                             annot_locus,
                             annot_protid,
                             annot_st,
                             annot_end,
                             annot_strand,
                             annot_seq)        
        outHandle.write(result_str)

def main(genome_files,
         bacteriocin_file,
         intergene_file,
         annotations_file,
         bacteriocinsOut,
         annotationsOut,
         intermediate,
         bac_evalue,
         num_threads,
         formatdb,
         bacteriocin_radius,
         verbose,
         keep_tmp):
    for gnome in genome_files:
        bacthr = BacteriocinHandler(gnome,
                                    intermediate,
                                    bac_evalue,
                                    num_threads,
                                    verbose,
                                    keep_tmp)
        bacteriocins = bacthr.getAlignedBacteriocins(bacteriocin_file,
                                                     bac_evalue,
                                                     num_threads,
                                                     formatdb)
        if verbose and bacteriocins == None: print "No bacteriocins found"
        if bacteriocins == None: continue
        if verbose: print "Bacteriocins found\n","\n".join(map(str,bacteriocins))
        if verbose: print "Number of original bacteriocins",len(bacteriocins)
        intergeneDict = identifyIntergenic(bacteriocins,intergene_file)
        writeBacteriocins(bacteriocins,intergeneDict,bacteriocinsOut)
        annots = [annot for annot in annotated_genes.AnnotatedGenes(annotations_file)]
        annots,bacteriocinNeighborhoods = filter.annotatedGenes(annots,bacteriocins,bacteriocin_radius)
        annot_bact_pairs = zip(annots,bacteriocinNeighborhoods)
        writeAnnotatedGenes(annot_bact_pairs, annotationsOut)
        
        #pickle.dump(intergeneDict,open("intergene.dict",'w'))

if __name__=="__main__":
    parser = argparse.ArgumentParser(description=\
        'Finds bacteriocins')
    parser.add_argument(\
        '--genome-files', type=str, nargs="+", required=False,
        help='FASTA files containing bacterial genomes')
    parser.add_argument(\
        '--genes', type=str,required=False,default=None,
        help='A FASTA file containing all of the target genes of interest')
    parser.add_argument(\
        '--intergenes', type=str, required=False,
        help='FASTA files containing intergenic regions')
    parser.add_argument(\
        '--annotated-genes', type=str, required=False,
        help='FASTA files containing annotated genetic regions')
    parser.add_argument(\
        '--bacteriocins', type=str, required=False,
        help='The bacteriocin proteins that are to be blasted')
    #parser.add_argument(\
    #    '--gene_radius', type=int, required=False, default=-1,
    #    help='The search radius around every specified gene. By default, this filter option is off')
    parser.add_argument(\
        '--bacteriocin-radius', type=int, required=False, default=5000,
        help='The search radius around every specified bacteriocin')
    parser.add_argument(\
        '--gene-evalue', type=float, required=False, default=0.00001,
        help='The evalue for gene hits')
    parser.add_argument(\
        '--bac-evalue', type=float, required=False, default=0.00001,
        help='The evalue for bacteriocin hits')
    #parser.add_argument(\
    #    '--num-threads', type=int, required=False, default=1,
    #    help='The evalue for bacteriocin hits')
    parser.add_argument(\
        '--intermediate', type=str, required=False,default='.',
        help='Directory for storing intermediate files')
    parser.add_argument(\
        '--output', type=str, required=False,
        help='The output file basename for filtered annotationed regions and bacteriocins')
    parser.add_argument(\
        '--formatdb', action='store_const', const=True, default=False,
        help='Indicates if formatdb should be run')
    parser.add_argument(\
        '--verbose', action='store_const', const=True, default=False,
        help='Messages for debugging')
    #parser.add_argument(\
    #    '--keep-tmp', action='store_const', const=True, default=False,
    #    help='Indicates if formatdb should be run')
    parser.add_argument(\
        '--test', action='store_const', const=True, default=False,
        help='Run unittests')

    blast.addArgs(parser)
    args = parser.parse_args()
    #outHandle.write("bacteriocin\tbacteriocin_location\torganism\tgene\tbacterciocin_sequence\n")
    if not args.test:
        bacteriocinsOut = open("%s.bacteriocins.txt"%(args.output),'w')
        #filteredOut = open("%s_filtered.txt"%(args.output),'w')
        annotationsOut  = open("%s.annotated.txt"%(args.output),'w')
        main(args.genome_files,
             args.bacteriocins,
             args.intergenes,
             args.annotated_genes,
             bacteriocinsOut,
             #filteredOut,
             annotationsOut,
             args.intermediate,
             args.bac_evalue,
             args.num_threads,
             args.formatdb,
             #args.gene_radius,
             args.bacteriocin_radius,
             args.verbose,
             args.keep_tmp)
    else:
        del sys.argv[1:]
        import unittest
        import test_modules
        import test_genbank
        import annotated_genes
        import intergene
        
        class TestPipeline(unittest.TestCase):
            def setUp(self):
                self.root = os.environ['BACFINDER_HOME']
                self.exampledir = "%s/example/Streptococcus_pyogenes"%self.root
                self.bacdir = "%s/bacteriocins"%self.root
                self.annotated_genes = "test_genes.fa"
                annotated_genes.go(self.exampledir,self.annotated_genes) 
                self.genome_files = test_modules.getFNA(self.exampledir)
                self.bacteriocins = "%s/bagel.fa"%self.bacdir
                self.genes = "%s/genes.fa"%self.bacdir
                self.intergenes = "test_intergenes.fa"
                intergene.go(self.root,self.intergenes)
                self.bacteriocinsOut = "test_out_bacteriocins.txt"
                #self.filteredOut,
                self.annotationsOut = "neighbor_genes.txt"
                self.intermediate = "intermediate"
                if not os.path.exists(self.intermediate):
                    os.mkdir(self.intermediate)
                self.gene_evalue = 0.000001
                self.bac_evalue = 0.000001
                self.num_threads = 10
                self.formatdb = True
                #self.gene_radius = 50000
                self.bacteriocin_radius = 50000
                self.verbose = True
                self.keep_tmp = False
                
            def tearDown(self):
                os.remove(self.annotated_genes)
                os.remove(self.intergenes)
                os.remove(self.bacteriocinsOut)
                os.remove(self.annotationsOut)
                shutil.rmtree(self.intermediate)
                pass
        
            def testrun(self):
                main(self.genome_files,
                     self.bacteriocins,
                     self.intergenes,
                     self.annotated_genes,
                     open(self.bacteriocinsOut,'w'),
                     open(self.annotationsOut,'w'),
                     self.intermediate,
                     self.bac_evalue,
                     self.num_threads,
                     self.formatdb,
                     #self.gene_radius,
                     self.bacteriocin_radius,
                     self.verbose,
                     self.keep_tmp)
                self.assertTrue(os.path.getsize(self.annotated_genes) > 0)
                self.assertTrue(os.path.getsize(self.intergenes) > 0)
                self.assertTrue(os.path.getsize(self.bacteriocinsOut) > 0)
                self.assertTrue(os.path.getsize(self.annotationsOut) > 0)
                
        class TestFilters(unittest.TestCase):
            def setUp(self):
                test_input = test_genbank.yeast
                self.test_file = "test.gbk"
                self.out_file = "out.fa"
                handle = open(self.test_file,'w')
                handle.write(test_input)
                handle.close()
                annotated_genes.parseAnnotations("NC_12345",self.test_file,open(self.out_file,'w'))
            def tearDown(self):
                os.remove(self.test_file)
                os.remove(self.out_file)
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
                filtered,hoods = filter.bacteriocins(bacteriocins,genes,radius)
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
                filtered,hoods = filter.bacteriocins(bacteriocins,genes,radius)
                #print '\n'.join(map(str,filtered))
                #print '\n'.join(map(str,hoods))
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
                                                sbjct_id = "gi|123|gb|NC_12345.1",
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
                                                sbjct_id = "gi|123|gb|NC_12345.1",
                                                sbjct="ACGTACGTTT",
                                                sbjct_start = 450,
                                                sbjct_end   = 460,
                                                strand = "-")]
                radius = 100
                annots = [A for A in annotated_genes.AnnotatedGenes(self.out_file)]
                filtered,hoods = filter.annotatedGenes(annots,bacteriocins,radius)
                self.assertEquals(1,len(filtered))
                self.assertEquals(1,len(hoods)) 
                self.assertTrue(bacteriocins[0] in hoods)
                
                     
        unittest.main()
