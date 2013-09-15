"""
A blast pipeline to blast bacteriocin proteins against intergenic regions


Word of caution: blastx is not thread safe :(
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

parser = argparse.ArgumentParser(description=\
    'Finds intergenic regions from genback file')
parser.add_argument(\
    '--genbank-path', type=str, required=True,
    help='The path of the genbank file')
parser.add_argument(\
    '--intergenes', type=str, required=True,
    help='The path of the fasta file containing the intergenes')
parser.add_argument(\
    '--bacteriocins', type=str, required=False,
    help='FASTA file containing candidate bacteriocin proteins to blast against')
parser.add_argument(\
    '--output-file', type=str, required=True,
    help='The output file containing the BLAST output')
parser.add_argument(\
    '--keep-tmp', action='store_const', const=True, default=False,
    help='Keeps temporary files such as blast database and blast output xml')
parser.add_argument(\
    '--radius', type=str, required=False, default=10000,
    help='The search radius around every specified gene')
parser.add_argument(\
    '--num-threads', type=int, required=False, default=1,
    help='The number of threads to be run on blast')

args = parser.parse_args()

import genbank


class Record():
    """
    Records information for each blast entry that is important for multiple alignment
    """
    def __init__(self,
                 record_number,     #Record number
                 description,       #Description of the sequence
                 expected_value,    #Expected value of the alignment
                 score,             #Score of the alignment
                 query,             #Portion of the adapter sequence
                 query_start,       #Beginning of the adapter alignment
                 query_end,         #End of the adapter alignment
                 sbjct,             #Portion of the subject sequence
                 sbjct_start,       #Beginning of the subject alignment
                 sbjct_end):        #End of the subject alignment
        self.record_number = record_number
        self.description = description
        self.score = score
        self.expected_value = expected_value
        self.query = query
        self.query_start = query_start
        self.query_end = query_end
        self.sbjct = sbjct
        self.sbjct_start = sbjct_start
        self.sbjct_end = sbjct_end

    def __str__(self):
        string = '****Alignment****\n' \
          +'Record%d\n  '%self.record_number \
          +'sequence%s\n: '+self.description \
          +'e value:%lf\n'%self.expected_value \
          +'Score:%lf\n'%self.score \
          +'Query:\t%s\t%s\t%s\n'%(self.query_start,self.query,self.query_end) \
          +'Sbjct:\t%s\t%s\t%s\n'%(self.sbjct_start,self.sbjct,self.sbjct_end)
        return string

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


class BLAST(object):
    def __init__(self,bacteriocins_file,intergenes):
        self.pid = os.getpid() #Use current pid to name temporary files
        self.protein_db = bacteriocins_file
        self.genomic_query = "genomicQuery.%d.fa"%self.pid
        self.blastxml = "%d.xml"%self.pid
        SeqIO.write(intergenes,self.genomic_query,"fasta")

    """
    Build database using intergenic regions
    """
    def buildDatabase(self):
        cmd="formatdb -i %s -p T -o T"%(self.protein_db)
        proc = subprocess.Popen(cmd,shell=True)
        proc.wait()

    """    Blast sequences    """
    def run(self,num_threads):
        outHandle = open(self.blastxml,'w')
        cmd="blastall -p blastx -d %s -i %s -m 7 -o %s -a %d"%(self.protein_db,self.genomic_query, self.blastxml, num_threads)
        proc = subprocess.Popen(cmd,shell=True)
        proc.wait()


    def cleanup(self):

        os.system("rm %s.*"%self.protein_db)
        os.remove(self.genomic_query)
        os.remove(self.blastxml)

    def parseBLAST(self):
        input_file = self.blastxml
        hits = []
        handle = open(input_file,'r')
        blast_hits = NCBIXML.parse(handle)
        blast_records = list(blast_hits)
        i = 0
        for record in blast_records:
            for alignment in record.alignments:
                for hsp in alignment.hsps:
                    if hsp.expect<0.04:
                        record=Record(record_number = i,
                                      description = alignment.title,
                                      expected_value = hsp.expect,
                                      score = hsp.score,
                                      query = hsp.query,
                                      query_start = hsp.query_start,
                                      query_end = hsp.query_end,
                                      sbjct = hsp.sbjct,
                                      sbjct_start = hsp.sbjct_start,
                                      sbjct_end = hsp.sbjct_end)
                        print record
                        hits.append(record)
        return hits

if __name__=="__main__":
    geneDict = genbank.GenBank(args.genbank_path)
    sagB_gene = geneDict.findGene("sagB")
    intergenes = getIntergenes(args.intergenes,sagB_gene,args.radius)
    blast_obj = BLAST(args.bacteriocins,intergenes)
    blast_obj.buildDatabase()
    blast_obj.run(args.num_threads)
    hits = blast_obj.parseBLAST()
    outHandle = open(args.output_file,'w')

    outHandle.write("\n".join( map( str, hits)))
    if not args.keep_tmp: blast_obj.cleanup()
