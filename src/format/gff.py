"""
Parses gff files and creates a fasta file complete with all of the genomic coordinates and sequences
"""


import os,sys,site
import numpy
base_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
for directory_name in os.listdir(base_path):
    site.addsitedir(os.path.join(base_path, directory_name))
import fasta
import argparse

from Bio.Seq import Seq
from Bio import SeqIO

class GFF():
    def __init__(self,gff_file,output_file,fasta_file,fasta_index,createIndex=False):
        self.gff_file = gff_file
        self.fasta = fasta_file
        self.fasta_index = fasta_index
        self.indexer = fasta.Indexer(self.fasta,self.fasta_index)
        if createIndex: self.indexer.index()
        self.indexer.load()
        self.output_file = output_file
    """Parses gff file and spits out a fasta file for all of the predicted orfs"""
    def parse(self,outhandle=None):
        if outhandle==None:
            outhandle = open(self.output_file,'w')
        with open(self.gff_file,'r') as handle:
            for ln in handle:
                if ln[0]=='#': continue
                ln = ln.rstrip()
                toks = ln.split('\t')
                
                species,type,st,end,_,strand,_,text = toks
                if type!="Genbank gene":continue
                st,end = map(int,[st,end])
                if strand=='+':
                    sequence = self.indexer.fetch(species,st,end)
                else:
                    sequence = self.indexer.reverse_fetch(species,st,end)
                outhandle.write(">%s|%s|%d|%d|%s\n%s\n"%(species,text,st,end,strand,sequence))
        
    
def go(gff_files,fasta,fasta_index,output_file,createIndex):
    outhandle = open(output_file,'w')
    for gff_file in gff_files:
        gff = GFF(gff_file,output_file,fasta,fasta_index,createIndex)
        gff.parse(outhandle)

if __name__=="__main__":
    parser = argparse.ArgumentParser(description=\
                                     'Parses gff files and creates a fasta file and index')
    parser.add_argument(\
        '--gff-files', type=str, nargs="+", required=False,
        help='Gene feature files')
    parser.add_argument(\
        '--fasta', type=str, required=False,
        help='Input fasta file')
    parser.add_argument(\
        '--fasta-index', type=str, required=False,
        help='Fasta index.  Can be empty')
    parser.add_argument(\
        '--output', type=str, required=False,
        help='Output basename of fasta file')
    parser.add_argument(\
        '--create-index', action='store_const', const=True, default=False,
        help='Create fasta index')
    parser.add_argument(\
        '--test', action='store_const', const=True, default=False,
        help='Run unittests')
    args = parser.parse_args()
    if not args.test:
        go(args.gff_file,args.output_file,args.fasta,args.fasta_index,args.create_index)
        
    else:
        del sys.argv[1:]
        import unittest
        
        class TestGFF(unittest.TestCase):
            def setUp(self):
                gff = '\n'.join([
                '##gff-version 3',
                '#!gff-spec-version 1.2',
                '#!processor NCBI annotwriter',
                '##sequence-region CP002987.1 1 4044777',
                '##species http://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=931626',
                '\t'.join(map(str,['CP002987.1','Genbank region',1,10000,'.','+','.','ID=id0;Name=ANONYMOUS;Dbxref=taxon:931626;Is_circular=true;'])),
                '\t'.join(map(str,['CP002987.1','Genbank gene',51,100,'.','+','.','ID=gene0;Name=dnaA;gbkey=Gene;gene=dnaA;locus_tag=Awo_c00010'])),
                '\t'.join(map(str,['CP002987.1','Genbank CDS',51,100,'.','+','0','ID=cds0;Name=AFA46815.1;Parent=gene0;Dbxref=NCBI_GP:AFA46815.1;'])),             
                '\t'.join(map(str,['CP002987.1','Genbank gene',51,100,'.','-','.','ID=gene0;Name=dnaA;gbkey=Gene;gene=dnaA;locus_tag=Awo_c00011'])),
                '\t'.join(map(str,['CP002987.1','Genbank CDS',51,100,'.','-','0','ID=cds0;Name=AFA46815.1;Parent=gene0;Dbxref=NCBI_GP:AFA46815.1;']))             
                ])
                seqs=['>CP002987.1',
                      'ACGTACGTAGACGTACGTAGACGTACGTAGACGTACGTAGACGTACGTAGACGTACGTAGACGTACGTAGACGTACGTAGACGTACGTAGACGTACGTAGACGTACGTAG\n'*100]
                self.gff = "test.gff"
                self.fasta = "test.fasta"
                self.faidx = "test.faidx"
                self.outfasta = "out.fasta"
                open(self.gff,'w').write(gff)
                open(self.fasta,'w').write('\n'.join(seqs))
            def tearDown(self):
                os.remove(self.fasta)
                os.remove(self.faidx)
                os.remove(self.outfasta)
                os.remove(self.gff)
            def test1(self):
                gff = GFF(self.gff,self.outfasta,self.fasta,self.faidx,True)
                gff.parse()
                result = ''.join(open(self.outfasta,'r').readlines())
                seqs = [str(s.seq) for s in SeqIO.parse(self.outfasta,"fasta")]
                headers = [s.id for s in SeqIO.parse(self.outfasta,"fasta")]
                self.assertEquals(seqs[0],'ACGTACGTAG'*5)
                self.assertEquals(seqs[1],'CTACGTACGT'*5)
                self.assertEquals(headers[0],'CP002987.1|ID=gene0;Name=dnaA;gbkey=Gene;gene=dnaA;locus_tag=Awo_c00010|51|100|+')
                self.assertEquals(headers[1],'CP002987.1|ID=gene0;Name=dnaA;gbkey=Gene;gene=dnaA;locus_tag=Awo_c00011|51|100|-')
        unittest.main()
        
        
        
        
        
        