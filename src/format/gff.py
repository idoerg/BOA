"""
Parses gff files and creates a fasta file complete with all of the genomic coordinates and sequences
"""


import os,sys,site
import numpy
base_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
for directory_name in os.listdir(base_path):
    site.addsitedir(os.path.join(base_path, directory_name))
import fasta
import faa
import argparse
import hmmer
from Bio.Seq import Seq
from Bio import SeqIO
import re
from collections import defaultdict
from bx.intervals import *

prot_reg = re.compile("Name=([A-Za-z0-9_]+.\d);")

class GFF():
    def __init__(self,gff_file,
                 output_file="",fasta_file="",fasta_index="",
                 createIndex=False,
                 sixFrame=False):
        self.sixFrame = sixFrame
        self.gff_file = gff_file
        self.fasta = fasta_file
        self.fasta_index = fasta_index
        self.indexer = fasta.Indexer(self.fasta,self.fasta_index)
        if createIndex: self.indexer.index()
        self.indexer.load()
        self.output_file = output_file
        self.inttrees = defaultdict(IntervalTree)
        self.proteins = {}
        #self.indextree()
    """ Clean up """
    def __del__(self):
        del self.inttrees
        """Parses gff file and spits out a fasta file for all of the predicted orfs"""
    def parse(self,outhandle=None):
        if outhandle==None:
            outhandle = open(self.output_file,'w')
        with open(self.gff_file,'r') as handle:
            for ln in handle:
                if ln[0]=='#': continue
                ln = ln.rstrip()
                toks = ln.split('\t')
                species,_,type,st,end,_,strand,_,text = toks
                if type!="gene":continue
                st,end = map(int,[st,end])
                if strand=='+':
                    sequence = self.indexer.fetch(species,st,end)
                else:
                    sequence = self.indexer.reverse_fetch(species,st,end)
                outhandle.write(">%s|%s|%d|%d|%s\n%s\n"%(species,text,st,end,strand,sequence))
   
    """ Create an interval tree"""
    def indextree(self):
        
        with open(self.gff_file,'r') as handle:
            for ln in handle:
                if ln[0]=='#': continue
                
                ln = ln.rstrip()
                toks = ln.split('\t')
               
                species,_,type,st,end,_,strand,_,text = toks
                if type=='CDS': 
                    st,end = map(int,[st,end])
                    
                    self.inttrees[(species,strand)].add( st,end,
                                                         (st,end,text) )
    """ Create dictionary indexed by protein id""" 
    def indexdb(self):
        self.proteins = {}
        with open(self.gff_file,'r') as handle:
            for ln in handle:
                if ln[0] == '#': continue
                ln = ln.rstrip()
                toks = ln.split('\t')
                species,_,type,st,end,_,strand,_,text = toks
                if type!='CDS': continue
                if prot_reg.findall(text)==[]:
                    print "No protid in",text
                    continue
                protid = prot_reg.findall(text)[0]
                self.proteins[protid] = (species,st,end,strand)
                
    
    """ Figures out which orfs overlap with HMMER hits """
    def call_orfs(self,hits,faaindex):
        records = []
        for hit in hits:
            protid,clrname,full_evalue,hmm_st,hmm_end,env_st,env_end,description=hit
            if protid not in self.proteins:
                print "%s is missing !"%protid
            else:
                species,st,end,strand = self.proteins[protid]
                env_st,env_end = map(int,[st,end])
                seq = faaindex[protid]
                seqinfo = (species,clrname,full_evalue,hmm_st,hmm_end,env_st,env_end,description,strand,protid) 
                records.append( (seqinfo, str(seq)) ) 
        return records
    
    """
    
    def call_orfs(self,hits):
        #Create interval tree from gff file
        
        newHits = []
        
        for hit in hits:
            acc,clrname,full_evalue,hmm_st,hmm_end,env_st,env_end,description=hit
            if self.sixFrame:
                curOrg = fasta.getName(acc)
                hitSt,curStrand  = self.indexer.sixframe_to_nucleotide(acc,env_st)
                hitEnd,curStrand = self.indexer.sixframe_to_nucleotide(acc,env_end)
            else:
                curOrg = acc
                
            orfs = self.inttrees[(curOrg,curStrand)].find(hitSt,hitEnd)
            
            if len(orfs)>0:
                orf_st,orf_end,_ = orfs[0]
                env_st,env_end = orf_st,orf_end
                newHits.append((acc,clrname,full_evalue,hmm_st,hmm_end,env_st,env_end,description))
        return newHits
    """
    
    """ Figures out which orfs overlap HMMER hits AND retreives orf from faa file
        Note: This assumes that coordinates are in genome coordinates
        
        Returns a list of sequence ids and their corresponding protein sequence
    """ 
    """
    def translate_orfs(self,hits,faaindex):
        
        records = []
        for hit in hits:
            acc,clrname,full_evalue,hmm_st,hmm_end,env_st,env_end,description=hit
            if self.sixFrame:
                curOrg = fasta.getName(acc)
            else:
                curOrg = acc
            #hitSt,curStrand  = self.indexer.sixframe_to_nucleotide(acc,env_st)
            #hitEnd,curStrand = self.indexer.sixframe_to_nucleotide(acc,env_end)
            hitSt,hitEnd = map(int,[env_st,env_end])
            curStrand = fasta.getStrand(acc)
            orfs = self.inttrees[(curOrg,curStrand)].find(hitSt,hitEnd)
            if len(orfs)>0:
                orf_st,orf_end,text = orfs[0]
                if prot_reg.findall(text)==[]:
                    print "No protid in",text
                    continue
                protid = prot_reg.findall(text)[0]
                
                seq = faaindex[protid]
                #print record
                seqinfo = map(str,[acc,clrname,full_evalue,env_st,env_end,protid])
                
                records.append( (seqinfo, str(seq)) ) 
        return records
    """
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
        go(args.gff_files,args.fasta,args.fasta_index,args.output,args.create_index)
        
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
                '\t'.join(map(str,['CP002987.1','Genbank\tregion',1,10000,'.','+','.','ID=id0;Name=ANONYMOUS;Dbxref=taxon:931626;Is_circular=true;'])),
                '\t'.join(map(str,['CP002987.1','Genbank\tgene',51,100,'.','+','.','ID=gene0;Name=dnaA;gbkey=Gene;gene=dnaA;locus_tag=Awo_c00010'])),
                '\t'.join(map(str,['CP002987.1','Genbank\tCDS',51,100,'.','+','0','ID=cds0;Name=AFA46815.1;Parent=gene0;Dbxref=NCBI_GP:AFA46815.1;'])),             
                '\t'.join(map(str,['CP002987.1','Genbank\tgene',51,100,'.','-','.','ID=gene0;Name=dnaA;gbkey=Gene;gene=dnaA;locus_tag=Awo_c00011'])),
                '\t'.join(map(str,['CP002987.1','Genbank\tCDS',51,100,'.','-','0','ID=cds0;Name=AFA46815.1;Parent=gene0;Dbxref=NCBI_GP:AFA46815.1;']))             
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
                gff.indextree()
                gff.parse() 
                result = ''.join(open(self.outfasta,'r').readlines())
                seqs = [str(s.seq) for s in SeqIO.parse(self.outfasta,"fasta")]
                headers = [s.id for s in SeqIO.parse(self.outfasta,"fasta")]
                self.assertEquals(seqs[0],'ACGTACGTAG'*5)
                self.assertEquals(seqs[1],'CTACGTACGT'*5)
                self.assertEquals(headers[0],'CP002987.1|ID=gene0;Name=dnaA;gbkey=Gene;gene=dnaA;locus_tag=Awo_c00010|51|100|+')
                self.assertEquals(headers[1],'CP002987.1|ID=gene0;Name=dnaA;gbkey=Gene;gene=dnaA;locus_tag=Awo_c00011|51|100|-')
        """
        class TestTranslateOrfs(unittest.TestCase):
            def setUp(self):
                indexes = [ '\t'.join(map(str,('CP002279.1_1',2294815, 185896721,60,61))),
                             '\t'.join(map(str,('CP002279.1_2',2294815, 188229850,60,61))),
                             '\t'.join(map(str,('CP002279.1_3',2294814, 190562979,60,61))),
                             '\t'.join(map(str,('CP002279.1_4',2294814, 192896107,60,61))),
                             '\t'.join(map(str,('CP002279.1_5',2294815, 195229235,60,61))),
                             '\t'.join(map(str,('CP002279.1_6',2294815, 197562364,60,61)))]
                
                gff = '\n'.join([
                '##gff-version 3',
                '#!gff-spec-version 1.2',
                '#!processor NCBI annotwriter',
                '##sequence-region CP002987.1 1 4044777',
                '##species http://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=931626',
                '\t'.join(map(str,['CP002279.1','Genbank\tregion',1,10000,'.','+','.','ID=id0;Name=ANONYMOUS;Dbxref=taxon:931626;Is_circular=true;'])),
                '\t'.join(map(str,['CP002279.1','Genbank\tgene',100,160,'.','+','.','ID=gene0;Name=dnaA;gbkey=Gene;gene=dnaA;locus_tag=Awo_c00010'])),
                '\t'.join(map(str,['CP002279.1','Genbank\tCDS',100,160,'.','+','0','ID=cds0;Name=AFA46815.1;Parent=gene0;Dbxref=NCBI_GP:AFA46815.1;'])),             
                '\t'.join(map(str,['CP002279.1','Genbank\tgene',1500,1800,'.','+','.','ID=gene0;Name=dnaA;gbkey=Gene;gene=dnaA;locus_tag=Awo_c00011'])),
                '\t'.join(map(str,['CP002279.1','Genbank\tCDS',1500,1800,'.','+','0','ID=cds0;Name=AFA46816.1;Parent=gene0;Dbxref=NCBI_GP:AFA46815.1;'])),             
                '\t'.join(map(str,['CP002279.1','Genbank\tgene',3200,3800,'.','+','.','ID=gene0;Name=dnaA;gbkey=Gene;gene=dnaA;locus_tag=Awo_c00011'])),
                '\t'.join(map(str,['CP002279.1','Genbank\tCDS',3200,3800,'.','+','0','ID=cds0;Name=AFA46817.1;Parent=gene0;Dbxref=NCBI_GP:AFA46815.1;']))             
                ])
                self.queries   = [('CP002279.1_1','toxin.fa.cluster2.fa',0,0,1,125,150,
                                    'Mesorhizobium opportunistum WSM2075, complete genome'),
                                   ('CP002279.1_1','transport.fa.cluster2.fa',0,0,1,570,600,
                                    'Mesorhizobium opportunistum WSM2075, complete genome'),
                                   ('CP002279.1_1','transport.fa.cluster2.fa',0,0,1,1220,1280,
                                    'Mesorhizobium opportunistum WSM2075, complete genome')] 
                seqs=['>CP002279.1',
                      'ACGTACGTAGACGTACGTAGACGTACGTAGACGTACGTAGACGTACGTAGACGTACGTAGACGTACGTAGACGTACGTAGACGTACGTAGACGTACGTAGACGTACGTAG\n'*100]
                prots=['>gi|093|gb|AFA46815.1|',
                       'TYVDVRRRTX'*2,
                       '>gi|093|gb|AFA46816.1|',
                       'TYVDVRRRTX'*10,
                       '>gi|093|gb|AFA46817.1|',
                       'TYVDVRRRTX'*20]
                      
                self.gff = "test.gff"
                self.fasta = "test.fasta"
                self.faidx = "test.faidx"
                self.outfasta = "out.fasta"
                open(self.gff,'w').write(gff)
                open(self.fasta,'w').write('\n'.join(seqs))
                print '\n'.join(indexes)
                open(self.faidx,'w').write('\n'.join(indexes))
                self.out = "prot_out.fasta"
                self.faa = "test.faa"
                self.faaidx = "test.faaidx"
                
                open(self.faa,'w').write('\n'.join(prots))
            def tearDown(self):
                os.remove(self.fasta)
                os.remove(self.faidx)
                os.remove(self.gff)
                #os.remove(self.out)
                
            def test(self):
                tmpfile = "tmp%d.faa"%(os.getpid())
                faa.reformat(self.faa,tmpfile)
                os.rename(tmpfile,self.faa)
                gff = GFF(self.gff,self.outfasta,self.fasta,self.faidx,False)
                gff.indextree()
                faaindex = fasta.Indexer(self.faa,self.faaidx)
                faaindex.index()
                faaindex.load()
                records = gff.translate_orfs(self.queries,faaindex)
                ids,seqs = zip(*records)
                
                self.assertGreater(len(records),0)
                self.assertEquals(ids[0],
                                  ['CP002279.1_1', 'toxin.fa.cluster2.fa', '0', '125', '150', 'AFA46815.1'])
                                   
                self.assertEquals(str(seqs[0]),'TYVDVRRRTX'*2)
        """
                    
        class TestIndexTree(unittest.TestCase):
            def setUp(self):
                indexes = [ '\t'.join(map(str,('CP002279.1_1',2294815, 185896721,60,61))),
                             '\t'.join(map(str,('CP002279.1_2',2294815, 188229850,60,61))),
                             '\t'.join(map(str,('CP002279.1_3',2294814, 190562979,60,61))),
                             '\t'.join(map(str,('CP002279.1_4',2294814, 192896107,60,61))),
                             '\t'.join(map(str,('CP002279.1_5',2294815, 195229235,60,61))),
                             '\t'.join(map(str,('CP002279.1_6',2294815, 197562364,60,61)))]
                
                gff = '\n'.join([
                '##gff-version 3',
                '#!gff-spec-version 1.2',
                '#!processor NCBI annotwriter',
                '##sequence-region CP002987.1 1 4044777',
                '##species http://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=931626',
                '\t'.join(map(str,['CP002279.1','Genbank\tregion',1,10000,'.','+','.','ID=id0;Name=ANONYMOUS;Dbxref=taxon:931626;Is_circular=true;'])),
                '\t'.join(map(str,['CP002279.1','Genbank\tgene',51,100,'.','+','.','ID=gene0;Name=dnaA;gbkey=Gene;gene=dnaA;locus_tag=Awo_c00010'])),
                '\t'.join(map(str,['CP002279.1','Genbank\tCDS',51,100,'.','+','0','ID=cds0;Name=AFA46815.1;Parent=gene0;Dbxref=NCBI_GP:AFA46815.1;'])),             
                '\t'.join(map(str,['CP002279.1','Genbank\tgene',1551,2000,'.','+','.','ID=gene0;Name=dnaA;gbkey=Gene;gene=dnaA;locus_tag=Awo_c00011'])),
                '\t'.join(map(str,['CP002279.1','Genbank\tCDS',1551,2000,'.','+','0','ID=cds0;Name=AFA46816.1;Parent=gene0;Dbxref=NCBI_GP:AFA46816.1;'])),             
                '\t'.join(map(str,['CP002279.1','Genbank\tgene',3551,4000,'.','+','.','ID=gene0;Name=dnaA;gbkey=Gene;gene=dnaA;locus_tag=Awo_c00011'])),
                '\t'.join(map(str,['CP002279.1','Genbank\tCDS',3551,4000,'.','+','0','ID=cds0;Name=AFA46817.1;Parent=gene0;Dbxref=NCBI_GP:AFA46817.1;']))             
                ])
                self.queries   = [('AFA46815.1','toxin.fa.cluster2.fa',0,0,1,55,150,
                                    'Mesorhizobium opportunistum WSM2075, complete genome'),
                                   ('AFA46816.1','transport.fa.cluster2.fa',0,0,1,570,600,
                                    'Mesorhizobium opportunistum WSM2075, complete genome'),
                                   ('AFA46817.1','transport.fa.cluster2.fa',0,0,1,1220,1280,
                                    'Mesorhizobium opportunistum WSM2075, complete genome')] 
                seqs=['>CP002279.1',
                      'ACGTACGTAGACGTACGTAGACGTACGTAGACGTACGTAGACGTACGTAGACGTACGTAGACGTACGTAGACGTACGTAGACGTACGTAGACGTACGTAGACGTACGTAG\n'*100]
                self.gff = "test.gff"
                self.fasta = "test.fasta"
                self.faidx = "test.faidx"
                self.outfasta = "out.fasta"
                open(self.gff,'w').write(gff)
                open(self.fasta,'w').write('\n'.join(seqs))
                print '\n'.join(indexes)
                open(self.faidx,'w').write('\n'.join(indexes))
                
                prots=['>gi|093|gb|AFA46815.1|',
                       'TYVDVRRRTX'*2,
                       '>gi|093|gb|AFA46816.1|',
                       'TYVDVRRRTX'*10,
                       '>gi|093|gb|AFA46817.1|',
                       'TYVDVRRRTX'*20]
                self.faa = "test.faa"
                self.faaidx = "test.faaidx"
                open(self.faa,'w').write('\n'.join(prots))
                self.maxDiff=1000
            def tearDown(self):
                os.remove(self.fasta)
                os.remove(self.faidx)
                os.remove(self.gff)
                pass
            def test(self):
                tmpfile = "tmp%d.faa"%(os.getpid())
                faa.reformat(self.faa,tmpfile)
                os.rename(tmpfile,self.faa)
                faaindex = fasta.Indexer(self.faa,self.faaidx)
                faaindex.index()
                faaindex.load()
                gff = GFF(self.gff,self.outfasta,self.fasta,self.faidx,False)
                #gff.indextree()
                gff.indexdb()
                hits = gff.call_orfs(self.queries,faaindex)
                print hits
                ids,seqs = zip(*hits)
                correct_queries = [('CP002279.1','toxin.fa.cluster2.fa',0,0,1,51,100,
                                    'Mesorhizobium opportunistum WSM2075, complete genome'),
                                   ('CP002279.1','transport.fa.cluster2.fa',0,0,1,1551,2000,
                                    'Mesorhizobium opportunistum WSM2075, complete genome'),
                                   ('CP002279.1','transport.fa.cluster2.fa',0,0,1,3551,4000,
                                    'Mesorhizobium opportunistum WSM2075, complete genome')] 
                self.assertItemsEqual(ids,correct_queries)
                                  
        unittest.main()
        
        
        
        
        
        