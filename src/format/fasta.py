from Bio.Seq import Seq
from Bio import SeqIO
import os,sys,site
base_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
for directory_name in os.listdir(base_path):
    site.addsitedir(os.path.join(base_path, directory_name))
import quorum
import subprocess
""" Remove duplicate entries"""
def removeDuplicates(items):
    uniqueDict = {tuple(x[-5:-1]):x for x in items}
    return uniqueDict.values()
""" Preprocess fasta file """
def preprocess(blastTab,fastaout):
    items = []
    with open(blastTab,'r') as handle:
        for ln in handle:
            ln = ln.rstrip()
            toks = ln.split('\t')
            assert len(toks)>=10
            toks = [tok.replace(' ','') for tok in toks] #remove white space
            items.append(tuple(toks))
    items = removeDuplicates(items)
    with open(fastaout,'w') as handle:
        for item in items:
            bacID,gi,bst,bend,bstrand,species,ast,aend,astrand,seq = item
            seqstr = ">%s|%s|%s|%s|%s|%s\n%s\n"%(bacID,gi,bst,bend,ast,aend,seq)
            handle.write(seqstr)

def reverse_complement(fastain,revfasta):
    out = open(revfasta,'w')
    for seq_record in SeqIO.parse(fastain,"fasta"):
        rev_seq = seq_record.reverse_complement()
        SeqIO.write(rev_seq,out,"fasta")
    out.close()

def remove_duplicates(fastain,fastaout):
    ids = set()
    out = open(fastaout,'w')
    for seq_record in SeqIO.parse(fastain,"fasta"):
        ID = seq_record.id
        if ID not in ids:
            ids.add(ID)
            SeqIO.write(seq_record,out,"fasta")
    out.close() 
""" Make human readable fasta format """
def format(seqin,width=60):
    seq = []
    j = 0
    for i in xrange(width,len(seqin),width):
        seq.append(seqin[j:i])
        j = i
    seq.append(seqin[j:])
    return '\n'.join(seq)
""" Flatten each fasta record into one tab delimited line """
def flatten(fastain,fastaout):
    out = open(fastaout,'w')
    for seq_record in SeqIO.parse(fastain,"fasta"):
        ID = seq_record.id
        seq = seq_record.seq
        out.write("%s\t%s\n"%(ID,seq))
    out.close() 

class Indexer():
    def __init__(self,fasta,fastaidx):
        self.fasta = fasta       #Input fasta
        self.fastaidx = fastaidx #Output fasta index
        self.faidx = {} #internal dict for storing byte offset values
        
    """ 
    Develop an index similar to samtools faidx
    speciesName: give the option to pick out only the accession id 
    """
    def index(self,accessionID=False):
        idx = []
        #name = name of sequence
        #seqLen = length of sequence without newline characters
        #lineLen = number of characters per line
        #byteLen = length of sequence, including newline characters
        #myByteoff = byte offset of sequence
        name, seqLen, byteoff, myByteoff, lineLen, byteLen = None, 0, 0, 0, 0, 0     
        index_out = open(self.fastaidx,'w')
        
        with open(self.fasta,'r') as handle:
            for ln in handle:
                lnlen = len(ln)
                
                if len(ln)==0: break
                if ln[0]==">":
                    #print ">",myByteoff,byteoff
                    if name is not None:
                        if accessionID: acc = name.split('|')[3]
                        else: acc = name 
                        index_out.write('\t'.join(map(str, [acc, seqLen, myByteoff, 
                                                            lineLen, byteLen])))
                        index_out.write('\n')
                        seqLen = 0
                    myByteoff = byteoff + lnlen
                    seqLen = 0
                    if ' ' in ln:
                        name = ln[1:ln.index(' ')].rstrip()
                    else:
                        name = ln[1:].rstrip()
                    byteoff+=lnlen
                else:
                    
                    byteLen = max(byteLen, len(ln))
                    ln = ln.rstrip()
                    lineLen = max(lineLen, len(ln))
                    seqLen += len(ln)
                    byteoff += lnlen
        if name is not None:
            if accessionID: acc = name.split('|')[3]
            else: acc = name 
            index_out.write('\t'.join(map(str, [acc, seqLen, myByteoff, 
                                                lineLen, byteLen])))
            index_out.write('\n')
        index_out.close()
        print "Index created"
    """ Load fasta index """
    def load(self):
        with open(self.fastaidx,'r') as handle:
            for line in handle:
                line=line.strip()
                cols=line.split('\t')
                
                chrom = cols[0]
                seqLen,byteOffset,lineLen,byteLen = map(int,cols[1:])
                self.faidx[chrom]=(seqLen,byteOffset,lineLen,byteLen)
    """ Retrieve a sequence based on fasta index """
    def fetch(self, defn, start, end):
        assert type(1)==type(start)
        assert type(1)==type(end)
        self.fasta_handle = open(self.fasta,'r')
        seq=""
        if not self.faidx.has_key(defn):
            raise ValueError('Chromosome %s not found in reference' % defn)
        seqLen,byteOffset,lineLen,byteLen=self.faidx[defn]
        start = start-1
        pos = byteOffset+start/lineLen*byteLen+start%lineLen
        self.fasta_handle.seek(pos)
        while len(seq)<end-start:
            line=self.fasta_handle.readline()
            line=line.rstrip() 
            seq=seq+line
        self.fasta_handle.close()
        return seq[:end-start]
    def complement(self,letter):
        basecomplement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A','R':'Y','Y':'R'}
        if letter in basecomplement:
            return basecomplement[letter]
        else:
            return letter
    """ Retrieve a sequence on the reverse strand based on fasta index """
    def reverse_fetch(self,defn,start,end):
        start,end = self.reverse(defn,end)+1,self.reverse(defn,start)+1
        assert start<=end,"Start=%d > End=%d"%(start,end)
        sequence = self.fetch(defn,start,end)
        sequence = sequence.upper()
        return ''.join(map(self.complement,sequence[::-1]))
        
    """ Returns corresponding coordinates on opposite strand """
    def reverse(self,seqname,pos):
        seqLen,byteOffset,lineLen,byteLen = self.faidx[seqname]
        assert type(seqLen)==type(int(0))
        assert type(pos)==type(int(0))
        return seqLen-pos
    def reverseOperon(self,operon):
        newOperon = []
        for gene in operon:
            s = strand(getFrame(gene[0]))
            acc,clrname,full_evalue,hmm_st,hmm_end,env_st,env_end, description = gene
            env_st,env_end = map(int,[env_st,env_end])
            env_st,env_end = self.reverse(acc,env_st),self.reverse(acc,env_end)    
            newOperon.append((acc,clrname,full_evalue,hmm_st,hmm_end,env_st,env_end, description))
        return newOperon
    
    """ Build a six-frame translation of entire fasta file"""
    def translate(self,transout,module=subprocess):
        cmd = "transeq -sequence %s -outseq %s -frame 6 -clean"%(self.fasta,transout)
        proc = module.Popen(cmd,shell=True)
        proc.wait()
        
    """ Translate six frame translated positions into nucleotides """
    def sixframe_to_nucleotide(self,seqname,aa_pos):
        frame = seqname.split('_')[-1]
        seqLen,byteOffset,lineLen,byteLen = self.faidx[seqname]
        nuc_pos = 0
        if frame=="1":
            return (aa_pos*3,"+")
        elif frame=="2":
            return (aa_pos*3+1,"+")
        elif frame=="3":
            return (aa_pos*3+2,"+")
        elif frame=="4":
            return (3*(seqLen-aa_pos),"-")
        elif frame=="5":
            return (3*(seqLen-aa_pos)+1,"-")
        elif frame=="6":
            return (3*(seqLen-aa_pos)+2,"-")
        else:
            raise Exception
""" Get frame number for EMBOSS transeq output"""
def getFrame(seqname):
    frame = seqname.split('_')[-1]
    return frame
def getName(seqname):
    name = seqname.split('_')[0]
    return name
def strand(frame):
    if int(frame)>=1 and int(frame)<=3:
        return '+'
    elif int(frame)>=4 and int(frame)<=6:
        return '-'
    else:
        raise Exception
     
"""Merge all fasta files together, create fasta index, six-frame translated
genomes and six-frame translation index"""
def go(rootdir,
       all_fasta,
       all_faidx,
       six_frame_genome,
       six_frame_faidx):
    outhandle = open(all_fasta,'w')
    for root, subFolders, files in os.walk(self.genome_dir):
        for fname in files:
            genome_files = []
            organism,ext = os.path.splitext(os.path.basename(fname))
            absfile=os.path.join(root,fname)
            if ext==".fna":
                shutil.copyfileobj(open(absfile),outhandle)
    
    indexer = Indexer(all_fasta,all_faidx)
    indexer.index()
    indexer.translate(six_frame_genomes)
    
    six_frame_index = fasta.Indexer(six_frame_genomes,
                                    six_frame_faidx) 
    six_frame_index.index()
    
if __name__=="__main__":
    parser = argparse.ArgumentParser(description=\
        'Finds bacteriocins and context genes')
    parser.add_argument(\
        '--root-dir', type=str, required=False,
        help='The output file basename for filtered annotationed regions and bacteriocins')
    parser.add_argument(\
        '--output', type=str, required=False,
        help='The output file basename for filtered annotationed regions and bacteriocins')
    parser.add_argument(\
        '--test', action='store_const', const=True, default=False,
        help='Run unittests')
    args = parser.parse_args()
    
    if not args.test:
        all_fasta = "%s_all.fna"%args.output
        all_faidx = "%s_all.fai"%args.output
        six_fasta = "%s_six.fna"%args.output
        six_faidx = "%s_six.fai"%args.output
        go(args.root_dir,
           all_fasta,
           all_faidx,
           six_fasta,
           six_faidx)
    else:    
        import unittest 
        class TestIndex(unittest.TestCase):
            def setUp(self):
                entries = ['>testseq10_1',
                           'AGCTACT',
                           '>testseq10_2',
                           'AGCTAGCT',
                           '>testseq40_2',
                           'AAGCTAGCT',
                           '>testseq40_5',
                           'AAGCTAGCT\n'*100
                           ]
                self.fasta = "test.fa"
                self.fastaidx = "test.fai"
                self.revfasta = "rev.fa"
                open(self.fasta,'w').write('\n'.join(entries))
            def tearDown(self):
                if os.path.exists(self.fasta):
                    os.remove(self.fasta)
                if os.path.exists(self.fastaidx):
                    os.remove(self.fastaidx)
                if os.path.exists(self.revfasta):
                    os.remove(self.revfasta)
            def testIndex(self):
                indexer = Indexer(self.fasta,self.fastaidx)
                indexer.index()
                indexer.load()
      
                seq = indexer.fetch("testseq10_1",1,4)
                self.assertEquals("AGCT",seq)
                seq = indexer.fetch("testseq40_5",1,13)
                self.assertEquals("AAGCTAGCTAAGC",seq)
                seq = indexer.fetch("testseq40_5",1,900)
                self.assertEquals("AAGCTAGCT"*100,seq)
            def testReverseFetch(self):
                indexer = Indexer(self.fasta,self.fastaidx)
                indexer.index()
                indexer.load()
                seq = indexer.reverse_fetch("testseq10_1",1,4)
                self.assertEquals("AGTA",seq)
                seq = indexer.reverse_fetch("testseq40_5",1,9)
                self.assertEquals("AGCTAGCTT",seq)
                seq = indexer.reverse_fetch("testseq40_5",1,900)
                self.assertEquals("AGCTAGCTT"*100,seq)
                
            def testTransform(self):
                indexer = Indexer(self.fasta,self.fastaidx)
                indexer.index()
                indexer.load()
                pos,_ = indexer.sixframe_to_nucleotide("testseq40_2", 5)
                self.assertEquals(pos,16)
                pos,_ = indexer.sixframe_to_nucleotide("testseq40_5", 5)
                self.assertEquals(pos,2686)
        class TestFasta(unittest.TestCase):
            def setUp(self):
                entries = ['>testseq1',
                           'AGCTACT',
                           '>testseq2',
                           'AGCTAGCT',
                           '>testseq2',
                           'AAGCTAGCT'
                           '>testseq3',
                           'AAGCTAGCT\n'*100
                           ]
                self.fasta = "test.fa"
                self.fastaidx = "test.fai"
                self.revfasta = "rev.fa"
                open(self.fasta,'w').write('\n'.join(entries))
            def tearDown(self):
                if os.path.exists(self.fasta):
                    os.remove(self.fasta)
                if os.path.exists(self.fastaidx):
                    os.remove(self.fastaidx)
                if os.path.exists(self.revfasta):
                    os.remove(self.revfasta)
            def test1(self):
                reverse_complement(self.fasta,self.revfasta)
                seqs = [s for s in SeqIO.parse(self.revfasta,"fasta")]
                self.assertEquals(str(seqs[0].seq),"AGTAGCT")
                self.assertEquals(str(seqs[1].seq),"AGCTAGCT")
            def test2(self):
                remove_duplicates(self.fasta,self.revfasta)
                seqs = [s for s in SeqIO.parse(self.revfasta,"fasta")]
                self.assertEquals(len(seqs),2)
            def test3(self):
                seq = "ACACGCGACGCAGCGACGCAGCAGCAGCAGCA"
                newseq = format(seq,5)
                self.assertEquals(newseq,
                                  '\n'.join([
                                  "ACACG",
                                  "CGACG",
                                  "CAGCG",
                                  "ACGCA",
                                  "GCAGC",
                                  "AGCAG",
                                  "CA"])
                                  )
    
        class TestTranslate(unittest.TestCase):
            def setUp(self):
              self.frames = ["MAIVMGRX",
                             "WPLXWAAX",
                             "GHCNGPLX",
                             "HGNYHAAS"[::-1],
                             "XPWQLPGS"[::-1],
                             "XAMTIPRQ"[::-1]
                             ]
              seq = "atggccattguaatgggccgctga"
              self.fasta = "test.fasta"
              self.sixframe = "six.fasta"
              open(self.fasta,'w').write(">test\n%s\n"%(seq))
              #open(self.sixframe,'w').write(">test_1\n%s\n"%(frame1)+\
              #                              ">test_2\n%s\n"%(frame2)+\
              #                              ">test_3\n%s\n"%(frame3)+\
              #                              ">test_4\n%s\n"%(frame4)+\
              #                              ">test_5\n%s\n"%(frame5)+\
              #                              ">test_6\n%s\n"%(frame6)
              #                              )
            def tearDown(self):
                os.remove(self.fasta)
                os.remove(self.sixframe)          
            def test1(self):
                indexer = Indexer(self.fasta,"")
                indexer.translate(self.sixframe)
                i = 0
                for record in SeqIO.parse(open(self.sixframe,'r'),'fasta'):
                    seq = record.seq
                    self.assertEquals(str(seq),self.frames[i])
                    i+=1
                
        unittest.main()
        
        
        
        
        