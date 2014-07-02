from Bio.Seq import Seq
from Bio import SeqIO
import os,sys,site

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
def format(seqin,width=60):
    seq = []
    j = 0
    for i in xrange(width,len(seqin),width):
        seq.append(seqin[j:i])
        j = i
    seq.append(seqin[j:])
    return '\n'.join(seq)
if __name__=="__main__":
    import unittest
    class TestFasta(unittest.TestCase):
        def setUp(self):
            entries = ['>testseq1',
                       'AGCTACT',
                       '>testseq2',
                       'AGCTAGCT',
                       '>testseq2',
                       'AAGCTAGCT'
                       ]
            self.fasta = "test.fa"
            self.revfasta = "rev.fa"
            open(self.fasta,'w').write('\n'.join(entries))
        def tearDown(self):
            if os.path.exists(self.fasta):
                os.remove(self.fasta)
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
    unittest.main()