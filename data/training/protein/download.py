from Bio import Entrez
from Bio import SeqIO, SeqFeature
from Bio.SeqRecord import SeqRecord
import re
Entrez.email = "mortonjt@miamioh.edu" 

fname = "../training_proteins.txt"

with open(fname,'r') as handle:
    for ln in handle:
        ln = ln.rstrip()
        if ln[0]=="#":continue
        toks = re.split('\s+',ln)
        protid = toks[0]

        handle = Entrez.efetch(db="nucleotide", id=protid, rettype="gbwithparts", retmode="text")
        fout = "%s.gbk"%protid
        open(fout,'w').write(str(handle.read()))
            
