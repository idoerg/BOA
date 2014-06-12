from Bio import Entrez

Entrez.email = "mortonjt@miamioh.edu" 

handle = Entrez.efetch(db="nucleotide", id="2664256", rettype="gb", retmode="text")
open("Y12234.gbk",'w').write(str(handle.read()))

handle = Entrez.efetch(db="nucleotide", id="19577293", rettype="gb", retmode="text")
open("AJ438950.gbk",'w').write(str(handle.read()))

handle = Entrez.efetch(db="nucleotide", id="239976949", rettype="gb", retmode="text")
open("GG688628.gbk",'w').write(str(handle.read()))

handle = Entrez.efetch(db="nucleotide", id="4731431", rettype="gb", retmode="text")
open("AJ438950.gbk",'w').write(str(handle.read()))

handle = Entrez.efetch(db="nucleotide", id="385829589", rettype="gb", retmode="text")
open("NC_017486.gbk",'w').write(str(handle.read()))

handle = Entrez.efetch(db="nucleotide", id="209558587", rettype="gb", retmode="text")
open("NC_011375.gbk",'w').write(str(handle.read()))

handle = Entrez.efetch(db="nucleotide", id="30018278", rettype="gb", retmode="text")
open("NC_004722.gbk",'w').write(str(handle.read()))
