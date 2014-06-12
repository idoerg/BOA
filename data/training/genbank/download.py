from Bio import Entrez
from Bio import SeqIO, SeqFeature
from Bio.SeqRecord import SeqRecord

Entrez.email = "mortonjt@miamioh.edu" 

GIs = ["2664256","19577293","239976949","4731431","385829589","209558587","30018278"]
files = ["Y12234.gbk","AJ438950.gbk","GG688628.gbk","AJ438950.gbk","NC_017486.gbk","NC_011375.gbk","NC_004722.gbk"]

for i in range(len(GIs)):
    handle = Entrez.efetch(db="nucleotide", id=GIs, rettype="gbwithparts", retmode="text")
    
    open(files[i],'w').write(str(handle.read()))
    seq_record = SeqIO.parse(open(files[i]),"genbank")
    #for feature in seq_record.features:
    #    try:
    #        protid = feature.qualifiers["protein_id"]
    #        pfile = "%s.gbk"%protid
    #        phandle = Entrez.efetch(db="nucleotide",
    #                                id=protid,rettype="gbwithparts",retmode="text")
    #        open(pfile,'w').write(str(phandle.read()))
    #    except Exception as e:
    #        continue
            