
"""
The purpose of this is to create phylogenetic trees through iTOL

Given 16SRNAs and a bunch of operons output the following

1. A fasttree for all of the species that have these operons
2. iTOL friendly format with 
   a. number of operons
   b. composition of functional genes
"""

import os,sys,site
base_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
for directory_name in os.listdir(base_path):
    site.addsitedir(os.path.join(base_path, directory_name))
import quorum
import subprocess
from Bio import SeqIO, SeqFeature
from Bio.SeqRecord import SeqRecord
from fasttree import *
from raxml import *
from collections import *
from muscle import Muscle
from mafft import MAFFT
from collections import defaultdict
acc_reg = re.compile("accession=(\S+).\d+|")
func_reg = re.compile("function=(\S+)|")
start_reg = re.compile("start=(\d+)")
end_reg = re.compile("end=(\d+)")

def formatName(name):
    name = name.replace("'",'')
    name = name.replace("\"",'')
    name = name.replace(",",'')
    name = name.replace("]",'')
    name = name.replace("[",'')
    name = name.replace("(",'')
    name = name.replace(")",'')
    name = name.replace("\\",'')
    name = name.replace("/",'')
    name = name.replace(";",'')
    name = name.replace(":",'')
    name = name.replace(">",'')
    name = name.replace("<",'')
    name = name.replace("+",'')
    name = name.replace("-",'')
    name = name.replace("=",'')
    name = name.replace("|",'')
    name = name.replace("{",'')
    name = name.replace("}",'')
    name = name.replace("?",'')
    name = name.replace("~",'')
    name = name.replace("`",'')
    name = name.replace("*",'')
    name = name.replace("&",'')
    name = name.replace("^",'')
    name = name.replace("%",'')
    name = name.replace("$",'')
    name = name.replace("#",'')
    name = name.replace("@",'')
    name = name.replace("!",'')
    name = name.replace(",",'')
    name = name.replace(" ",'_')
    name = '_'.join(name.split("_")[:2])
    return name

class iTOL():
    
    def __init__(self,operonFile,allrrna,rrnaFile=None,alignFile=None,treeFile=None):
        
        self.operonFile = operonFile
        self.allrrna = allrrna
        basename = os.path.basename(operonFile)
        #Store RRNA fasta file
        if rrnaFile==None:
            self.rrnaFile = "%s.rrna"%basename
        else:
            self.rrnaFile = rrnaFile
            #For multiple alignment
        if alignFile == None:
            self.alignFile = "%s.align"%basename
        else:
            self.alignFile = alignFile    
        #For fasttree
        if treeFile == None:
            self.treeFile = "%s.tree"%basename
        else:
            self.treeFile = treeFile
        self.rrna_dict = SeqIO.to_dict(SeqIO.parse(open(self.allrrna,'r'), "fasta"))
        
        
    def setOperonFile(self,operonFile):
        self.operonFile = operonFile
            
    def cleanUp(self):
        if os.path.exists(self.treeFile):
            os.remove(self.treeFile)
        if os.path.exists(self.rrnaFile):
            os.remove(self.rrnaFile)
        if os.path.exists(self.alignFile):
            os.remove(self.alignFile)
            
    """ Extract 16S rnas """
#     def getRRNAs(self):
#         rrnas = []
#         seen = set()
#         rrna_dict = SeqIO.to_dict(SeqIO.parse(open(self.allrrna,'r'), "fasta"))
#         with open(self.operonFile,'r') as handle:
#             for ln in handle:
#                 if ln[0]=="-": continue
#                 ln = ln.rstrip()
#                 toks = ln.split('|')
#                 acc,clrname,full_evalue,hmm_st,hmm_end,env_st,env_end,description=toks
#                 description = formatName(description)
#                 accession = acc.split('.')[0]
#                 if accession in rrna_dict:
#                     name = description.split(',')[0]
#                     name = formatName(name)
#                     if name not in seen:
#                         record = rrna_dict[accession]
#                         record.id = name
#                         rrnas.append(record)
#                         seen.add(name)
#                 else:
#                     print "Accession %s is missing"%accession
#         SeqIO.write(rrnas, open(self.rrnaFile,'w'), "fasta")
#    
    def getRRNAs(self):
        rrnas = []
        seen = set()
        for record in SeqIO.parse(open(self.operonFile,'r'),'fasta'):
            acc = acc_reg.findall(record.id)[0]
                    
            accession = acc.split('.')[0]
            if accession in self.rrna_dict:
                record = self.rrna_dict[accession]
                acc,taxa = record.description.split("\t")
                #record.id = formatName(taxa)
                record.id = taxa
                if taxa not in seen:
                    rrnas.append(record)
                    seen.add(taxa)
            else:
                print "Accession %s is missing"%accession
        SeqIO.write(rrnas, open(self.rrnaFile,'w'), "fasta")
        
    """ Build phylogenetic tree """
    def buildTree(self,module=subprocess,TREE=UnAlignedRaxml,MSA=Muscle,iters=4,threads=8,hours=12):
        ft = TREE(self.rrnaFile,self.treeFile)
        ft.align(module=module,MSA=MSA,iters=iters,threads=threads,hours=hours) #Run multiple sequence alignment and spit out aligned fasta file
        proc=ft.run(module=module) #Run fasttree on multiple alignment and spit out newick tree
        proc.wait()
        self.treeFile=ft.treeFile
        ft.cleanUp() #Clean up!
        
    """ Pick top k biggest operons """
    def sizeFilter(self,filterout,k=100):
        buf = []
        outhandle = open(filterout,'w')
        lengths = []
        with open(self.operonFile,'r') as handle:
            for ln in handle:
                if ln[0]=='-':
                    lengths.append(len(buf))
                    buf = []                    
                else:
                    buf.append(ln)
        lengths = sorted(lengths,reverse=True)
        if k>len(lengths): threshold = lengths[-1]
        else:  threshold = lengths[k]
        with open(self.operonFile,'r') as handle:
            for ln in handle:
                if ln[0]=='-':
                    if len(buf)>threshold:
                        for line in buf:
                            outhandle.write(line)
                        outhandle.write("----------\n")
                    buf = []                    
                else:
                    buf.append(ln)
        outhandle.close()
        
    """ Spit out iTOL file for functions distribution """        
    def functionDistribution(self,itolout):
        outhandle = open(itolout,'w')
        outhandle.write("LABELS\timmunity\tmodifier\tregulator\ttoxin\ttransport\n")
        outhandle.write("COLORS\t#0000ff\t#00ff00\t#ff0000\t#ff00ff\t#ff8c00\n")
        seen = set()    
        func_counts = defaultdict(Counter)
        for record in SeqIO.parse(open(self.operonFile,'r'),'fasta'):
            acc = acc_reg.findall(record.id.split("|")[0])[0]
            acc = acc.split(".")[0]
            curCluster = record.id.split("|")[7]
            #if acc not in self.rrna_dict:
            #    print "Missing accession",acc,record.id
            #    continue
            #taxa_record = self.rrna_dict[acc]
            env_st  = start_reg.findall(record.id)[0]
            env_end = end_reg.findall(record.id)[0]
            #taxa = taxa_record.description.split("\t")[1]
            function = func_reg.findall(record.id.split("|")[1])[0]
            taxa = acc
            if (taxa,env_st,env_end) in seen: continue
            if taxa not in func_counts:
                func_counts[taxa]=Counter({'immunity':0,'modifier':0,'regulator':0,'toxin':0,'transport':0,})
            func_counts[taxa][function] += 1
        for k,cnts in func_counts.items():
            functions = cnts.items()
            functions = sorted(functions,key=lambda x:x[0])
            functions,counts = zip(*functions)
            outhandle.write("%s\t%s\n"%(k,'\t'.join(map(str,counts))))
        outhandle.close()

    """ Spit out iTOL file for operon distribution """        
    def operonDistribution(self,itolout):
        outhandle = open(itolout,'w')
        outhandle.write("LABELS\tOperons\n")
        outhandle.write("COLORS\t#0000ff\n")
        operon_counts = defaultdict(set)
        for record in SeqIO.parse(open(self.operonFile,'r'),'fasta'):
            acc = acc_reg.findall(record.id.split("|")[0])[0]
            curCluster = record.id.split("|")[7]
            acc = acc.split(".")[0]
            taxa = acc
            #if acc not in self.rrna_dict:
            #    print "Missing accession",acc,record.id
            #    continue
            #taxa_record = self.rrna_dict[acc]
            env_st  = start_reg.findall(record.id)[0]
            env_end = end_reg.findall(record.id)[0]
            #taxa = taxa_record.description.split("\t")[1]
            cluster = record.id.split("|")[7]
            operon_counts[taxa].add(cluster)
        for k,cnts in operon_counts.items():
            outhandle.write("%s\t%d\n"%(k,len(operon_counts[k])))
        outhandle.close()

    """ Spit out graphlan file for functions distribution """        
    def functionRings(self,itolout):
        outhandle = open(itolout,'w')
        #outhandle.write("LABELS\timmunity\tmodifier\tregulator\ttoxin\ttransport\n")
        #outhandle.write("COLORS\t#0000ff\t#00ff00\t#ff0000\t#ff00ff\t#ff8c00\n")
        seen = set()
        colors = {"immunity":"#0000ff","modifier":"#00ff00","regulator":"#ff0000","toxin":"#ff00ff","transport":"#ff8c00"}
        func_counts = defaultdict(Counter)
        operon_counts = defaultdict(set)
        for record in SeqIO.parse(open(self.operonFile,'r'),'fasta'):
            acc = acc_reg.findall(record.id.split("|")[0])[0]
            acc = acc.split(".")[0]
            curCluster = record.id.split("|")[7]
            #if acc not in self.rrna_dict:
            #    print "Missing accession",acc,record.id
            #    continue
            #taxa_record = self.rrna_dict[acc]
            env_st  = start_reg.findall(record.id)[0]
            env_end = end_reg.findall(record.id)[0]
            #taxa = taxa_record.description.split("\t")[1]
            function = func_reg.findall(record.id.split("|")[1])[0]
            taxa = acc
            if (taxa,env_st,env_end) in seen: continue
            if taxa not in func_counts:
                func_counts[taxa]=Counter({'immunity':0,'modifier':0,'regulator':0,'toxin':0,'transport':0,})
            func_counts[taxa][function] += 1
            cluster = record.id.split("|")[7]
            operon_counts[taxa].add(cluster)
            
        for k,cnts in func_counts.items():
            functions = cnts.items()
            functions = sorted(functions,key=lambda x:x[0])
            #functions,counts = zip(*functions)
            level=1
            for function,counts in functions:
                outhandle.write("%s\tring_height\t%d\t%lf\n"%(k,level,counts/10.0))
                outhandle.write("%s\tring_color\t%d\t%s\n"%(k,level,colors[function]))
                level+=1
            outhandle.write("%s\tring_height\t%d\t%lf\n"%(k,level,len(operon_counts[k])/10.0))

            #outhandle.write("%s\t%s\n"%(k,'\t'.join(map(str,counts))))
            accs = [k]*5
            
        outhandle.close()
        
if __name__=="__main__":
    import unittest
    
    class TestiTOL(unittest.TestCase):
        def setUp(self):            
            self.rrnaFile="test.fa"
            self.operonFile = "test_operons.txt"
            #self.accSeqs= "acc.fa"
            rrnas   = [">CP003901\tbug1",
                       "AACGAACGCTGGCGGCATGCCTAACACATGCAAGTCGAACGAGACCTTCGGGTCTAGTGGCGCACGGGTGCGTAACGCGTGGGAATCTGCCCTTGGGTACGGAATAACAGTTAGAAATGACTGCTAATACCGTATAATGACTTCGGTCCAAAGATTTATCGCCCAGGGATGAGCCCGCGTAGGATTAGCTTGTTGGTGAGGTAAAGGCTCACCAAGGCGACGATCCTTAGCTGGTCTGAGAGGATGATCAGCCACACTGGGACTGAGACATGGCCCAGACTCCTACGGGAGGCAGCAGTGGGGAATATTGGACAATGGGCGAAAGCCTGATCCAGCAATGCCGCGTGAGTGATGAAGGCCTTAGGGTTGTAAAGCTCTTTTACCCGGGATGATAATGACAGTACCGGGAGAATAAGCCCCGGCTAACTCCGTGCCAGCAGCCGCGGTAATACGGAGGGGGCTAGCGTTGTTCGGAATTACTGGGCGTAAAGCGCACGTAGGCGGCTTTGTAAGTTAGAGGTGAAAGCCCGGGGCTCAACTCCGGAATTGCCTTTAAGACTGCATCGCTAGAATTGTGGAGAGGTGAGTGGAATTCCGAGTGTAGAGGTGAAATTCGTAGATATTCGGAAGAACACCAGTGGCGAAGGCGACTCACTGGACACATATTGACGCTGAGGTGCGAAAGCGTGGGGAGCAAACAGGATTAGATACCCTGGTAGTCCACGCCGTAAACGATGATGACTAGCTGTCGGGGCGCTTAGCGTTTCGGTGGCGCAGCTAACGCGTTAAGTCATCCGCCTGGGGAGTACGGCCGCAAGGTTAAAACTCAAAGAAATTGACGGGGGCCTGCACAAGCGGTGGAGCATGTGGTTTAATTCGAAGCAACGCGCAGAACCTTACCAGCGTTTGACATGGTAGGACGGTTTCCAGAGATGGATTCCTACCCTTACGGGACCTACACACAGGTGCTGCATGGCTGTCGTCAGCTCGTGTCGTGAGATGTTGGGTTAAGTCCCGCAACGAGCGCAACCCTCGTCTTTGGTTGCTACCATTTAGTTGAGCACTCTAAAAAAACTGCCGGTGATAAGCCGGAGGAAGGTGGGGATGACGTCAAGTCCTCATAGCCCTTACGCGCTGGGCTACACACGTGCTACAATGGCGGTGACAGAGGGCAGCAAACCCGCGAGGGTGAGCTAATCTCCAAAAGCCGTCTCAGTTCGGATTGTTCTCTGCAACTCGAGAGCATGAAGGCGGAATCGCTAGTAATCGCGGATCAGCACGCCGCGGTGAATACGTTCCCAGGCCTTGTACACACCGCCCGTCACATCACGAAAGTCGGTTGCACTAGAAGTCGGTGGGCTAACCCGCAAGGGAGGCAGCCGCCTAAAGTGTGATCGGTAATTGGGGTG",
                       ">CP003902\tbug2",
                       "AGAGTTTGATCCTGGCTCAGAATGAACGCTGGCGGCGTGCCTAACACATGCAAGTCGTACGAGAAATCCCGAGCTTGCTTGGGAAAGTAAAGTGGCGCACGGGTGAGTAACGCGTGGGTAACCCACCCCCGAATTCGGGATAACTCCGCGAAAGCGGTGCTAATACCGGATAAGACCCCTACCGCTTCGGCGGCAGAGGTAAAAGCTGACCTCTCCATGGAAGTTAGCGTTTGGGGACGGGCCCGCGTCCTATCAGCTTGTTGGTGGGGTAACAGCCCACCAAGGCAACGACGGGTAACTGGTCTGAGAGGATGATCAGTCACACTGGAACTGGAACACGGTCCAGACTCCTACGGGAGGCAGCAGTGAGGAATTTTGCGCAATGGGCGAAAGCCTGACGCAGCAACGCCGCGTGGGTGAAGAAGGCTTTCGGGTCGTAAAGCCCTGTCAGGTGGGAAGAAACCTTTCCGGTACTAATAATGCCGGAAATTGACGGTACCACCAAAGGAAGCACCGGCCAACTCCGTGCCAGCAGCCGCGGTAATACGGAGGGTGCAAGCGTTGTTCGGAATTATGGGGCGTAAAGAGCGTGTGGGCGGTTAGGAAAGTCAGATGTGAAAGCCCTGGGCTCAACCCAGGAAGTGCATTTGAAACTGCCTAACTTGAGTACGGGAGAGGAAGGGGGAATTCCCGGTGTAGAGGTGAAATTCGTAGATATCGGGAGGAATACCGGTGGCGAAGGCGCCCTTCTGGACCGATACTGACGCTGAGACGCGAAAGCGTGGGGAGCAAACAGGATTAGATACCCTGGTAGTCCACGCCGTAAACGATGAGCACTAGGTGTAGCGGGTATTGACCCCTGCTGTGCCGTAGCTAACGCATTAAGTGCTCCGCCTGGGGATTACGGTCGCAAGACTAAAACTCAAAGGAATTGACGGGGGCCCGCACAAGCGGTGGAGCATGTGGTTTAATTCGACGCAACGCGAAGAACCTTACCTGGGCTTGACATCCCCGGACAGCCCTGGAAACAGGGTCTCCCACTTCGGTGGGCTGGGTGACAGGTGCTGCATGGCTGTCGTCAGCTCGTGTCGTGAGATGTTGGGTTAAGTCCCGCAACGAGCGCAACCCCTGCCTTTAGTTGCCATCATTTAGCTGGGCACTCTAAAGGGACTGCCGGTGTTAAACCGGAGGAAGGTGGGGACGACGTCAAGTCCTCATGGCCTTTATGCCCAGGGCTACACACGTGCTACAATGGGCGGTACAAAGGGCAGCGACATCGTGAGGTGAAGCAAATCCCAAAAAACCGCTCTCAGTTCGGATCGGAGTCTGCAACTCGACTTCGTGAAGGTGGAATCACTAGTAATCGTGGATCAGCATGCCACGGTGAATACGTTCCCGGGCCTTGTACA",
                       ">CP002079\tbug3",
                       "GCTGGCGGCGTGCCTAACACATGTAAGTCGAACGGGACTGGGGGCAACTCCAGTTCAGTGGCAGACGGGTGCGTAACACGTGAGCAACTTGTCCGACGGCGGGGGATAGCCGGCCCAACGGCCGGGTAATACCGCGTACGCTCGTTTAGGGACATCCCTGAATGAGGAAAGCCGTAAGGCACCGACGGAGAGGCTCGCGGCCTATCAGCTAGTTGGCGGGGTAACGGCCCACCAAGGCGACGACGGGTAGCTGGTCTGAGAGGATGGCCAGCCACATTGGGACTGAGACACGGCCCAGACTCCTACGGGAGGCAGCAGTGGGGAATCTTGCGCAATGGCCGCAAGGCTGACGCAGCGACGCCGCGTGTGGGATGACGGCCTTCGGGTTGTAAACCACTGTCGGGAGGAACGAATACTCGGCTAGTCCGAGGGTGACGGTACCTCCAAAGGAAGCACCGGCTAACTCCGTGCCAGCAGCCGCGGTAATACGGAGGGTGCGAGCGTTGTCCGGAATCACTGGGCGTAAAGGGCGCGTAGGTGGCCCGTTAAGTGGCTGGTGAAATCCCGGGGCTCAACTCCGGGGCTGCCGGTCAGACTGGCGAGCTAGAGCACGGTAGGGGCAGATGGAATTCCCGGTGTAGCGGTGGAATGCGTAGATATCGGGAAGAATACCAGTGGCGAAGGCGTTCTGCTGGGCCGTTGCTGACACTGAGGCGCGACAGCGTGGGGAGCAAACAGGATTAGATACCCTGGTAGTCCACGCCGTAAACGATGGACACTAGACGTCGGGGGGAGCGACCCTCCCGGTGTCGTCGCTAACGCAGTAAGTGTCCCGCCTGGGGAGTACGGCCGCAAGGCTGAAACTCAAAGGAATTGACGGGGGCCCGCACAAGCGGTGGAGCATGTGGTTTAATTCGAAGCAACGCGAAGAACCTTACCTGGGCTTGACATGCTGGTGCAAGCCGGTGGAAACATCGGCCCCTCTTCGGAGCGCCAGCACAGGTGCTGCATGGCTGTCGTCAGCTCGTGTCGTGAGATGTTGGGTTAAGTCCCGCAACGAGCGCAACTCTCGCTCCCAGTTGCCAGCGGTTCGGCCGGGGACTCTGGGGGGACTGCCGGCGTTAAGCCGGAGGAAGGTGGGGACGACGTCAAGTCATCATGGCCCTTACGTCCAGGGCGACACACGTGCTACAATGCCTGGTACAGCGCGTCGCGAACTCGCAAGAGGGAGCCAATCGCCAAAAGCCGGGCTAAGTTCGGATTGTCGTCTGCAACTCGACGGCATGAAGCCGGAATCGCTAGTAATCGCGGATCAGCCACGCCGCGGTGAATACGTTCCCGGGCCTTGTACACACCGCCCGTCACGCCATGGAAGCCGGAGGGACCCGAAACCGGTGGGCCAACCGCAAGGGGGCAGCCGTCTAAGGT",
                       ">CP002089\tbug4",
                       "AGAGTTTGATCATGGCTCAGGATGAACGCTAGCGGCAGGCCTAACACATGCAAGTCGAGGGGTAGAGGCTTTCGGGCCTTGAGACCGGCGCACGGGTGCGTAACGCGTATGCAATCTGCCTTGTACTAAGGGATAGCCCAGAGAAATTTGGATTAATACCTTATAGTATATAGATGTGGCATCACATTTCTATTAAAGATTTATCGGTACAAGATGAGCATGCGTCCCATTAGCTAGTTGGTATGGTAACGGCATACCAAGGCAATGATGGGTAGGGGTCCTGAGAGGGAGATCCCCCACACTGGTACTGAGACACGGACCAGACTCCTACGGGAGGCAGCAGTGAGGAATATTGGTCAATGGGCGCAAGCCTGAACCAGCCATGCCGCGTGCAGGATGACGGTCCTATGGATTGTAAACTGCTTTTGTACGGGAAGAAACACTCCTACGTGTAGGGGCTTGACGGTACCGTAAGAATAAGGATCGGCTAACTCCGTGCCAGCAGCCGCGGTAATACGGAGGATCCAAGCGTTATCCGGAATCATTGGGTTTAAAGGGTCCGTAGGCGGTTTTATAAGTCAGTGGTGAAATCCGGCAGCTCAACTGTCGAACTGCCATTGATACTGTAGAACTTGAATTACTGTGAAGTAACTAGAATATGTAGTGTAGCGGTGAAATGCTTAGATATTACATGGAATACCAATTGCGAAGGCAGGTTACTAACAGTATATTGACGCTGATGGACGAAAGCGTGGGGAGCGAACAGGATTAGATACCCTGGTAGTCCACGCCGTAAACGATGGATACTAGCTGTTTGGCAGCAATGCTGAGTGGCTAAGCGAAAGTGTTAAGTATCCCACCTGGGGAGTACGAACGCAAGTTTGAAACTCAAAGGAATTGACGGGGGCCCGCACAAGCGGTGGAGCATGTGGTTTAATTCGATGATACGCGAGGAACCTTACCAGGGCTTAAATGTAGAGTGACAGGACTGGAAACAGTTTTTTCTTCGGACACTTTACAAGGTGCTGCATGGTTGTCGTCAGCTCGTGCCGTGAGGTGTCAGGTTAAGTCCTATAACGAGCGCAACCCCTGTTGTTAGTTGCCAGCGAGTAATGTCGGGAACTCTAACAAGACTGCCGGTGCAAACCGTGAGGAAGGTGGGGATGACGTCAAATCATCACGGCCCTTACGTCCTGGGCTACACACGTGCTACAATGGCCGGTACAGAGAGCAGCCACCTCGCGAGGGGGAGCGAATCTATAAAGCCGGTCACAGTTCGGATTGGAGTCTGCAACCCGACTCCATGAAGCTGGAATCGCTAGTAATCGGATATCAGCCATGATCCGGTGAATACGTTCCCGGGCCTTGTACACACCGCCCGTCAAGCCATGGAAGCTGGGGGTACCTGAAGTCGGTGACCGCAAGGAGCTGCCTAGGGTAAAACTGGTAACTGGGGCTAAGTCGTACAAGGTAGCCGTA",
                       ">CP002987",
                       "CCTAATGCATGCAAGTCGAACGCAGCAGGCGTGCCTGGCTGCGTGGCGAACGGCTGACGAACACGTGGGTGACCTGCCCCGGAGTGGGGGATACCCCGTCGAAAGACGGGACAATCACGCATACGCTCTTTGGAGGAAAGCCATCCGGCGCTCTGGGAGGGGCCTGCGGCCCATCAGGTAGTTGGTGTGGTAACGGCGCACCAAGCCAATGACGGGTACCCGGTCTGAGAGGACGACCGGCCAGACTGGAACTGCGACACGGCCCAGACTCCTACGGGAGGCAGCAGCAAGGAATTTTCCCCAATGGGCGCAAGCCTGAGGCAGCAACGCCGCGTGCGGGATGACGGACTTCGGGTTGTAAACCGCTTTTCGGGGGGACAACCCTGACGGTACCCCCGGAACAAGCCCCGGCTAACTCTGTGCCAGCAGCCGCGGTAAGACAGAGGGGGCAAGCGTTGTCCGGAGTCACTGGGCGTAAAGCGCGCGCAGGCGGCTGCCTAAGTGTCGTGTGAAAGCCCCCGGCTCAACCGGGGGAGGCCATGGCAAACTGGGTGGCTCGAGCGGCGGAGAGGTCCCTCGAATTGCCGGTGTAGCGGTGAAATGCGTAGAGATCGGCAGGAAGACCAAGGGGGAAGCCAGGGGGCTGGCCGCCGGCTGACGCTGAGGCGCGACAGCGTGGGGAGCAAACCGGATTAGATACCCGGGTAGTCCACGCCGTAAACGATGACCACTCGGCGTGTGGCGACTATTAACGTCGCGGCGCGCCCTAGCTCACGCGATAAGTGGTCCGCCTGGGAACTACGAGCGCAAGCTTAAAACTCAAAGGAATTGACGGGGGCCCGCACAAGCAGCGGAGCGTGTGGTTTAATTCGACGCAACCCGCAGAACCTTACCCAGACTGGACATGACGGTGCAGACGGCGGAAACGTCGTCGCCTGCGAGGGTCCGTCACAGGTGCTGCATGGCTGTCGTCAGCTCGTGTCGTGAGATGTTGGGTTAAGTCCCGCAACGAGCGCAACCCCTGCGGTTAGTTACCCGTGTCTAACCGGACTGCCCTTCGGGGAGGAAGGCGGGGATGACGTCAAGTCCGCATGGCCCTTACGTCTGGGGCGACACACACGCTACAATGGCGCCGACAATGCGTCGCTCCCGCGCAAGCGGATGCTAATCGCCAAACGGCGCCCCAGTGCAGATCGGGGGCTGCAACTCGCCCCCGTGAAGGCGGAGTTGCTAGTAACCGCGTATCAGCCATGGCGCGGTGAATACGTACCCGGGCCTTGTACACACCGCCCGTCACGTCATGGAGTTGTCAATGCCTGAAGTCCGCCAGCTAACC"
                       ]
            open(self.rrnaFile,'w').write('\n'.join(rrnas))
            operons = [ ">accession=CP003901.1|function=modifier|start=558785|end=559735|strand=+|score=745.5|protein_id=AFV37692.1|cluster_1105|-",
                        "MSFFTKEQQPKENCPPITVEKARQLFEFNTNHLSLSDYHHQTVLKTSKQLVAQHLMPNAT",
                        "DNLSQHFLMNYKANNNYLGFQASIVDFFTDSAVANFSSSYVYESQEKIIRLPKPTKISTA",
                        "LSTCIIKRRSHRQFSDRQMPLQDLSNILYYACGVSSQASIRDGASDKITLRNCASGGGLY",
                        "PIHLVFYARNISKLIDGFYEYLPYQHALRCYRHSSEENVRDFAEYGAINAENCNIIIIYV",
                        "YHYMKNTRKYGNQATAYAFIESGEIAQNIQLTATALTYGSIDIGGYNKEYLQELLDLDGL",
                        "GEHVIHMTLVGTKESQ",
                        ">accession=CP003901.1|function=immunity|start=562143|end=562814|strand=+|score=513.3|protein_id=AFV37695.1|cluster_1105|-",
                        "MPLSIQCLNLCFLLVTFCPSIPMQAIFGKEDSGYAFNLIGFLRATLIYDILALVSIYVLS",
                        "PQITLSLESIDSKTFFMGLVFCVLIVLIELVFLHGLRCWQKKQWLPATFSFVGTTNDWSK",
                        "IGYPLLLALFEETIYRFLWFNILAFQWHLPTIIVLIVTSFCYALNHLLMGKSIFYAKLVT",
                        "GIIYGSIYMLTSQLWLVVIMHVGGNLLVECLSHLQTKKKKEVT",
                        ">accession=CP003901.1|function=modifier|start=559732|end=560790|strand=+|score=829.3|protein_id=AFV37693.1|cluster_1105|-",
                        "MKYQLNSNVRVVTFQDTFCFRKGIWDFNEAILDLTQEPQNLKIVYQEIVSQLVKGVAIDT",
                        "ETYENTLEPETFANLMEVISGLYYNDMLMLEDDYNLEENVMKILMGNFRFMAQAGHMVSN",
                        "DPVLFISDSTYVNESAELLAEQLHLQLQVASDDLKMLIQQTDVSSRLDALEHRRNMNCLS",
                        "DALCDYQSIIICQERLNIMMLRHLNEVSVAMKKQLVIGFVDGPFLHTCTLNPPHSADFDS",
                        "LERRVLARLQDSTLYQHFANQVLPATQDVSQAYLPLLNVLMNLVVSEAFIIAQTGSSKFE",
                        "GRLLSIYLPTLEIQVQDILKMSNSQTQGALAKLKYEDQQISTREIVKKLLDE",
                        ">accession=CP003901.1|function=toxin|start=558402|end=558563|strand=+|score=75.0|protein_id=AFV37691.1|cluster_1105|-",
                        "MLKFTSNILATSVAETTQVAPGGCCCCCTTCCFSIATGSGNSQGGSGSYTPGK",
                        ">accession=CP003901.1|function=modifier|start=560810|end=562168|strand=+|score=945.3|protein_id=AFV37694.1|cluster_1105|-",
                        "MLYYYPSFNHIFDELKSLSGNRTGILNQSQVPVCNHPHDVYLKSITGQMPDYHKQFIGEL",
                        "SQVSYHIIGYGSYYEEALIKYLGESIERYATVIAGDLLSDRIVYASYNELKLLHKVMPLE",
                        "YLQVFTQEQIALSCDLQMMMCDKMVTENDVLGWVKCPMFFEDAEMYVPAQMLCVGYKTNE",
                        "TVGERRIIPGFSTGTASHKTLEAAMCNSLIEYIQIDSMMLSWYTKKPCPKIIVDDPDIEV",
                        "ILEEARLGKDSLYDIIPIDMTVGEDNPLYTFGIILKNKYDEGPYLLFGVQAGLDPKHALL",
                        "RGIMEASAISYSYYYNLLYQKASLANIECEEPLFLDLDSNVFYYAHPKDQDHKWKAFEPL",
                        "ISGEVLLSDLEDHSGKDKKEDLKTLLAYAKKVSPNAVFLDITPPEALEKGWYVTRVLMPE",
                        "LLEMCIPAFPFANHPRMRQFGGVTNAFVHPMP",
                        ">accession=CP003901.1|function=modifier|start=562811|end=563494|strand=+|score=472.7|protein_id=AFV37696.1|cluster_1105|-",
                        "MMLLVLLSFLGLTTLWLTSYRRRCLARWHLQWLVALSYQDFLDVLLSLFQFVVIILVLFF",
                        "YSATINLGEVLTFLTQTSWHWQILCYLVLYLMAIIEMTLLVLILIFDVLLQKDSRLLFKK",
                        "ITWLPFRQDKPVLSILLLGTVILVDSLFYLGLLILLGEQSLTSLATLILGYGIVKACRYS",
                        "GWLNVLLAFCLFTIVGLWAVTATLLYGWLVGAMILTVTYLMISCKEY",
                        ">accession=CP003901.1|function=immunity|start=563517|end=564440|strand=+|score=126.2|protein_id=AFV37697.1|cluster_1105|-",
                        "MSFVQLTNVVKSYKNGKKAVNDVSLSIEAGNIYGLLGPNGAGKSTLINLILGLIPLSSGK",
                        "ITVLGQSQKTIRKISSQIGYVPQDIAVYPDLTAYENVELFGSLYGLKGAQLKKQVLKSLE",
                        "FVGLHSQAKQFPSQFSGGMKRRLNIACALVHSPKLIIFDEPTVGIDPQSRNHILESIRLL",
                        "NKEGATVIYTTHYMEEVEALCDYIFIMDHGQVIEEGPKFELEKRYVANLANQIIVTLTDS",
                        "RHLELADKPDWSLIEDGEKLMLKIDNSDMTSVVHQLTQANITFSEIRHNHLNLEEIFLHL",
                        "TGKKLRD",
                        ">accession=CP003901.1|function=transport|start=564449|end=565576|strand=+|score=859.6|protein_id=AFV37698.1|cluster_1105|-",
                        "MVLFHLIKKESLQIFRNRTALLMMVIFPILMIVILSFAFKSSFNTATTVPKLTIRYQLEG",
                        "EKTDYQKNFLAFLKVLNQKLHLETKPSNSLEKDRQRVSEGALTAVLEVKKNQTIKVITNN",
                        "INQQNADLINMLVKNYVDNAKTYDSIAALYPQQLNHIRKRSVDYVKVSSIQTSKGMTSAD",
                        "YYAISMFTMITFYSMMSAMNLVLSDRQQRITNRIHLTGVSPSFLVFGKLIGAMLATTVQL",
                        "SLLYIFTRFVLRVNWGTNEWMLIGITASLVYLSVAIGIGLGISIKNEAFLTVASNTIIPI",
                        "FAFLGGSYVPLTTLHSSIINQLSNISPIKWVNDSLFYLIFGGQYNPIPVTLIVNISIGTI",
                        "FIILALIGMRKQVTT",
                        ">accession=CP003901.1|function=transport|start=565573|end=566691|strand=+|score=737.1|protein_id=AFV37699.1|cluster_1105|-",
                        "MICFIKTLFVKIKRKKTSYVTFFLMPILTTLLALSLSFSNNNQAKIGILDKDNSQISKQF",
                        "IAQLKQNKKYDIFTKIKKEHIDHYLQDKSLEAVLTIDKGFSDKVLQGKSQKLNIRSIANS",
                        "EITEWVKAQTNYLLENYNIIGDVALGNEDTFNRILQKNQQLNYDVKQVTLTDRSRSKAVS",
                        "STTTGFLLILMLGSTSVIYSGILADKSSQLYHRLMLSNLSRFRYMLSYVCVGFVAFTIQI",
                        "VIMLSLLKVFNISFFVPTSLLLIIFFLFSLLAIGFGLLIGAITQNSQQSSQLANLIVMPT",
                        "SMLAGCLWPLSITPSYMQAIGKLLPQNWVLSAIAIFQSGGTLSQAWPYLLALMGTALALI",
                        "SFSSLLLKPTKL",
                        ">accession=CP003902.1|function=transport|start=1592456|end=1593193|strand=-|score=65.8|protein_id=AFV38739.1|cluster_1107|-",
                        "MKQLVLKDVCKKYPNQLNYALDHINLIVEKGEFVAVMGRSGSGKTTLLNVTSIIDKIDSG",
                        "NIYCADKEISTFSDKEATNFRKNDIGFVFQDYMLLDSLTIRENISVALSLKNIDSSKIDD",
                        "LINSYAKRFNLYEQLKKYPYQLSGGQRQRVSIIRAIIKEPEIIFADEPTGALDLKSSEET",
                        "MMILSEINKTEKVTILMVTHDVLSASYADRVILLKDGKLHMEIDKKDCGESFYDVISQAL",
                        "SDRGE",
                        ">accession=CP003902.1|function=toxin|start=1595433|end=1595579|strand=-|score=133.7|protein_id=AFV38742.1|cluster_1107|-",
                        "MKNSKDILTNAIEEVSEKELMEVAGGKKGSGWFATITDDCPNSVFVCC",
                        ">accession=CP003902.1|function=transport|start=1593190|end=1593651|strand=-|score=243.6|protein_id=AFV38740.1|cluster_1107|-",
                        "MYNILLGRESVTEKQVIEICKRVSLYEDIRSMPMKFHTPLFRDNPSLSGGQKQRISLARE",
                        "LVTTPRILVLDEPTSALDVKTERIIQKNVEALHCTRVLVTHRLNTVEKADKILIMDNGKI",
                        "IDYGNHHYLYKNNKDYCDLYDSYMNKYQEEEVK",
                        ">accession=CP003902.1|function=regulator|start=1588945|end=1590507|strand=-|score=1075.5|protein_id=AFV38737.1|cluster_1107|-",
                        "MIRIKNITKSKFFGTAIILLQQLIALLILVYNRENLSLLFSEKVALVMTLIDTAFIWLAT",
                        "ILRQKQGDIFKRIISIISLTIWQYLVSVLTRGTPLFLSSGLQIILLYCYTVEITNLILYG",
                        "HKQFKDKLDKGLLIIIFVSITSLFVNRILFNFLFLMVFTILHLYPLLVIVLYYRSFRQQI",
                        "SVVRRSLILFSLLLLVILGSELYGEMLDVNQAFNNLGWYLFPLIMSVIYYFKTIHDKLSF",
                        "VTQRWLGDYKARIELLFLLLILCWIVLMKVLVKDFLLFFIVVDASTLFCLVVISCIFYYL",
                        "ENSKQNFDYENRRLNYFMKSEENMRVEFSNYLHDDVLQNIIAIKNLFSLENSNITHGFIV",
                        "NELNDLVSGIRDEMDTYHPIVPANQTMKENIQSLFDDIVKSRKSNTLLYFNCSDNMVVPS",
                        "PYGDIVYRFIKELINNAIKYGDGKDIRLSLTIQSDIIIIEESNQVVEKVHSINYGRGLKS",
                        "FQETLAAFDGDLELQMDTKQFTIRILLPIDWKLCYEDFIN",
                        ">accession=CP003902.1|function=transport|start=1590547|end=1592454|strand=-|score=1379.6|protein_id=AFV38738.1|cluster_1107|-",
                        "MIWSITKSNIKKNFSLYRIYFLATIGLLSIFIAFLNFISDKIITEKIGDSGQALVIANGS",
                        "LIFLIVFLVVFLIYFNNFFVKKRSQELGVLAILGFSKRELTKLLTLENLVILVLSYLVSL",
                        "LLGPTLYFLAVLAITHLLNLTMEVQWFITVNEIIESLGILVVVFLINVITNGLIISKQSL",
                        "IEFVNFSRKAEKKIKIRKVRAIIAITALLLSYILCLATVFSSTRNMLLSIGMVPVSLLII",
                        "VLVVLGTVFTIRYGLAFVVSLLKENKKRLYRPLSNIIYPKFNYRIATKNKLLTVLGGLLT",
                        "VTVSVAGMMVMLYAYSLNGIERLTPSAIEYNVESENGQVNVTTILENDQVSLVDVGLLRL",
                        "NTIPEVTITDSGQTIPYFDIINYSDYKELMKAQGRTNSIEGSKSLPLLINYYPTEISLGK",
                        "TFNLGNAYDVTVKQVSTNNVFSFSTSVTTLVVSDKLYAKLSSRFPEKEMTIRTFNGTSIR",
                        "SSEAFYNQFSMVPDVISSYSKEHTVKTANIATYIFITFLSILFIICTGSILYFTSLIEIM",
                        "ENKEEYGYLSKLGYSKKMIHRILRYETGILFLIPVFIGIVNGGMLLIYYKYLFMDTLVAG",
                        "NIIMLSLLLCLLFFLIIYGTFYVLTLRLVTSIIKN",
                        ">accession=CP003902.1|function=regulator|start=1588359|end=1588964|strand=-|score=426.0|protein_id=AFV38736.1|cluster_1107|-",
                        "MKILLIDDHRLFAKSIQLLFQQYDEVDVIDTITSHFNDVTIDLSKYDIILLDINLANISK",
                        "ENGLEIAKELIQSTPHLKVVMLTGYVKSIYRERAKKVGAYGFVDKNIDPKQLISILKKVD",
                        "SGKKYFEQIESQDYVESLTDQEIAILNLSKKGFSIKEIEETLQISRRTVFNHLTHIYSKL",
                        "LVNNKQEAIYKAEQLGYFMDF",
                        ">accession=CP003902.1|function=modifier|start=1593725|end=1595350|strand=-|score=101.5|protein_id=AFV38741.1|cluster_1107|-",
                        "MIKRNELKLEYQQHQFYEYFNSIFDYDILDTSIRETSEVLRKKTHKLFYSEFENQLFETI",
                        "MFLSMKTLVLDINHFSKEIENKSEAYEQYIQQIREENGINHFFDRYPYLLKQINKEVGLI",
                        "EESYSLLFDRFLEDLSEIKSCFNISEPLSNVAFSLGDSHSKKQTVVKIAFKEKSVYYKPK",
                        "SYHSHSILLELTSLLKSSNIPSFSLPKSLVKADYCWQLGVAYTSSNKDEVAKIYFKYGVL",
                        "AAFSEIFSITDLHMENVIVSGGDLYLIDVETFFQRKLNVQNQNFEGITVDTYQRIYETSL",
                        "SNGLFPVQFEKNSAPNVSGISGKGGKRKKGKYELINKNRGDMKLVKVDYFQEDGFNIPTL",
                        "NGKVVEPLDYANEIISGFRECYIFLLSQRSKIKEIVEGFPELKSRVPFRNTSDYGKFLQA",
                        "STNPKYLFSEKKRKNLFSILYETKHIEHFIVDNEIKDLMNGDIPYFSMDTRGNVYNSVGT",
                        "LIGNLGDTTSLFDSITILNDERLKFTCELLEIVLKKPIKHWEREKGKSYQFLSISSEHYL"]
                    
 
            open(self.operonFile,'w').write('\n'.join(operons))
            self.itol = None
            
        def tearDown(self):
            os.remove(self.rrnaFile)
            os.remove(self.operonFile)
            #if not self.itol==None:
            #    self.itol.cleanUp()
        def testGetRRNAS(self):
            self.itol = iTOL(self.operonFile,self.rrnaFile)
            self.itol.getRRNAs() 
            self.assertTrue(os.path.getsize(self.itol.rrnaFile)>0)
            seq_records = list(SeqIO.parse(open(self.itol.rrnaFile,'r'), "fasta"))
            self.assertEquals(len(seq_records),2)
        def testBuildTree(self):
            self.itol = iTOL(self.operonFile,self.rrnaFile)
            self.itol.getRRNAs()
            self.itol.buildTree(TREE=UnAlignedFastTree)
            self.assertTrue(os.path.getsize(self.itol.treeFile)>0)

        def testBootstrap(self):
            
            pass
            
    unittest.main()
    
    
    
