"""
Wrapper class for genbank files
"""

import Bio
from Bio import SeqIO, SeqFeature
from Bio.SeqRecord import SeqRecord
from Bio import Entrez
from nltk import PorterStemmer
import sys
import os,site
import argparse
import glob
import sqlite3
from collections import defaultdict
base_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
for directory_name in os.listdir(base_path):
    site.addsitedir(os.path.join(base_path, directory_name))
import text
import cPickle
import gzip
class GenBankTable(object):
    def __init__(self,genbank_files=[]):
        self.genbank_files = genbank_files
        self.locusDict = defaultdict(list)
        self.proteinDict = {}
    """ Create dictionary for loci and corresponding proteins"""
    def buildLocusTable(self):
        for genbank_file in self.genbank_files:
            seq_record = SeqIO.parse(open(genbank_file), "genbank").next()
            for feature in seq_record.features:        
                try:
                    if "locus_tag" in feature.qualifiers:
                        locus = feature.qualifiers["locus_tag"][0]
                    elif "gene" in feature.qualifiers:
                        locus = feature.qualifiers["gene"][0]
                    else:
                        continue
                    proteinID = feature.qualifiers["protein_id"][0]
                    self.locusDict[locus].append(proteinID)                    
                except KeyError as k:
                    continue
    """ Create dictionary for proteins and corresponding text"""
    def buildProteinTable(self):
        for genbank_file in self.genbank_files:
            seq_record = SeqIO.parse(open(genbank_file), "genbank").next()
            for feature in seq_record.features:
                try:
                    proteinID = feature.qualifiers["protein_id"][0]
                    note = ""
                    if "note" in feature.qualifiers:
                        note+= text.formatText(feature.qualifiers["note"][0])
                    if "function" in feature.qualifiers:
                        note+= text.formatText(feature.qualifiers["function"][0])
                    if "product" in feature.qualifiers:
                        note+= text.formatText(feature.qualifiers["product"][0])
                    self.proteinDict[proteinID] = note
                except KeyError as k:
                    continue
    def getProteinText(self,protid):
        return self.proteinDict[protid]
    def getLocusText(self,locus_tag):
        note = ""
        protein_ids = self.locusDict[locus_tag]
        for prot in protein_ids:
            note+=self.getProteinText(prot)
        return note
    """ Dump object file into pickle file"""
    def dump(self,outfile):
        cPickle.dump( (self.locusDict,self.proteinDict,self.genbank_files) ,gzip.GzipFile(outfile,'wb'))
    """ Load object from pickle file """
    def load(self,infile):
        self.locusDict,self.proteinDict,self.genbank_files = cPickle.load(gzip.GzipFile(infile,'rb'))
        
def entrezProteinDescription(protid):
    handle = Entrez.efetch(db="nucleotide", 
                           id=protid, rettype="gb", retmode="text")
    description = ""
    seq_record = SeqIO.read(handle, "genbank")
    
    for feature in seq_record.features:
        description+=" "+getDescription(feature)
    return description
    
""" Look for descriptive fields """
def getDescription(feature):
    description = ''
    try:
        description += " "+feature.qualifiers["note"][0]
    except KeyError as k:
        pass
    try:
        description += " "+feature.qualifiers["function"][0]
    except KeyError as k:
        pass
    try:
        description += " "+feature.qualifiers["product"][0]
    except KeyError as k:
        pass
    try:
        protid = feature.qualifiers["protein_id"][0]
        description+=" "+ entrezProteinDescription(protid)
    except KeyError as k:
        pass
    return description

def go(root_dir,output_file):
    genbank_files = []
    outHandle = open(output_file,'w')
    for root, subFolders, files in os.walk(root_dir):
        for fname in files:
            genome_files = []
            organism,ext = os.path.splitext(os.path.basename(fname))
            if ext==".gbk":
                absfile=os.path.join(root,fname)
                genbank_files.append(absfile)
    gbk = GenBankTable(genbank_files)
    gbk.buildLocusTable()
    gbk.buildProteinTable()
    gbk.dump(output_file)
    
if __name__=="__main__":
    
    parser = argparse.ArgumentParser(description=\
        'Finds intergenic regions from genback files')
    parser.add_argument(\
        '--root-dir', type=str, required=False,
        help='Genbank files containing bacterial genomes')
    parser.add_argument(\
        '--output', type=str, required=False,
        help='The output file containing pickle file ')
    parser.add_argument(\
        '--test', action='store_const', const=True, default=False,
        help='Run unittests')
    args = parser.parse_args()
 
    if not args.test:
        go(args.root_dir,args.output)
    else:
        del sys.argv[1:]
        import unittest
        class TestBuildPickle(unittest.TestCase):
            def setUp(self):
                self.root = os.environ['BACFINDER_HOME']
                self.genome_dirs = [ 'Acetobacterium',
                                     'Butyrivibrio_fibrisolvens',
                                     'Catenulispora_acidiphila',
                                     'Saccharopolyspora_erythraea',
                                     'Staphylococcus_aureus',
                                     'Streptococcus_pyogenes',
                                     'Streptomyces_avermitilis']
                self.exampledir   = '%s/example'%self.root
                self.genome_dirs  = ["%s/%s"%(self.exampledir,g) for g in self.genome_dirs]
                print self.genome_dirs
                self.genbank_files = []
                for gdir in self.genome_dirs:
                    for file in os.listdir(gdir):
                        if file.endswith(".gbk"):
                            self.genbank_files.append("%s/%s"%(gdir,file))
                print self.genbank_files
                self.pickle = "test.pickle"
            def tearDown(self):
                os.remove(self.pickle) 
            def test1(self):
                gbk = GenBankTable(self.genbank_files)
                gbk.buildLocusTable()
                gbk.buildProteinTable()
                
                gbk2 = GenBankTable(self.genbank_files)
                gbk.dump(self.pickle)
                self.assertTrue(os.path.getsize(self.pickle) > 0)
                gbk2.load(self.pickle)
                self.assertTrue('BAC70490.1' in gbk2.proteinDict)
        unittest.main()
        
        
        
        
        
        
        
        
        
        
        
        
        