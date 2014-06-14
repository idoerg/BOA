"""
Wrapper class for genbank files
"""

import Bio
from Bio import SeqIO, SeqFeature
from Bio.SeqRecord import SeqRecord
from Bio import Entrez
from nltk import PorterStemmer
import sys
import os
import argparse
import glob
import sqlite3

class GenBank(object):
    def __init__(self,fname,dicttype="gene"):
        self.fname = fname
        """
        if dicttype=="protein":
            self.dict = self.buildProteinDictionary(fname)
        else:
            self.dict = self.buildGeneDictionary(fname)
        """

    def buildProteinDictionary(self,fname):
        proteins = dict()
        print fname
        seq_record = SeqIO.parse(open(fname), "genbank").next()
        for feature in seq_record.features:
            if feature.type == 'CDS':
                try:
                    protein = feature.qualifiers["protein_id"][0]
                    proteins[ protein ]=feature
                except KeyError:
                    continue
        return proteins

    def buildGeneDictionary(self,fname):
        genes = dict()
        seq_record = SeqIO.parse(open(fname), "genbank").next()
        for feature in seq_record.features:
            if feature.type == 'CDS':
                try:
                    gene = feature.qualifiers["gene"][0]
                    genes[ gene ]=feature
                except KeyError:
                    continue
        return genes

    def findGene(self,gname):
        try:
            return self.dict[gname]
        except KeyError:
            print >> sys.stderr,"No such gene"

    def findProtein(self,pname):
        try:
            return self.dict[pname]
        except KeyError:
            print >> sys.stderr,"No such protein"


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

"""Format the text to aid naive bayes"""
def formatText(text):
    text = text.lower()
    text = text.replace('.',' ')
    text = text.replace('\\',' ')
    text = text.replace('/',' ')
    text = text.replace('\"',' ')
    text = text.replace('\'',' ')
    text = text.replace(':',' ')
    text = text.replace(';',' ')
    text = text.replace('(',' ')
    text = text.replace(')',' ')
    porter = PorterStemmer()
    return ' '.join([porter.stem(word) for word in text.split(' ')])
""" Loads database into memory"""
def loadTextDB(dbin):
    proteinDict = {}
    with open(dbin,'r') as handle:
        for ln in handle:
            ln = ln.rstrip()
            toks = ln.split('\t')
            proteinID = toks[0]
            text = ' '.join(toks[1:])
            proteinDict[proteinID] = text
    return proteinDict

""" Creates database whose primary key is the protein ID """
def buildProteinTextDB(genbank_files,dbout):
    self.db = splite3.connect(dbout)
    cursor = self.db.cursor()
    cursor.execute('''
            CREATE TABLE protein_text(protein_id TEXT , note TEXT)
            ''')
    db.commit()
    for genbank_file in genbank_files:
        seq_record = SeqIO.parse(open(genbank_file), "genbank").next()
        for feature in seq_record.features:
            try:
                proteinID = feature.qualifiers["proteinID"][0]
                note = ""
                if "note" in feature.qualifiers:
                    note+= formatText(feature.qualifiers["note"][0])
                if "function" in feature.qualifiers:
                    note+= formatText(feature.qualifiers["function"][0])
                if "product" in feature.qualifiers:
                    note+= formatText(feature.qualifiers["function"][0])
                 
                cursor.execute('''INSERT INTO protein_text(protein_id,note)
                                    VALUES(?,?)''',(proteinID,note))
            except KeyError as k:
                continue
        db.commit()
        
    """
    proteinDict = {}        
    for genbank_file in genbank_files:
        seq_record = SeqIO.parse(open(genbank_file), "genbank").next()
        for feature in seq_record.features:
            try:
                proteinID = feature.qualifiers["proteinID"][0]
                note = formatText(feature.qualifiers["note"][0])
            except KeyError as k:
                continue
                
            if proteinID not in proteinDict:
                proteinDict[proteinID] = note
            else:
                proteinDict[proteinID]= "%s %s"%(proteinDict[proteinID],note)
    """
    
""" Creates database whose primary key is the locus tag"""
def buildGeneTextDB(genbank_files,dbout):
    db = splite3.connect(dbout)
    cursor = db.cursor()
    cursor.execute('''
        CREATE TABLE coding_regions(locus_tag TEXT,
                                    protein_id TEXT)
    '''
    )
    self.db.commit()
    
    for genbank_file in genbank_files:
        seq_record = SeqIO.parse(open(genbank_file), "genbank").next()
        for feature in seq_record.features:        
            
            try:
                if "locus_tag" in feature.qualifiers:
                    locus = feature.qualifiers["locus_tag"][0]
                elif "gene" in feature.qualifiers:
                    locus = features.qualifiers["gene"][0]
                else:
                    continue
                proteinID,sequence = "",""
                if "protein_id" in feature.qualifiers:
                    proteinID = feature.qualifiers["protein_id"][0]
                    
                cursor.execute('''INSERT INTO coding_regions(locus_tag,protein_id)
                              VALUES(?,?,?,?,?,?)''',(locus,seq,st,end,note,protid))
                
            except KeyError as k:
                continue
            
            """
            if proteinID not in proteinDict:
                proteinDict[proteinID] = note
            else:
                proteinDict[proteinID]= "%s %s"%(proteinDict[proteinID],note)
            """
                
    proteinIDs,textBlobs = proteinDict.keys(),proteinDict.values()
    handle = open(dbout,'w')
    textFeatures = zip(proteinIDs,textBlobs)
    
    handle.write( '\n'.join(["%s\t%s"%f for f in textFeatures]) )
    handle.close()
    
if __name__=="__main__":
    
    parser = argparse.ArgumentParser(description=\
        'Finds intergenic regions from genback files')
    parser.add_argument(\
        '--genome-files', type=str, nargs="+", required=False,
        help='FASTA files containing bacterial genomes')
    parser.add_argument(\
        '--output-db', type=str, required=False,
        help='The output file containing the tab-delimited output')
    parser.add_argument(\
        '--test', action='store_const', const=True, default=False,
        help='Run unittests')
    args = parser.parse_args()
    if not args.test:

        buildTextDB(args.genome_files,args.output_db)

    else:
        del sys.argv[1:]
        import unittest
        class TestBuildDB(unittest.TestCase):
            def setUp(self):
                self.genome_dirs = [ 'Acetobacterium',
                                     'Butyrivibrio_fibrisolvens',
                                     'Catenulispora_acidiphila',
                                     'Saccharopolyspora_erythraea',
                                     'Staphylococcus_aureus',
                                     'Streptococcus_pyogenes',
                                     'Streptomyces_avermitilis']
                self.exampledir   = '/media/HD/Documents/Jamie/MiamiBio/Bacfinder/example'
                self.genome_dirs  = ["%s/%s"%(self.exampledir,g) for g in self.genome_dirs]
                print self.genome_dirs
                self.genome_files = []
                for gdir in self.genome_dirs:
                    for file in os.listdir(gdir):
                        if file.endswith(".gbk"):
                            self.genome_files.append("%s/%s"%(gdir,file))
                self.outdb = "testdb"
            def testdb1(self):
                buildTextDB(self.genome_files,self.outdb)
                self.assertTrue(os.path.getsize(self.outdb) > 0)
        unittest.main()
        
        
        
        
        
        
        
        
        
        