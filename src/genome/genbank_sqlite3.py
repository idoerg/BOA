
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

base_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
for directory_name in os.listdir(base_path):
    site.addsitedir(os.path.join(base_path, directory_name))
import text

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

""" Creates database whose primary key is the protein ID 
    that maps to annotation text"""
def buildProteinTable(dbout):
    db = sqlite3.connect(dbout)
    cursor = db.cursor()
    cursor.execute('''CREATE TABLE protein_text(protein_id TEXT , 
                                                note TEXT)''')
    db.commit()
""" Inserts all protein entries in genbank files into database"""
def insertGenbankProteins(genbank_files,dbout):
    db = sqlite3.connect(dbout)
    cursor = db.cursor()
    
    for genbank_file in genbank_files:
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
                 
                cursor.execute('''INSERT INTO protein_text(protein_id,note)
                                    VALUES(?,?)''',(proteinID,note))
            except KeyError as k:
                continue
        db.commit()
    db.close()
""" Inserts a single query """
def insertProteinQuery(proteinID,note,dbfile):
    db = sqlite3.connect(dbfile)
    cursor = db.cursor()
    cursor.execute('''INSERT INTO protein_text(protein_id,note)
                        VALUES(?,?)''',(proteinID,note))
    db.commit()
    db.close()
""" Returns rows with (proteinID, note) """
def proteinQuery(protID,dbfile):
    db = sqlite3.connect(dbfile)
    cursor = db.cursor()
    cursor.execute('''SELECT * FROM protein_text WHERE protein_id=?''',(protID,))
    #cursor.execute('''SELECT * FROM protein_text''')
    rows = cursor.fetchall()
    db.close()
    return rows

""" Returns all entries in table """
def proteinQueryAll(dbfile):
    db = sqlite3.connect(dbfile)
    cursor = db.cursor()
    cursor.execute('''SELECT * FROM protein_text''')
    #cursor.execute('''SELECT * FROM protein_text''')
    rows = cursor.fetchall()
    db.close()
    return rows
 
""" Creates database whose primary key is the locus tag
    that maps to protein IDs"""
def buildLocusTable(dbout):
    db = sqlite3.connect(dbout)
    cursor = db.cursor()
    cursor.execute('''CREATE TABLE loci(locus_tag TEXT,
                                       protein_id TEXT)''')
    db.commit()
    db.close()
    
""" Inserts all entries in genbank files into database"""
def insertGenbankLoci(genbank_files,dbout):
    db = sqlite3.connect(dbout)
    cursor = db.cursor()
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
                proteinID = feature.qualifiers["protein_id"][0]
                cursor.execute('''INSERT INTO loci(locus_tag,protein_id)
                              VALUES(?,?)''',(locus,proteinID))
            except KeyError as k:
                continue
            
        db.commit()
    db.close()

""" Inserts a single query """
def insertLocusQuery(locus,proteinID,dbfile):
    db = sqlite3.connect(dbfile)
    cursor = db.cursor()
    cursor.execute('''INSERT INTO loci(locus_tag,protein_id)
                  VALUES(?,?)''',(locus,proteinID))
    db.commit()
    db.close()

""" Returns all entries in table"""    
def locusQueryAll(dbfile):
    db = sqlite3.connect(dbfile)
    cursor = db.cursor()
    cursor.execute('''SELECT * FROM loci''')
    #cursor.execute('''SELECT * FROM loci''')
    rows = cursor.fetchall()
    db.close()
    return rows
    
""" Returns rows of (locus_tag,protein_id)"""    
def locusQuery(locus,dbfile):
    db = sqlite3.connect(dbfile)
    cursor = db.cursor()
    cursor.execute('''SELECT * FROM loci WHERE locus_tag=?''',(locus,))
    #cursor.execute('''SELECT * FROM loci''')
    rows = cursor.fetchall()
    db.close()
    return rows

def buildTextDB(genbank_files,db):
    insertGenbankLoci(genbank_files,db)
    insertGenbankProteins(genbank_files,dbout)
    
if __name__=="__main__":
    
    parser = argparse.ArgumentParser(description=\
        'Finds intergenic regions from genback files')
    parser.add_argument(\
        '--genome-files', type=str, nargs="+", required=False,
        help='Genbank files containing bacterial genomes')
    parser.add_argument(\
        '--output', type=str, required=False,
        help='The output file containing either pickle file or sqlite3 database')
    parser.add_argument(\
        '--test', action='store_const', const=True, default=False,
        help='Run unittests')
    args = parser.parse_args()
 
    if not args.test:
        buildTextDB(args.genome_files,args.output)
    

    else:
        del sys.argv[1:]
        import unittest
        class TestBuildDB(unittest.TestCase):
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
                self.genome_files = []
                for gdir in self.genome_dirs:
                    for file in os.listdir(gdir):
                        if file.endswith(".gbk"):
                            self.genome_files.append("%s/%s"%(gdir,file))
                self.db = "testdb"
                if os.path.exists(self.db):
                    os.remove(self.db)
                
            def tearDown(self):
                #os.remove(self.locusdb)
                #os.remove(self.protdb)
                os.remove(self.db)
                pass
            def testdb1(self):
                buildLocusTable(self.db)
                insertGenbankLoci(self.genome_files,self.db)
                self.assertTrue(os.path.getsize(self.db) > 0)
                entries = locusQuery("Spy49_0568",self.db)
                self.assertTrue(len(entries) > 0)
                    
            def testdb2(self):
                buildProteinTable(self.db)
                insertGenbankProteins(self.genome_files,self.db)
                self.assertTrue(os.path.getsize(self.db) > 0)
                    
                entries = proteinQuery("NP_828769.1",self.db)
                #entries = proteinQuery("ACI60894.1",self.db)
                self.assertTrue(len(entries) > 0)
             
        unittest.main()
        
        
        
        
