"""
Obtains training dataset

This assumes that all of the training genbank files 
are organized in the same directory (trainingDir)
"""

from collections import defaultdict
import re
import argparse
import os, site, sys
import itertools
from Bio import SeqIO
from Bio import Entrez
import random


base_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
for directory_name in os.listdir(base_path):
    site.addsitedir(os.path.join(base_path, directory_name))
import genbank
import text

""" Stores labels as keys and descriptions as values """
accession_reg = re.compile("([A-Z]+_?\d+)")
gi_reg = re.compile("GI:\s+(\d+)")
word_reg = re.compile("[a-z]+")
class Labels(object):
    def __init__(self,trainingDir,labelfile):
        self.labels = defaultdict( dict )
        self.seqs = []
        self.parseLabels(trainingDir,labelfile)
    def getKeys(self):#returns filename basically
        return self.labels.keys()
    def numSpecies(self):
        return len(self.labels)
    def filterWords(self,words):
        #Remove all single letter words
        return [text.formatText(w) for w in words if len(w)>1]
    def decomposeDocument(self,text):
        words = [self.filterWords(word_reg.findall(t.lower()))
                 for t in text]
        return words
    """ Get training sequences sorted by label"""
    def getTrainingSequences(self):
        sorted(self.seqs,key=lambda x:x[2])
        return self.seqs
    
    """ Get descriptions associated with training sequences"""
    def getTrainingText(self):
        L =  self.labels.values()
        K = [l.items() for l in L]
        trainset = list(itertools.chain(*K))
        #print trainset
        loci,labels = zip(*trainset)
        labs,text = zip(*labels)
        #print labs,list(text)
        #Decompose text into words
        words = self.decomposeDocument(list(text))
        return zip(words,labs) 
    def __str__(self):
        s = ""
        for k,v in self.labels.iteritems():
            s+= "Organism:%s\n%s\n"%(k,str(v))
        return s
    """ Parses training text file """
    def parseLabels(self,directory,labelfile):
        with open(labelfile,'r') as handle:            
            for ln in handle:
                ln = ln.rstrip()
                if ln[0]=='#':
                    #if "GI:" in ln:
                    if "Organism:" in ln:
                        key = accession_reg.findall(ln)[0]
                        if key not in self.labels:
                            self.labels[key] = {}
                        else:
                            print "Already in labels!!!"                        
                else:
                    toks = re.split("\s+",ln)
                    proteinID,label = toks[0],toks[1]
                    #print proteinID,label
                    #Each label has a dictionary
                    #proteinID:(label, description)       
                    note = self.getProteinGenbankDescription(directory,proteinID)  
                    seq  = self.getProteinGenbankSequence(directory,proteinID)  
                    self.labels[key][proteinID]=(label,note)
                    self.seqs.append( (key,proteinID,label,seq) )
                    #print key,self.labels[key]

    """ Parses a genbank file for a particular protein and returns its description"""
    def getProteinGenbankDescription(self,directory,proteinID):
        fname = "%s/%s.gbk"%(directory,proteinID)
        seq_record = SeqIO.parse(open(fname), "genbank").next()
        note = ""
        for feature in seq_record.features: 
            try:
                if "note" in feature.qualifiers:
                    note+= text.formatText(feature.qualifiers["note"][0])
                if "function" in feature.qualifiers:
                    note+= text.formatText(feature.qualifiers["function"][0])
                if "product" in feature.qualifiers:
                    note+= text.formatText(feature.qualifiers["product"][0])    
            except KeyError as k:
                continue
        return note 
    """ Parses a genbank file for a particular protein and returns its description"""
    def getProteinGenbankSequence(self,directory,proteinID):
        fname = "%s/%s.gbk"%(directory,proteinID)
        seq_record = SeqIO.parse(open(fname), "genbank").next()
        return str(seq_record.seq)
            
        
    """Parses genbank record"""
    def parseRecord(self,seq_record):
        acc = seq_record.annotations['accessions'][0]
        #print acc,self.labels[acc].keys()
        for feature in seq_record.features:
            try:                    
                locus = None                   
                if "locus_tag" in feature.qualifiers:
                    locus = feature.qualifiers["locus_tag"][0]
                    
                elif "gene" in feature.qualifiers:
                    locus = feature.qualifiers["gene"][0]
                else:
                    continue
                if locus !=None:
                    if locus in self.labels[acc]:
                        label,description = self.labels[acc][locus]
                        description = genbank.getDescription(feature)
                        self.labels[acc][locus]=(label,description)
                    elif (random.random()<0.0007): #sample 7/1000 for null entries
                        label,description = "null",''
                        print "Null label"
                        description = genbank.getDescription(feature)
                        self.labels[acc][locus]=(label,description)
                    #print locus,label,description
                
            except KeyError as k:
                print "Exception",k                    
                continue
    
    "Parses genbank file"
    """
    def parseGenbank(self,handle):
        try:
            seq_record = SeqIO.read(handle, "genbank")
            self.parseRecord(seq_record)  
        except Exception as e:
            print "Exception at",e
    """
""" Create a labels object """
"""
def setup(rootdir,labelFile):
    print rootdir
    labs = Labels(labelFile)
    #outHandle = open(out,'w')
    for root, subFolders, files in os.walk(rootdir):
        for fname in files:
            genome_files = []
            organism,ext = os.path.splitext(os.path.basename(fname))
            if ext==".gbk":
                absfile=os.path.join(root,fname)
                print absfile
                labs.parseGenbank(absfile)
    return labs
"""
def setup(labelFile):
    labs = Labels(labelFile)
    Entrez.email = "mortonjt@miamioh.edu"
    for key in labs.getKeys():
        print key
        handle = Entrez.efetch(db="nucleotide", 
                               id=key, rettype="gbwithparts", retmode="text")
        #if "NC_011375" in key: 
        #    print(handle.read())
        #print(handle.read())
        labs.parseGenbank(handle)
    
    return labs



if __name__=="__main__":
    parser = argparse.ArgumentParser(description=\
        'Builds the training dataset')
    parser.add_argument(\
        '--training-labels', type=str, required=False,
        help='Text file containing training labels')    
    parser.add_argument(\
        '--root-dir', type=str,required=False,default="",
        help='Root directory of all of the files of interest')
    parser.add_argument(\
        '--test', action='store_const', const=True, default=False,
        help='Run the unittests')
    args = parser.parse_args()
    if not args.test:
        labels = args.training_labels
        #root,labels,out = args.root_dir,args.training_labels,args.output_file
        setup(labels)
    else:
        del sys.argv[1:]
        import unittest
        class TestParser(unittest.TestCase):
            def setUp(self):
                self.root = os.environ['BACFINDER_HOME']
                self.genbankDir = "%s/example/Streptococcus_pyogenes"%self.root
                self.genbankFile = "%s/example/Streptococcus_pyogenes/NC_011375.gbk"%self.root
                self.test_file = "test_labels.txt"
                self.trainDir = "%s/data/training/protein"%self.root
                string = "\n".join(["#Organism: Y12234.1 (as-48A-D1) and AJ438950.1 (as-48E - H), Enterococcus faecalis subsp. liquefaciens plasmid submitted as separate sequences)",
                            "#Reference: http://jb.asm.org/content/190/1/240.full, http://aem.asm.org/content/69/2/1229.full.pdf",
                            "#locus_tag label name",
                            "CAA72917.1 toxin",
                            "CAA72918.1 modifier",
                            "CAA72919.1 transport",
                            '#Organism: AF061787.1, Escherichia coli plasmid pTUC100',
                            '#Reference: http://www.ncbi.nlm.nih.gov/pmc/articles/PMC93700',
                            '#locus_tag label name',
                            'AAD28494.1 toxin',
                            'AAD28495.1 modifier',
                            'AAD28496.1 modifier',
                            'AAD28497.1 transport'])
                open(self.test_file,'w').write("%s\n"%string)
            def tearDown(self):
                os.remove(self.test_file)
            def testParse(self):
                labs = Labels(self.trainDir,self.test_file)
                self.assertEquals(labs.numSpecies(),2)
                trainingText = labs.getTrainingText()
                print '\n'.join(map(str,trainingText))
        unittest.main()
