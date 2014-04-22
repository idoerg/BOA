"""
Obtains training dataset
"""

from collections import defaultdict
import re
import argparse
import os, site, sys
import itertools
from Bio import SeqIO

""" Stores labels as keys and descriptions as values """
accession_reg = re.compile("([A-Z]+_?\d+.\d)")
class Labels(object):
    def __init__(self,labelfile):
        self.labels = defaultdict( dict )
        self.parseLabels(labelfile)
    def numSpecies(self):
        return len(self.labels)
    def getTrainingText(self):
        L =  self.labels.values()
        print L
        K = [l.items() for l in L]
        print K
        return list(itertools.chain(*K))
    def parseLabels(self,labelfile):
        with open(labelfile,'r') as handle:            
            headerCnt = 0
            for ln in handle:
                ln = ln.rstrip()
                if ln[0]=='#' and headerCnt==0:
                    acc = accession_reg.findall(ln)
                    key = accession_reg.findall(ln)[0]
                    if key not in self.labels:
                        self.labels[key] = {}
                    else:
                        print "Already in labels!!!"
                    headerCnt+=1
                elif ln[0]=="#": continue
                else:
                    headerCnt=0
                    toks = re.split("\S+",ln)
                    locus,label = toks[0],toks[1]
                    #Each label has a dictionary
                    #locus:(label, description)                    
                    self.labels[key]= {locus:(label,"")}    
    """Parses genbank record"""
    def parseRecord(self,seq_record):
            acc = seq_record.annotations['accessions'][0]
            for feature in seq_record.features:
                try:                         
                    locus = feature.qualifiers["locus_tag"][0]
                    note = feature.qualifiers["note"][0]
                    if locus=="Spy49_0568": print locus
                    if locus in self.labels[acc]:
                        label,description = self.labels[acc][locus]
                        description += note
                        self.labels[acc]={[locus]:(label,description)}
                    else:
                        description = note
                        self.labels[acc]={[locus]:(label,description)}
                    print locus,label,description
                except KeyError as k:
                    #print "Exception",k                    
                    continue
    "Parses genbank file"
    def parseGenbank(self,genbank_file):
        try:
            seq_record = SeqIO.read(open(genbank_file), "genbank")
            self.parseRecord(seq_record)  
        except Exception as e:
            print "Exception at",genbank_file,e

""" Create a labels object """
def setup(rootdir,labels):
    print rootdir
    labs = Labels(labels)
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
        root,labels,out = args.root_dir,args.training_labels,args.output_file
        setup(root,labels,out)
    else:
        del sys.argv[1:]
        import unittest
        class TestParser(unittest.TestCase):
            def setUp(self):
                self.genbankDir = "../example/Streptococcus_pyogenes"
                self.genbankFile = "../example/Streptococcus_pyogenes/NC_011375.gbk"
                self.test_file = "test_labels.txt"
                string = "\n".join(["#Organism: Y12234.1 (as-48A-D1) and AJ438950.1 (as-48E - H), Enterococcus faecalis subsp. liquefaciens plasmid submitted as separate sequences)",
                                    "#Reference: http://jb.asm.org/content/190/1/240.full, http://aem.asm.org/content/69/2/1229.full.pdf",
                                    "#locus_tag label name",
                                    "as-48 toxin",
                                    "as-48B modifier",
                                    "as-48C transport",
                                    '#Organism: AF061787.1, Escherichia coli plasmid pTUC100',
                                    '#Reference: http://www.ncbi.nlm.nih.gov/pmc/articles/PMC93700',
                                    '#locus_tag label name',
                                    'mcjA toxin',
                                    'mcjB modifier',
                                    'mcjC modifier',
                                    'mcjD transport'])
                open(self.test_file,'w').write("%s\n"%string)
            def tearDown(self):
                os.remove(self.test_file)
            def testParse(self):
                labs = Labels(self.test_file)
                self.assertEquals(labs.numSpecies(),2)
            def testText(self):
                labs = setup(self.genbankDir,self.test_file)
                trainingText = labs.getTrainingText()
                print trainingText
        unittest.main()
