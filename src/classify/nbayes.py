import os,sys
import nltk
import itertools
import argparse
from collections import defaultdict
import random
import Bio
import re
from Bio import SeqIO, SeqFeature
from Bio.SeqRecord import SeqRecord
import training

class NBayes(object):
    def __init__(self,labelFile):
        self.classifier = None
        self.labelFile = labelFile
        self.labels = training.setup(labelFile)
        self.train()
        #self.loci = dict() #classifcation of loci
        #self.uploadTable()
        
    """
    def uploadTable(self):
        with open(self.tableFile,'r') as tableIn:
            for ln in tableIn:
                ln = ln.rstrip()
                if ln[0]=="#": continue
                toks = ln.split('\t')
                assert len(toks)==2
                locus,category = toks
                self.loci[locus] = category
    """
    def gene_features(self,gene_annotations,all_words):
        gene_words = set(gene_annotations)
        features = {}
        for word in all_words:
            features['contains(%s)'%word] = (word in gene_words)
        return features
    def train(self):
        trainingText = self.labels.getTrainingText()
        random.shuffle(trainingText)
        text,labs = zip(*trainingText)
        all_words = list(itertools.chain(*text))
        #all_words = re.split("\S+"," ".join(map(str,text)))
        all_words = nltk.FreqDist(w.lower() for w in all_words).keys()[:2000]
        feature_sets = [(self.gene_features(d,all_words),c) for (d,c) in trainingText]
        self.classifier = nltk.NaiveBayesClassifier.train(feature_sets)    
        
    """ Make sure that the algorithm works on training data """
    def crossvalidation(self):
        trainingText = self.labels.getTrainingText()
        random.shuffle(trainingText)
        text,labs = zip(*trainingText)
        all_words = list(itertools.chain(*text))
        #all_words = re.split("\S+"," ".join(map(str,text)))
        all_words = nltk.FreqDist(w.lower() for w in all_words).keys()[:2000]
        feature_sets = [(self.gene_features(d,all_words),c) for (d,c) in trainingText]
        p = nltk.classify.accuracy(self.classifier,feature_sets)
        return p
    
    """Classifies proteins based on its text"""
    def classify(self):
        pass
def go():
    pass

if __name__=="__main__":
    parser = argparse.ArgumentParser(description=\
        'A naive bayes classifier that attempts to categorize context genes')
    parser.add_argument(\
        '--training-labels', type=str, required=False,
        help='A training data set to serve as a template for categorizing context genes')
    parser.add_argument(\
        '--genbank-files', type=str, nargs="+", required=False,
        help='Genbank files containing annotations of bacterial genes')
    parser.add_argument(\
        '--test', action='store_const', const=True, default=False,
        help='Run unittests')
    args = parser.parse_args()
    
    if not args.test:
        go()
    else:
        del sys.argv[1:]
        import unittest
        
        class TestTraining1(unittest.TestCase):
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
            def testText(self):
                nb = NBayes(self.test_file)
                p = nb.crossvalidation()
                print "Accuracy:",p
            def test1(self):
                #Labs = training.setup(self.genbankDir,self.labelFile)
                nb = NBayes(self.test_file)
                nb.classifier.show_most_informative_features()
                
        class TestTraining2(unittest.TestCase):
            def setUp(self):
                #self.genbankDir = "../example/Streptococcus_pyogenes"
                #self.genbankFile = "../example/Streptococcus_pyogenes/NC_011375.gbk"
                self.labelFile = "../data/training/training.txt"                
                #Obtain training labels
            def test1(self):
                #Labs = training.setup(self.genbankDir,self.labelFile)
                nb = NBayes(self.labelFile)
                nb.classifier.show_most_informative_features()
            def test2(self):
                #Obtain training labels
                #Labs = training.setup(self.genbankDir,self.labelFile)
                nb = NBayes(self.labelFile)
                p = nb.crossvalidation()
                print "Accuracy:",p
                self.assertTrue(p>0.5)
                
        unittest.main()
