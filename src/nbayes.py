import os,sys
import nltk
import argparse
from collections import defaultdict
import random
import Bio
import re
from Bio import SeqIO, SeqFeature
from Bio.SeqRecord import SeqRecord
import training

class NBayes(object):
    def __init__(self,tableFile):
        self.tableFile = tableFile
        self.loci = dict() #classifcation of loci

        #self.uploadTable()
        #self.classifier = None
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
    
    """Version 1: Only using manually curated genbank files"""
    """
    def train(self,genbank_files):
        trainingText = []
        all_words = []
        for genbank in genbank_files:
            seq_record = SeqIO.parse(open(genbank), "genbank").next()
            for feature in seq_record.features:
                try:
                    locus = feature.qualifiers["locus_tag"][0]
                    note = feature.qualifiers["note"][0]
                    words = [x for x in note.split(' ') if x!='']
                    if locus in self.loci:
                        category = self.loci[locus]
                    else:
                        category = "Null"
                    all_words+=words
                    trainingText.append( (words,category) )
                except KeyError as k:
                    continue
        random.shuffle(trainingText)
        all_words = nltk.FreqDist(w.lower() for w in all_words).keys()[:2000]
        feature_sets = [(self.gene_features(d,all_words),c) for (d,c) in trainingText]
        self.classifier = nltk.NaiveBayesClassifier.train(feature_sets)
    """
    """ Version 2: Using manually curated label objects """
    def train(self,labels):
        trainingText = labels.getTrainingText()
        random.shuffle(trainingText)
        labs,text = zip(*trainingText)
        all_words = re.split("\S+"," ".join(map(str,text)))
        all_words = nltk.FreqDist(w.lower() for w in all_words).keys()[:2000]
        feature_sets = [(self.gene_features(d,all_words),c) for (d,c) in trainingText]
        self.classifier = nltk.NaiveBayesClassifier.train(feature_sets)    
    
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
        class TestTraining(unittest.TestCase):
            def setUp(self):
                self.genbankDir = "../example/Streptococcus_pyogenes"
                self.genbankFile = "../example/Streptococcus_pyogenes/NC_011375.gbk"
                self.labelFile = "../data/training/training.txt"                
            def test1(self):
                #Obtain training labels
                Labs = training.setup(self.genbankDir,self.labelFile)
                nb = NBayes(self.labelFile)
                nb.train(Labs)
                nb.classifier.show_most_informative_features(n=20)
        unittest.main()
