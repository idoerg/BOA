""" 
Random forests classifier
"""
from sklearn.ensemble import RandomForestClassifier
import nltk
from nltk.classify.scikitlearn import SklearnClassifier
import os,sys
import itertools
import argparse
from collections import defaultdict
import random
import Bio
import re
from Bio import SeqIO, SeqFeature
from Bio.SeqRecord import SeqRecord
from Bio import Entrez
import training
import genbank
import cPickle
import gzip
import copy
import text_classifier

word_reg = re.compile("[a-z]+")
class RForests(text_classifier.TextClassifier):
    def __init__(self,trainDir,labelFile,numTrees=10):
        self.classifier = None
        self.labelFile = labelFile
        self.trainingDir = trainDir
        self.labels = None
        self.all_words = None
        self.numTrees = numTrees
        self.classifier = SklearnClassifier(RandomForestClassifier(
                                            n_estimators=self.numTrees),sparse=False)
        #self.labels = training.setup(labelFile)
        #self.train()
    
    def train(self):
        feature_sets = self.getFeatures()
        self.classifier.train(feature_sets)
        
    """ Determines training error"""
    def trainingError(self):
        feature_sets = self.getFeatures()
        p = nltk.classify.accuracy(self.classifier,feature_sets)
        return p
        
    """ Make sure that the algorithm works on training data using a k fold 
        cross validation scheme """
    def kfoldCrossValidation(self,k):
        feature_sets = self.getFeatures()
        error = 0
        for i in range(k):
            self.classifier = SklearnClassifier(RandomForestClassifier(
                                                n_estimators=self.numTrees),sparse=False)
            n = len(feature_sets)/k
            train_set,test_set = feature_sets[:n*i],feature_sets[n*i:]
            test_set1 = feature_sets[:n*i]
            train_set   = feature_sets[n*i:n*(i+1)]
            test_set2 = feature_sets[i+1:]
            test_set = test_set1+test_set2
            self.classifier.train(feature_sets)
            p = nltk.classify.accuracy(self.classifier,test_set)
        return p
    """ Make sure that the algorithm works on training data using a leave one out 
        cross validation scheme """
    def leave1OutCrossValidation(self):
        error = 0
        feature_sets = self.getFeatures()
        N = len(feature_sets)
        for i in range(N):
            self.classifier = SklearnClassifier(RandomForestClassifier(
                                                n_estimators=self.numTrees),sparse=False)
            train_set1,test_set,train_set2 = feature_sets[:i],feature_sets[i],feature_sets[i+1:]
            train_set = train_set1+train_set2
            test_set = [test_set]
            self.classifier.train(feature_sets)
            p = nltk.classify.accuracy(self.classifier,test_set)
            error+=p
        return error/N
            
    """ Construct a learning curve to see if there is overfitting"""
    def learningCurve(self,numTrials=4):
        accuracies = []
        feature_sets = self.getFeatures()
        for k in xrange(1,len(feature_sets)-1):
            total = 0
            for i in xrange(numTrials):
                self.classifier = SklearnClassifier(RandomForestClassifier(
                                                    n_estimators=self.numTrees),sparse=False)
                random.shuffle(feature_sets)
                train_set,test_set = feature_sets[:k],feature_sets[k:]
                self.classifier.train(feature_sets)
                p = nltk.classify.accuracy(self.classifier,test_set)
                total+=p
            accuracies.append(total/numTrials)
        return accuracies
    
    """ Train on only k features and return training labels and predicted labels """
    def testClassify(self,k):
        feature_sets = self.getFeatures()
        random.shuffle(feature_sets)
        self.classifier = SklearnClassifier(RandomForestClassifier(
                                            n_estimators=self.numTrees),sparse=False)
        
        self.classifier.train(feature_sets[k:])
        features,ref_labels = zip(*feature_sets[:k])
        pred_labels = self.classifier.prob_classify_many(features)   
        return ref_labels,pred_labels
    
    """ nltk confusion matrix """
    def confusionMatrix(self,ref,test):
        ref.sort(key=lambda x: x[0])
        test.sort(key=lambda x: x[0])
        _,ref_labels = zip(*ref)
        _,test_labels = zip(*test)
        cm = ConfusionMatrix(ref_labels, test_labels)
        return cm
    
    """ Classifies proteins based on its text """
    def classify(self,db,fastain):
        proIDs,features,labels = [],[],[]
        prevFeatureset = ''
        prevText = ''
        for seq_record in SeqIO.parse(fastain, "fasta"):
            title = seq_record.id
            toks = title.split("|")
            proteinID = toks[5]
            query_rows = genbank.proteinQuery(proteinID,db)
            ids,text = zip(*query_rows)
            text = ''.join(map(str,text))
            if text=='': 
                label = ['na']
            else:
                text = word_reg.findall(text)
                featureset = self.gene_features(text)
                assert text!=prevText
                assert featureset!=prevFeatureset
                prevFeatureset = featureset
                prevText = text
                label = self.classifier.classify_many([featureset])    
            
            proIDs.append(proteinID)  
            labels+=label
        return zip(proIDs,labels)


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
                self.root = os.environ['BACFINDER_HOME']
                self.genbankDir = "%s/example/Streptococcus_pyogenes"%self.root
                self.trainDir = "%s/data/training/protein"%self.root
                self.genbankFile = "%s/example/Streptococcus_pyogenes/NC_011375.gbk"%self.root
                self.test_file = "test_labels.txt"
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
            def testText(self):
                nb = RForests(self.trainDir,self.test_file)
                nb.train()
                p = nb.trainingError()
                print "Accuracy:",p
            def test1(self):
                #Labs = training.setup(self.genbankDir,self.labelFile)
                nb = RForests(self.trainDir,self.test_file)
                nb.train()
                #nb.classifier.show_most_informative_features()
                
        class TestTraining2(unittest.TestCase):
            def setUp(self):
                #self.genbankDir = "../example/Streptococcus_pyogenes"
                #self.genbankFile = "../example/Streptococcus_pyogenes/NC_011375.gbk"
                self.root = os.environ['BACFINDER_HOME']
                self.trainDir = "%s/data/training/protein"%self.root
                self.labelFile = "%s/data/training/training_proteins.txt"%self.root
                self.zip = "test_serial.zip"                
                #Obtain training labels
            def test1(self):
                #Labs = training.setup(self.genbankDir,self.labelFile)
                nb = RForests(self.trainDir,self.labelFile)
                nb.train()
                original = copy.deepcopy(nb)
                #nb.classifier.show_most_informative_features()
                nb.dump(self.zip)
                nb.load(self.zip)
                self.assertEquals(nb.labelFile,original.labelFile)
                
            def test2(self):
                #Obtain training labels
                #Labs = training.setup(self.genbankDir,self.labelFile)
                nb = RForests(self.trainDir,self.labelFile)
                nb.train()
                p = nb.trainingError()
                print "Accuracy:",p
                self.assertTrue(p>0.5)
        class TestClassify(unittest.TestCase):
             def setUp(self):
                #self.genbankDir = "../example/Streptococcus_pyogenes"
                #self.genbankFile = "../example/Streptococcus_pyogenes/NC_011375.gbk"
                self.root = os.environ['BACFINDER_HOME']
                self.trainDir = "%s/data/training/protein"%self.root
                self.labelFile = "%s/data/training/training_proteins.txt"%self.root
                self.zip = "test_serial.zip"                
                #Obtain training labels
             def test1(self):
                #Labs = training.setup(self.genbankDir,self.labelFile)
                nb = RForests(self.trainDir,self.labelFile)
                ref,pred = nb.testClassify(30)
                self.assertTrue(len(ref)>0)
                self.assertTrue(len(pred)>0)
                
        unittest.main()

