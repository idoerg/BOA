"""
Trys to assign a function to a gene based on its annotated text
"""
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
from Bio import Entrez
import training
import genbank
import cPickle
import gzip
import copy
word_reg = re.compile("[a-z]+")
class NBayes(object):
    def __init__(self,trainDir,labelFile):
        self.classifier = None
        self.labelFile = labelFile
        self.trainingDir = trainDir
        self.labels = None
        self.all_words = None
        #self.labels = training.setup(labelFile)
        #self.train()
        
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
    def gene_features(self,gene_annotations):
        gene_words = set(gene_annotations)
        features = {}
        for word in self.all_words:
            features['contains(%s)'%word] = (word in gene_words)
        return features
    def train(self):
        self.labels = training.Labels(self.trainingDir,self.labelFile)
        trainingText = self.labels.getTrainingText()
        random.shuffle(trainingText)
        text,labs = zip(*trainingText)
        self.all_words = list(itertools.chain(*text))
        #all_words = re.split("\S+"," ".join(map(str,text)))
        self.all_words = nltk.FreqDist(w.lower() for w in self.all_words).keys()[:2000]
        feature_sets = [(self.gene_features(d),c) for (d,c) in trainingText]
        self.classifier = nltk.NaiveBayesClassifier.train(feature_sets)    
        
    """ Make sure that the algorithm works on training data """
    def twofoldcrossvalidation(self):
        self.labels = training.Labels(self.trainingDir,self.labelFile)
        trainingText = self.labels.getTrainingText()
        random.shuffle(trainingText)
        text,labs = zip(*trainingText)
        self.all_words = list(itertools.chain(*text))
        #all_words = re.split("\S+"," ".join(map(str,text)))
        self.all_words = nltk.FreqDist(w.lower() for w in self.all_words).keys()[:2000]
        feature_sets = [(self.gene_features(d),c) for (d,c) in trainingText]
        train_set,test_set = feature_sets[:70],feature_sets[70:]
        self.classifier = nltk.NaiveBayesClassifier.train(train_set)
        p = nltk.classify.accuracy(self.classifier,test_set)
        return p
    def crossvalidation(self):
        self.labels = training.Labels(self.trainingDir,self.labelFile)
        trainingText = self.labels.getTrainingText()
        random.shuffle(trainingText)
        text,labs = zip(*trainingText)
        self.all_words = list(itertools.chain(*text))
        #all_words = re.split("\S+"," ".join(map(str,text)))
        self.all_words = nltk.FreqDist(w.lower() for w in self.all_words).keys()[:2000]
        feature_sets = [(self.gene_features(d),c) for (d,c) in trainingText]
        p = nltk.classify.accuracy(self.classifier,feature_sets)
        return p
    def leaveOneOutCrossValidation(self):
        self.labels = training.Labels(self.trainingDir,self.labelFile)
        trainingText = self.labels.getTrainingText()
        random.shuffle(trainingText)
        text,labs = zip(*trainingText)
        self.all_words = list(itertools.chain(*text))
        #all_words = re.split("\S+"," ".join(map(str,text)))
        self.all_words = nltk.FreqDist(w.lower() for w in self.all_words).keys()[:2000]
        feature_sets = [(self.gene_features(d),c) for (d,c) in trainingText]
        error = 0
        N = len(feature_sets)
        for i in range(N):
            train_set1,test_set,train_set2 = feature_sets[:i],feature_sets[i],feature_sets[i+1:]
            train_set = train_set1+train_set2
            test_set = [test_set]
            self.classifier = nltk.NaiveBayesClassifier.train(train_set)
            p = nltk.classify.accuracy(self.classifier,test_set)
            error+=p
        return error/N
    
        
    """ Classifies proteins based on its text """
    def classify(self,db,fastain):
        proIDs,features,labels = [],[],[]
        #Entrez.email = "mortonjt@miamioh.edu"
        prevFeatureset = ''
        prevText = ''
        for seq_record in SeqIO.parse(fastain, "fasta"):
            title = seq_record.id
            toks = title.split("|")
            proteinID = toks[5]
            #text = genbank.entrezProteinDescription(proteinID)
            query_rows = genbank.proteinQuery(proteinID,db)
            ids,text = zip(*query_rows)
            text = ''.join(map(str,text))
            if text=='': 
                label = 'na'
            else:
                text = word_reg.findall(text)
                featureset = self.gene_features(text)
                assert text!=prevText
                assert featureset!=prevFeatureset
                prevFeatureset = featureset
                prevText = text
                label = self.classifier.classify(featureset)    
                pd = self.classifier.prob_classify(featureset)
            
            proIDs.append(proteinID)  
            labels.append(label)
            features.append(text)
        #print features
        #labels = self.classifier.batch_classify(features)
        return zip(proIDs,labels)
        
    """ Dump object into pickle file """
    def dump(self,outfile):
        print "Dumped pickle file"
        cPickle.dump( (self.classifier,self.all_words,self.labels) ,gzip.GzipFile(outfile,'wb'))
    """ Load object from pickle file """
    def load(self,infile):
        self.classifier,self.all_words,self.labels = cPickle.load(gzip.GzipFile(infile,'rb'))
        
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
                nb = NBayes(self.trainDir,self.test_file)
                nb.train()
                p = nb.crossvalidation()
                print "Accuracy:",p
            def test1(self):
                #Labs = training.setup(self.genbankDir,self.labelFile)
                nb = NBayes(self.trainDir,self.test_file)
                nb.train()
                nb.classifier.show_most_informative_features()
                
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
                nb = NBayes(self.trainDir,self.labelFile)
                nb.train()
                original = copy.deepcopy(nb)
                nb.classifier.show_most_informative_features()
                nb.dump(self.zip)
                nb.load(self.zip)
                self.assertEquals(nb.labelFile,original.labelFile)
                
            def test2(self):
                #Obtain training labels
                #Labs = training.setup(self.genbankDir,self.labelFile)
                nb = NBayes(self.trainDir,self.labelFile)
                nb.train()
                p = nb.crossvalidation()
                print "Accuracy:",p
                self.assertTrue(p>0.5)
                
        unittest.main()
