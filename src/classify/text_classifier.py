"""
This is the abstract class for all text classifiers
"""
import abc
import training
import genbank
import cPickle
import gzip
import copy
import random
import itertools
import re
from collections import defaultdict

class TextClassifier(object):
    __metaclass__ = abc.ABCMeta
    def __init__(selfself,trainDir,labelFile):
        self.classifier = None
        self.labelFile = labelFile
        self.trainingDir = trainDir
        self.labels = None
        self.all_words = None
    def gene_features(self,gene_annotations):
        gene_words = set(gene_annotations)
        features = {}
        for word in self.all_words:
            #features['contains(%s)'%word] = (word in gene_words)
            features[word] = (word in gene_words)
        return features
    def getFeatures(self):
        self.labels = training.Labels(self.trainingDir,self.labelFile)
        trainingText = self.labels.getTrainingText()
        random.shuffle(trainingText)
        text,labs = zip(*trainingText)
        self.all_words = list(set(itertools.chain(*text)))
        feature_sets = [(self.gene_features(d),c) for (d,c) in trainingText]
        return feature_sets
    @abc.abstractmethod 
    def train(self): pass
    @abc.abstractmethod 
    def confusionMatrix(self,ref,test): pass
    @abc.abstractmethod
    def trainingError(self): pass
    @abc.abstractmethod
    def leave1OutCrossValidation(self,k): pass
    @abc.abstractmethod
    def kfoldCrossValidation(self,k): pass
    @abc.abstractmethod
    def testClassify(self,k): pass
    @abc.abstractmethod
    def learningCurve(self): pass
    
    """ Reads classifier output"""
    def readOutput(self,fname):
        protids,functions = [],[]
        with open(fname,'r') as handle:
            for ln in handle:
                ln = ln.rstrip()
                toks = ln.split('\t')
                protid,function = toks[0],toks[1]
                protids.append(protid)
                functions.append(function)
        return zip(protids,functions)
    """ Dump object file into pickle file"""
    def dump(self,outfile):
        print "Dumped pickle file"
        cPickle.dump( (self.classifier,self.all_words,self.labels) ,gzip.GzipFile(outfile,'wb'))
    """ Load object from pickle file """
    def load(self,infile):
        self.classifier,self.all_words,self.labels = cPickle.load(gzip.GzipFile(infile,'rb'))
        
