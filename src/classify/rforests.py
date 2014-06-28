""" 
Random forests classifier
Note: Need to fix unittest.  Need to check the version of sklearn
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
import genbank_sqlite3
import cPickle
import gzip
import copy
import text_classifier

word_reg = re.compile("[a-z]+")
class RForests(text_classifier.TextClassifier):
    def __init__(self,trainDir,labelFile,numTrees=10,numJobs=1):
        self.classifier = None
        self.labelFile = labelFile
        self.trainingDir = trainDir
        self.labels = None
        self.all_words = None
        self.numTrees = numTrees
        self.numJobs = numJobs
        self.classifier = SklearnClassifier(RandomForestClassifier(
                                            n_estimators=self.numTrees,
                                            n_jobs=numJobs),sparse=False)
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
                                                    n_estimators=self.numTrees),
                                                    sparse=False)
                random.shuffle(feature_sets)
                train_set,test_set = feature_sets[:k],feature_sets[k:]
                self.classifier.train(train_set)
                p = nltk.classify.accuracy(self.classifier,test_set)
                print len(train_set),len(test_set),p
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
        pred_labels = self.classifier.batch_classify(features)   
        return ref_labels,pred_labels
    
    """ nltk confusion matrix """
    def confusionMatrix(self,ref,test):
        ref.sort(key=lambda x: x[0])
        test.sort(key=lambda x: x[0])
        _,ref_labels = zip(*ref)
        _,test_labels = zip(*test)
        cm = ConfusionMatrix(ref_labels, test_labels)
        return cm

    def prob_classify(self,db,fastain):
        proIDs,pds,labels = [],[],[]
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
                pd = None
            else:
                text = word_reg.findall(text)
                
            
                featureset = self.gene_features(text)
                assert text!=prevText
                assert featureset!=prevFeatureset
                prevFeatureset = featureset
                prevText = text
                label = self.classifier.batch_classify(featureset)    
                pd = self.classifier.prob_classify([featureset])[0]
                    
            proIDs.append(proteinID)  
            pds.append(pd)
            labels+=label
        return proIDs,labels,pds

    def classifyPickle(self,pickle,fastain):
        proIDs,features,labels = [],[],[]
        prevFeatureset = ''
        prevText = ''
        gbkTable = genbank.GenBankTable()
        gbkTable.load(pickle)
        for seq_record in SeqIO.parse(fastain, "fasta"):
            title = seq_record.id
            toks = title.split("|")
            locus_tag = toks[5]
            text = gbkTable.getLocusText(locus_tag)
            if text=='': 
                label = 'na'
            else:
                text = word_reg.findall(text)
                featureset = self.gene_features(text)
                #assert text!=prevText
                #assert featureset!=prevFeatureset
                prevFeatureset = featureset
                prevText = text
                label = self.classifier.classify(featureset)    
                #print label,text
            proIDs.append(locus_tag)  
            labels.append(label)
        return zip(proIDs,labels)
        
    """ Classifies proteins based on its text from sqlite3 database"""
    def classifyDB(self,db,fastain):
        proIDs,features,labels = [],[],[]
        prevFeatureset = ''
        prevText = ''
        for seq_record in SeqIO.parse(fastain, "fasta"):
            title = seq_record.id
            toks = title.split("|")
            locus_tag = toks[5]
            locus_rows = genbank_sqlite3.locusQuery(locus_tag,db)
            protein_rows = []
            for row in locus_rows:
                locus,proteinID = row
                query_rows = genbank_sqlite3.proteinQuery(proteinID,db)
                protein_rows+=query_rows
            #print len(protein_rows),locus_tag
            if len(protein_rows)==0:
                label = 'na'
            else:
                ids,text = zip(*protein_rows)
                text = ''.join(map(str,text))
                if text=='': 
                    label = 'na'
                else:
                    text = word_reg.findall(text)
                    featureset = self.gene_features(text)
                    #assert text!=prevText
                    #assert featureset!=prevFeatureset
                    prevFeatureset = featureset
                    prevText = text
                    label = self.classifier.classify(featureset)    
                    #print label,text
            proIDs.append(locus_tag)  
            labels.append(label)
        return zip(proIDs,labels)

    def classify(self,dbin,fastain,type='sqlite3'):
        if type=='sqlite3':
            return self.classifyDB(dbin,fastain)
        else:
            return self.classifyPickle(dbin,fastain)

if __name__=="__main__":
    parser = argparse.ArgumentParser(description=\
        'A naive bayes classifier that attempts to categorize context genes')
    parser.add_argument(\
        '--training-labels', type=str, required=False,
        help='A training data set to serve as a template for categorizing context genes')
    parser.add_argument(\
        '--training-directory', type=str, required=False,
        help='A directory containing all of the genbank files required for training')
    parser.add_argument(\
        '--text-database', type=str, required=False,
        help='Path to database containing text annotations')
    parser.add_argument(\
        '--database-type', type=str, required=False,default="sqlite3",
        help='Type of database either "sqlite3" or "pickle"')
    parser.add_argument(\
        '--fasta', type=str, required=False,
        help='Fasta file of sequences to be classified')
    parser.add_argument(\
        '--num-trees', type=int, required=False,default=100,
        help='Number of trees to construct')
    parser.add_argument(\
        '--num-jobs', type=int, required=False,default=3,
        help='Number of jobs to run in parallel')
    parser.add_argument(\
        '--output', type=str, required=False,
        help='Output file containing classifications')
    parser.add_argument(\
        '--test', action='store_const', const=True, default=False,
        help='Run unittests')
    parser.add_argument(\
        '--profile', action='store_const', const=True, default=False,
        help='Run speed tests')
    args = parser.parse_args()
    
    if not args.test and not args.profile:
        classifier = RForests(args.training_directory,
                              args.training_labels,
                              args.num_trees,
                              args.num_jobs)
        classifier.train()
        print "Trained"
        sets = classifier.classify(args.text_database,args.fasta,args.database_type)
        print "Classified"
        titles,labels = zip(*sets)
        open(args.output,'w').write('\n'.join(["%s\t%s"%x for x in sets])+"\n")
    else:
        del sys.argv[1:]
        import cProfile
        import unittest
        import test_fasta
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
                self.labelFile = "%s/data/training/training_proteins4.txt"%self.root
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
                self.labelFile = "%s/data/training/training_proteins4.txt"%self.root
                self.db = "%s/db/bacteria_database"%(self.root,)
                self.zip = "test_serial.zip"
                self.pickle = "test_pickle.zip"
                self.seqs = test_fasta.cluster
                self.fasta = "test.fa"
                open(self.fasta,'w').write(self.seqs)
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
                gbk = genbank.GenBankTable(self.genbank_files)
                gbk.buildLocusTable()
                gbk.buildProteinTable()
                gbk.dump(self.pickle)
                
            def tearDown(self):
                os.remove(self.fasta)
                os.remove(self.pickle)
                os.remove(self.zip)
                pass
            def test1(self):
                #Labs = training.setup(self.genbankDir,self.labelFile)
                nb = RForests(self.trainDir,self.labelFile)
                ref,pred = nb.testClassify(30)
                self.assertTrue(len(ref)>0)
                self.assertTrue(len(pred)>0)
            """
            def test2(self):
                nb = RForests(self.trainDir,self.labelFile)
                nb.train()
                classes = nb.classify(self.db,self.fasta)
                self.assertTrue(len(classes)>0)
            """
            def test3(self):
                nb = RForests(self.trainDir,self.labelFile)
                nb.train()
                classes = nb.classify(self.pickle,self.fasta,type="pickle")
                self.assertTrue(len(classes)>0)
        if args.test:
            unittest.main()
        else:
            cProfile.run('unittest.main()')
            
