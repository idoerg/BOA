"""
Simulator used to help train the naive bayes model

There are two stages of simulation
1) Simulation of text
Output: Text database
2) Simulation of biological sequences
Output: Fasta file
"""
import argparse
import random
import string
import os,site,sys
from collections import defaultdict
from collections import Counter
base_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
for directory_name in os.listdir(base_path):
    site.addsitedir(os.path.join(base_path, directory_name))
import training
import genbank
class TextSeqSim(object):
    """ This assumes that all of the genbank files associated with
        the training labels have already been downloaded"""
    def __init__(self,genbankDir,trainingFile):
        self.trainingFile = trainingFile
        self.textsim = TextSim(genbankDir,trainingFile)
        self.seqsim = SeqSim()    
        self.textsim.textBuild()
        self.seqsim.seqBuild()
        self.functions = ['toxin','modifier','regulator','transport','immunity','null','na']
    """ Generate unique locus tags"""
    def genLocusTags(self,numEntries):
        loci = set()
        for i in range(numEntries):
            locus = ''.join([random.choice(string.ascii_uppercase) 
                             for i in range(10)])
            while locus in loci:
                locus = ''.join([random.choice(string.ascii_uppercase) 
                                 for i in range(10)])
                print "Locus",locus
            loci.add(locus)    
        loci = list(loci)
        return loci
    """ Generate unique protein IDs"""
    def genProteinIDs(self,numEntries):
        proteinIDs = set()
        for i in range(numEntries):
            proteinID = ''.join([random.choice(string.ascii_uppercase) 
                                 for i in range(10)]+['.1'])
            while proteinID in proteinIDs:
                proteinID = ''.join([random.choice(string.ascii_uppercase) 
                                 for i in range(10)]+['.1'])
                print "Protein",proteinID
            proteinIDs.add(proteinID)    
        proteinIDs = list(proteinIDs)
        return proteinIDs
    """ Generate simulated text databases and fasta files"""
    def simulationOutput(self,numEntries,outputFunction,outputDB,outputFASTA):
        loci    = self.genLocusTags(numEntries)
        protIDs = self.genProteinIDs(numEntries)
        fastahandle = open(outputFASTA,'w')
        genbank.buildLocusTable(outputDB)
        genbank.buildProteinTable(outputDB)
        funcHandle = open(outputFunction,'w')
        for i in range(numEntries):
            locus,protID = loci[i],protIDs[i]
            function = random.choice(self.functions)
            numWords = 20
            numMutations = 30
            words = self.textsim.textSample(function,numWords)
            sequence = self.seqsim.seqSample(function,numMutations)
            note = " ".join(words)
            genbank.insertProteinQuery(protID,note,outputDB)
            genbank.insertLocusQuery(locus,protID,outputDB)
            seqid = "%s|_|_|_|%s"%(locus,protID)
            funcHandle.write("%s\t%s\t%s\n"%(locus,protID,function) )
            fastahandle.write(">%s\n%s\n"%(seqid,sequence))
        fastahandle.close()
        funcHandle.close()
    """Read simulated output"""
    def readOutput(self,fname):
        protids,functions = [],[]
        with open(fname,'r') as handle:
            for ln in handle:
                ln = ln.rstrip()
                toks = ln.split('\t')
                protid,function = toks[1],toks[2]
                protids.append(protid)
                functions.append(function)
        return protids,functions
    
class TextSim(object):
    def __init__(self,genbankDir,trainingFile):
        self.trainingFile = trainingFile
        self.genbankDir = genbankDir
        self.distribution = defaultdict(list)
    """ Build a distribution of words based on genbank files and training labels"""
    def textBuild(self):
        #nb = nbayes.NBayes(self.genbankDir,self.trainingFile)
        #nb.train()
        labs = training.Labels(self.genbankDir,self.trainingFile)
        values =  labs.getTrainingText()
        text,functions = zip(*values)
        for entry in zip(functions,text):
            function,words = entry
            self.distribution[function]+=words
    """ Given a function, sample a bunch of words"""
    def textSample(self,function,numWords):
        #pd = nb.classifier.prob
        words = []
        if function not in self.distribution:
            return ['']
        for i in range(numWords):
            word = random.choice(self.distribution[function])
            words.append(word)
        return words
        
class SeqSim(object):
    def __init__(self):
        self.distribution = defaultdict(list)
        self.amino_acids = ['A','R','N',"D",'C',
                            'Q','E','G',"H",'I',
                            'L',"K","M",'F','P',
                            'S','T','W','V','Y']
        self.functions = ['toxin','modifier','regulator','transport','immunity','null','na']
    """ Builds a set of related biological sequences for each training label"""
    def seqBuild(self):
        length = 200
        for function in self.functions:
            #length = random.randint(50,400)
            seq = [random.choice(self.amino_acids) for i in range(0,length)]
            self.distribution[function] = seq 
    """ Grabs sequence and applies mutations using uniform distribution"""   
    def seqSample(self,function,numMutations):
        seq = self.distribution[function]
        pos = random.randint(0,len(seq)-1)
        for i in range(len(seq)):
            samples = random.sample(self.amino_acids,2)
            aa = seq[pos]
            if samples[0]==aa: seq[pos] = samples[1]
            else:              seq[pos] = samples[0]
        return ''.join(seq)
        

if __name__=="__main__":
    parser = argparse.ArgumentParser(description=\
        'Simulates annotation text and sequences for testing naive bayes')
    parser.add_argument(\
        '--genome-files', type=str, nargs="+", required=False,
        help='FASTA files containing bacterial genomes')
    parser.add_argument(\
        '--output-db', type=str, required=False,
        help='The output file containing the sqlite3 output')
    parser.add_argument(\
        '--output-fasta', type=str, required=False,
        help='The output file containing the sqlite3 output')
    parser.add_argument(\
        '--test', action='store_const', const=True, default=False,
        help='Run unittests')
    args = parser.parse_args()
    if not args.test:
        
        pass
    else:
        del sys.argv[1:]
        import unittest
        class TestSimulation1(unittest.TestCase):
            def setUp(self):
                #self.genbankDir = "../example/Streptococcus_pyogenes"
                #self.genbankFile = "../example/Streptococcus_pyogenes/NC_011375.gbk"
                self.root = os.environ['BACFINDER_HOME']
                self.trainDir = "%s/data/training/protein"%self.root
                self.labelFile = "%s/data/training/training_proteins.txt"%self.root
                self.genbank_files    = []
                for file in os.listdir(self.trainDir):
                        if file.endswith(".gbk"):
                            self.genbank_files.append("%s/%s"%(self.trainDir,file)) 
                self.db = "testdb"
                if os.path.exists(self.db):
                    os.remove(self.db)
                
                self.fasta = "test.fa"
                self.function_check="correct_functions.txt"
            def tearDown(self):
                os.remove(self.fasta)
                os.remove(testdb)
                os.remove(self.function_check)
            def testTextSim(self):
                testsim = TextSim(self.trainDir,self.labelFile)
                testsim.textBuild()
                words = testsim.textSample("modifier",100)
                self.assertTrue(len(words)>0)
                words = testsim.textSample("toxin",100)
                self.assertTrue(len(words)>0)
                self.assertTrue("bacteriocin" in words)
                self.assertTrue("precursor" in words)
            def testSeqSim(self):
                seqsim = SeqSim()
                seqsim.seqBuild()
                toxin1   = seqsim.seqSample("toxin",10)
                toxin2   = seqsim.seqSample("toxin",10)
                modifier = seqsim.seqSample("modifier",10)
                self.assertFalse(toxin1==toxin2)
                self.assertFalse(toxin1==modifier)
            def testSim(self):
                nsim = TextSeqSim(self.trainDir,self.labelFile)
                nsim.simulationOutput(10,self.function_check,self.db,self.fasta)
                rows = genbank.proteinQueryAll(self.db)
                print '\n'.join(map(str,rows))
                self.assertEquals(10,len(rows))
                rows = genbank.locusQueryAll(self.db)
                print '\n'.join(map(str,rows))
                self.assertEquals(10,len(rows))
                self.assertTrue(os.path.getsize(self.fasta)>0)                
        unittest.main()
