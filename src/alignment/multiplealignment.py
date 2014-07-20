
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

class MultipleAlignment(object):
    __metaclass__ = abc.ABCMeta
    def __init__(self,input_file,output_file):
        self.input = input_file #input fasta file
        self.output = output_file
        basename,_ = os.path.splitext(input_file)
        #output from clustalw
        self.aln = "%s.aln"%basename
    
    @abc.abstractmethod 
    def wait(self): pass
    @abc.abstractmethod 
    def outputFASTA(self): pass
    @abc.abstractmethod 
    def outputSTO(self): pass
    @abc.abstractmethod 
    def erase_files(self): pass
    
