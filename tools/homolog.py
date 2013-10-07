#in this version i have removed the seq, as it serves no purpose for what i am doing now
import math

class Homolog:
    """This is a class that will hold the values that i wish to store about homologs"""
    #def __init__(self, Accession, Organism, Locus, Gene, Predicted_gene, Synonyms, Eval, Percent_ident, Bits_score, GC, Start, Stop, Strand, Product_type, HGT_candidate = {'likelyhood':'not_eval', 'method': 'none' , 'eval_score': 0.0, 'eval_thresh': 0.0}): #, Seq):
    def __init__(self, Accession, Organism, Locus, Gene, Predicted_gene, Synonyms, Eval, Percent_ident, Bits_score, GC, Start, Stop, Strand, Product_type, Alignment_length, Method, Source_accession, Source_common, Source_locus, Source_start, HGT): #, Seq):
        "This will initialize Homolog, the assumed format is all strings, but i am making an  exception for HGT_candidate here since it is an ad-hoc improvement right now"
        # i want to rename accession to nc
        self.__accession = str(Accession)
        # i want to rename organism to common (since it is the common name)
        self.__organism = str(Organism)
        self.__locus = str(Locus)
        self.__gene = str(Gene)
        self.__predicted_gene = str(Predicted_gene)
        self.__synonyms = Synonyms
        self.__eval = float(Eval)
        self.__percent_ident = float(Percent_ident)
        # i do not think that i will ever use this (bits_score) or need more than 640k or ram!
        self.__bits_score = float(Bits_score)
        self.__gc = float(GC)
        self.__start = int(Start)
        self.__stop = int(Stop)
        self.__strand = int(Strand)
        self.__product_type = str(Product_type)
        self.__alignment_length = int(Alignment_length)
        # Method currently has two values, {exact, BLAST} no idea if more should be enacted or if i should make BLAST read blast
        self.__method = str(Method)
        self.__source_accession = str(Source_accession)
        self.__source_common = str(Source_common)
        # The following two lines have a unique locator for the source gene, which i will likely need later
        self.__source_locus = str(Source_locus)
        self.__source_start = int(Source_start)
        self.__hgt = str(HGT)
        # i think that this variable should have 4 parameters: going to just dictionary this by name
        # 1) likelyhood :'not_eval' - no evaluation performed. 'probable' - evaluation suggests HGT likely. 
        #    'improbable' - evaluation suggests HGT is unlikely. 'disagreement' - more than one method used, no consensus (for later)
        # 2) method: used to determine which method used in evaluation
        # 3) eval_score
        # 4) eval_thresh
        # I would like this value to be retained as a Dict internally, but i would like to have the ability to read it in as a string which is first ',' separated then ':' separated
        #if str(type(HGT_candidate)) == "<type 'dict'>":
        #    self.__hgt_candidate = HGT_candidate
        #else: # txt is assumed... piss off otherwise, you deserve the catastrophic errors the program will throw as a result :)
        #    self.__hgt_candidate = {}
        #    tmp = HGT_candidate.replace(' ', '').split(',')
        #    for i in tmp[:2]:
        #        self.__hgt_candidate.update({i.split(':')[0]:i.split(':')[1]})
        #    for i in tmp[2:]:
        #        self.__hgt_candidate.update({i.split(':')[0]: float(i.split(':')[1])})
                
        #self.__seq = Seq
    
    @classmethod
    def from_file(cls, line):
        a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u = line.strip().split('\t')
        return Homolog(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u)
        
        
    
    
        
    def accession(self):
        """Return the accession number of the organism where the homolog is detected. The return type is string."""
        return self.__accession
    
    def organism(self):
        """Return the organism's common name of the organism where the homolog is detected. The return type is string."""
        return self.__organism
        
    def locus(self):
        """Return the locus of the homolog detected. The return type is string."""
        return self.__locus
        
    def gene(self):
        return self.__gene
        
    def predicted_gene(self):
        return self.__predicted_gene
        
    def synonyms(self): # i should just standardize this outright.... ugh do this after i a test the new version a bit
        return self.__synonyms
    
    def e_val(self):
        return self.__eval
        
    def percent_ident(self):
        return self.__percent_ident
    
    def bits_score(self):
        return self.__bits_score
        
    def gc(self):
        return self.__gc
    
    def start(self):
        return self.__start
            
    def stop(self):
        return self.__stop
        
    def strand(self):
        return self.__strand
        
    def product_type(self):
        return self.__product_type
    
    # TO DO: add back when i get HGT shit together
    #def hgt_candidate(self):
    #    return self.__hgt_candidate
        
    def alignment_length(self):
        return self.__alignment_length
    
    def method(self):
        return self.__method
        
    def source_accession(self):
        return self.__source_accession
        
    def source_common(self):
        return self.__source_common

    def source_locus(self):
        return self.__source_locus
        
    def source_start(self):
        return self.__source_start
    
    def hgt(self):
        return self.__hgt
        
    def nucleotide_length(self):
        return int(math.fabs(self.__start - self.__stop))
        
    def amino_acid_length(self):
        return int(math.fabs(self.__start - self.__stop))/3
    
    
    # TO DO: HGT has to be dealt with in an appropriate way. I do not have a way i want to specify it currently.
        
    #def hgt_candidate_str(self):
    #    ret_list = []
    #    order_list = ['likelyhood', 'method', 'eval_score', 'eval_thresh']
    #    for att in order_list:
    #        ret_list.append(att + ':' + str(self.__hgt_candidate[att]))
    #    return ','.join(ret_list)

    # i have not changed most of these functions to accomodate HGT candidate yet, i would like to define this part more carefully.
    
    def ret_str(self, delim = '\t'):
        #return delim.join([self.accession(), self.organism(), self.locus(), self.gene(), self.predicted_gene(), self.synonyms(), str(self.e_val()), str(self.percent_ident()),  str(self.bits_score()), str(self.gc()), str(self.start()), str(self.stop()), str(self.strand()), self.product_type(), self.hgt_candidate_str()])
        return delim.join([str(i) for i in[self.accession(), self.organism(), self.locus(), self.gene(), self.predicted_gene(), self.synonyms(), str(self.e_val()), str(self.percent_ident()),  str(self.bits_score()), str(self.gc()), str(self.start()), str(self.stop()), str(self.strand()), self.product_type(), self.alignment_length(), self.method(), self.source_accession(), self.source_common(), self.source_locus(), self.source_start(), self.hgt()]])

    
    # TO DO: overload the print function, so it is not Print.  see http://stackoverflow.com/questions/550470/overload-print-python for details here. 
    def Print(self):
        #length = len(self.seq) - 1
        print self.ret_str()
    
    # I think that i will be removing this but will have to do the testing on it.  ugh, hate screwing with this part of my code    
    def ReturnVals(self): # I DO NOT LIKE THIS... 
        return self.accession(), self.organism(), self.locus(), self.gene(), self.predicted_gene(), ':'.join(self.synonyms), self.e_val(), self.percent_ident(), self.bits_score(), self.gc(), self.start(), self.stop(), self.strand(), self.product_type(), self.hgt_candidate_str()
    
    
    # i have no idea what the hell the point of this even is.... i have something above that is capable of dealing with this functionality already.
    def ReturnHomologStr(self, delim = '\t'):
        result = delim.join([self.accession(), self.organism(), self.locus(), self.gene(), self.predicted_gene(), self.synonyms(), str(self.e_val()), str(self.percent_ident()),  str(self.bits_score()), str(self.gc()), str(self.start()), str(self.stop()), str(self.strand()), self.product_type()])
        return result
        
    
