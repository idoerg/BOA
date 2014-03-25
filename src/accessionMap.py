
"""
Convert accession ids to green gene ids
"""
class AccessionGG(object):
    def __init__(self,input_file):
        self.input = input_file
        self.ggtable = dict()
        self.buildTable()
    def buildTable(self):
        with open(self.input,'r') as handle:
            for ln in handle:
                if ln[0]=='#': continue
                ln = ln.rstrip()
                toks = ln.split('\t')
                gg,db,acc = toks #green genes, database type, accession
                #print gg,db,acc
                self.ggtable[(db,acc)] = gg
    def lookupGenbank(self,acc):
        try:
            return self.ggtable[("Genbank",acc)]
        except:
            print "Accession ID %s not found in table"%acc
            return None
            #raise

"""
Convert accession ids to green gene ids
"""
class GGAccession(object):
    def __init__(self,input_file):
        self.input = input_file
        self.ggtable = dict()
        self.buildTable()
    def buildTable(self):
        with open(self.input,'r') as handle:
            for ln in handle:
                if ln[0]=='#': continue
                toks = ln.split('\t')
                gg,db,acc = toks #green genes, database type, accession
                self.ggtable[(db,gg)] = accession
    def lookupGenbank(self,gg):
        try:
            return self.ggtable[("Genbank",gg)]
        except:
            print "Accession ID %s not found in table"%acc
            return None
            #raise


                
