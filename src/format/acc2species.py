"""
Converts accesion id into a species name
"""

class AccessionToSpecies(object):
    def __init__(self,tableFile):
        self.tableFile = tableFile
        self.accessions = dict()
        self.buildTable()
    def buildTable(self):
        with open(self.tableFile,'r') as handle:
            for ln in handle:
                ln = ln.rstrip()
                if ln[0]=="#": continue
                toks = ln.split('\t')
                accID,seqType = toks[0],toks[1]
                species = '_'.join(toks[2:])
                self.accessions[accID] = (seqType,species)
    def lookUp(self,accID):
           return self.accessions[accID]
       

            