"""
This module conducts a majority vote on each cluster to determine the function
There should only be 7 possible function classes: toxins, modifiers, regulators, transporters, immunity, null, NA (not assigned)
"""
from collections import Counter
base_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
for directory_name in os.listdir(base_path):
    site.addsitedir(os.path.join(base_path, directory_name))


class MVote(object):
    def __init__(self,titles,labels,clusters):
        self.titles = titles
        self.labels = labels
        self.clusters = clusters
        self.labeldict = {self.titles[i]:self.labels[i] 
                          for i in range(len(self.titles))}
    def count(self):
        
        for cluster in self.clusters:
            funcCounts = Counter() #Counts the number of assigned functions
            continue
        pass
    def vote(self):
        pass
    def correct(self):
        pass
    
if __name__=="__main__":
    import unittest
    
    
    
    unittest.main()
    