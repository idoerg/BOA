"""
Literally a collection of intervals. 
This data structure makes sure that no two intervals overlap each other
"""
import bisect

class Intervals(object)

    def __init__(self):
        self.intervals = []

    def append(self,refid,start,end):
        """
        TODO: Will need do the following
        1) Determine where it belongs in the list of intervals
        2) Merge intervals if they overlap
        """
        self.intervals.append( (start,end,refid) )
        
if __name__=="__main__":
    import unittest
    unittest.main()
