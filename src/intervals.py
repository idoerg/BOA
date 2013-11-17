"""
Literally a collection of intervals. 
This data structure makes sure that no two intervals overlap each other
"""
import bisect

class Intervals(object):

    def __init__(self):
        self.intervals = []

    def append(self,x):
        self.intervals.append(x)
        #bisect.insort_left(self.intervals, x )
        
    def __contains__(self,x):
        for interval in self.intervals:
            qstart,qend = x[0],x[1]
            start,end = interval[0],interval[1]
            if qstart>=start and qstart<=end or qend>=start and qend<=end:
                return True
        return False

    def search(self,x):
        for interval in self.intervals:
            qstart,qend = x[0],x[1]
            start,end = interval[0],interval[1]
            if qstart>=start and qstart<=end or qend>=start and qend<=end:
                return interval
        return None

    def __str__(self):
        return "\n".join( map( str,self.intervals))
if __name__=="__main__":
    import unittest
    class TestIntervals(unittest.TestCase):
        def test1(self):
            ints = Intervals()
            ints.append( ( 1,10) )
            self.assertTrue((5,6) in ints)
            self.assertTrue((1,10) in ints)

        def test2(self):
            ints = Intervals()
            ints.append( ( 1,10) )
            ints.append( ( 21,30) )
            self.assertTrue((25,26) in ints)
            self.assertTrue((21,28) in ints)


    unittest.main()
