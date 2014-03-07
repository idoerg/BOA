"""
Literally a collection of intervals. 
This data structure makes sure that no two intervals overlap each other
"""
import bisect

    

class Intervals(object):

    def __init__(self):
        self.intervals = []
        self.starts= []
        self.ends= []
        self.strands = []
        
    def append(self,x):
        self.intervals.append(x)

        # self.intervals.sort(key=lambda tup: tup[2])
        # self.intervals.sort(key=lambda tup: tup[1])
        # self.intervals.sort(key=lambda tup: tup[0])
        # self.starts,self.ends,self.strands = zip(*self.intervals)
        #bisect.insort_left(self.intervals, x )
        
    def setIntervals(self,intervals):
        self.intervals = intervals
        #self.intervals.sort(key=lambda tup: tup[2])
        #self.intervals.sort(key=lambda tup: tup[1])
        #self.intervals.sort(key=lambda tup: tup[0])
        #self.starts,self.ends,self.strands = zip(*self.intervals)
        
    def __contains__(self,x):
        qstart,qend = x[0],x[1]
        # l_index = bisect.bisect_left(self.starts,qstart)
        # r_index = bisect.bisect_left(self.ends,qend)
        # print l_index,r_index
        # start,end = self.starts[l_index],self.ends[r_index]
        # if l_index == len(self.starts):
        #     return False
        # if l_index==r_index:
        #     return True
        # if ((qstart>=start and qstart<=end) or
        #     (qend>=start and qend<=end) or
        #     (qstart<=start and qend>=end)):
        #     return True
        # return False
        """TODO: This is linear time.
           I will eventually have to implement an interval tree to speed it up
        """
        for interval in self.intervals: 
            start,end = interval[0],interval[1]
            if ((qstart>=start and qstart<=end) or
                (qend>=start and qend<=end) or
                (qstart<=start and qend>=end)):
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

#Obtains original format of the interval
def reformat(interval,radius):
    start_bound,end_bound,refid,gene = interval
    start,end = start_bound+radius,end_bound-radius
    return (start,end,refid,gene)
    
if __name__=="__main__":
    import unittest
    class TestIntervals(unittest.TestCase):
        def test1(self):
            ints = Intervals()
            ints.append( ( 1,10,'+') )
            self.assertTrue((5,6,'+') in ints)
            self.assertTrue((1,10,'+') in ints)

        def test2(self):
            ints = Intervals()
            ints.append( ( 1,10,'+') )
            ints.append( ( 21,30,'+') )
            self.assertTrue((25,26,'+') in ints)
            self.assertTrue((21,28,'+') in ints)

        def test3(self):
            ints = Intervals()
            ints.append( ( 1,10,'+') )
            ints.append( ( 21,30,'+') )
            self.assertTrue((25,26,'-') in ints)
            self.assertTrue((21,28,'-') in ints)

    unittest.main()
