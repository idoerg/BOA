"""
Given a newick tree, this computes a hypermetric tree
"""

import tpg
import argparse
import os,sys,site

class Leaf(object):
    def __init__(self, label):
        self.label = label
    def __str__(self):
        return self.label
    def str2(self):#another string representation
        return self.label
    
class Internal(object):
    def __init__(self, label, subtrees):
        self.label = label
        self.subtrees = subtrees

    def __str__(self):
        return '(%s)%s' % (','.join(str(st) for st in self.subtrees),self.label)
    def str2(self):#Another string representation
        return '(%s)%s' % (','.join(str(st.str2()) for st in self.subtrees),self.label)

class Branch(object):
    def __init__(self, subtree, length):
        self.subtree = subtree
        self.length = length
        self.depth = 0
    def getDepth(self):
        return self.depth
    def setDepth(self,depth):
        self.depth = depth
    def setLength(self,length):
        self.length = length
    def __str__(self):
        return '%s:%s' % (str(self.subtree),str(self.length))
    def str2(self):#Another string representation
        return '%s:%s' % (str(self.subtree.str2()),str(self.depth))
        
class HyperTree(object):
    def __init__(self,treeFile,totLen=1,pctLen=0.5):
        treeStr = open(treeFile).readline().rstrip()
        treeParser = Parser()
        self.tree = treeParser(treeStr)
        self.pctSplit = pctLen    # for calculating branch length
        self.totalLength = totLen #total length of tree
        #self.tree  = tree
        self.setDepths(self.tree,1)
        self.depth = self.maxDepth(self.tree,1)
        self.normalizeLengths(self.tree,self.totalLength)
         
        """ Sets depths for all nodes """
    def setDepths(self,tree,depth):
        if type(tree)==Leaf:
            return
        for br in tree.subtrees:
            br.setDepth(depth)
            self.setDepths(br.subtree,depth+1)
    """ Sets max depth of tree """
    def maxDepth(self,tree,maxD):
        if type(tree)==Leaf:
            return
        for br in tree.subtrees:
            if br.getDepth()>maxD:
                maxD = br.getDepth()
            else:
                d = self.maxDepth(br.subtree,maxD)
                if d > maxD:
                    maxD = d
        return maxD
    
    def __str__(self):#Another string representation
        return "%s;"%str(self.tree)
    def str2(self):#Another string representation
        return "%s;"%str(self.tree.str2())
    
    """ Normalize lengths"""
    def normalizeLengths(self,tree,lenRemaining):
        for br in tree.subtrees:
            if type(br.subtree)==Leaf:
                br.setLength(lenRemaining)
            else:           
                brLen = (1-self.pctSplit)*lenRemaining     
                newlen = self.pctSplit*lenRemaining #length of branch
                br.setLength(brLen)
                self.normalizeLengths(br.subtree,newlen)
        
class Parser(tpg.Parser):
    r"""

    separator space	'\s+' 							;

    token label		'[^,:;() \t\n]+' ;

    START/tree ->
    	Subtree/tree ';'
    ;

    Subtree/tree ->
    	Internal/tree
      | Leaf/tree
    ;

    Leaf/leaf ->
    	Name/name				$ leaf = Leaf(name)
    ;

    Internal/tree ->
    	'\(' BranchList/st '\)' Name/n		$ tree = Internal(n,st)
    ;

    BranchList/blist ->
        Branch/b (',' BranchList/bl $ blist = bl+[b] $ )+
      | Branch/branch				$ blist = [branch]
    ;

    Branch/branch ->
    	Subtree/st Length/l			$ branch = Branch(st,l)
    ;

    Name/name ->
    	label/name
      | Empty					$ name = ""
    ;

    Length/length ->
    	':' Number/length
      | Empty					$ length = 0
    ;

    Number/n ->
        label/l					$ n = float(l)
    ;

    Empty -> ;


    """
if __name__=="__main__":
    parser = argparse.ArgumentParser(description=\
        'Parses a newick tree and transforms it into a hypermetric tree')
    parser.add_argument(\
        '--tree', type=str, required=False,
        help='A training data set to serve as a template for categorizing context genes')
    parser.add_argument(\
        '--length', type=int, required=False,
        help='Total length of the tree')
    parser.add_argument(\
        '--pct', type=float, required=False,
        help='Percentage of branch lengths dedicated to internal nodes')
    parser.add_argument(\
        '--test', action='store_const', const=True, default=False,
        help='Run the unittests')
    args = parser.parse_args()
    if not args.test:
        #tree = open(args.tree).readline().rstrip()
        #treeParser = Parser()
        #t=treeParser(tree)
        h = HyperTree(args.tree,args.length,args.pct)
        print h
    else:
        del sys.argv[1:]
        import unittest
        class TestTree(unittest.TestCase):
            def setUp(self):
                self.file = "test.tree"
                tree = "(A:0.1,B:0.2,(C:0.3,D:0.4):0.5);"
                open(self.file,'w').write(str(tree))
            def tearDown(self):
                os.remove(self.file)
            def test1(self):
                #Just a test case for tinkering
                tree = "(A:0.1,B:0.2,(C:0.3,D:0.4):0.5);"
                treeParser = Parser()
                t = treeParser(tree)
                print dir(t)
                print str(t.subtrees)
                print str(t.subtrees[0])
                print str(t.subtrees[1])
                print t.subtrees[2].subtree
                print str(t)
            def testDepth(self):
                tree = "(A:0.1,B:0.2,(C:0.3,D:0.4):0.5);"
                #treeParser = Parser()
                #t = treeParser(tree)
                h = HyperTree(self.file)
                self.assertEquals(h.str2(),
                                  "((D:2,C:2):1,B:1,A:1);")
            def testMaxDepth(self):
                tree = "(A:0.1,B:0.2,(C:0.3,D:0.4):0.5);"
                #treeParser = Parser()
                #t = treeParser(tree)
                h = HyperTree(self.file)
                self.assertEquals(h.depth,2)
            def testNormalizeDepth(self):
                #tree = "(A:0.1,B:0.2,(C:0.3,D:0.4):0.5);"
                #treeParser = Parser()
                #t = treeParser(tree)
                h = HyperTree(self.file)
                self.assertEquals(str(h),
                                  "((D:0.5,C:0.5):0.5,B:1,A:1);")
                
        unittest.main()
    
    
    
    
    
    