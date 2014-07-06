
"""

Creates a colored graph from HMMER output
Edges are created when two genes are within a specified radius of each other.
Colors are assigned according to function (toxin, modifier, immunity, transport, regulator)

"""
import numpy
import networkx as nx

class CliqueFilter():
    def __init__(self,radius=50000):
        self.radius = radius
        self.functions = ["toxin","modifier","immunity","transport","regulator"]
    def createGraph(self,hits):
        self.graph = nx.Graph()
        handle = open("graph_error.txt",'w')
        for i in xrange(len(hits)):
            for j in xrange(0,i):
                hiti = hits[i]
                hitj = hits[j]
                
                ienv_st,ienv_end = int(hiti[5]),int(hiti[6])
                jenv_st,jenv_end = int(hitj[5]),int(hitj[6])
                midi = (ienv_st+ienv_end)/2
                midj = (jenv_st+jenv_end)/2
                if abs(midi-midj)<self.radius:
                    nodei = "|".join(map(str,hiti))
                    nodej = "|".join(map(str,hitj))
                    self.graph.add_edge(nodei,nodej)
    """
    The output will be cliques that contain a full coloring (aka contain all functions) 
    """
    def fullcolorFilter(self):
        clique_gen = nx.find_cliques(self.graph)
        clusters = [] #Context gene clusters 
        
        for clique in clique_gen:
            
            if len(clique)<len(self.functions): continue
            functions = set()
            for node in clique:
                toks = node.split("|")
                query = toks[1]
                func = query.split(".")[0]
                functions.add(func)
            
            if functions.issuperset(set(self.functions)):
                clusters.append(clique)
        return clusters
"""
Locate all context gene clusters
"""
def findContextGeneClusters(hits,radius=50000):
    err_handle = open('error.log','w')
    prevGenome = None
    buf,clusters = [],[]
    cfilter = CliqueFilter(radius)
    for hit in hits:
        if prevGenome == None:      
            prevGenome = hit[-1]
        elif prevGenome == hit[-1]: 
            buf.append(hit)
        else: 
            print >>err_handle,prevGenome,'\n'
            print >>err_handle,'Buffer'
            print >>err_handle,"\n".join(map(str,buf)),'\n'
            
            cfilter.createGraph(buf)    
            
            cliques = cfilter.fullcolorFilter()
            print >>err_handle,'Cliques'
            print >>err_handle,"\n".join(map(str,cliques)),'\n'
            clusters+= cliques
            buf = [hit]
            prevGenome = hit[-1]
    return clusters
"""
If a hit overlaps with a previous hit,
throw the hit away
"""
def collapseOverlaps(hits):
    newHits = []
    prevHit = None
    for hit in hits:
        if prevHit == None:
            prevHit = hit
            newHits.append(hit)
        else:
            prevSt,prevEnd = prevHit[5],prevHit[6]
            hitSt,hitEnd = hit[5],hit[6]
            """ prev |===============|
                hit    |=========|"""
            if prevSt<=hitSt and prevEnd>=hitEnd:
                continue
                """ prev    |==========|
                    hit  |===============|"""
            elif prevSt>=hitSt and prevEnd<=hitEnd:
                continue
                """ prev    |==========|
                    hit         |===============|"""
            elif hitSt>=prevSt and hitSt<=prevEnd:
                continue
                """ prev           |==========|
                    hit   |===============|        """
            elif prevSt>=hitSt and prevSt<=hitEnd:
                continue
            else:
                newHits.append(hit)
            prevHit = hit
    return newHits

if __name__=="__main__":
    import unittest
    class TestCase(unittest.TestCase):
        def setUp(self):
            self.queries   = [('CP002279.1_3','toxin.fa.cluster2.fa',0,0,100,25000,25100,
                               'Mesorhizobium opportunistum WSM2075, complete genome'),
                              ('CP002279.1_3','transport.fa.cluster2.fa',0,0,100,25200,25300,
                               'Mesorhizobium opportunistum WSM2075, complete genome'),
                              ('CP002279.1_3','modifier.fa.cluster2.fa',0,0,100,25400,25500,
                               'Mesorhizobium opportunistum WSM2075, complete genome'),
                              ('CP002279.1_3','regulator.fa.cluster2.fa',0,0,100,25600,25700,
                               'Mesorhizobium opportunistum WSM2075, complete genome'),
                              ('CP002279.1_3','immunity.fa.cluster2.fa',0,0,100,25800,26900,
                               'Mesorhizobium opportunistum WSM2075, complete genome'),
                              ('CP002279.1_3','transport.fa.cluster2.fa',0,0,100,35127,35356,
                               'Mesorhizobium opportunistum WSM2075, complete genome'),
                              ('CP002279.1_3','transport.fa.cluster2.fa',0,0,100,35127,35456,
                               'Mesorhizobium opportunistum WSM2075, complete genome'),
                              ('CP002279.1_3','transport.fa.cluster2.fa',0,0,100,45127,45356,
                               'Mesorhizobium opportunistum WSM2075, complete genome')]              

            self.maxDiff=1000
        def test1(self):
            cfilter = CliqueFilter(radius=2000)
            cfilter.createGraph(self.queries)
            clusters = cfilter.fullcolorFilter()
            self.assertTrue(len(clusters)>0)
            for cluster in clusters:
                print cluster
                self.assertTrue(len(cluster)<=5)
                if len(cluster)==5:
                    answer = ["|".join(map(str,s)) for s in self.queries[:5]]
                    self.assertItemsEqual(answer,cluster)
        def test2(self):
            reduced = collapseOverlaps(self.queries)
            print "Reduced",reduced
            self.assertItemsEqual(reduced,
                                  [('CP002279.1_3','toxin.fa.cluster2.fa',0,0,100,25000,25100,
                                   'Mesorhizobium opportunistum WSM2075, complete genome'),
                                  ('CP002279.1_3','transport.fa.cluster2.fa',0,0,100,25200,25300,
                                   'Mesorhizobium opportunistum WSM2075, complete genome'),
                                  ('CP002279.1_3','modifier.fa.cluster2.fa',0,0,100,25400,25500,
                                   'Mesorhizobium opportunistum WSM2075, complete genome'),
                                  ('CP002279.1_3','regulator.fa.cluster2.fa',0,0,100,25600,25700,
                                   'Mesorhizobium opportunistum WSM2075, complete genome'),
                                  ('CP002279.1_3','immunity.fa.cluster2.fa',0,0,100,25800,26900,
                                   'Mesorhizobium opportunistum WSM2075, complete genome'),
                                  ('CP002279.1_3','transport.fa.cluster2.fa',0,0,100,35127,35356,
                                   'Mesorhizobium opportunistum WSM2075, complete genome'),
                                  ('CP002279.1_3','transport.fa.cluster2.fa',0,0,100,45127,45356,
                                   'Mesorhizobium opportunistum WSM2075, complete genome')]     )
    unittest.main()








