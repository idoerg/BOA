
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
        for i in xrange(len(hits)):
            for j in xrange(0,i):
                hiti = hits[i]
                hitj = hits[j]
                
                ienv_st,ienv_end = int(hiti[7]),int(hiti[8])
                jenv_st,jenv_end = int(hitj[7]),int(hitj[8])
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
                query = toks[2]
                func = query.split(".")[0]
                function.add(func)
            if functions.issuperset(set(self.functions)):
                clusters.append(clique)
        return clusters
    
if __name__=="__main__":
    import unittest
    class TestCase(unittest.TestCase):
        def setUp(self):
            self.queries   = [('CP002279.1_3','toxin.fa.cluster2.fa',0,9.7e-13,1.7e-12,0,100,25000,25100,
                               'Mesorhizobium opportunistum WSM2075, complete genome'),
                              ('CP002279.1_3','transport.fa.cluster2.fa',0,9.7e-13,1.7e-12,0,100,25200,25300,
                               'Mesorhizobium opportunistum WSM2075, complete genome'),
                              ('CP002279.1_3','modifier.fa.cluster2.fa',0,9.7e-13,1.7e-12,0,100,25400,25500,
                               'Mesorhizobium opportunistum WSM2075, complete genome'),
                              ('CP002279.1_3','regulator.fa.cluster2.fa',0,9.7e-13,1.7e-12,0,100,25600,25700,
                               'Mesorhizobium opportunistum WSM2075, complete genome'),
                              ('CP002279.1_3','immunity.fa.cluster2.fa',0,9.7e-13,1.7e-12,0,100,25800,26900,
                               'Mesorhizobium opportunistum WSM2075, complete genome'),
                              ('CP002279.1_3','transport.fa.cluster2.fa',0,9.7e-13,1.7e-12,0,100,35127,35356,
                               'Mesorhizobium opportunistum WSM2075, complete genome'),
                              ('CP002279.1_3','transport.fa.cluster2.fa',0,9.7e-13,1.7e-12,0,100,45127,45356,
                               'Mesorhizobium opportunistum WSM2075, complete genome')]              

        
        def test1(self):
            cfilter = CliqueFilter(radius=1000)
            cfilter.createGraph(self.queries)
            clusters = cfilter.fullcolorFilter()
            for cluster in clusters:
                self.assertTrue(len(cluster)<=5)
                if len(cluster)==5:
                     self.assertItemsEqual(cluster,map(str,self.queries[:5]))
    unittest.main()








