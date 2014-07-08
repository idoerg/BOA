
"""

Creates a colored graph from HMMER output
Edges are created when two genes are within a specified radius of each other.
Colors are assigned according to function (toxin, modifier, immunity, transport, regulator)

"""
import os,sys,site
import numpy
import networkx as nx
base_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
for directory_name in os.listdir(base_path):
    site.addsitedir(os.path.join(base_path, directory_name))
import fasta
class CliqueFilter():
    def __init__(self,fasta_index,radius=50000):
        self.faidx = fasta.Indexer("",fasta_index)
        self.faidx.load()
        self.radius = radius
    
    def createGraph(self,hits):
        self.graph = nx.Graph()
        handle = open("graph_error.txt",'w')
        for i in xrange(len(hits)):
            for j in xrange(0,i):
                hiti = hits[i]
                hitj = hits[j]
                #acc,clrname,full_evalue,hmm_st,hmm_end,env_st,env_end,description=toks
                iname,ienv_st,ienv_end = hiti[0],int(hiti[5]),int(hiti[6])
                jname,jenv_st,jenv_end = hitj[0],int(hitj[5]),int(hitj[6])
                # Translate the six-frame translated coordinates
                # into nucleotide reference coordinates
                ienv_st,istrand = self.faidx.sixframe_to_nucleotide(iname,ienv_st)
                ienv_end,istrand= self.faidx.sixframe_to_nucleotide(iname,ienv_end)
                jenv_st,jstrand = self.faidx.sixframe_to_nucleotide(jname,jenv_st)
                jenv_end,jstrand= self.faidx.sixframe_to_nucleotide(jname,jenv_end)
                midi = (ienv_st+ienv_end)/2
                midj = (jenv_st+jenv_end)/2
                print "midpoints",midi,midj,"strands",istrand,jstrand
                if abs(midi-midj)<self.radius and istrand==jstrand:
                    print "Create edge"
                    #Record genome coordinates of operons
                    iacc,iclrname,ifull_evalue,ihmm_st,ihmm_end,_,_,idescription=hiti
                    jacc,jclrname,jfull_evalue,jhmm_st,jhmm_end,_,_,jdescrjptjon=hitj
                    nodei = "|".join(map(str,[iacc,iclrname,ifull_evalue,ihmm_st,ihmm_end,ienv_st,jenv_end,idescription]))
                    nodej = "|".join(map(str,[jacc,jclrname,jfull_evalue,jhmm_st,jhmm_end,jenv_st,jenv_end,jdescrjptjon]))
                    print nodei
                    print nodej
                    self.graph.add_edge(nodei,nodej)
        print "graph edges",'\n'.join(map(str,self.graph.edges()))
    """    
    The output will be cliques that contain all of the functions specified
    """
    def filter(self,keyfunctions = ["toxin","modifier","immunity","transport","regulator"] ):
        clique_gen = nx.find_cliques(self.graph)
        #print "cliques",'\n'.join(map(str,list(clique_gen)))
        clusters = [] #Context gene clusters 
        for clique in clique_gen:
            print "Len",len(clique)
            #if len(clique)<len(keyfunctions): continue
            functions = set()
            for node in clique:
                toks = node.split("|")
                query = toks[1]
                func = query.split(".")[0]
                functions.add(func)
            print "functions",functions
            if functions.issuperset(set(keyfunctions)):
                clusters.append(clique)
        return clusters

        
"""
Locate all context gene clusters
"""
def findContextGeneClusters(hits,faidx,radius=50000, functions = ["toxin","modifier","immunity","transport","regulator"]):
    err_handle = open('error.log','w')
    prevGenome = None
    buf,clusters = [],[]
    cfilter = CliqueFilter(faidx,radius)
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
            
            cliques = cfilter.filter(functions)
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
def collapseOverlaps(hits,fasta_index):
    newHits = []
    prevHit = None
    faidx = fasta.Indexer("",fasta_index)
    faidx.load()
    print "Before filtering",len(hits)
    for hit in hits:
        if prevHit == None:
            prevHit = hit
            newHits.append(hit)
        else:
            prevName,hitName = prevHit[0],hit[0]
            prevSt,prevEnd = prevHit[5],prevHit[6]
            hitSt,hitEnd = hit[5],hit[6]
           
            prevSt,prevStrand = faidx.sixframe_to_nucleotide(prevName,prevSt)
            prevEnd,prevStrand= faidx.sixframe_to_nucleotide(prevName,prevEnd)
            hitSt,hitStrand  = faidx.sixframe_to_nucleotide(hitName,hitSt)
            hitEnd,hitStrand = faidx.sixframe_to_nucleotide(hitName,hitEnd) 
            prevName = fasta.getName(prevName)
            hitName = fasta.getName(hitName)
            if prevName!=hitName: 
                newHits.append(hit)
            else:
                if hitSt>hitEnd:  hitSt,hitEnd = hitSt,hitEnd
                if prevSt>prevEnd:  hitSt,hitEnd = hitSt,hitEnd
                if hitStrand==prevStrand: 
                    if prevSt>hitEnd: newHits.append(hit)
                    elif prevEnd<hitSt: newHits.append(hit)
                    else: continue
                else:
                    newHits.append(hit)
            prevHit = hit 
    print "After filtering",len(newHits)
    return newHits

if __name__=="__main__":
    import unittest
    class TestCase(unittest.TestCase):
        def setUp(self):
            indexes = [ '\t'.join(map(str,('CP002279.1_1',2294815, 185896721,60,61))),
                        '\t'.join(map(str,('CP002279.1_2',2294815, 188229850,60,61))),
                        '\t'.join(map(str,('CP002279.1_3',2294814, 190562979,60,61))),
                        '\t'.join(map(str,('CP002279.1_4',2294814, 192896107,60,61))),
                        '\t'.join(map(str,('CP002279.1_5',2294815, 195229235,60,61))),
                        '\t'.join(map(str,('CP002279.1_6',2294815, 197562364,60,61)))]
            self.testfai = "test.fai"
            open(self.testfai,'w').write('\n'.join(indexes))
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
                              ('CP002279.1_4','immunity.fa.cluster2.fa',0,0,100, 740038, 740138,
                               'Mesorhizobium opportunistum WSM2075, complete genome'), 
                              ('CP002279.1_3','transport.fa.cluster2.fa',0,0,100,35127,35356,
                               'Mesorhizobium opportunistum WSM2075, complete genome'),
                              ('CP002279.1_3','transport.fa.cluster2.fa',0,0,100,35127,35456,
                               'Mesorhizobium opportunistum WSM2075, complete genome'),
                              ('CP002279.1_2','transport.fa.cluster2.fa',0,0,100,35127,35356,
                               'Mesorhizobium opportunistum WSM2075, complete genome'),
                              ('CP002279.1_4','transport.fa.cluster2.fa',0,0,100,35127,35456,
                               'Mesorhizobium opportunistum WSM2075, complete genome'),
                              ('CP002279.1_3','transport.fa.cluster2.fa',0,0,100,45127,45356,
                               'Mesorhizobium opportunistum WSM2075, complete genome')]   
            self.faidx = fasta.Indexer("",self.testfai)           
            self.faidx.load()
            self.maxDiff=10000
        def tearDown(self):
            if os.path.exists(self.testfai):
                os.remove(self.testfai)
        def test1(self):
            cfilter = CliqueFilter(self.testfai,radius=10000)
            cfilter.createGraph(self.queries)
            clusters = cfilter.filter()
            self.assertTrue(len(clusters)>0)
            for cluster in clusters:
               
                self.assertTrue(len(cluster)<=5)
                if len(cluster)==5:
                    answer = ["|".join(map(str,s)) for s in self.queries[:5]]
                    self.assertItemsEqual(answer,cluster)
        def test2(self):
            #all_hits=sorted(self.queries,key=lambda x: x[6])        
            #all_hits=sorted(all_hits,key=lambda x: x[5])
            #Sort by genome name
            #all_hits=sorted(all_hits,key=lambda x: x[-1])  
            reduced = collapseOverlaps(self.queries,self.testfai)
           
            self.assertItemsEqual(reduced,
                                  [   ('CP002279.1_3','toxin.fa.cluster2.fa',0,0,100,25000,25100,
                                      'Mesorhizobium opportunistum WSM2075, complete genome'),
                                      ('CP002279.1_3','transport.fa.cluster2.fa',0,0,100,25200,25300,
                                       'Mesorhizobium opportunistum WSM2075, complete genome'),
                                      ('CP002279.1_3','modifier.fa.cluster2.fa',0,0,100,25400,25500,
                                       'Mesorhizobium opportunistum WSM2075, complete genome'),
                                      ('CP002279.1_3','regulator.fa.cluster2.fa',0,0,100,25600,25700,
                                       'Mesorhizobium opportunistum WSM2075, complete genome'),
                                      ('CP002279.1_3','immunity.fa.cluster2.fa',0,0,100,25800,26900,
                                       'Mesorhizobium opportunistum WSM2075, complete genome'),
                                      ('CP002279.1_4','immunity.fa.cluster2.fa',0,0,100, 740038, 740138,
                                       'Mesorhizobium opportunistum WSM2075, complete genome'), 
                                      ('CP002279.1_3','transport.fa.cluster2.fa',0,0,100,35127,35356,
                                       'Mesorhizobium opportunistum WSM2075, complete genome'),
                                      ('CP002279.1_4','transport.fa.cluster2.fa',0,0,100,35127,35456,
                                       'Mesorhizobium opportunistum WSM2075, complete genome'),
                                      ('CP002279.1_3','transport.fa.cluster2.fa',0,0,100,45127,45356,
                                       'Mesorhizobium opportunistum WSM2075, complete genome')])  
                                      
                                  
    unittest.main()








