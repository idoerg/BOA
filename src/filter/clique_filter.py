
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
from bx.intervals import *
import matplotlib.pyplot as plt
import interval_filter
import itertools

class CliqueFilter():
    def __init__(self,fasta_index="",radius=50000):
        self.faidx = fasta.Indexer("",fasta_index)
        if fasta_index!="":
            self.faidx.load()
        self.radius = radius
    
    def createGraph(self,hits,backtrans=True):
        print "Creating graph %d hits"%(len(hits))
        self.graph = nx.Graph()
        handle = open("graph_error.txt",'w')
        for i in xrange(len(hits)):
            for j in xrange(0,i):
                #print i,j
                hiti = hits[i]
                hitj = hits[j]
                #acc,clrname,full_evalue,hmm_st,hmm_end,env_st,env_end,description=toks
                iname,ienv_st,ienv_end = hiti[0],int(hiti[5]),int(hiti[6])
                jname,jenv_st,jenv_end = hitj[0],int(hitj[5]),int(hitj[6])
                # Translate the six-frame translated coordinates
                # into nucleotide reference coordinates
                if backtrans:
                    inenv_st,istrand  = self.faidx.sixframe_to_nucleotide(iname,ienv_st)
                    inenv_end,istrand = self.faidx.sixframe_to_nucleotide(iname,ienv_end)
                    jnenv_st,jstrand  = self.faidx.sixframe_to_nucleotide(jname,jenv_st)
                    jnenv_end,jstrand = self.faidx.sixframe_to_nucleotide(jname,jenv_end)
                else:
                    inenv_st,inenv_end = ienv_st,ienv_end
                    jnenv_st,jnenv_end = jenv_st,jenv_end
                    istrand = fasta.strand(fasta.getFrame(iname))
                    jstrand = fasta.strand(fasta.getFrame(jname))
                    
                assert inenv_st>=0,"Less than 0, transformed:%d original:%d"%(inenv_st,ienv_st)
                assert jnenv_st>=0,"Less than 0, transformed:%d original:%d"%(jnenv_st,jenv_st)
                midi = (inenv_st+inenv_end)/2
                midj = (jnenv_st+jnenv_end)/2
                
                if abs(midi-midj)<self.radius and istrand==jstrand:
                    
                    #Record genome coordinates of operons
                    iacc,iclrname,ifull_evalue,ihmm_st,ihmm_end,_,_,idescription=hiti
                    jacc,jclrname,jfull_evalue,jhmm_st,jhmm_end,_,_,jdescription=hitj
                    nodei = "|".join(map(str,[iacc,iclrname,ifull_evalue,ihmm_st,ihmm_end,inenv_st,inenv_end,idescription]))
                    nodej = "|".join(map(str,[jacc,jclrname,jfull_evalue,jhmm_st,jhmm_end,jnenv_st,jnenv_end,jdescription]))
                    
                    self.graph.add_edge(nodei,nodej)
                    
        #nx.draw(self.graph)
        #plt.show()
    """    
    The output will be cliques that contain all of the functions specified
    """
    def filter(self,merge=False,keyfunctions = ["toxin","modifier","immunity","transport","regulator"] ):
        clique_gen = nx.find_cliques(self.graph)
        #print "cliques",'\n'.join(map(str,list(clique_gen)))
        clusters = [] #Context gene clusters 
        for clique in clique_gen:
            
            if len(clique)<len(keyfunctions): continue
            functions = set()
            for node in clique:
                toks = node.split("|")
                query = toks[1]
                func = query.split(".")[0]
                functions.add(func)
            
            if functions.issuperset(set(keyfunctions)):
                clusters.append(clique)
        #Merge overlapping cliques
        if merge:
            return self.merge(clusters)
        else:
            return clusters

    """ Returns the start and end coordinates of the entire clique"""
    def envelop(self,clique):
        starts,ends = [],[]
        for node in clique:
            toks = node.split("|")
            st,end = map(int,toks[5:7])
            starts.append(st)
            ends.append(end)
        return min(starts),max(ends)

    """Sort clusters by start/end position"""
    def sort(self,clusters):
        tups = []
        for clique in clusters:
            st,end = self.envelop(clique)
            tups.append((st,end,clique))
        tups = sorted(tups,key=lambda x:x[0])
        tups = sorted(tups,key=lambda x:x[1])
        return zip(*tups)[2]

    def overlap(self,st1,end1,st2,end2):
        assert st1<end1
        assert st2<end2
        if end1<st2:
            return False
        elif st1>end2:
            return False
        else:
            return True

    """ Merge overlapping cliques together """
    def merge(self,clusters):
        clique_intervals = IntervalTree()
        clusters = self.sort(clusters)
        newClusters = []
        #print clusters
        for i in xrange(len(clusters)):
            j = i+1
            icluster = clusters[i]
            if j<len(clusters):   
                jcluster = clusters[j]
                ist,iend = self.envelop(icluster)
                jstart,jend = self.envelop(jcluster)
                overlaps = [icluster]
                while self.overlap(ist,iend,jstart,jend):
                    overlaps.append(jcluster)
                    j+=1
                    if j==len(clusters): break
                    
                    jcluster = clusters[j]
                    jstart,jend = self.envelop(jcluster)
            else:
                overlaps=[icluster]
            group = list(set(itertools.chain(*overlaps)))
            #print '\n'.join(map(str,set(group)))
            newClusters.append(group)
        return newClusters
    
"""
Locate all context gene clusters
"""
def findContextGeneClusters(hits,faidx,radius=50000,backtrans=True, functions = ["toxin","modifier","immunity","transport","regulator"]):
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
            
            cfilter.createGraph(buf,backtrans)    
            
            cliques = cfilter.filter(functions)
            print >>err_handle,'Cliques'
            print >>err_handle,"\n".join(map(str,cliques)),'\n'
            clusters+= cliques
            buf = [hit]
            prevGenome = hit[-1]
    
    cfilter.createGraph(buf,backtrans)    
    cliques = cfilter.filter(functions)
    print >>err_handle,'Cliques'
    print >>err_handle,"\n".join(map(str,cliques)),'\n'
    clusters+= cliques
            
    return clusters



def go(input,faidx,radius,functions):
    cfilter = CliqueFilter(faidx,radius)
    hits = []
    with open(input,'r') as handle:
        for ln in handle:
            ln = ln.rstrip()
            toks = ln.split('|')
            acc,clrname,full_evalue,hmm_st,hmm_end,env_st,env_end,description=toks
            hmm_st,hmm_end,env_st,env_end = map(int,[hmm_st,hmm_end,env_st,env_end])
            full_evalue = float(full_evalue)
            hits.append((acc,clrname,full_evalue,hmm_st,hmm_end,env_st,env_end,description))
    cfilter.createGraph(hits)
    cliques = cfilter.filter(functions)    
    #with open(outfile,'w') as outhandle:
    for clique in cliques:
        if len(clique)>0:
            for gene in clique:
                sys.stdout.write("%s\n"%gene)
            sys.stdout.write('----------\n')
        
    
if __name__=="__main__":
    
    import argparse
    parser = argparse.ArgumentParser(description=\
                                     'Finds operons by looking for cliques. Reads from stdin and writes to stdout')
    parser.add_argument(\
        '--input', type=str, required=False,
        help='Input hmmer hits (in amino acid coordinates)')
    parser.add_argument(\
        '--faidx', type=str, required=False,
        help='The fasta index file')
    parser.add_argument(\
        '--clique-radius', type=int, required=False, default=30000,
        help='The maximum radius of the operons')  
    parser.add_argument(\
        '--functions', type=str, nargs="+", default=None, required=False,
        help='The list of functions to look for (e.g. toxin, modifier, transport, immunity, regulator)')
    parser.add_argument(\
        '--output', type=str, required=False,
        help='Output predicted cliques (converted to genomic coordinates)')
    parser.add_argument(\
        '--test', action='store_const', const=True, default=False,
        help='Run unittests')
    args = parser.parse_args()
    
    if not args.test:
        functions = args.functions
        if args.functions==None:
            functions = ["toxin","modifier","immunity","transport","regulator"]
        go(args.input,args.faidx,args.clique_radius,functions)
        
    else:
        del sys.argv[1:]    
        import unittest
        class TestMerge1(unittest.TestCase):
            def setUp(self):
                indexes = [
                            '\t'.join(map(str,('HE577328.1_4',    588676,  8720859786,      60,      61))),
                            '\t'.join(map(str,('HE577328.1_5',    588676,  8721458351,      60,      61))),
                            '\t'.join(map(str,('HE577328.1_6',    588676,  8722056916,      60,      61))),
                            '\t'.join(map(str,('HE577330.1_1',    259600,  8722655481,      60,      61))),
                            '\t'.join(map(str,('HE577330.1_2',    259599,  8722919485,      60,      61))),
                            '\t'.join(map(str,('HE577330.1_3',    259599,  8723183488,      60,      61))),
                            '\t'.join(map(str,('HE577330.1_4',    259599,  8723447491,      60,      61))),
                            '\t'.join(map(str,('HE577330.1_5',    259600,  8723711494,      60,      61))),
                            '\t'.join(map(str,('HE577330.1_6',    259599,  8723975498,      60,      61)))]
    
                self.testfai = "test.fai"
                open(self.testfai,'w').write('\n'.join(indexes))
                self.queries = [
                ('HE577328.1_5','transport.fa.cluster10.fa',3.4e-172,15,215,200000,200300,'Azospirillum brasilense Sp245 plasmid AZOBR_p1 complete genome'),
                ('HE577328.1_5','transport.fa.cluster2.fa' ,5.6e-14,462,545,200500,200600,'Azospirillum brasilense Sp245 plasmid AZOBR_p1 complete genome'),
                ('HE577328.1_5','transport.fa.cluster10.fa',3.4e-172,15,215,201000,201400,'Azospirillum brasilense Sp245 plasmid AZOBR_p1 complete genome'),
                ('HE577328.1_5','transport.fa.cluster2.fa' ,5.6e-14,462,545,201500,201800,'Azospirillum brasilense Sp245 plasmid AZOBR_p1 complete genome'),
                ('HE577328.1_5','toxin.fa.cluster105.fa'   ,8.2e-16,37, 120,202000,202400,'Azospirillum brasilense Sp245 plasmid AZOBR_p1 complete genome'),
                ('HE577328.1_5','toxin.fa.cluster190.fa'   ,1.2e-14,3,   96,312000,312400,'Azospirillum brasilense Sp245 plasmid AZOBR_p1 complete genome'),
                ('HE577328.1_4','transport.fa.cluster11.fa',1e-154, 6,  205,312600,312800,'Azospirillum brasilense Sp245 plasmid AZOBR_p1 complete genome'),
                ('HE577328.1_5','toxin.fa.cluster195.fa'   ,1.8e-10,3,   86,321000,321000,'Azospirillum brasilense Sp245 plasmid AZOBR_p1 complete genome'),
                ('HE577328.1_4','transport.fa.cluster2.fa' ,4e-36,462,  537,325200,325400,'Azospirillum brasilense Sp245 plasmid AZOBR_p1 complete genome')]
                self.faidx = fasta.Indexer("",self.testfai)           
                self.faidx.load()
                self.maxDiff=10000
            def test1(self):
                cfilter = CliqueFilter(self.testfai,radius=10000)
                cfilter.createGraph(self.queries,backtrans=False)
                clusters = cfilter.filter(keyfunctions=["toxin","transport"])
                self.assertEquals(len(clusters),3)
                print '\n'.join(map(str,clusters[0]))
                print '\n'
                print '\n'.join(map(str,clusters[1]))
                print '\n'
                print '\n'.join(map(str,clusters[2]))
                print '\n'

                clusters = cfilter.merge(clusters)
                print '\n'.join(map(str,clusters[1]))
                self.assertEquals(set(clusters[0]),
                    set([
                    '|'.join(map(str,('HE577328.1_5','transport.fa.cluster10.fa',3.4e-172,15,215,200000,200300,'Azospirillum brasilense Sp245 plasmid AZOBR_p1 complete genome'))),
                    '|'.join(map(str,('HE577328.1_5','transport.fa.cluster2.fa' ,5.6e-14,462,545,200500,200600,'Azospirillum brasilense Sp245 plasmid AZOBR_p1 complete genome'))),
                    '|'.join(map(str,('HE577328.1_5','transport.fa.cluster10.fa',3.4e-172,15,215,201000,201400,'Azospirillum brasilense Sp245 plasmid AZOBR_p1 complete genome'))),
                    '|'.join(map(str,('HE577328.1_5','transport.fa.cluster2.fa' ,5.6e-14,462,545,201500,201800,'Azospirillum brasilense Sp245 plasmid AZOBR_p1 complete genome'))),
                    '|'.join(map(str,('HE577328.1_5','toxin.fa.cluster105.fa'   ,8.2e-16,37, 120,202000,202400,'Azospirillum brasilense Sp245 plasmid AZOBR_p1 complete genome')))]
                    ))
                self.assertEquals(set(clusters[1]),
                    set([
                    '|'.join(map(str,(('HE577328.1_5','toxin.fa.cluster190.fa'   ,1.2e-14,3,   96,312000,312400,'Azospirillum brasilense Sp245 plasmid AZOBR_p1 complete genome')))),
                    '|'.join(map(str,(('HE577328.1_4','transport.fa.cluster11.fa',1e-154, 6,  205,312600,312800,'Azospirillum brasilense Sp245 plasmid AZOBR_p1 complete genome')))),
                    '|'.join(map(str,(('HE577328.1_5','toxin.fa.cluster195.fa'   ,1.8e-10,3,   86,321000,321000,'Azospirillum brasilense Sp245 plasmid AZOBR_p1 complete genome')))),
                    '|'.join(map(str,(('HE577328.1_4','transport.fa.cluster2.fa' ,4e-36,462,  537,325200,325400,'Azospirillum brasilense Sp245 plasmid AZOBR_p1 complete genome'))))]
                    ))

        class TestCase1(unittest.TestCase):
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
                 correct = [  ('CP002279.1_3','toxin.fa.cluster2.fa',0,0,100,
                                self.faidx.sixframe_to_nucleotide('CP002279.1_3',25000)[0],
                                self.faidx.sixframe_to_nucleotide('CP002279.1_3',25100)[0],
                                'Mesorhizobium opportunistum WSM2075, complete genome'),
                               ('CP002279.1_3','transport.fa.cluster2.fa',0,0,100,
                                self.faidx.sixframe_to_nucleotide('CP002279.1_3',25200)[0],
                                self.faidx.sixframe_to_nucleotide('CP002279.1_3',25300)[0],
                                'Mesorhizobium opportunistum WSM2075, complete genome'),
                               ('CP002279.1_3','modifier.fa.cluster2.fa',0,0,100,
                                self.faidx.sixframe_to_nucleotide('CP002279.1_3',25400)[0],
                                self.faidx.sixframe_to_nucleotide('CP002279.1_3',25500)[0],
                                'Mesorhizobium opportunistum WSM2075, complete genome'),
                               ('CP002279.1_3','regulator.fa.cluster2.fa',0,0,100,
                                self.faidx.sixframe_to_nucleotide('CP002279.1_3',25600)[0],
                                self.faidx.sixframe_to_nucleotide('CP002279.1_3',25700)[0],
                                'Mesorhizobium opportunistum WSM2075, complete genome'),
                               ('CP002279.1_3','immunity.fa.cluster2.fa',0,0,100,
                                self.faidx.sixframe_to_nucleotide('CP002279.1_3',25800)[0],
                                self.faidx.sixframe_to_nucleotide('CP002279.1_3',26900)[0],
                                'Mesorhizobium opportunistum WSM2075, complete genome')
                            ]
                 for cluster in clusters:
                    
                     self.assertTrue(len(cluster)<=5)
                     if len(cluster)==5:
                         answer = ["|".join(map(str,s)) for s in correct]
                         self.assertItemsEqual(answer,cluster)
             def test2(self):
                 #all_hits=sorted(self.queries,key=lambda x: x[6])        
                 #all_hits=sorted(all_hits,key=lambda x: x[5])
                 #Sort by genome name
                 #all_hits=sorted(all_hits,key=lambda x: x[-1])  
                 reduced = interval_filter.overlaps(self.queries,self.testfai)
                
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

                
        class TestCase2(unittest.TestCase):
            def setUp(self):
                indexes = [
                            '\t'.join(map(str,('HE577328.1_4',    588676,  8720859786,      60,      61))),
                            '\t'.join(map(str,('HE577328.1_5',    588676,  8721458351,      60,      61))),
                            '\t'.join(map(str,('HE577328.1_6',    588676,  8722056916,      60,      61))),
                            '\t'.join(map(str,('HE577330.1_1',    259600,  8722655481,      60,      61))),
                            '\t'.join(map(str,('HE577330.1_2',    259599,  8722919485,      60,      61))),
                            '\t'.join(map(str,('HE577330.1_3',    259599,  8723183488,      60,      61))),
                            '\t'.join(map(str,('HE577330.1_4',    259599,  8723447491,      60,      61))),
                            '\t'.join(map(str,('HE577330.1_5',    259600,  8723711494,      60,      61))),
                            '\t'.join(map(str,('HE577330.1_6',    259599,  8723975498,      60,      61)))]
    
                self.testfai = "test.fai"
                open(self.testfai,'w').write('\n'.join(indexes))
                self.queries = [
                ('HE577328.1_4','transport.fa.cluster12.fa',1e-181,23,200,347781,347181,'Azospirillum brasilense Sp245 plasmid AZOBR_p1 complete genome    '),
                ('HE577328.1_4','transport.fa.cluster14.fa',3.4e-76,20,78,347451,347175,'Azospirillum brasilense Sp245 plasmid AZOBR_p1 complete genome    '),
                ('HE577328.1_4','transport.fa.cluster11.fa',1e-154,85,192,347604,347178,'Azospirillum brasilense Sp245 plasmid AZOBR_p1 complete genome    '),
                ('HE577328.1_4','transport.fa.cluster4.fa' ,8.9e-53,35,178,347769,347223,'Azospirillum brasilense Sp245 plasmid AZOBR_p1 complete genome    '),
                ('HE577328.1_4','transport.fa.cluster10.fa',3.2e-181,16,205,347799,347148,'Azospirillum brasilense Sp245 plasmid AZOBR_p1 complete genome  '),
                ('HE577328.1_4','transport.fa.cluster13.fa',5.6e-47,7,169,347802,347211,'Azospirillum brasilense Sp245 plasmid AZOBR_p1 complete genome    '),
                ('HE577328.1_4','transport.fa.cluster0.fa' ,1.6e-56,485,668,347814,347130,'Azospirillum brasilense Sp245 plasmid AZOBR_p1 complete genome   '),
                ('HE577328.1_4','immunity.fa.cluster2.fa'  ,1.1e-111,20,195,347790,347172,'Azospirillum brasilense Sp245 plasmid AZOBR_p1 complete genome    '),
                ('HE577328.1_5','toxin.fa.cluster9.fa'     ,8.7e-09,26,117,226738,226396,'Azospirillum brasilense Sp245 plasmid AZOBR_p1 complete genome  	 '),
                ('HE577328.1_4','transport.fa.cluster10.fa',3.2e-181,15,215,287472,286815,'Azospirillum brasilense Sp245 plasmid AZOBR_p1 complete genome  '),
                ('HE577328.1_6','immunity.fa.cluster2.fa'  ,1.8e-62,105,198,257786,257279,'Azospirillum brasilense Sp245 plasmid AZOBR_p1 complete genome    '),
                ('HE577328.1_5','toxin.fa.cluster207.fa'   ,1.7e-12,3,87,226681,226423,'Azospirillum brasilense Sp245 plasmid AZOBR_p1 complete genome  	 '),
                ('HE577328.1_5','toxin.fa.cluster196.fa'   ,6.2e-13,3,87,226684,226426,'Azospirillum brasilense Sp245 plasmid AZOBR_p1 complete genome  	    '),
                ('HE577328.1_5','toxin.fa.cluster187.fa'   ,4.6e-10,3,84,226681,226423,'Azospirillum brasilense Sp245 plasmid AZOBR_p1 complete genome  	    '),
                ('HE577328.1_5','toxin.fa.cluster211.fa'   ,1.4e-11,3,85,226681,226429,'Azospirillum brasilense Sp245 plasmid AZOBR_p1 complete genome  	    '),
                ('HE577328.1_6','transport.fa.cluster12.fa',6.5e-58,79,217,257777,257273,'Azospirillum brasilense Sp245 plasmid AZOBR_p1 complete genome   '),
                ('HE577328.1_6','transport.fa.cluster13.fa',7.3e-12,24,57,257144,256991,'Azospirillum brasilense Sp245 plasmid AZOBR_p1 complete genome    '),
                ('HE577328.1_4','transport.fa.cluster13.fa',5.6e-47,18,189,287454,286866,'Azospirillum brasilense Sp245 plasmid AZOBR_p1 complete genome   '),
                ('HE577328.1_5','transport.fa.cluster12.fa',2.2e-138,98,198,256909,256534,'Azospirillum brasilense Sp245 plasmid AZOBR_p1 complete genome  '),
                ('HE577328.1_6','modifier.fa.cluster16.fa' ,2.4e-32,3,170,266213,265649,'Azospirillum brasilense Sp245 plasmid AZOBR_p1 complete genome     '),
                ('HE577328.1_5','toxin.fa.cluster186.fa'   ,1.4e-21,1,98,226726,226426,'Azospirillum brasilense Sp245 plasmid AZOBR_p1 complete genome  	    '),
                ('HE577328.1_6','regulator.fa.cluster2.fa' ,0.0002,38,102,303671,303326,'Azospirillum brasilense Sp245 plasmid AZOBR_p1 complete genome     '),
                ('HE577328.1_4','immunity.fa.cluster2.fa'  ,1.1e-111,8,198,287478,286842,'Azospirillum brasilense Sp245 plasmid AZOBR_p1 complete genome     '),
                ('HE577328.1_6','transport.fa.cluster10.fa',3.2e-75,46,214,257813,257255,'Azospirillum brasilense Sp245 plasmid AZOBR_p1 complete genome   '),
                ('HE577328.1_4','transport.fa.cluster12.fa',1e-181,17,211,287460,286836,'Azospirillum brasilense Sp245 plasmid AZOBR_p1 complete genome    '),
                ('HE577328.1_5','transport.fa.cluster13.fa',5.6e-44,16,191,234625,234019,'Azospirillum brasilense Sp245 plasmid AZOBR_p1 complete genome   '),
                ('HE577328.1_5','toxin.fa.cluster109.fa'   ,9.2e-18,47,141,226774,226378,'Azospirillum brasilense Sp245 plasmid AZOBR_p1 complete genome	    '),
                ('HE577328.1_5','transport.fa.cluster14.fa',3.7e-62,20,104,234265,233956,'Azospirillum brasilense Sp245 plasmid AZOBR_p1 complete genome   '),
                ('HE577328.1_5','transport.fa.cluster10.fa',3.4e-172,73,213,256975,256504,'Azospirillum brasilense Sp245 plasmid AZOBR_p1 complete genome  '),
                ('HE577328.1_6','transport.fa.cluster11.fa',1.2e-40,70,204,257726,257279,'Azospirillum brasilense Sp245 plasmid AZOBR_p1 complete genome   '),
                ('HE577328.1_4','transport.fa.cluster4.fa' ,8.9e-53,133,207,287085,286809,'Azospirillum brasilense Sp245 plasmid AZOBR_p1 complete genome   '),
                ('HE577328.1_5','immunity.fa.cluster2.fa'  ,7.1e-120,9,199,234622,233995,'Azospirillum brasilense Sp245 plasmid AZOBR_p1 complete genome     '),
                ('HE577328.1_5','toxin.fa.cluster111.fa'   ,9.2e-15,43,138,226747,226384,'Azospirillum brasilense Sp245 plasmid AZOBR_p1 complete genome	    '),
                ('HE577328.1_5','toxin.fa.cluster203.fa'   ,4.6e-10,3,86,226681,226426,'Azospirillum brasilense Sp245 plasmid AZOBR_p1 complete genome  	    '),
                ('HE577328.1_4','transport.fa.cluster10.fa',3.2e-181,12,44,257907,257757,'Azospirillum brasilense Sp245 plasmid AZOBR_p1 complete genome   '),
                ('HE577328.1_4','transport.fa.cluster14.fa',3.4e-76,20,98,287109,286809,'Azospirillum brasilense Sp245 plasmid AZOBR_p1 complete genome    '),
                ('HE577328.1_5','immunity.fa.cluster2.fa'  ,7.1e-120,110,197,256966,256531,'Azospirillum brasilense Sp245 plasmid AZOBR_p1 complete genome   '),
                ('HE577328.1_5','transport.fa.cluster10.fa',3.4e-172,15,215,234622,233968,'Azospirillum brasilense Sp245 plasmid AZOBR_p1 complete genome  '),
                ('HE577328.1_5','transport.fa.cluster2.fa' ,5.6e-14,462,545,234247,233923,'Azospirillum brasilense Sp245 plasmid AZOBR_p1 complete genome   '),
                ('HE577328.1_5','toxin.fa.cluster105.fa'   ,8.2e-16,37,120,226753,226423,'Azospirillum brasilense Sp245 plasmid AZOBR_p1 complete genome	    '),
                ('HE577328.1_5','toxin.fa.cluster190.fa'   ,1.2e-14,3,96,226681,226399,'Azospirillum brasilense Sp245 plasmid AZOBR_p1 complete genome  	    '),
                ('HE577328.1_4','transport.fa.cluster11.fa',1e-154,6,205,287445,286842,'Azospirillum brasilense Sp245 plasmid AZOBR_p1 complete genome     '),
                ('HE577328.1_5','toxin.fa.cluster195.fa'   ,1.8e-10,3,86,226681,226408,'Azospirillum brasilense Sp245 plasmid AZOBR_p1 complete genome  	    '),
                ('HE577328.1_4','transport.fa.cluster2.fa' ,4e-36,462,537,287091,286809,'Azospirillum brasilense Sp245 plasmid AZOBR_p1 complete genome     '),
                ('HE577328.1_5','transport.fa.cluster12.fa',2.2e-138,17,217,234610,233989,'Azospirillum brasilense Sp245 plasmid AZOBR_p1 complete genome  '),
                ('HE577328.1_5','transport.fa.cluster0.fa' ,1.9e-18,605,688,234580,233950,'Azospirillum brasilense Sp245 plasmid AZOBR_p1 complete genome   '),
                ('HE577328.1_5','transport.fa.cluster11.fa',5.1e-103,5,204,234595,233995,'Azospirillum brasilense Sp245 plasmid AZOBR_p1 complete genome   '),
                ('HE577328.1_4','transport.fa.cluster3.fa' ,9.8e-31,408,493,287100,286803,'Azospirillum brasilense Sp245 plasmid AZOBR_p1 complete genome   '),
                ('HE577328.1_5','transport.fa.cluster3.fa' ,6.6e-12,408,440,234241,234106,'Azospirillum brasilense Sp245 plasmid AZOBR_p1 complete genome   '),
                ('HE577328.1_5','regulator.fa.cluster2.fa' ,8e-48,2,222,302848,302155,'Azospirillum brasilense Sp245 plasmid AZOBR_p1 complete genome 	    '),
                ('HE577328.1_5','toxin.fa.cluster209.fa'   ,7.4e-14,3,85,226681,226426,'Azospirillum brasilense Sp245 plasmid AZOBR_p1 complete genome  	    '),
                ('HE577328.1_4','transport.fa.cluster0.fa' ,1.6e-56,487,683,287439,286791,'Azospirillum brasilense Sp245 plasmid AZOBR_p1 complete genome   '),
                ('HE577328.1_5','toxin.fa.cluster205.fa'   ,1.2e-13,3,87,226681,226432,'Azospirillum brasilense Sp245 plasmid AZOBR_p1 complete genome  	    '),
                ('HE577328.1_5','transport.fa.cluster4.fa' ,3.6e-106,133,208,234241,233953,'Azospirillum brasilense Sp245 plasmid AZOBR_p1 complete genome  '),
                ('HE577328.1_6','modifier.fa.cluster15.fa' ,5.1e-49,2,190,266210,265637,'Azospirillum brasilense Sp245 plasmid AZOBR_p1 complete genome     '),
                ('HE577328.1_5','toxin.fa.cluster107.fa'   ,3.5e-16,5,87,226690,226426,'Azospirillum brasilense Sp245 plasmid AZOBR_p1 complete genome        ')]
                self.faidx = fasta.Indexer("",self.testfai)           
                self.faidx.load()
                self.maxDiff=10000
            def test(self):
                cfilter = CliqueFilter(self.testfai,radius=50000)
                cfilter.createGraph(self.queries,backtrans=False)
                clusters = cfilter.filter()
                self.assertTrue(len(clusters)==0)
                for cluster in clusters:
                    for i in xrange(len(cluster)):
                        for j in xrange(i):
                            igene,jgene = cluster[i],cluster[j]
                            iacc,iclrname,ifull_evalue,ihmm_st,ihmm_end,inenv_st,inenv_end,idescription = igene.split('|')
                            jacc,jclrname,jfull_evalue,jhmm_st,jhmm_end,jnenv_st,jnenv_end,jdescrjptjon = jgene.split('|')
                            inenv_st,inenv_end = map(int,[inenv_st,inenv_end])
                            jnenv_st,jnenv_end = map(int,[jnenv_st,jnenv_end])
                            midi = (inenv_st+inenv_end)/2
                            midj = (jnenv_st+jnenv_end)/2
                            self.assertLessEqual(abs(midi-midj), 50000,"midi: %d midj: %d"%(midi,midj))
                print "Clusters",len(clusters)
                pass
            def test2(self):
                cfilter = CliqueFilter(self.testfai,radius=100000)
                cfilter.createGraph(self.queries,backtrans=False)
                clusters = cfilter.filter()
               
                for cluster in clusters:
                    for i in xrange(len(cluster)):
                        for j in xrange(i):
                            igene,jgene = cluster[i],cluster[j]
                            iacc,iclrname,ifull_evalue,ihmm_st,ihmm_end,inenv_st,inenv_end,idescription = igene.split('|')
                            jacc,jclrname,jfull_evalue,jhmm_st,jhmm_end,jnenv_st,jnenv_end,jdescrjptjon = jgene.split('|')
                            inenv_st,inenv_end = map(int,[inenv_st,inenv_end])
                            jnenv_st,jnenv_end = map(int,[jnenv_st,jnenv_end])
                            assert type(inenv_st)==type(0)
                            assert type(inenv_end)==type(0)

                            midi = (inenv_st+inenv_end)/2
                            midj = (jnenv_st+jnenv_end)/2
                            self.assertLessEqual(abs(midi-midj), 100000,"midi: %d midj: %d"%(midi,midj))
                print "Clusters",len(clusters)
                pass
                
        unittest.main()















