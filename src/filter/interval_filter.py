"""
Performs filtering methods

Todo:
1) Set intersection of species/translation sequences
2) Set intersection of 50 kb regions
3) Set intersections of all context genes present
4) Set intersection of key context genes present
5) Use e-values for tie breakers
"""

from bx.intervals import *
from collections import defaultdict
import fasta
import re
"""
Filters out all query sequences that don't fall in one of the intervals
Note: returns a list of all of the overlapping intervals
intervalDict: A dictionary of intervals indexed by (organismID,strand)
queries: A list of query objects with
        (start,end,orgid,strand,query_obj)
"""
def spatialFilter(queries,intervalDict,radius):
    filtered = []
    geneNeighborhoods = []
    for query in queries:
        header, query_obj = query
        start,end,orgid,strand = header
        #print "Org id",orgid,orgid not in intervalDict
        if orgid not in intervalDict: continue
        nearTargets = intervalDict[orgid].find( start,end )
        if len(nearTargets)>0:
            for geneobj in nearTargets:
                filtered.append(query_obj)
                geneNeighborhoods.append(geneobj)
    return filtered,geneNeighborhoods
               
"""
Filter out all annotations that are not within the radius of a bacteriocin
"""
#annot_reg = re.compile("([A-Z0-9_]+_[0-9]+).")
annot_reg = re.compile("([A-Z][A-Z0-9_]+).")
def annotatedGenes(annots,bacteriocins,radius):
    intervalDict = defaultdict( IntervalTree )
    for bact in bacteriocins:
        start,end,orgid,strand = bact.sbjct_start,bact.sbjct_end,bact.sbjct_id,bact.strand
        #print start,end,orgid,strand
        orgid = orgid.split('|')[3]
        orgid = orgid.split('.')[0]
        #orgid = annot_reg.findall(orgid)[0]
        
        stBound,endBound = start-radius,end+radius
        intervalDict[orgid].add( stBound,endBound,bact )                                                 
    headers = []
    for a in annots: 
        start,end,orgid,strand = a[0],a[1],a[2],a[3];
        headers.append( (start,end,orgid,strand) )
    queries = zip(headers,annots) 
    
    filtered,annotNeighborhoods = spatialFilter(queries,intervalDict,radius)
    return filtered,annotNeighborhoods

"""
Filters out bacteriocins not contained in a gene neighborhood
bacteriocins: list of bacteriocins blast hits
genes: list of target genes blast hits 
radius: radius around target gene
filtered: list of target genes blast hits that make it
          through the filtering criteria imposed by the radius
geneNeighborhoods: a list of intervals with a radius surrounding the
          target genes
header = (orgid,start,end,strand)
"""
def bacteriocins(bacteriocins,genes,radius):
    if radius<0: radius = 50000000 #largest bacterial genome
    intervalDict = defaultdict( IntervalTree )
    for gene in genes:    
        start,end,refid,orgid,strand = gene.sbjct_start,gene.sbjct_end,gene.query_id,gene.sbjct_id,gene.strand
        stBound,endBound = start-radius,end+radius
        intervalDict[orgid].add( stBound,endBound,gene )
    headers = []
    for b in bacteriocins: 
        headers.append( (b.sbjct_start,b.sbjct_end,b.sbjct_id,b.strand) )
    queries = zip(headers,bacteriocins) 
    filtered,geneNeighborhoods = spatialFilter(queries,intervalDict,radius)
    return filtered,geneNeighborhoods

"""
Collapses hmmer hits using Interval trees
"""
def overlaps(hits,fasta_index,backtrans=True):
    tree = IntervalTree()
    faidx = fasta.Indexer('',fasta_index)
    faidx.load()
    prevOrg,curOrg = None,None
    prevStrand,curStrand = None,None
    newHits = []
    print "Before",len(hits)
    for hit in hits:
        acc,clrname,full_evalue,hmm_st,hmm_end,env_st,env_end,description=hit
        curOrg = fasta.getName(acc)
        if backtrans:
            hitSt,curStrand  = faidx.sixframe_to_nucleotide(acc,env_st)
            hitEnd,curStrand = faidx.sixframe_to_nucleotide(acc,env_end) 
        if prevOrg == None:
            prevOrg = curOrg
            prevStrand = curStrand
            tree.add(hitSt,hitEnd,hit)
            newHits.append(hit)            
        elif prevOrg!=curOrg or prevStrand!=curStrand:
            tree = IntervalTree()
            tree.add(hitSt,hitEnd,hit)
            prevOrg = curOrg
            prevStrand = curStrand
            newHits.append(hit)
        else:
            overlaps = tree.find(hitSt,hitEnd)
            if len(overlaps)==0:
                tree.add(hitSt,hitEnd,hit)
                newHits.append(hit)
    print "After",len(newHits)
    
    return newHits


