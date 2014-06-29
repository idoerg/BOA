"""
Commands to obtain preprocessed data
cat classify | tr '|' '\t' | sort -k17 | sort -k16 | sort -k9 > classified.txt

1  bacID                bacteriocin ID  		   
2  bac_org  		    organism    			   
3  bac_st   		    bacteriocin start   	   
4  bac_end  		    bacteriocin end 		   
5  bac_strand   	    bacteriocin strand  	   
6  annot_org		    annotated gene organism    
7  annot_locus  	    annotated gene locus       
8  annot_st 		    annotated gene start       
9  annot_end		    annotated gene end  	   
10  annot_strand	    annotated gene strand      
11 context_org  	    context gene organism      
12 context_id   	    context gene protein id    
13 context_label	    context gene class label   
"""

from collections import defaultdict
import os
cgenes = []
root = os.environ['BACFINDER_HOME']
fname = "%s/workspace/quorum_data/sag.txt"%root
out = "%s/workspace/quorum_data/predicted_sag.txt"%root
outhandle=open(out,'w')
with open(fname,'r') as handle:
    for ln in handle:
        ln = ln.rstrip()
        toks = ln.split('\t')
        
        bacID,bac_org,bac_st,bac_end,bac_strand=toks[0],toks[4],int(toks[5]),int(toks[6]),toks[7]
        annot_org,annot_locus,annot_st,annot_end,annot_strand=toks[8],toks[9],int(toks[11]),int(toks[12]),toks[13]
        context_org,context_id,context_label = toks[14],toks[15],toks[16]
        cgenes.append( (bacID,bac_org,bac_st,bac_end,bac_strand,
                        annot_org,annot_locus,annot_st,annot_end,annot_strand,
                        context_org,context_id,context_label) )
        
""" Group all entries by annotated gene organism """
annot_dict = defaultdict(list)
for entry in cgenes:
    annot_org = entry[5]
    annot_dict[annot_org].append(entry)

""" Drop all duplicate context genes """
for key in annot_dict.keys():
    block = annot_dict[key]
    newblock = []
    context_ids = set()
    
    for b in block:
        cid = b[11]
        if cid not in context_ids:
            context_ids.add(cid)
            newblock.append(b)
    annot_dict[key] = newblock

""" 
Drop all organisms that only have one context gene
and only transport genes
"""
new_annot_dict = defaultdict()
for key in annot_dict.keys():
    block = annot_dict[key]
    labels = zip(*block)[12]
    if len(block)>1 and len(set(labels))>1:
        new_annot_dict[key] = block
        print >>outhandle,"\n".join(["\t".join(map(str,b)) for b in block])
        print >>outhandle
    
    
    
    
    
    
    
    