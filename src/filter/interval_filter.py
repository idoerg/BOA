A:wq

A
A
A
A
A
A
A
A
A
A
A
A
A
A
A
A
A
A
A
A
A
A
A
A
A
A
A
A
A
A
A
A
A
A

    faidx = fasta.Indexer('',fasta_index)
    faidx.load()
    prevOrg,curOrg = None,None
    prevStrand,curStrand = None,None
    newHits = []
    #print "Before",len(hits)
    for hit in hits:
        acc,clrname,full_evalue,hmm_st,hmm_end,env_st,env_end,description=hit
        curOrg = fasta.getName(acc)
        if backtrans:
            hitSt,curStrand  = faidx.sixframe_to_nucleotide(acc,env_st)
<<<<<<< HEAD
            hitEnd,curStrand = faidx.sixframe_to_nucleotide(acc,env_end) 
        else:
            hitSt,hitEnd = env_st,env_end
            
=======
            hitEnd,curStrand = faidx.sixframe_to_nucleotide(acc,env_end)
        else:
            hitSt,hitEnd = map(int,[env_st,env_end])
>>>>>>> 85866840a7fc60c97a830a097b05cb14cdedd925
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
    #print "After",len(newHits)
    
    return newHits


