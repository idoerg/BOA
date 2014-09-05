
"""
Filters out all HMMER hits that have a score below a certain threshold
"""
def filter(hits,threshold):
    goodHits = []
    for hit in hits:
        seqid,seq = hit
        score = seqid[2]
        if score>threshold:
            goodHits.append(hit)
    return goodHits
    

