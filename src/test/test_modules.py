"""
A set of test functions for unittests
"""
import os, site

def createFasta(fname,refid,refseq):
    fastaH = open(fname,'w')
    fastaH.write(">%s\n%s\n"%(refid,refseq))
    fastaH.close()
    
def createFastaIdx(fname,refid,refseq):
    fastaH = open(fname,'w')
    fastaIdx = open(fname+".fai",'w')
    fastaH.write(">%s\n%s\n"%(refid,refseq))
    fastaIdx.write("%s\t%d\t%d\t%d\t%d\n"%(refid,len(refseq),len(refid)+2,len(refseq),len(refseq)+1))
    fastaH.close()
    fastaIdx.close()
    
""" Get all fna files within all a directory and all of its subdirs"""
def getFNA(root_dir):
    fnafiles = []
    for root, subFolders, files in os.walk(root_dir):
        for fname in files:
            genome_files = []
            organism,ext = os.path.splitext(os.path.basename(fname))
            if ext==".fna":
                absfile=os.path.join(root,fname)
                fnafiles.append(absfile)
    return fnafiles
