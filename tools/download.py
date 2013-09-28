import os
import argparse
import time
import sys
import subprocess
import re
import os.path

from BeautifulSoup import BeautifulSoup

parser = argparse.ArgumentParser(description='Automatically download all genomes from a URL')

parser.add_argument(\
    "-i", "--URL", dest="URL", metavar="FILE", required=True,
    help="web url of root directory containing all genomes to be downloaded")

parser.add_argument(\
    "-o", "--outfolder", dest="outfolder", metavar="FOLDER", required=False,default='.',
    help="Folder where results will be stored.")

args = parser.parse_args()

exts = {'TXT','DIR'} #file extensions


"""
Parses HTML file
"""
def parseHTML(fname):
    lines = open(fname,'r').readlines()
    tree = BeautifulSoup(''.join(lines))
    tags = tree('a',href=re.compile(r"^[a-zA-Z0-9_]+"))
    links = [str(tag['href']).lstrip() for tag in tags]
    return links

"""
Retrieves file from html site
"""
def wgetRetrieve(url,fname,outdir):
    print >> sys.stderr,"Out dir",outdir,"fname",fname
    cmd = "wget %s -O %s/%s; "%(url,outdir,fname)
    down_proc = subprocess.Popen(cmd,shell=True)
    down_proc.wait()

def isDir(string):
    return string[-1]=="/"

"""
Recurses through online directory
"""
def recursiveWget(url,fname,outfolder):
    wgetRetrieve(url,fname,outfolder)
    indexFile="%s/%s"%(outfolder,fname)

    links = parseHTML(indexFile)
    
    for link in links:
        newURL="%s%s"%(url,link)
        out = "%s/%s"%(outfolder,link)
        try:
            if isDir(link):
                os.mkdir(out)
                recursiveWget(newURL,"index.html",out)
            else:
                out = os.path.dirname(out)
                wgetRetrieve(newURL,link,out)
        except OSError:
            continue
            
        
if __name__=="__main__":
    fname = args.URL.split('/')[-2]
    
    recursiveWget(args.URL,fname,args.outfolder)
