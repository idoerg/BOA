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
    print links[0]
    return links

"""
Retrieves file from html site
"""
def wgetRetrieve(url,outdir):
    fname = url.split("/")[-1]
    print "fname",fname
    cmd = "wget %s; mv %s %s"%(url,fname,outdir)
    down_proc = subprocess.Popen(cmd,shell=True)
    down_proc.wait()
    if fname=="":
        return None
    return fname

def isDir(string):
    return string[-1]=="/"

"""
Recurses through online directory
"""
def recursiveWget(url,outfolder):
    fname = wgetRetrieve(url,outfolder)
    if os.path.exists("index.html"):        
        indexFile="%sindex.html"%(outfolder)
        os.remove("index.html")
    else:
        indexFile="%s%s"%(outfolder,fname)

    print "index file",indexFile
    links = parseHTML(indexFile)
    
    for link in links:
        newURL="%s%s"%(url,link)
        out = "%s%s"%(outfolder,link)
        try:
            if isDir(link):
                print link,out,outfolder,newURL
                os.mkdir(out)
                recursiveWget(newURL,out)
            else:
                wgetRetrieve(newURL,out)
        except OSError:
            continue
            
        
if __name__=="__main__":
    recursiveWget(args.URL,args.outfolder)
