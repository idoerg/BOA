import argparse

import os
import sys
import site
import re
import shutil
base_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

site.addsitedir(os.path.join(base_path, 'src'))
base_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
base_path="%s/src"%base_path
for directory_name in os.listdir(base_path):
    site.addsitedir(os.path.join(base_path, directory_name))

from accessionMap import AccessionGI
acc_reg = re.compile("accession=(\S+)|")
from Bio import SeqIO

def accession_set(accFile):
    accs = set()
    with open(accFile) as handle:
        for ln in handle:
            ln = ln.rstrip()
            tok = ln.split("\t")
            acc = tok[0]
            acc = acc.split(".")[0]
            accs.add(acc)
    return accs

def go(rootdir,acc_set,outdir):
    
    for root, subFolders, files in os.walk(rootdir):
        for fname in files:
            genome_files = []
            acc,ext = os.path.splitext(os.path.basename(fname))
            absfile=os.path.join(root,fname)
            directory = os.path.dirname(absfile)
            taxa = os.path.basename(directory)
            if ext==".faa":
                print acc
                if acc in acc_set:
                    shutil.copy(absfile,outdir)
                    

if __name__=="__main__":

    parser = argparse.ArgumentParser(description=\
        'Finds bacteriocins and context genes')
    parser.add_argument(\
        '--root-dir',type=str, required=True,
        help='Root directory')
    parser.add_argument(\
        '--accession-file',type=str, required=True,
        help='Root directory')
    parser.add_argument(\
        '--output-dir',type=str, required=True,
        help='Output directory')
    args = parser.parse_args()
    acc_set = accession_set(args.accession_file)
    go(args.root_dir,acc_set,args.output_dir)
