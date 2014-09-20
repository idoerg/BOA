import os, sys, site
import argparse
import re
import Bio
from Bio import SeqIO, SeqFeature
from Bio.SeqRecord import SeqRecord

def parse(organism,genbank_file,out):
    seq_set = set()
    try:
        seq_record = SeqIO.read(open(genbank_file), "genbank")
        gi = seq_record.annotations['gi']
        acc = seq_record.id
        out.write(">%s\t%s\n"%(gi,acc)
    except Exception as e:
        print "Exception at",organism,genbank_file,e


def go(rootdir,out):
    outHandle = open(out,'w')
    for root, subFolders, files in os.walk(rootdir):
        for fname in files:
            genome_files = []
            organism,ext = os.path.splitext(os.path.basename(fname))
            if ext==".gbk":
                absfile=os.path.join(root,fname)
                parse(organism,absfile,outHandle)

if __name__=="__main__":
    parser = argparse.ArgumentParser(description=\
        'Creates accession gi database for bacterial genomes')
    parser.add_argument(\
        '--root-dir', type=str,required=False,default="",
        help='Root directory of all of the files of interest')
    parser.add_argument(\
        '--output-file', type=str, required=False,
        help='The output file containing the fasta output')
    args = parser.parse_args()
    root,out = args.root_dir,args.output_file
    go(root,out)

