"""
Finds all 16S RNA in bacterial genomes and stores it in a database
"""
import os, sys, site
import argparse
import re
import Bio
from Bio import SeqIO, SeqFeature
from Bio.SeqRecord import SeqRecord

loc_reg = re.compile("(\d+):>?(\d+)\S\((\S)\)")

def parse16SRNA(organism,genbank_file,out):
    seq_set = set()
    try:
        seq_record = SeqIO.read(open(genbank_file), "genbank")  
        for feature in seq_record.features:
            #print feature
            if feature.type == 'rRNA':
                try: 
                    product = feature.qualifiers["product"][0]
                    if product.lower()=='16s ribosomal rna' and organism not in seq_set:
                        st,end,strand = loc_reg.findall(str(feature.location))[0]
                        st,end = int(st),int(end)
                        if strand=="+":
                            seq = seq_record.seq
                            seq = seq[st:end]
                            if len(seq)>0:
                                out.write(">%s\n%s\n"%(organism,str(seq)))
                                seq_set.add(organism)
                        else:
                            seq = seq_record.reverse_complement().seq
                            seq = seq[st:end]
                            print st,end
                            if len(seq)>0:
                                out.write(">%s\n%s\n"%(organism,str(seq)))
                                seq_set.add(organism)
                except KeyError as k:                    
                    continue
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
                parse16SRNA(organism,absfile,outHandle)
    pass
if __name__=="__main__":
    parser = argparse.ArgumentParser(description=\
        'Creates a 16SRNA database for bacterial genomes')
    parser.add_argument(\
        '--root-dir', type=str,required=False,default="",
        help='Root directory of all of the files of interest')
    parser.add_argument(\
        '--output-file', type=str, required=False,
        help='The output file containing the tab-delimited output')
    parser.add_argument(\
        '--test', action='store_const', const=True, default=False,
        help='Run the unittests')
    args = parser.parse_args()
    if not args.test:
        root,out = args.root_dir,args.output_file
        go(root,out)
    else:
        del sys.argv[1:]
        import unittest
        unittest.main()
