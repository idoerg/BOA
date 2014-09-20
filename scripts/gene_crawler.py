"""
A simple genbank crawler that identifies all genbank records that 
identify bacteriocins, lanbiotics or ABC transporters
"""

import os, sys, site
import argparse
from Bio import SeqIO
base_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
site.addsitedir(os.path.join(base_path, 'src'))

from  acc2species import AccessionToSpecies

def checkWords(note,words):
    for word in words:
        if word in note:
            return True
    return False

def filterGeneAnnotations(accID,genbank_file,out,words,accTable):
    seqType,speciesName = accTable.lookUp(accID)
    out.write("%s\n"%speciesName)
    seqs = []
    try:
        seq_record = SeqIO.parse(open(genbank_file,'r'), "genbank").next()
        for feature in seq_record.features:
            try: 
                note = feature.qualifiers["note"][0]
                misc = feature.qualifiers["note"][0]
                if checkWords(note,words):
                    out.write("Note: %s\n"%str(note))
                    out.write("Position: %s\n"%feature.location)                
                    continue
            except KeyError as k:
                continue
    except Exception as e:
        print "Parsing error ? ",accID,e
    out.write('\n')
    
def go(root,out,words,accTable):
    outHandle = open(out,'w')
    for root, subFolders, files in os.walk(root):
        for fname in files:
            genome_files = []
            organism,ext = os.path.splitext(os.path.basename(fname))
            if ext==".gbk":
                absfile=os.path.join(root,fname)
                filterGeneAnnotations(organism,absfile,outHandle,words,accTable)
                
if __name__=="__main__":
    parser = argparse.ArgumentParser(description=\
        'Finds intergenic regions from genback files')
    parser.add_argument(\
        '--test', action='store_const', const=True, default=False,
        help='Run the unittests')
    parser.add_argument(\
        '--accession-table', type=str, required=False,
        help='A table that maps accession ids to species')
    parser.add_argument(\
        '--root-dir', type=str,required=False,default="",
        help='Root directory of all of the files of interest')
    parser.add_argument(\
        '--output-file', type=str, required=False,
        help='The output file containing the tab-delimited output')
    args = parser.parse_args()
    words = ['bacteriocin','lantibiotic']
    if not args.test:
        accTable = AccessionToSpecies(args.accession_table)
        go(args.root_dir,args.output_file,words,accTable)
    
