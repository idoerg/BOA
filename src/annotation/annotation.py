"""
Base class for annotations
"""
import Bio
from Bio import SeqIO, SeqFeature
from Bio.SeqRecord import SeqRecord
import genbank_annotation
import gff_annotation

"""
A container to process fasta objects and obtain information for annotations
"""
class AnnotatedGenes(object):
    def __init__(self,infile):
        self.infile = infile
        self.iterator = SeqIO.parse(infile,"fasta")
    def __iter__(self):
        return self
    def next(self):
        try:
            record = next(self.iterator)
            toks = record.description.split('\t')
            #print toks
            assert len(toks)==7
            index,orgid,start,end,strand,locus,protid = toks
            start,end = int(start),int(end)
            sequence = str(record.seq)
            #Ignore the protein id. Locus tag is more general
            return start,end,orgid,strand,locus,protid,sequence
        except StopIteration as s:
            raise StopIteration()

def go(root_dir,output_file):
    outHandle = open(output_file,'w')
    for root, subFolders, files in os.walk(root_dir):
        for fname in files:
            genome_files = []
            organism,ext = os.path.splitext(os.path.basename(fname))
            if ext==".gbk":
                absfile=os.path.join(root,fname)
                genbank_annotations.parse(organism,absfile,outHandle)
                
if __name__=="__main__":
    parser = argparse.ArgumentParser(description=\
        'Finds intergenic regions from genback files')
    parser.add_argument(\
        '--root-dir', type=str,required=False,default="",
        help='Root directory of all of the files of interest')
    parser.add_argument(\
        '--output-file', type=str, required=False,
        help='The output file containing the tab-delimited output')
    args = parser.parse_args()
    go(args.root_dir,args.output_file)
