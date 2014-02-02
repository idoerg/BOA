#!/usr/bin/env python
import sys
import Bio
from Bio import SeqIO, SeqFeature
from Bio.SeqRecord import SeqRecord
import os

import site
import argparse
import string
import numpy
import re

base_path = os.path.dirname(os.path.abspath(__file__))
for directory_name in ['test']:
    site.addsitedir(os.path.join(base_path, directory_name))
import test_genbank

# Copyright(C) 2009 Iddo Friedberg & Ian MC Fleming
# Released under Biopython license. http://www.biopython.org/DIST/LICENSE
# Do not remove this comment
def get_interregions(genbank_path,output_file,intergene_length=1,write_flag = 'w'):
    try:
        seq_record = SeqIO.parse(open(genbank_path), "genbank").next()
        cds_list_plus = []
        cds_list_minus = []
        intergenic_records = []
        # Loop over the genome file, get the CDS features on each of the strands
        for feature in seq_record.features:
            if feature.type == 'CDS':
                mystart = feature.location._start.position
                myend = feature.location._end.position
                if feature.strand == -1:
                    cds_list_minus.append((mystart,myend,-1))
                elif feature.strand == 1:
                    cds_list_plus.append((mystart,myend,1))
                else:
                    sys.stderr.write("No strand indicated %d-%d. Assuming +\n" %
                                      (mystart, myend))
                    cds_list_plus.append((mystart,myend,1))

        for i,pospair in enumerate(cds_list_plus[1:]):
            # Compare current start position to previous end position
            last_end = cds_list_plus[i][1]
            this_start = pospair[0]
            strand = pospair[2]
            if this_start - last_end >= intergene_length:
                intergene_seq = seq_record.seq[last_end:this_start]
                strand_string = "+"
                intergenic_records.append(
                      SeqRecord(intergene_seq,id="%s-ign-%d:%d-%d%s" % (seq_record.name,i,last_end+1,this_start,strand_string),
                      description="%s %d-%d %s" % (seq_record.name, last_end+1,
                                                            this_start,strand_string)))
        for i,pospair in enumerate(cds_list_minus[1:]):
            last_end = cds_list_minus[i][1]
            this_start = pospair[0]
            strand = pospair[2]
            if this_start - last_end >= intergene_length:
                intergene_seq = seq_record.seq[last_end:this_start]
                strand_string = "-"
                intergenic_records.append(
                      SeqRecord(intergene_seq,id="%s-ign-%d:%d-%d%s" % (seq_record.name,i,last_end+1,this_start,strand_string),
                      description="%s %d-%d %s" % (seq_record.name, last_end+1,
                                                            this_start,strand_string)))
        SeqIO.write(intergenic_records, open(output_file,write_flag), "fasta")
    except Exception as e:
        print "Error at",genbank_path,e

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=\
                                         'Finds intergenic regions from genback file')
    parser.add_argument(\
        '--genbank-path', type=str, required=False,
        help='The path of the genbank file')
    parser.add_argument(\
        '--intergene-length', type=int, default=1,required=False,
        help='The path of the genbank file')
    parser.add_argument(\
        '--output-file', type=str, default='.',required=False,
        help='The fasta file containing the intergenic regions')
    parser.add_argument(\
        '--test', action='store_const', const=True, default=False,
        help='Run the unittests')

    args = parser.parse_args()

    if not args.test:
        get_interregions(args.genbank_path,args.output_file,args.intergene_length)
    else:
        del sys.argv[1:]
        import unittest
        print "Testing ..."
        class TestIntegene(unittest.TestCase):
            def setUp(self):
                test_input = test_genbank.yeast
                self.test_file = "test.gbk"
                self.out_file = "out.fa"
                handle = open(self.test_file,'w')
                handle.write(test_input)
                handle.close()
            def tearDown(self):
                os.remove(self.test_file)

            def test1(self):
                get_interregions(self.test_file,self.out_file)
                reg = re.compile("([A-Z]+\d+)-ign-\d+:(\d+)-(\d+)(\S)")
                records = [r for r in SeqIO.parse(open(self.out_file,'r'),"fasta")]
                ids = [r.id for r in records]
                #Each record has (accession,start,end,strand)
                descriptions = [reg.findall(ids[i])[0] for i in range(len(ids))]
                accession,start,end,strand = descriptions[0]
                self.assertEquals(int(start),207)
                self.assertEquals(int(end),686)
    unittest.main()
