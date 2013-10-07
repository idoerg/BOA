#!/usr/bin/python

from multiprocessing import Pool
import time
import os
import sys
import argparse
from homolog2 import *
from collections import defaultdict

# Copyright(C) 2013 David Ream
# Released under Biopython license. http://www.biopython.org/DIST/LICENSE
# Do not remove this comment


# This exists to  make the main function easier to read. It contains code to run the argument parser, and does nothing else.
def parser_code():

    parser = argparse.ArgumentParser(description="Parse the results of a BLAST search and organize the results by specific operons. The program will save the results in a folder designated by the user or the default './blast_parse/'.")

    parser.add_argument("-i", "--infolder", dest="infolder", default='./blast_result/', metavar="FOLDER",
                help="A folder that contains all BLAST results in tabular form that you are interested in.")
                 
    parser.add_argument("-o", "--outfolder", dest="outfolder", metavar="FOLDER", default='./blast_parse/',
                help="Folder where the parsed BLAST results will be stored. Default is the folder './blast_parse/'.")
    
    parser.add_argument("-f", "--filter", dest="filter", metavar="FILE", default='',
                help="File restrictiong which accession numbers this script will process. If no file is provided, filtering is not performed.")
          
    parser.add_argument("-n", "--num_proc", dest="num_proc", metavar="INT", default = os.sysconf("SC_NPROCESSORS_CONF"), type=int,
                help="Number of processors that you want this script to run on. The default is every CPU that the system has.")

    # This program will utterly ignore the genes, but will use the named operons as a filter for what operon(s) we investigate.
    parser.add_argument("-p", "--operon_file", dest="operon_file", metavar="FILE", default='',
                help="File which contains operon information for use in custom operon queries. The file format is opern_name followed by the constituent gene names, tab delineated.")
                             
    return parser.parse_args()
    

def check_options(parsed_args):

    if os.path.isdir(parsed_args.infolder):
        infolder = parsed_args.infolder
    else:
        print "The folder %s does not exist." % parsed_args.infolder
        sys.exit()
    
    # if the directory that the user specifies does not exist, then the program makes it for them. 
    if not os.path.isdir(parsed_args.outfolder):
        os.makedirs(parsed_args.outfolder)
    outfolder = parsed_args.outfolder
    
    # Check the filter file
    if parsed_args.filter == '' or os.path.exists(parsed_args.filter):
        filter_file = parsed_args.filter
    else:
        print "The file %s does not exist." % parsed_args.filter
        sys.exit()
    
    # section of code that deals determining the number of CPU cores that will be used by the program
    if parsed_args.num_proc > os.sysconf("SC_NPROCESSORS_CONF"):
        num_proc = os.sysconf("SC_NPROCESSORS_CONF")
    elif parsed_args.num_proc < 1:
        num_proc = 1
    else:
        num_proc = parsed_args.num_proc
            
    # Check the operon file
    if parsed_args.operon_file == '' or os.path.exists(parsed_args.operon_file):
        operon_file = parsed_args.operon_file
    else:
        print "The file %s does not exist." % parsed_args.operon_file
        sys.exit()
        
    return infolder, outfolder, operon_file, filter_file, num_proc


#this function will return all of the files that are in a directory. os.walk is recursive traversal.
def returnRecursiveDirFiles(root_dir):
    result = []
    for path, dir_name, flist in os.walk(root_dir):
        for f in flist:
            fname = os.path.join(path, f)
            if os.path.isfile(fname):
                result.append(fname)
    return result

def raw_parse(arg_tuple):
    infile, out_folder = arg_tuple
    outfile = out_folder + infile.split('/')[-1].replace('.', '_raw_data.')
    result = [i.strip() for i in open(infile).readlines() if i[0] != '#']
    handle = open(outfile, 'w')
    handle.write('\n'.join(result))
    handle.close()
    # note that in the event that we wish to remove the intermediate files here, htat we could simple return 'result'

def parallel_raw_parse(outfolder, infolder, num_proc):
    blast_result_list = returnRecursiveDirFiles(infolder)
    
    # This folder is not necessary later (after validation), however it is very useful to me to see the raw data, as parsed out
    # by each operon. So I am adding this location for raw data, and will take out if I need to.
    raw_parse_folder = outfolder + 'organism_raw_info/'
    if not os.path.isdir(raw_parse_folder):
        os.makedirs(raw_parse_folder)
    raw_data_tuple_list = [(i, raw_parse_folder) for i in blast_result_list]
    
    pool = Pool(processes = num_proc)
    pool.map(raw_parse, raw_data_tuple_list)


# The filter file here i think should be either a vacant value (such as '') or a user defined
# value.  I do not think by default it should be given, like I have made the default behavior.    
def parallel_blast_parse_dict(in_folder, out_folder, num_proc, filter_file, operon_dict):
    print len(operon_dict)
    operon_out_folder = out_folder
    if not os.path.isdir(operon_out_folder):
        os.makedirs(operon_out_folder)
    do_filter = not filter_file == ''
    print 'do_filter', do_filter
    nc_list = []
    if do_filter:
        nc_list = [i.strip() for i in open(filter_file).readlines()]
    
    result = {}
    for fname in returnRecursiveDirFiles(in_folder):
        for line in [i.strip() for i in open(fname).readlines()]:
            hlog = blast_hit_to_homolog(line)
            accession = hlog.accession()
            predicted_gene = hlog.predicted_gene()
            store_entry = True
            # quickly see if the gene is part of the set of operons that we are looking at
            try:
                operon = operon_dict[predicted_gene]
            except:
                operon = ''
                store_entry = False
                
            # If we are filtering the organism list, check if the organism is part of the set that we are interested in
            if do_filter:
                if accession not in nc_list: # 
                    stare_entry = False
            
            if store_entry: # check to see if the cirteria for retaining an entry has been met
                if operon in result.keys():    # check if the operon is in the result dict already, if not make a new entry in the else clause
                    if accession in result[operon].keys(): # Check if the organism has been added to the operon
                        result[operon][accession].append(hlog.ret_str())
                    else: # if the organims is not part of the operon, add it
                        result[operon].update({accession:[hlog.ret_str()]})
                else: # add the operon to the result
                    result.update({operon: {accession: [hlog.ret_str()]}})
    #print sorted(result.keys())
    print len(result['atpIBEFHAGDC'].keys())
    
    for op in result.keys():
        outfile = out_folder + op +'.txt'
        handle = open(outfile, 'w')
        for org in sorted(result[op].keys()):
            handle.write('\n'.join(result[op][org]) + '\n')
        handle.close()
                    

# This function is a parser for converting a BLAST hit into a Homolog class object. 
# To clarify what is going on here: BLAST hits run from the -m 9 tabular with comment lines.  I filter out the comment lines, but makes everything 
# more readable in my opinion.  
# Fields: Query id, Subject id, % identity, alignment length, mismatches, gap openings, q. start, q. end, s. start, s. end, e-value, bit score

def blast_hit_to_homolog(hit):
    field = hit.split('\t')
    query = field[0].split('|')
    source_accession = query[0]
    source_common = query[1]
    source_locus = query[2]
    predicted_gene = query[3]
    source_start = query[4]
    source_start = query[5]
    source_strand = query[6]
    product_type = query[7]
    synonyms = query[8]
    synonym_list = synonyms.split(':')
    source_gc = query[9]
    try:
        #accession, locus, gene, start, stop, strand, organism, gc, hgt = field[1].split('|')[1:]
        #accession, locus, gene, start, stop, strand, organism, gc, hgt = field[1].split('|')
        accession, organism, locus, gene, start, stop, strand, gc = field[1].split('|')
        hgt = 'N/A'
    except:
        print 'error', field[1].split('|')
        #accession
    percent_ident = field[2]
    
    alignment_length = field[3]
    mismatches = field[4]
    gap_openings = field[5]
    query_start = field[6]
    query_end = field[7]
    subject_start = field[8]
    subject_end = field[9]
    e_val = field[10]
    bit_score = field[11]
    alignment_length = field[3]
    #operon_name = gene_to_operon_dict[predicted_gene]
    method = 'BLAST' # might have to change this to the version of blast that is used... not sure there.
    if source_accession == accession and source_locus == locus: # eliminates problem of homologous genes within a single operon
        method = "exact"
    
        
    return Homolog(accession, organism, locus, gene, predicted_gene, synonyms, e_val, percent_ident, bit_score, gc, start, stop, strand, product_type, alignment_length, method, source_accession, source_common, source_locus, source_start, hgt)


# i have to figure out a better name for this function. The jist of what i am doing here is as follows:
# First, i will provide a file name that contains all of the hits for every organism that we are interested in.
# Then it sorts this homolog list first by organism, then by locus. (By the required input, the files already have 
# been screened for both eval cutoff and operon membership. The function will then return a dictionary for that
# operon. The intention is that this data structure will then be used to find the best hit for the locus out of the
# many redundant hits, however this functionality will be handled another function that i have yet to write/test.
def return_operon_list(fname):
    operon = fname.split('/')[-1].split('.')[0]
    hlog_list = [Homolog.from_file(i.strip()) for i in open(fname).readlines()]
    result_dict = {}
    for hlog in hlog_list:
        accession = hlog.accession()
        locus = hlog.locus()
        if accession not in result_dict.keys():
            result_dict.update({accession: {}})
        if locus not in result_dict[accession].keys():
            result_dict[accession].update({locus: [hlog]})
        else:
            result_dict[accession][locus].append(hlog)
    #print result_dict[accession]
    
    return result_dict
        

# might  not use this: will see
def parallel_return_operon_list(infolder, outfolder, num_proc):
    pool = Pool(processes = num_proc)
    organism_dict_for_recovery = dict(pool.map(parallel_operon_fasta, genome_of_interest_list))


# This function will take the organism-locus dict (per operon file) and determine the best homolog.
# 
def best_homolog_list(operon_dict, outfile):
    result = []
    
    for org in sorted(operon_dict.keys()):
        for locus in sorted(operon_dict[org].keys()):
            hlog_list = operon_dict[org][locus][1:]
            best_hit = operon_dict[org][locus][0]
            gene_count = defaultdict(int) # i am goign to use this, to see if a locus has more than one predicted gene, and the count ratio
            gene_count[best_hit.predicted_gene()] +=1
            for hlog in hlog_list:
                gene_count[hlog.predicted_gene()] +=1
                if best_hit.e_val() > hlog.e_val():
                    best_hit = hlog
            print gene_count.keys()
            result.append(best_hit)
    handle = open(outfile, 'w')
    handle.write('\n'.join([i.ret_str() for i in result]))
    handle.close()
            

  
def main():
    
    start = time.time()
    
    parsed_args = parser_code()
    
    infolder, outfolder, operon_file, filter_file, num_proc = check_options(parsed_args)

    parallel_raw_parse(outfolder, infolder, num_proc)
    
    operon_file = './regulonDB/operon_name_and_genes_prot_only.txt'
    operon_dict = {}
    
    for line in [i.strip().split('\t') for i in open(operon_file).readlines()]:
        operon = line[0]
        for gene_entry in line[1:]:
            gene_name, gene_type = gene_entry.split(':')
            operon_dict.update({gene_name: operon})
    #print operon_dict
    
    parallel_blast_parse_dict('./blast_parse/organism_raw_info/', './blast_parse/filtered_homologs/', num_proc, './genbank_pathway_lists/nc_filter_file.txt', operon_dict)
    
    operon_dict = return_operon_list('./blast_parse/filtered_homologs/atpIBEFHAGDC.txt')
    
    best_homolog_list(operon_dict, './blast_parse/processed_operon_files/atpIBEFHAGDC.txt')
    
    print time.time() - start

    
if __name__ == '__main__':
    main()
