#!/usr/bin/python

from multiprocessing import Pool
import time
import os
import sys
import argparse

# Copyright(C) 2013 David Ream
# Released under Biopython license. http://www.biopython.org/DIST/LICENSE
# Do not remove this comment

# This exists to  make the main function easier to read. It contains code to run the argument parser, and does nothing else.
def parser_code():

    parser = argparse.ArgumentParser(description='Run BLAST on a set of searchable database using a query that is specified by the user. Results will be stored by the accession number of the database. Currently this script will only accept protein queries, but I will update to automatically run on all types of genes, as most of the information needed for this behavior exists.')

    parser.add_argument("-i", "--infolder", dest="infolder", metavar="FOLDER", default='./db/',
                help="Folder containing all BLAST searchable databases to be used by the program.")
                 
    parser.add_argument("-o", "--outfolder", dest="outfolder", metavar="FOLDER", default='./blast_result/',
                help="Folder where the BLAST results will be stored.")
    
    parser.add_argument("-f", "--filter", dest="filter", metavar="FILE", default='',
                help="File restrictiong which accession numbers this script will process. If no file is provided, filtering is not performed.")
                
    parser.add_argument("-n", "--num_proc", dest="num_proc", metavar="INT", default = os.sysconf("SC_NPROCESSORS_CONF"), type=int,
                help="Number of processors that you want this script to run on. The default is every CPU that the system has.")
    
    # Fix this option, ultimately it should be a folder, that reads in a operon file with the name of operon(s) with a list of gene names and types.
    # The program will then take two files, protein and nucleic acid queries and run them. (This would require that there are two seperate, and complementary
    # blast databases within this folder. My desire is that it would be two subfolders, 'protein/' and 'rna/' which house these sets of data.
    parser.add_argument("-q", "--query", dest="query", default='./operon_query_files/protein_matches.fa', metavar="FILE",
                help="A file that contains the BLAST query for every gene of interest in the dataset.")
                
    parser.add_argument("-e", "--eval", dest="eval", default='1e-10', metavar="FLOAT", type=float,
                help="eval for the BLAST search.")
                
    return parser.parse_args()


def check_options(parsed_args):
    # section of code that checks the infolder entry    
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
        num_proc = int(parsed_args.num_proc)
    
    # Check the query file
    if parsed_args.query == '' or os.path.exists(parsed_args.query):
        query_file = parsed_args.query
    else:
        print "The file %s does not exist." % parsed_args.query
        sys.exit()
    
    e_val = parsed_args.eval
        
    return infolder, outfolder, filter_file, num_proc, query_file, e_val
 
        
#this function will return all of the files that are in a directory. os.walk is recursive traversal.
def returnRecursiveDirFiles(root_dir):
    result = []
    for path, dir_name, flist in os.walk(root_dir):
        for f in flist:
            fname = os.path.join(path, f)
            if os.path.isfile(fname):
                result.append(fname)
    return result


# This code right now only deals with protein, but I will add functionality later for nucleotides. 
# Just moving the project along here, but this is a critical flaw moving forward.
def do_parallel_blast(arg_tuple):
    db, query_file, blast_result_folder, num_processors, eval_threshold = arg_tuple
    out_file = "%s%s_prot.txt" % (blast_result_folder, db.split('/')[-1].split('.')[0])
    cmd = "blastall -p tblastn -a %i -i %s -d %s -e %s -o %s -m 9" % (num_processors, query_file, db, eval_threshold, out_file)
    print cmd
    os.system( cmd )


#def parallel_blast(infile, query, folder, num_proc, e_val = '1e-10'):
def parallel_blast(infolder, outfolder, filter_file, num_proc, query_file, e_val):
    # you kinda have to trust me here, but having blast run on as many threads per CPU as you have total processors is fastest
    # I have no idea why this is... ugh.
    
    unfiltered_db_list = [i for i in returnRecursiveDirFiles(infolder) if i.split('/')[-1].split('.')[-1] == 'ffc']

    if filter_file == '':
        db_list = unfiltered_db_list
    else:
        filter_list = [i.strip() for i in open(filter_file).readlines()]
        db_list = [i for i in unfiltered_db_list if i.split('/')[-1].split('.')[0] in filter_list]
    
    #print len(unfiltered_db_list), len(db_list)
    
    blast_arg_list = [(i, query_file, outfolder, num_proc, e_val) for i in db_list]
    pool = Pool(processes = num_proc)
    pool.map(do_parallel_blast, blast_arg_list)


def main():
    
    start = time.time()
    
    '''parser = argparse.ArgumentParser(description="Conduct a BLAST search over a list of  BLAST searchable databases using a common query file. The program will save the results in a folder designated by the user or the default './blast_result/'.")
                
    parser.add_argument("-i", "--infile", dest="infile", default='', metavar="FILE", required=True,
                help="A file that contains the path to every organism database that you are interested in.")
    
    parser.add_argument("-f", "--folder", dest="folder", metavar="FOLDER", default='./blast_result/',
                help="Folder where the BLAST results will be stored. Default is the folder './blast_result/'.")
    
    parser.add_argument("-q", "--query", dest="query", default='./blast_database_list.txt', metavar="FILE",
                help="A file that contains the BLAST query for every gene of interest in the dataset.")
    
    parser.add_argument("-n", "--num_proc", dest="num_proc", metavar="INT", default = os.sysconf("SC_NPROCESSORS_CONF"), type=int,
                help="Number of processors that you want this script to run on. The default is every CPU that the system has.")

    parser.add_argument("-e", "--eval", dest="eval", default='1e-10', metavar="FLOAT",
                help="eval for the BLAST search.")
    
    parsed_args = parser.parse_args()'''
    
    parsed_args = parser_code()
    
    infolder, outfolder, filter_file, num_proc, query_file, e_val = check_options(parsed_args)
    
    print infolder, outfolder, filter_file, num_proc, query_file, e_val
    
    #parallel_blast(infile, query, folder, num_proc, e_val)
    parallel_blast(infolder, outfolder, filter_file, num_proc, query_file, e_val)

    print time.time() - start

    # ./blast_script.py -i blast_database_list.txt -f ./blast_result/ -q ./fasta_generation/protein_matches.fa 
    
if __name__ == '__main__':
    main()

