Installation

Make sure that you have the following installed
1) biopython version 1.62
2) Beautiful Soup (only if you want to do online retrieval from the NCBI ftp site)


To start using this software, first do the following

1) Place your bacteriocin proteins in FASTA format in the bacteriocins folder

2) Place your genbank file in the genome folder.  This can be found on the NCBI website 

Then go into the drivers folder, and configure the following variables in the run.sh script

GENOME
BACTERIOCIN

Before running program, make sure that you have Biopython and blastall installed on your machine

To run the software execute run.sh in the drivers folder as follows

sh run.sh


