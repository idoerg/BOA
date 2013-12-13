Installation
============

Make sure that you have the following installed
1) biopython
2) Beautiful Soup (only if you want to do online retrieval from the NCBI ftp site)
3) blastall

Getting Started
===============
To run the run.sh in the drivers folder modify the four following parameters

GENOME

BACDIR

BACTERCIOCIN_FILE

GENES

Then the run script by executing the command 

sh run.sh

To run the python script in the src folder, use the something along the lines of the following command 

python bacteriocin.py --genes=genes.fa --genome-files=all.fna --bacteriocins=bacteriocin.fa --intermediate=tmp_dir --bac-evalue=1e-5 --num-threads=7 --radius=50000 --output-file=blast.out


