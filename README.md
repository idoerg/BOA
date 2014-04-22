Installation
============

Make sure that you have the following installed

1.  Python2.7 or greater 
2.  biopython 
3.  matplotlib
4.  numpy
5.  panda
6.  bx-python
7.  blastall
8.  hmmer
9.  cdhit
10. clustalw
11. ftputil

If you are using ubuntu, most of these packages can be installed using the following command
```
sudo apt-get install python python-biopython python-matplotlib python-panda python-numpy python-ftputil clustalw cd-hit hmmer
```
bx-python can be installed through the following link: https://bitbucket.org/james_taylor/bx-python/wiki/Home


Getting Started
===============

First set PATH=$PATH:Bacfinder/src:Bacfinder/scripts

Before running the pipeline, two databases must be setup, the annotated genes database, and the intergenes database.  This is assuming that you already have a set of genbank files and fasta genomes contained within a root directory.

To create the annotated genes database, execute the following command:
```
python annotated_genes.py --root-dir=< root directory of genbank files > 
                          --output-file=< output file of annotated regions >
```
To create the intergenes database, execute the following command:
```
python intergene.py --root-dir=<root directory genbank files> 
                    --output-file=< output file of intergenic regions >
```
Both of these scripts are hierarchy independent.  It will find all of the genbank files within the root-directory, regardless of how the root directory is organized.

Now you are ready to run the blast pipeline.  To run the blast pipeline, run the following command
```
python bacteriocin.py --genome-files=< Fasta files of genomes >
                      --bacteriocins=< known bacteriocins fasta > 
                      --annotated-genes=< annotated genes database >  
                      --intergenes=< intergenes database > 
                      --intermediate=< A folder to store extra files > 
                      --output=<basename of output file>  
```
The output option is the basename for two different files.  If you output option is test, then the files you expect to see are test.annotated.txt and test.bacteriocins.

test.annotated.txt contains the list of annotated genes within the a radius around all of the blasted bacteriocins.  This search radius can be specified in the bacteriocin.py script.
The format of test.annotated.txt is a tab-delimited format with column headers specified as follows

1.  bacteriocin name 
2.  ncbi id of anchor gene
3.  blast bacteriocin start
4.  blast bacteriocin end
5.  blast bacteriocin strand
6.  accession id of whole genome
7.  anchor gene start 
8.  anchor gene end
9.  anchor gene strand
10. sequence of bacteriocin


test.bacteriocin.txt contains the list of bacteriocins aligned against all of the bacterial genomes provided.
The format of test.bacteriocins.txt is a tab-delimited format with column headers specified as follows

1. bacteriocin name 
2. ncbi id of species blasted against
3. bacteriocin start
4. bacteriocin end
5. bacteriocin strand
6. overlaps intergene or gene
7. blasted bacteriocin sequence

If you want to visualize the results from the blast pipeline, run the following command
```
python analyze.py --accession-table=< a map between accession and species >
                  --bacteriocins=< blasted bacteriocins in tab format >
                  --anchor-genes=< overlapping anchor genes in tab format >
```
The accession table can be found under the data folder
