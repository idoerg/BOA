SRC=../src
GENOME=/home/jamie/Documents/Bacfinder/example/Streptococcus_pyogenes
BACTERIOCINS=../bacteriocins
GENBANK=$GENOME/NC_011375.gbk
INTERGENES=$GENOME/`basename $GENBANK .gbk`_ign.fasta

python $SRC/intergene.py --genbank-path=$GENBANK --output-dir=$GENOME 
python $SRC/bacteriocin.py \
    --genes=$BACTERIOCINS/genes.fa \
    --genbank-files=$GENOME/NC_011375.gbk \
    --intergenes=$INTERGENES \
    --bacteriocins=$BACTERIOCINS/bacteriocins.fa \
    --output-file=blast_results.txt 
    

