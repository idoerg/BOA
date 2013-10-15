SRC=../src
GENOME=../example/Streptococcus_pyogenes
BACTERIOCINS=../bacteriocins
GENBANK=$GENOME/NC_011375.gbk
INTERMEDIATE=intermediate
python $SRC/bacteriocin.py \
    --genes=$BACTERIOCINS/genes.fa \
    --genbank-files=$GENOME/NC_011375.gbk \
    --bacteriocins=$BACTERIOCINS/bacteriocins.fa \
    --intermediate=$INTERMEDIATE \
    --radius=50000 \
    --output-file=blast_results.txt 

    

