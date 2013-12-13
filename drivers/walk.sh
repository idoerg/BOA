SRC=../src
BACTERIOCINS=../bacteriocins
INTERMEDIATE=/data/tmp
#ROOT_DIR=/home/jamie/Documents/Miami_bio/Bacfinder/Bacfinder/example
ROOT_DIR=/data/genomes/Bacterial
python $SRC/walk.py \
    --genes=$BACTERIOCINS/genes.fa \
    --root-dir=$ROOT_DIR \
    --bacteriocins=$BACTERIOCINS/bacteriocins.fa \
    --intermediate=$INTERMEDIATE \
    --radius=50000 \
    --bac-evalue=1e-5 \
    --gene-evalue=1e-5 \
    --keep-tmp \
    --verbose \
    --output-file=blast_results.txt 

    

