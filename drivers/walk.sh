SRC=../src
GENOME=../example/Streptococcus_pyogenes
BACTERIOCINS=../bacteriocins
INTERMEDIATE=intermediate
ROOT_DIR=/home/jamie/Documents/Miami_bio/Bacfinder/Bacfinder/example
python $SRC/walk.py \
    --genes=$BACTERIOCINS/genes.fa \
    --root-dir=$ROOT_DIR \
    --bacteriocins=$BACTERIOCINS/bacteriocins.fa \
    --intermediate=$INTERMEDIATE \
    --radius=50000 \
    --keep-tmp \
    --verbose \
    --output-file=blast_results.txt 

    

