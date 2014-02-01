#File containing all of the bacterial genomes
GENOME=/home/jamie/Documents/Bacfinder/example/Streptococcus_pyogenes/NC_011375.fna
INTERGENES=/data/genomes/databases/intergenic.fa #File containing all of the bacterial genomes
BACDIR=../bacteriocins                           #Folder containing all of the bacteriocins/genes
BACTERCIOCIN_FILE=bacteriocins.fa                #Location of bacteriocin file (where known bacteriocins are stored)
GENES=genes.fa                                   #Location of gene file (where sagB is stored)

BACTERIOCINS=$BACDIR/$BACTERCIOCIN_FILE      
TARGET_GENES=$BACDIR/$GENES                      #Location of gene file (where sagB is stored)

INTERMEDIATE=intermediate                        #Intermediate directory
BLASTED=$INTERMEDIATE/blast_results.txt          #BLAST results (filtering outside of sagB neighborhoods)
ALIGN=$INTERMEDIATE/aligned.fa                   #Poor mans multiple alignment
MULTIALIGN=$INTERMEDIATE/aligned.faa             #---
HMMFILE=$INTERMEDIATE/bac.hmm                    #HMMER hmmbuild file
RESULT=out.txt                                   #HMMER output
					        
TRANSLATED=$INTERMEDIATE/sixpack.fa              #Temporary file where all of the translated proteins are

SRC=../src
UTIL=$SRC/util
WRAPPERS=$SRC/wrappers

if [ ! -d $INTERMEDIATE ]:
then
    echo "Creating Intermediate file"
    mkdir $INTERMEDIATE
fi

python $SRC/bacteriocin.py \
    --genes=$TARGET_GENES \
    --genome-files=$GENOME \
    --bacteriocins=$BACTERIOCINS \
    --intergenes=$INTERGENES \
    --intermediate=$INTERMEDIATE \
    --bac-evalue=1e-5 \
    --num-threads=7 \
    --radius=50000 \
    --verbose \
    --formatdb \
    --output-file=$BLASTED

   

# cat $BLASTED | tail -n+2 | awk '{print ">"$2,"\n"$9}' > $ALIGN
# LEN=`cat $ALIGN | awk 'NR%2==0' | wc -L`
# echo "longest sequence $LEN"
# cat $ALIGN | python $UTIL/formatFaa.py --length=$LEN > $MULTIALIGN 

# csplit -n 1 -f $INTERMEDIATE/genome $GENOME "%^>%" "/^>/" "{*}" -s
# FILES=$INTERMEDIATE/genome*
# python $WRAPPERS/sixpack.py --num-threads=8 --genome-dir=$INTERMEDIATE --basename="genome" --intermediate=$INTERMEDIATE

# for FILE in $INTERMEDIATE/genome*
# do
#     sixpack -sequence $FILE -outseq $FILE.six.faa -outfile $INTERMEDIATE/sixpack_visual
# done

# cat $INTERMEDIATE/genome*.six.faa > $TRANSLATED
# hmmbuild --informat afa $HMMFILE $MULTIALIGN
# hmmsearch $HMMFILE $TRANSLATED > $RESULT

# SRC=../src
# #GENOME=../example/Brachyspira_pilosicoli
# #GENOME=../example/all
# GENOME=/data/genomes/Bacterial/all
# BACTERIOCINS=../bacteriocins
# INTERMEDIATE=intermediate
# #INTERMEDIATE=/data/tmp
# python $SRC/bacteriocin.py \
#     --genes=$BACTERIOCINS/genes.fa \
#     --genome-files=$GENOME/*.fna \
#     --bacteriocins=$BACTERIOCINS/bacteriocins.fa \
#     --intermediate=$INTERMEDIATE \
#     --bac-evalue=1e-2 \
#     --num-threads=7 \
#     --verbose \
#     --radius=50000 \
#     --output-file=blast_results.txt 

    

