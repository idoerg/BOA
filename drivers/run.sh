##User specified parameters
GENOME=/data/genomes/Bacterial/all/all.fna       #File containing all of the bacterial genomes
INTERGENES=/data/genomes/databases/intergenic.fa #File containing all of the intergenic regions in bacterial genomes
ANNOTATIONS=/data/tmp/annotations.fa             #Location of the annotations

##TEST files
#ANNOTATIONS=/data/tmp/test_annotations.fa
#INTERGENES=../example/Streptococcus_pyogenes/NC_011375_intergene.fa
#GENOME=../example/Streptococcus_pyogenes/NC_011375.fna
BACDIR=../bacteriocins                       #Folder containing all of the bacteriocins and genes

BACTERCIOCIN_FILE=bagel.fa                   #Location of bacteriocin file (where known bacteriocins are stored)
#BACTERCIOCIN_FILE=shaun.fa
GENES=genes.fa                               #Location of gene file (where sagB is stored)
##########################################################################################

BACTERIOCINS=$BACDIR/$BACTERCIOCIN_FILE      
TARGET_GENES=$BACDIR/$GENES                  #Location of gene file (where sagB is stored)
INTERMEDIATE=/tmp/intermediate               #Intermediate directory
#INTERMEDIATE=intermediate                   #Intermediate directory
BLASTED=$INTERMEDIATE/all_bacteria           #BLAST results (filtering outside of sagB neighborhoods)
ALIGN=$INTERMEDIATE/aligned.fa               #Poor mans multiple alignment
MULTIALIGN=$INTERMEDIATE/aligned.faa         #---
HMMFILE=$INTERMEDIATE/bac.hmm                #HMMER hmmbuild file
RESULT=out.txt                               #HMMER output

TRANSLATED=$INTERMEDIATE/sixpack.fa          #Temporary file where all of the translated proteins are
BLASTFA=$INTERMEDIATE/blast.fa
CLUSTER=$INTERMEDIATE/cluster 

SRC=../src
UTIL=$SRC/util
WRAPPERS=$SRC/wrappers

if [ ! -d $INTERMEDIATE ]:
then
    echo "Creating Intermediate file"
    mkdir $INTERMEDIATE
fi
#sh generate_annotations.sh
python $SRC/bacteriocin.py \
    --genes=$TARGET_GENES \
    --genome-files=$GENOME \
    --bacteriocins=$BACTERIOCINS \
    --intergenes=$INTERGENES \
    --annotations=$ANNOTATIONS \
    --intermediate=$INTERMEDIATE \
    --bac-evalue=1e-5 \
    --num-threads=7 \
    --bacteriocin-radius=50000 \
    --output=$BLASTED

# cat $BLASTED | awk '{print ">"$2"\n"$10}' > $BLASTFA
# cdhit -i $BLASTFA -o $CLUSTER -c 0.7
   



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