SRC=../src
GENOME=/data/genomes/Bacterial/all/all.fna   #File containing all of the bacterial genomes
BACDIR=../bacteriocins                       #Folder containing all of the bacteriocins
INTERMEDIATE=intermediate
BACTERIOCINS=$BACDIR/bacteriocins.fa         #Contains bacteriocins
TARGET_GENES=$BACDIR/genes.fa
BLASTED=$INTERMEDIATE/blast_results.txt
MULTIALIGN=$INTERMEDIATE/aligned.faa
ALIGN=$INTERMEDIATE/aligned.fa
HMMFILE=$INTERMEDIATE/bac.hmm
RESULT=out.txt
TRANSLATED=$INTERMEDIATE/sixpack.fa
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
    --intermediate=$INTERMEDIATE \
    --bac-evalue=1e-5 \
    --num-threads=7 \
    --radius=50000 \
    --verbose \
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