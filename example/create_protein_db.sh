cat `find . -name "*.fna"` > all.fna
transeq all.fna -frame=6 -clean -outseq all_trans.fna