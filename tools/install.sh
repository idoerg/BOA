wget ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.2.28+-x64-linux.tar.gz
tar -zxvf ncbi-blast-2.2.28+-x64-linux.tar.gz
echo "export BLAST=$PWD/ncbi-blast-2.2.28+/bin" >> ~/.bashrc