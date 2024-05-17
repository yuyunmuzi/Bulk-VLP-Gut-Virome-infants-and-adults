# Purpose: Use CoverM to calculate coverage of contigs

infile=$1

bwa index vNGS_checkv95.fa
mkdir -p 08_coverm
coverm contig \
    --mapper bwa-mem \
    --min-read-aligned-percent 95 \
    --min-read-percent-identity 90 \
    --min-covered-fraction 75 \
    --methods rpkm \
    --output-file 08_coverm/${infile}.coverm \
    --reference vNGS_checkv95.fa \
    --threads 10 \
    --coupled ../../05_Removed/${infile}_1.fastq ../../05_Removed/${infile}_2.fastq

sed -i '1d' ${infile}.coverm
sed -i "s/^/${infile}\t/g" ${infile}.coverm
cat ${infile}.coverm >> all.coverm.txt





