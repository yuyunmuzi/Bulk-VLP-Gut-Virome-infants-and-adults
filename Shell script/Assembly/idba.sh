infile=$1
mNGS_path=$2

mkdir -p 02_trimmed 05_Removed 04_Abundance 06_Assembly 07_CD-HIT

TRIMMO_JAR_FILE='/mnt/raid2/tools/04_software/trimmomatic-0.39/trimmomatic-0.39.jar'
TRIMMO_ADAPTOR_FILE_PE='/mnt/raid2/tools/04_software/trimmomatic-0.39/adapters/TruSeq3-PE.fa'
TRIMMO_ADAPTOR_FILE_SE='/mnt/raid2/tools/04_software/trimmomatic-0.39/adapters/TruSeq3-SE.fa'

R1=${path}/${infile}.R1.fastq.gz
R2=${path}/${infile}.R2.fastq.gz
if [ ! -s 02_trimmed/${infile}_clean.1.fq ];then
    java -jar $TRIMMO_JAR_FILE PE -threads 16 ${R1} ${R2} 02_trimmed/${infile}_clean.1.fq 02_trimmed/${infile}_clean_unpaired.1.fq 02_trimmed/${infile}_clean.2.fq 02_trimmed/${infile}_clean_unpaired.2.fq ILLUMINACLIP:$TRIMMO_ADAPTOR_FILE_PE:2:15:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:15:30 MINLEN:50
fi

bowtie2 -p 10 --un-conc 05_Removed/${infile}_%.fastq --no-unal -k 20 -x /mnt/raid5/sunchuqing/Human_Gut_Phage/ref/hg38_ref -1 02_trimmed/${infile}_clean.1.fq -2 02_trimmed/${infile}_clean.2.fq >log
bowtie2 -p 16 --al-conc 05_Removed/${infile}_%_prophage.fastq --end-to-end -x /mnt/raid8/sunchuqing/Database/HGYG_prophage -1 05_Removed/${infile}_1.fastq -2 05_Removed/${infile}_2.fastq -S 04_Abundance/${infile}_NGS_HGYG_prophage.sam   >log
bowtie2 -p 16 --un-conc 05_Removed/${infile}_%_phage.fastq --end-to-end -x /mnt/raid8/sunchuqing/Database/HGYG -1 05_Removed/${infile}_1.fastq -2 05_Removed/${infile}_2.fastq -S 04_Abundance/${infile}_NGS_HGYG.sam  >log
cat 05_Removed/${infile}_1_prophage.fastq 05_Removed/${infile}_1_phage.fastq > 05_Removed/${infile}_virome_bft_1.fastq
cat 05_Removed/${infile}_2_prophage.fastq 05_Removed/${infile}_2_phage.fastq > 05_Removed/${infile}_virome_bft_2.fastq

R1=`ls 05_Removed/${infile}_virome_bft_1.fastq`
R2=`ls 05_Removed/${infile}_virome_bft_2.fastq`
java -jar $TRIMMO_JAR_FILE PE -threads 16 ${R1} ${R2} 05_Removed/${infile}_virome_1.fastq 02_trimmed/${infile}_clean_unpaired.vir.1.fq 05_Removed/${infile}_virome_2.fastq 02_trimmed/${infile}_clean_unpaired.vir.2.fq ILLUMINACLIP:$TRIMMO_ADAPTOR_FILE_PE:2:15:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:15:30 MINLEN:50


# IDBA
/mnt/raid5/sunchuqing/Softwares/idba/bin/fq2fa --merge --filter 05_Removed/${infile}_virome_1.fastq 05_Removed/${infile}_virome_2.fastq 05_Removed/${infile}.fa
/mnt/raid5/sunchuqing/Softwares/idba/bin/idba_ud -r 05_Removed/${infile}.fa --maxk 120 --step 10 -o 06_Assembly/${infile} --num_threads 16 --min_contig 1000

###cd-hit
cat 06_Assembly/${infile}/contig.fa >06_Assembly/${infile}.fasta
/mnt/raid9/liyun/Software/cd-hit-v4.6.8-2017-1208/cd-hit-est -i 06_Assembly/${infile}.fasta -o 07_CD-HIT/${infile}.fa -c 0.95 -n 5 -T 15 -M 16000 >07_CD-HIT/${infile}.log

