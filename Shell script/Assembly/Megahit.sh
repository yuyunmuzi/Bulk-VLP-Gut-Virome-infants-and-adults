infile=$1
mNGS_path=$2

mkdir -p 02_trimmed 05_Removed

TRIMMO_JAR_FILE='/mnt/raid2/tools/04_software/trimmomatic-0.39/trimmomatic-0.39.jar'
TRIMMO_ADAPTOR_FILE_PE='/mnt/raid2/tools/04_software/trimmomatic-0.39/adapters/TruSeq3-PE.fa'
TRIMMO_ADAPTOR_FILE_SE='/mnt/raid2/tools/04_software/trimmomatic-0.39/adapters/TruSeq3-SE.fa'

R1=${path}/${infile}.R1.fastq.gz
R2=${path}/${infile}.R2.fastq.gz
if [ ! -s 02_trimmed/${infile}_clean.1.fq ];then
    java -jar $TRIMMO_JAR_FILE PE -threads 16 ${R1} ${R2} 02_trimmed/${infile}_clean.1.fq 02_trimmed/${infile}_clean_unpaired.1.fq 02_trimmed/${infile}_clean.2.fq 02_trimmed/${infile}_clean_unpaired.2.fq ILLUMINACLIP:$TRIMMO_ADAPTOR_FILE_PE:2:15:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:15:30 MINLEN:50
fi

bowtie2 -p 10 --un-conc 05_Removed/${infile}_%.fastq --no-unal -k 20 -x /mnt/raid5/sunchuqing/Human_Gut_Phage/ref/hg38_ref -1 02_trimmed/${infile}_clean.1.fq -2 02_trimmed/${infile}_clean.2.fq >log


mkdir -p 06_megahit
megahit -t 10  --min-contig-len 1000 \
-1 05_Removed/${infile}_1.fastq -2 05_Removed/${infile}_2.fastq \
-o ./06_megahit/${infile} 


