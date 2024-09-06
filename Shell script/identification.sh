infile=$1
seq="./06_megahit/${infile}/final.contigs.fa"

##1.virfinder >0.9 ; P < 0.01
mkdir  -p 02_virfinder_result 04_Positive
if [ ! -s "04_Positive/${infile}.virfiner.genome" ];then

/usr/bin/Rscript  /mnt/raid5/tools/Buffalo_gut/VirFinder.R $seq 02_virfinder_result/${infile}.csv
awk 'BEGIN {FS=","} $3>0.9 && $4<0.01 {print $1}' 02_virfinder_result/${infile}.csv|awk 'BEGIN {FS=" "} {print $1}' > 04_Positive/${infile}.virfiner.genome
fi

##2.virsorter2 
mkdir -p 01_virsorter_result
virsorter run \
        -w ./01_virsorter_result/${infile}  \
        -i ${seq} \
        --db-dir /mnt/raid8/sunchuqing/Softwares/Virsorterdb\
        -j 10 --rm-tmpdir\
        --tmpdir ./all_virsorter_result/temp --min-score 0.5 --keep-original-seq
cat ./01_virsorter_result/${infile}/final-viral-score.tsv | awk 'BEGIN {FS="|"} {print $1}'|sort|uniq  > 04_Positive/${infile}.VirSorter.genome

##3.VIBRANT
  mkdir -p 07_vibrant_result
    VIBRANT_run.py \
        -folder 07_vibrant_result/$infile \
        -i $seq \
        -t 10 \
        -no_plot \
        -virome # used when input is a VLP but not a metagenome 
awk 'BEGIN {FS="\t"}  {print $1}'  07_vibrant_result/$infile/VIBRANT_final.contigs/VIBRANT_results_final.contigs/VIBRANT_summary_results_final.contigs.tsv |awk 'BEGIN {FS=" "} {print $1}' > 04_Positive/$infile.vibrant.genome

##4.合并
mkdir -p 05_phage_Positive

cat 04_Positive/${infile}.virfiner.genome  04_Positive/${infile}.VirSorter.genome 04_Positive/${infile}.vibrant.genome|sort |uniq -c|sort -n|awk 'BEGIN {FS=" "} $1>=1 {print $2}'|sed 's/_length/ length/g'|awk 'BEGIN {FS=" "} {print $1}'  >05_phage_Positive/${infile}.phage.genome
grep -w -A 1 -Ff 05_phage_Positive/${infile}.phage.genome $seq |grep -v -e "--" |awk '/^>/&&NR>1{print "";}{ printf "%s",/^>/ ? $0" ":$0 }' |awk '{print $1","length($NF)}'|sed 's/>//g' |awk 'BEGIN {FS=","} $2>3000 {print $1}'|sort|uniq > 05_phage_Positive/${infile}_fullfill2
grep -w -A 1 -Ff 05_phage_Positive/${infile}_fullfill2 $seq |grep -v -e "--" > 05_phage_Positive/${infile}.phage.3k.genome.fa

sed  "s/>/>${infile}_/g" 05_phage_Positive/${infile}.phage.3k.genome.fa >>vNGS_all.fa

export CHECKVDB=/mnt/raid9/liyun/Software/checkv-db-v1.0

checkv end_to_end vNGS_all.fa  vNGS_all_checkv.fa -t 10

#(1) with higher number of viral genes than host genes, 
#(2) longer than 3 kbp, 
#(3) the times the viral sequences is represented in the contig less than/equal to one(kmer_freq ≤ 1), and 
#(4) without warning “>1 viral region detected” and “contig >1.5× longer than expected genome length”

cat quality_summary.tsv |awk 'BEGIN {FS="\t"} $6>$7 && $2>3000 && $13<=1 ' > vNGS_all_.check1.list
grep -v ">1 viral region detected" vNGS_all_.check1.list >VLP150_all_checkv.fa.list1
grep -v "contig >1.5x longer than expected genome length" VLP150_all_checkv.fa.list1 |awk 'BEGIN {FS="\t"} {print $1}' > VLP150_all_checkv.fa.genomelist
grep -w -A 1 -Ff ./vNGS_all_checkv.fa/VLP150_all_checkv.fa.genomelist vNGS_all.fa |grep -v -e "--"> vNGS_checkv.fa

/mnt/raid9/liyun/Software/cd-hit-v4.6.8-2017-1208/cd-hit-est -i vNGS_checkv.fa -o vNGS_checkv95.fa -c 0.95 -M 16000 -aS 0.8 -d 0 -T 20
