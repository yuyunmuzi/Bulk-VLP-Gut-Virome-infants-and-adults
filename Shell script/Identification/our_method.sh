infile=$1
seq="07_CD-HIT/${infile}.fa"

##1.virfinder 
mkdir  -p 02_virfinder_result 04_Positive
if [ ! -s "04_Positive/${infile}.virfiner.genome" ];then
    mkdir -p 03_VirFinder
    /user/bin/Rscript  /mnt/raid5/sunchuqing/Buffalo_gut/VirFinder.R $seq 03_VirFinder/${infile}.csv
    #VirFinder输出文件
    awk 'BEGIN {FS=","} $3>0.6 {print $1}' 03_VirFinder/${infile}.csv|awk 'BEGIN {FS=" "} {print $1}' > 04_Positive/${infile}.virfiner.genome
fi

##2.virsorter2 
mkdir -p 01_virsorter_result
virsorter run \
        -w ./01_virsorter_result/${infile} \
        -i $seq\
        --db-dir /mnt/raid8/sunchuqing/Softwares/Virsorterdb\
        -j 16 --rm-tmpdir\
        --tmpdir ./01_virsorter_result/${infile}/temp --min-score 0.7
cat ./01_virsorter_result/${infile}/final-viral-score.tsv | awk 'BEGIN {FS="|"} {print $1}'|sort|uniq  > 04_Positive/${infile}.VirSorter.genome

mkdir -p 05_phage_Positive
cat 04_Positive/${infile}.virfiner.genome  04_Positive/${infile}.VirSorter.genome|sort |uniq -c|sort -n|awk 'BEGIN {FS=" "} $1>=2 {print $2}'|sed 's/_length/ length/g'|awk 'BEGIN {FS=" "} {print $1}'  >05_phage_Positive/${infile}.phage.genome
grep -w -A 1 -Ff 05_phage_Positive/${infile}.phage.genome $seq |grep -v -e "--" |awk '/^>/&&NR>1{print "";}{ printf "%s",/^>/ ? $0" ":$0 }' |awk '{print $1","length($NF)}'|sed 's/>//g' |awk 'BEGIN {FS=","} $2>5000 {print $1}'|sort|uniq > 05_phage_Positive/${infile}_fullfill2
grep -w -A 1 -Ff 05_phage_Positive/${infile}_fullfill2 $seq |grep -v -e "--" > 05_phage_Positive/${infile}.phage.5k.genome.fa

sed  "s/>/>${infile}_/g" ./05_phage_Positive/${infile}.phage.5k.genome.fa  >>mNGS_all.fa

cd-hit-est -i mNGS_all.fa  -o mNGS_all95.fa -c 0.95 -M 16000 -aS 0.8 -d 0 -t 20
