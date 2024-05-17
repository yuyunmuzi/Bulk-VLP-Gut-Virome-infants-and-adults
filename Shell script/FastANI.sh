fa1=/mnt/raid9/liyun/2024_NC_infant/vNGS/Phage_identity_megahit/vNGS_checkv95.fa
fa2=/mnt/raid9/liyun/2024_NC_infant/mNGS/Phage_identity_megahit_noremove/mNGS_checkv95.fa
awk '/^>/{s=++num}{print > "mNGS/mNGS_"s".fa"}' $fa1
awk '/^>/{s=++num}{print > "vNGS/vNGS_"s".fa"}' $fa2

###list
ls -d -1 "$PWD/"*  >../mcontiglist.txt
ls -d -1 "$PWD/"*  >../vcontiglist.txt

###fastANI设置ANI（Average Nucleotide Identity）>95,AF (Alignment Fraction)>0.85
fastANI --ql mcontiglist.txt --rl vcontiglist.txt  --fragLen 100 --minFraction 0.85  -t 10 -o ani_output85.csv

awk  'BEGIN {FS=" "} $3>95  {print $1}' ani_output.csv|sort|uniq >m_speciesp
awk  'BEGIN {FS=" "} $3>95  {print $2}' ani_output.csv|sort|uniq >v_speciesp

while read line;do
head -n1 $line|sed 's/>//g'|cut -d " " -f 1 
done < m_speciesp >m_speciesp.contig

while read line;do
head -n1 $line|sed 's/>//g'|cut -d " " -f 1 
done < v_speciesp >v_speciesp.contig
