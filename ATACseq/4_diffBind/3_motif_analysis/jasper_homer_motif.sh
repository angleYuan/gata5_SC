### This script takes the fasta file of an close or open DAR $1, used the corresponding shared peaks as backgound
# $1 has to be absolute directory
# $2 is the output directory

pre=`basename $1 .fa`
dir=`dirname $1`
sh=`basename $1 | rev | cut -d'.' -f3- | rev`
#num=`wc -l $1| cut -d' ' -f1`


echo $dir/$sh.shared.fa; echo $1
mkdir -p $2/$pre.compShare
cd $2/$pre.compShare

homer2 known -i $1 -b $dir/$sh.shared.fa -m /home/xuefei/mdwilson/xuefei/reference/motifs/JASPAR2020_CORE_vertebrates_non-redundant_pfms_homer.txt -o enrich.txt

#findMotifs.pl $1 fasta $2/$pre.compShare -fasta $dir/$sh.shared.fa -find /home/xuefei/mdwilson/xuefei/reference/motifs/JASPAR2020_CORE_vertebrates_non-redundant_pfms_jaspar.txt







