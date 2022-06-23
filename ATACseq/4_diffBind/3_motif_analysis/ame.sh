### This script takes the fasta file of an close or open DAR $1, used the corresponding shared peaks as backgound

# $2 is the output directory

pre=`basename $1 .fa`
dir=`dirname $1`
sh=`basename $1 | rev | cut -d'.' -f3- | rev`
#num=`wc -l $1| cut -d' ' -f1`


#echo $dir/$sh.shared.bed; echo $3/$sh.shared.$num.$4.bed


ame --oc $2/$pre.compShare --control $dir/$sh.shared.fa $1 ~/mdwilson/xuefei/reference/motifs/JASPAR2020_CORE_vertebrates_non_redundant_pfms_meme.txt






