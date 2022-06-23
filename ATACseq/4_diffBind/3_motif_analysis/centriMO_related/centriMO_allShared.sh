### This script takes the fasta file of an close or open DAR.250_extened ($1), used the corresponding shared peaks as backgound to conduct ame analysis

# $1 is the output from single centriMo analysis
# $2 is the output directory

pre=`basename $1 .fa`
dir1=`dirname $1`
dir=`dirname $dir1`
sh=`basename $1 | rev | cut -d'.' -f4- | rev`
#num=`wc -l $1| cut -d' ' -f1`


#echo $dir/$sh.shared.bed; echo $3/$sh.shared.$num.$4.bed

echo $dir/$sh.shared*/$sh.shared*.fa

centrimo --oc $2/$pre.compShare --neg $dir/$sh.shared*/$sh.shared*.fa $1 ~/mdwilson/xuefei/reference/motifs/JASPAR2020_CORE_vertebrates_non_redundant_pfms_meme.txt




