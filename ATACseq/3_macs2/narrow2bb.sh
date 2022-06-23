# This script convert a narrowPeak file (or other bed format) into a bigBed for display
# $1 is the narrowPeak file
# $2 is the ChromInfo file
# $3 is the output directory

pre=`basename $1 | cut -d'_' -f1`; echo $pre

awk -F'\t' '{OFT="\t";}{print $1,$2,$3, $1":"$2"-"$3}' $1 > $3/$pre.bed
bedToBigBed -type=bed4 $3/$pre.bed $2 $3/$pre.bb
rm $3/$pre.bed
