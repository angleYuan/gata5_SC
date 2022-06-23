### This script takes a bed file $1, get its midPoint
### then extend to certain base pairs ($2) to each end 
## get the fasta sequence and then submit to centriMO for motif enrichment analysis
## v2 version keeps the fasta file for other usages
## v2 version update JASPAR2020

# $3 is the output directory

pre=`basename $1 .bed`
#dir=`dirname $1`
#sh=`basename $1 | rev | cut -d'.' -f3- | rev`
#num=`wc -l $1| cut -d' ' -f1`

#echo $dir/$sh.shared.bed; echo $2/$sh.shared.$num.$3.bed

#bedtools sample -n $num -seed $3 -i $dir/$sh.shared.bed > $2/$sh.shared.$num.$3.bed

Rscript /hpf/largeprojects/mdwilson/xuefei/projects/2018_ATACseq/merged/zv11/local_scripts/getMidPoint.r $1 $3/$pre.mid.bed
check=`head -3 $3/$pre.mid.bed`
echo $3/$pre.mid.bed
echo $check
slopBed -b $2 -i $3/$pre.mid.bed -g ~/mdwilson/genomes/drer_zv11/ChromInfo.txt > $3/$pre.$2.bed
check2=`head -3 $3/$pre.$2.bed`
echo "$3/$pre.$2.bed"
echo $check2

fastaFromBed -fi ~/mdwilson/genomes/drer_zv11/Danio_rerio.Zv11.GRCz11.dna.chromosomes.fa -bed $3/$pre.$2.bed -fo $3/$pre.$2.fa
#fastaFromBed -fi ~/mdwilson/genomes/drer_zv11/Danio_rerio.Zv11.GRCz11.dna.chromosomes.fa -bed $2/$sh.shared.$num.$3.bed -fo $2/$sh.shared.$num.$3.fa

centrimo --oc $3/$pre.$2.single $3/$pre.$2.fa ~/mdwilson/xuefei/reference/motifs/JASPAR2020_CORE_vertebrates_non_redundant_pfms_meme.txt

rm $3/$pre.mid.bed
rm $3/$pre.$2.bed
mv $3/$pre.$2.fa $3/$pre.$2.single/
echo "done"
