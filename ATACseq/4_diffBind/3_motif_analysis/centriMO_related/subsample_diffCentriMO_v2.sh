### This script takes an close or open DAR $1, count how many lines existed
## then extend to certain base pairs ($2) to each end 
## then subsampled the same number*$5 of regions from the corresponding shared peaks (increase the background region numbers to reduce noise when the the real input only contains a small number of regions) (V2 specific features)
## get the fasta for both and then submit to centriMO for motif enrichment analysis
## (V2 specific features)keep the final input and background fasta file in case needed for future analysis


# $3 is the output directory
# $4 is the seed for sampling
# $5 is the multiply factor for backgound regions

pre=`basename $1 .bed`
dir=`dirname $1`
sh=`basename $1 | rev | cut -d'.' -f3- | rev`
num=`wc -l $1| cut -d' ' -f1`
let bnum=$5*$num
echo $bnum

#echo $dir/$sh.shared.bed; echo $3/$sh.shared.$num.$4.bed

Rscript /hpf/largeprojects/mdwilson/xuefei/projects/2018_ATACseq/merged/zv11/local_scripts/getMidPoint.r $1 $3/$pre.mid.bed
slopBed -b $2 -i $3/$pre.mid.bed -g ~/mdwilson/genomes/drer_zv11/ChromInfo.txt > $3/$pre.$2.bed
echo $3/$pre.$2.bed

bedtools sample -n $bnum -seed $4 -i $dir/$sh.shared.bed > $3/$sh.shared.$bnum.$4.bed
echo $3/$sh.shared.$bnum.$4.bed
Rscript /hpf/largeprojects/mdwilson/xuefei/projects/2018_ATACseq/merged/zv11/local_scripts/getMidPoint.r $3/$sh.shared.$bnum.$4.bed $3/$sh.shared.$bnum.$4.mid.bed
echo $3/$sh.shared.$bnum.$4.mid.bed
slopBed -b $2 -i $3/$sh.shared.$bnum.$4.mid.bed -g ~/mdwilson/genomes/drer_zv11/ChromInfo.txt > $3/$sh.shared.$bnum.$4.$2.bed
echo $3/$sh.shared.$bnum.$4.$2.bed

fastaFromBed -fi ~/mdwilson/genomes/drer_zv11/Danio_rerio.Zv11.GRCz11.dna.chromosomes.fa -bed $3/$pre.$2.bed -fo $3/$pre.$2.fa
fastaFromBed -fi ~/mdwilson/genomes/drer_zv11/Danio_rerio.Zv11.GRCz11.dna.chromosomes.fa -bed $3/$sh.shared.$bnum.$4.$2.bed -fo $3/$sh.shared.$bnum.$4.$2.fa

centrimo --oc $3/$pre.compShare.s$5_$4 --neg $3/$sh.shared.$bnum.$4.$2.fa $3/$pre.$2.fa ~/mdwilson/xuefei/reference/motifs/motif_databases.12.18/JASPAR/JASPAR2018_CORE_vertebrates_non-redundant.meme

rm $3/$pre.mid.bed
rm $3/$pre.$2.bed
mv $3/$pre.$2.fa $3/$pre.compShare.s$5_$4
rm $3/$sh.shared.$bnum.$4.bed
rm $3/$sh.shared.$bnum.$4.mid.bed
rm $3/$sh.shared.$bnum.$4.$2.bed
mv $3/$sh.shared.$bnum.$4.$2.fa $3/$pre.compShare.s$5_$4






