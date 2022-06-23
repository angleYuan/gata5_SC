## this script takes a bam file ($1), filter it for q30, index the filtered bam file, remove the bedGraph, bed, and another bam file generated in the alignment pipline

pre=`basename $1 .bam`

samtools view -b -q 30 $1 > $pre.filter.bam; echo "q30 filter done"
samtools index $pre.filter.bam; echo "index done"
samtools idxstats $pre.filter.bam > $pre.filter.idxstats; echo "idxtsats done"
cut -f1 $pre.filter.idxstats | grep -v chrM | grep -v "*" | xargs samtools view -b $pre.filter.bam > $pre.filter.NoChrM.bam; echo "NoChrM done"
samtools index $pre.filter.NoChrM.bam ; echo "index done"

rm $pre.bed.gz
rm $pre.bedGraph
rm $pre.R1.trimmed.bam
rm $pre.trimmed.bam
