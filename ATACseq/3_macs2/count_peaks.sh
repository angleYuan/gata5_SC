# $1 is the narrowPeak file
# $2 is the output name
## double check if you want the $2 file empty or not

pre=`basename $1 | cut -d'_' -f1`
n_peaks=`wc -l $1 | cut -d' ' -f1`
echo -e "$pre\t$n_peaks" >> $2