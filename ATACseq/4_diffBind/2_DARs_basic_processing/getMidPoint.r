#!/usr/bin/Rscript
# this script takes a bed file, get the midpoint of each region and output 
# intervals with just one base (midpoint)
# agrs[1] is the bed file
# args[2] is the output directory + file name

args <- commandArgs(trailingOnly = T)

df <- read.table(args[1], sep="\t", quote="", stringsAsFactors = F)
df$mid <- round((df$V3-df$V2)/2, digits = 0)
df$start <- df$V2+df$mid
df$end <- sapply(df$start, function(x) x+1)


df.out <- data.frame(chr=df$V1,
                     st=df$start,
                     en=df$end,
                     coor=df$V4
                     #cne=df$V5,
                     #type=df$V6,
                     #over=df$V7
					 )

write.table(df.out, file=args[2], sep = "\t", quote = F, row.names = F, col.names = F)