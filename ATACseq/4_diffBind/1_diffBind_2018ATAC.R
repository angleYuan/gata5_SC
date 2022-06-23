#### This file contains the main DiffBind analysis for identifying differentially accessible peaks 
#### besides GFP+ VS GFP- comparison and Gata5/6 KD VS WT comparison, some other samples that were not published in the Science Advance paper were also included here. 

library(DiffBind)
library(VennDiagram)
library(UpSetR)

dodiffAnalysis <- function(path, summits=250, compareOn, minMembers=2,
                           method=DBA_ALL_METHODS) {
  diff <- dba(sampleSheet = path)
  #plot(diff)
  diff <- dba.count(diff,summits = summits)
  
  diff <- dba.contrast(diff, categories = compareOn, minMembers = 2)
  diff <- dba.analyze(diff, method = DBA_ALL_METHODS, bSubControl = F)
  return(diff)
}


p.l <- list.files(path = ".", pattern = ".csv", full.names = F)
names(p.l) <- lapply(p.l, function(x) unlist(strsplit(x, "[.]"))[1])
compare.l <- list("DBA_TISSUE","DBA_TISSUE","DBA_CONDITION", "DBA_CONDITION", "DBA_CONDITION", "DBA_CONDITION","DBA_TISSUE")
names(compare.l) <- names(p.l)

gfpPosi_wtVSg56KD <- dodiffAnalysis(p.l[["gfpPosi_wtVSg56KD"]], compareOn = DBA_CONDITION)
gfpPosi_wtVSaplnrKD <- dodiffAnalysis(p.l[["gfpPosi_wtVSaplnrKD"]], compareOn = DBA_CONDITION)
wt_gfpPosiVSNega <- dodiffAnalysis(p.l[["wt_gfpPosiVSNega"]], compareOn = DBA_TISSUE)
gfpNega_wtVSg56KD <- dodiffAnalysis(p.l[["gfpNega_wtVSg56KD"]], compareOn = DBA_CONDITION)
gfpNega_wtVSaplnrKD <- dodiffAnalysis(p.l[["gfpNega_wtVSaplnrKD"]], compareOn = DBA_CONDITION)
aplnrKD_gfpPosiVSNega <- dodiffAnalysis(p.l[["aplnrKD_gfpPosiVSNega"]], compareOn = DBA_TISSUE)
g56KD_gfpPosiVSNega <- dodiffAnalysis(p.l[["g56KD_gfpPosiVSNega"]], compareOn = DBA_TISSUE)

dba.center.ls <- list("gfpPosi_wtVSg56KD"=gfpPosi_wtVSg56KD,
                      "gfpPosi_wtVSaplnrKD"=gfpPosi_wtVSaplnrKD,
                      "wt_gfpPosiVSNega"=wt_gfpPosiVSNega,
                      "gfpNega_wtVSg56KD"=gfpNega_wtVSg56KD,
                      "gfpNega_wtVSaplnrKD"=gfpNega_wtVSaplnrKD,
                      "g56KD_gfpPosiVSNega"=g56KD_gfpPosiVSNega,
                      "aplnrKD_gfpPosiVSNega"=aplnrKD_gfpPosiVSNega)

fixNegaBugs <- function(vec) {
  n_v <- sapply(vec, function(x) {
    if (x < 0) {
      x <- 0
      return(x)
    } else return(x)
  })
  return(n_v)
}
t <- df.all
t$End <- fixNegaBugs(t$End)
t[t$End == 0,]

diff2report <- function(diffA, name, method, DataType = DBA_DATA_FRAME, th=0.05, fold=0) {
  options(scipen = 999)
  
  df.all <- dba.report(diffA, method = method, DataType = DBA_DATA_FRAME, th = 1)
  df.all$Start <- fixNegaBugs(df.all$Start); df.all$End <- fixNegaBugs(df.all$End)
  #print("df.all")
  df <- dba.report(diffA, method = method, DataType = DBA_DATA_FRAME, th = th, fold = fold)
  df$Start <- fixNegaBugs(df$Start); df$End <- fixNegaBugs(df$End)
  #print("df")
  #df.notDB <- dba.report(diffA, method = method, DataType = DBA_DATA_FRAME, th = th, fold = fold, bNotDB = T)
  #df.notDB <- as.data.frame(df.notDB$peaks)
  #print(head(df.notDB)); print(head(df))
  
  if (is.null(df)) {
    print(paste0(name, " has no records"))
  } else {
    f.name <- paste0("./output_table/", name, ".", method, ".pV_", th, ".fold_",fold)
    DB.name <- paste0(f.name, ".DB.txt")
    nDB.name <- paste0(f.name, ".NotDB.txt")
    all.name <- paste0(f.name, ".all.txt")
    b.close.name <- paste0(f.name,".close.bed")
    b.open.name <- paste0(f.name,".open.bed")
    b.shared.name <- paste0(f.name, ".shared.bed")
    b.all.name <- paste0(f.name, ".all.bed")
    #print(f.name);print(b.name)
    
    df$id <- paste0(df$Chr,":", df$Start, "-", df$End)
    rownames(df) <- df$id
    df <- df[order(df$Fold, decreasing = T), ]
    df.close <- df[df$Fold > 0, c("Chr", "Start", "End", "id")]  
    df.open <- df[df$Fold < 0, c("Chr", "Start", "End", "id")]
    
    df.all$id <- paste0(df.all$Chr, ":", df.all$Start, "-", df.all$End)
    rownames(df.all) <- df.all$id
    df.all.bed <- df.all[, c("Chr", "Start", "End", "id")]
    
    df.shared <- df.all[setdiff(rownames(df.all), rownames(df)), ]
    
    df.shared.bed <- df.shared[, c("Chr", "Start", "End", "id")]
    
    write.table(df, file = DB.name, sep = '\t', 
                col.names=NA, row.names=T, quote = F)
    write.table(df.shared, file = nDB.name, sep = '\t', 
                col.names=NA, row.names=T, quote = F)
    write.table(df.all, file = all.name, sep = '\t', 
                col.names=NA, row.names=T, quote = F)
    
    write.table(df.close, file = b.close.name, sep = '\t', col.names=F, row.names=F, quote = F)
    write.table(df.open, file = b.open.name, sep = '\t', col.names=F, row.names=F, quote = F)
    write.table(df.shared.bed, file = b.shared.name, sep = '\t', col.names=F, row.names=F, quote = F)
    write.table(df.all.bed, file = b.all.name, sep = '\t', col.names=F, row.names=F, quote = F)
    print(name)
  }
  
}

## test
diff2report(diffA = dba.center.ls[["gfpPosi_wtVSg56KD"]], name = "gfpPosi_wtVSg56KD", 
            method = DBA_EDGER, th = 0.05, fold = 0)

lapply(names(dba.center.ls), 
       function(x) diff2report(diffA = dba.center.ls[[x]], name = x, 
                               method = DBA_EDGER, th = 0.05, fold = 0))
lapply(names(dba.center.ls), 
       function(x) diff2report(diffA = dba.center.ls[[x]], name = x, 
                               method = DBA_EDGER, th = 0.05, fold = 1))
lapply(names(dba.center.ls), 
       function(x) diff2report(diffA = dba.center.ls[[x]], name = x, 
                               method = DBA_EDGER, th = 0.05, fold = 1.5))

lapply(names(dba.center.ls), 
       function(x) diff2report(diffA = dba.center.ls[[x]], name = x, 
                               method = DBA_DESEQ2, th = 0.05, fold = 0))
lapply(names(dba.center.ls), 
       function(x) diff2report(diffA = dba.center.ls[[x]], name = x, 
                               method = DBA_DESEQ2, th = 0.1, fold = 0))
############################## have venn comparison of different results #################################
bedToVennlist <- function(filename){
  df <- try(read.table(filename, sep='\t', header=F))
  if(class(df)=='try-error'){
    df=NULL
  } else {
    df$V4 <- paste0(df$V1,":",df$V2,"-", df$V3)
    return(df$V4)
  }
  #df <- read.table(filename, sep='\t', header=F)
  #df$V4 <- paste(df$V1, df$V2, df$V3, sep = '_')
  #return(df$V4)
}

makeUpsetPlots_oneContrast <- function(comparison, openOrclose, path="./output_table/",
                                       output="./upset_plot/", w=7, h=5, text.scale = 1.4) {
  pa <- paste0("^", comparison, ".*", openOrclose, ".bed$")
  print(pa)
  in.ls <- list.files(path = path, pattern = pa , full.names = T)
  n.ls <- list.files(path = path, pattern = pa , full.names = F)
  names(in.ls) <- lapply(n.ls, function(x) {
    v <- unlist(strsplit(x, '[.]'))
    v <- v[-length(v)]
    v <- v[-length(v)]
    v <- v[-1]
    n <- paste(v, collapse = ".")
    return(n)
  })
  
  pre <- paste0(unlist(strsplit(n.ls[[1]], '[.]'))[1], ".",
                unlist(strsplit(n.ls[[1]], '[.]'))[length(unlist(strsplit(n.ls[[1]], '[.]')))-1])
  dir.create(output, showWarnings = F)
  f.name <- paste0(output, "/", pre, ".stats_thres.pdf")
  
  in_bed.ls <- lapply(in.ls, bedToVennlist)
  names(in_bed.ls) <- names(in.ls)
  in_bed.ls <- Filter(Negate(is.null),in_bed.ls)
  print("list finished")
  
  if (is.null(in_bed.ls)) {
    print("all lists are empty")
  } else {
    pdf(f.name, w,h, useDingbats = F)
    upset(fromList(in_bed.ls), order.by = "freq", number.angles = 45, 
          text.scale = text.scale, keep.order = T)
    dev.off()
  }
  
  
  
  
  #print(in.ls)
}
makeUpsetPlots_oneContrast(comparison = "wt_gfpPosiVSNega", openOrclose = "open")
sapply(names(dba.center.ls), function(x) makeUpsetPlots_oneContrast(comparison = x, openOrclose = "open"))
sapply(names(dba.center.ls), function(x) makeUpsetPlots_oneContrast(comparison = x, openOrclose = "close"))

################################ compare overlaps between different DARs ###############################
makeUpsetPlots <- function(pa="u.bed",path, w=12, h=5, text.scale = 1.3, 
                           f.n="upset.pdf", sets=NULL) {
  print(pa)
  in.ls <- list.files(path = path, pattern = pa , full.names = T)
  n.ls <- list.files(path = path, pattern = pa , full.names = F)
  names(in.ls) <- lapply(n.ls, function(x) {
    v <- unlist(strsplit(x, '[.]'))
    n <- paste0(v[1],".",v[6])
    return(n)
  })
  
  #pa <- unlist(strsplit(pa, '[.]'))[1]
  #dir.create(output, showWarnings = F)
  f.name <- paste0(path, "/", f.n)
  
  in_bed.ls <- lapply(in.ls, bedToVennlist)
  names(in_bed.ls) <- names(in.ls)
  in_bed.ls <- Filter(Negate(is.null),in_bed.ls)
  print("list finished")
  print(names(in_bed.ls))
  
  if (!is.null(sets)){
    sets <- names(in_bed.ls)
  }
  
  if (is.null(in_bed.ls)) {
    print("all lists are empty")
  } else {
    pdf(f.name, w,h, useDingbats = F)
    upset(fromList(in_bed.ls), nintersects = NA, nsets = length(in_bed.ls), 
          sets=sets, order.by = "freq", number.angles = 45, 
          text.scale = text.scale, keep.order = T)
    dev.off()
  }
  
  
  
  
  #print(in.ls)
}

makeUpsetPlots(path = "./overlap_comparison/sf6_all_interested/")
makeUpsetPlots(path = "./overlap_comparison/sf6_all_interested/", sets = c(2:5),
               f.n = "gfpPosi.upset.pdf")

in.ls <- list.files(path = "./overlap_comparison/sf6_all_interested/", pattern = "u.bed$" , full.names = T)
n.ls <- list.files(path = "./overlap_comparison/sf6_all_interested/", pattern = "u.bed$" , full.names = F)
names(in.ls) <- lapply(n.ls, function(x) {
  v <- unlist(strsplit(x, '[.]'))
  n <- paste0(v[1],".", v[6])
  return(n)
})

in_bed.ls <- lapply(in.ls, bedToVennlist)
names(in_bed.ls) <- names(in.ls)
in_bed.ls <- Filter(Negate(is.null),in_bed.ls)

pdf("./overlap_comparison/sf6_all_interested/close.upset2.pdf", 9,5, useDingbats = F)
upset(fromList(in_bed.ls), nsets = length(in_bed.ls), order.by = "freq", number.angles = 45, group.by = "sets",
      text.scale = 1.4, keep.order = T, sets = c("gfpPosi_wtVSaplnrKD.close", "gfpPosi_wtVSg56KD.close", "wt_gfpPosiVSNega.close"))
dev.off()

pdf("./overlap_comparison/sf6_all_interested/close_gfpNega.upset.pdf", 9,5, useDingbats = F)
upset(fromList(in_bed.ls), nsets = length(in_bed.ls), order.by = "freq", number.angles = 45, group.by = "sets",
      text.scale = 1.4, keep.order = T, sets = c("gfpPosi_wtVSaplnrKD.close", "gfpPosi_wtVSg56KD.close", "wt_gfpPosiVSNega.open"))
dev.off()

pdf("./overlap_comparison/sf6_all_interested/open.upset2.pdf", 9,5, useDingbats = F)
upset(fromList(in_bed.ls), nsets = length(in_bed.ls), order.by = "freq", number.angles = 45, group.by = "sets",
      text.scale = 1.4, keep.order = T, sets = c("gfpPosi_wtVSaplnrKD.open", "gfpPosi_wtVSg56KD.open", "wt_gfpPosiVSNega.open"))
dev.off()

pdf("./overlap_comparison/sf6_all_interested/open.upset.pdf", 9,5, useDingbats = F)
upset(fromList(in_bed.ls), nsets = length(in_bed.ls), order.by = "freq", number.angles = 45, group.by = "sets",
      text.scale = 1.4, keep.order = T, sets = c("gfpPosi_wtVSaplnrKD.open", "gfpPosi_wtVSg56KD.open", "wt_gfpPosiVSNega.close"))
dev.off()

pdf("./overlap_comparison/sf6_all_interested/open_all.upset.pdf", 9,5, useDingbats = F)
upset(fromList(in_bed.ls), nsets = length(in_bed.ls), order.by = "freq", number.angles = 45, group.by = "sets",
      text.scale = 1.4, keep.order = T, sets = c("gfpPosi_wtVSaplnrKD.open", "gfpPosi_wtVSg56KD.open", "wt_gfpPosiVSNega.close", "wt_gfpPosiVSNega.open"))
dev.off()


library(VennDiagram)
venn.diagram(list(in_bed.ls[["gfpPosi_wtVSaplnrKD.close"]], 
                  in_bed.ls[["gfpPosi_wtVSg56KD.close"]],
                  in_bed.ls[["wt_gfpPosiVSNega.close"]]),
             filename = "overlap_comparison/test.tiff", fill = c("red", "green", "blue"),
             alpha = c(0.5, 0.5), cex = 2,cat.fontface = 4,lty =2, fontfamily =3)
library(venneuler)
################################# generate bed file withs colors to show on UCSC #####################
agrs <- commandArgs(trailingOnly = F)
bedForUCSC <- function(pre, fold=1, out_dir="./bed4UCSC/") {
  all.path <- paste0("./shift_corrected/", pre, ".edgeRGLM.pV_0.05.fold_", fold, ".all.bed")
  close.path <- paste0("./shift_corrected/", pre, ".edgeRGLM.pV_0.05.fold_", fold, ".close.bed")
  open.path <- paste0("./shift_corrected/", pre, ".edgeRGLM.pV_0.05.fold_", fold, ".open.bed")
  print(c(all.path, close.path, open.path))
  
  dir.create(out_dir, showWarnings = F)
  f.name <- paste0(out_dir, pre, ".all.9c.bed")
  
  bed.all <- read.table(all.path, sep="\t",quote = "")
  rownames(bed.all) <- bed.all$V4
  bed.close <- read.table(close.path, sep = "\t", quote = "")
  bed.open <- read.table(open.path,sep = "\t", quote = "")
  bed.all$V5 <- 1000
  bed.all$V6 <- "."
  bed.all$V7 <- bed.all$V2; bed.all$V8 <- bed.all$V3
  bed.all$V9 <- "153,153,153"
  bed.all[rownames(bed.all) %in% bed.close$V4, "V9"] <- "166,97,26"
  bed.all[rownames(bed.all) %in% bed.open$V4, "V9"] <- "1,133,113"
  
  write.table(bed.all, file=f.name, sep="\t", col.names = F, row.names = F, quote = F)
  
}

#bedForUCSC(pre="gfpPosi_wtVSg56KD")

sapply(names(dba.center.ls)[1:3], function(x) bedForUCSC(pre = x))


