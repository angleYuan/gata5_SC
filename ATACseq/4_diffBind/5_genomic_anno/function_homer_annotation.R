########## function 1: to use customized setting for promoter annotation and get the distance to TSS information ################
#### this function can takes a very long time if there are many regions, better to run it 
extractAnnoDist <- function(name, pref) {
  
  
  anno <- paste0(pref, name)
  
  print(anno)
  
  anno.df <- read.table(file = anno, sep = "\t", quote = "", header = T )
  #prom.df <- read.table(file = prom, sep = "\t", quote = "", header = F)
  
  colnames(anno.df)[1] <- "ID"
  anno.df <- subset(anno.df, select=c(ID, Chr, Start, End, Annotation, Distance.to.TSS))
  
  rownames(anno.df) <- anno.df$ID
  anno.df$distance_type_1 <- NA
  anno.df$distance_type_2 <- NA
  
  anno.df$New_anno <- sapply(anno.df$Annotation, function(x) {
    category <- unlist(strsplit(as.character(x), "[()]"))[1]
    category <- trimws(category)
    
    
    return(category)
  })
  
  for(i in rownames(anno.df)) {
    if(anno.df[i, "Distance.to.TSS"] < 1000 & anno.df[i, "Distance.to.TSS"] > -3000) {
      anno.df[i, "New_anno"] <- "promoter"
    }
    
    if(anno.df[i, "Distance.to.TSS"] < 3000 & anno.df[i, "Distance.to.TSS"] >= 0) {
      anno.df[i, "distance_type_1"] <- "0_+3kb"
      anno.df[i, "distance_type_2"] <- "<3kb"
    } else if(anno.df[i, "Distance.to.TSS"] <10000 & anno.df[i, "Distance.to.TSS"] >= 3000){
      anno.df[i, "distance_type_1"] <- "+3_+10kb"
      anno.df[i, "distance_type_2"] <- "3-10kb"
    } else if(anno.df[i, "Distance.to.TSS"] < 50000 & anno.df[i, "Distance.to.TSS"] >= 10000) {
      anno.df[i, "distance_type_1"] <- "+10_+50kb"
      anno.df[i, "distance_type_2"] <- "10-50kb"
    } else if(anno.df[i, "Distance.to.TSS"] < 100000 & anno.df[i, "Distance.to.TSS"] >= 50000) {
      anno.df[i, "distance_type_1"] <- "+50_+100kb"
      anno.df[i, "distance_type_2"] <- "50-100kb"
    } else if(anno.df[i, "Distance.to.TSS"] >= 100000) {
      anno.df[i, "distance_type_1"] <- ">+100kb"
      anno.df[i, "distance_type_2"] <- ">100kb"
    } else if(anno.df[i, "Distance.to.TSS"] < 0 & anno.df[i, "Distance.to.TSS"] >= -3000) {
      anno.df[i, "distance_type_1"] <- "-3_0kb"
      anno.df[i, "distance_type_2"] <- "<3kb"
    } else if(anno.df[i, "Distance.to.TSS"] < -3000 & anno.df[i, "Distance.to.TSS"] >= -10000) {
      anno.df[i, "distance_type_1"] <- "-10_-3kb"
      anno.df[i, "distance_type_2"] <- "3-10kb"
    } else if(anno.df[i, "Distance.to.TSS"] < -10000 & anno.df[i, "Distance.to.TSS"] >= -50000) {
      anno.df[i, "distance_type_1"] <- "-50_-10kb"
      anno.df[i, "distance_type_2"] <- "10-50kb"
    } else if(anno.df[i, "Distance.to.TSS"] < -50000 & anno.df[i, "Distance.to.TSS"] >= -100000) {
      anno.df[i, "distance_type_1"] <- "-100_-50kb"
      anno.df[i, "distance_type_2"] <- "50-100kb"
    } else if(anno.df[i, "Distance.to.TSS"] < -100000) {
      anno.df[i, "distance_type_1"] <- "<-100kb"
      anno.df[i, "distance_type_2"] <- ">100kb"
    }
    
    
  }
  
  anno.df$distance_type_1 <- factor(anno.df$distance_type_1, 
                                    levels = c("<-100kb", "-100_-50kb", "-50_-10kb", 
                                               "-10_-3kb", "-3_0kb",  
                                               "0_+3kb", "+3_+10kb",
                                               "+10_+50kb", "+50_+100kb", ">+100kb"))
  
  anno.df$distance_type_2 <- factor(anno.df$distance_type_2, 
                                    levels = c(">100kb", "50-100kb", "10-50kb", "3-10kb", "<3kb"))
  
  
  
  #print(head(anno.df))
  #print(head(prom.df))
  #print(as.data.frame(table(anno.df$New_anno)))
  #print(dim(anno.df))
  #df <- as.data.frame(table(anno.df$New_anno))
  return(anno.df)
}
########## function 2: combine the list into a dataframe for plotting #######################################
ls2df <- function(ls) {
  ls <- lapply(names(ls), function(x) {
    df <- ls[[x]]
    colnames(df)[2] <- x
    return(df)
    
  })
  
  #merge.df <- do.call(function(...) merge(..., by="Var1", all=T), ls)
  merge.df <- Reduce(function(...) merge(..., by="Var1", all=T), ls)
  merge.df[is.na(merge.df)] <- 0
  rownames(merge.df) <- merge.df$Var1
  merge.df$Var1 <- NULL
  merge.df["sum", ] <- c(sum(as.numeric(merge.df$all)), sum(as.numeric(merge.df[, "close"])),
                         sum(as.numeric(merge.df[, "open"])), sum(as.numeric(merge.df[, "shared"])))
  
  for (f in colnames(merge.df)[1:4]) {
    new_f <- paste0(f, ".percent")
    merge.df[, new_f] <- sapply(merge.df[, f], 
                                function(x) as.numeric(x)/as.numeric(merge.df["sum", f]))
    
  }
  
  print(merge.df)
  return(merge.df)
  #return(merge.df[-dim(merge.df)[1], 5:8])
}

########## function 3: customized color for ggplotting ###################################
choose_color <- function(ncol) {
  if(ncol==7) {
    colorsToUse <- rev(c("#a6cee3","#1f78b4", "#b2df8a", "#e31a1c", "#fdbf6f",
                         "#ffff99", "#b15928"))
  }
  if(ncol==8) {
    colorsToUse <- rev(c("#a6cee3","#1f78b4", "#b2df8a", "#e31a1c", "#fdbf6f",
                         "#6a3d9a", "#ffff99", "#b15928"))
    
  }
  if(ncol==5){
    colorsToUse <- c("#F0F9E8","#BAE4BC","#7BCCC4","#43A2CA","#0868AC")
  }
  if(ncol==10){
    colorsToUse <- c("#A50026","#D73027", "#F46D43", "#FDAE61", "#FEE090",
                     "#E0F3F8","#ABD9E9","#74ADD1","#4575B4","#313695")
  }
  #print(colorsToUse)
  return(colorsToUse)
}

########## function 4:  plot the percentage results for genomic annotation ################
df2plot <- function(df, w=5, h=3.5, fname="", cat_level) {
  library(reshape2)
  library(ggplot2)
  df <- df[-dim(df)[1], 5:8]
  m.df <- melt(t(df), measure.vars=rownames(df))
  m.df$Var2 <- factor(m.df$Var2, levels = cat_level)
  m.df$Var1 <- factor(m.df$Var1, levels = rev(levels(m.df$Var1)))
  print(head(m.df))
  
  num <- length(cat_level)
  
  p <- ggplot(m.df)+
    geom_bar(aes(x=Var1, y=value, fill=Var2),stat="identity")+
    scale_fill_manual(values = choose_color(num))+
    guides(fill=guide_legend(reverse = T))+coord_flip()+theme_bw()+
    ylab("percentage")+xlab(NULL)+labs(fill="Feature")
  
  pdf(fname, w, h, useDingbats = F)
  print(p)
  dev.off()
  
  
  return(m.df)
  
  
}

########## function 5: plot the percentage result for distance to TSS bind ###############
df2plotTSS <- function(df, w=5, h=3.5, fname="", cat_level) {
  library(reshape2)
  library(ggplot2)
  df <- df[-dim(df)[1], 5:8]
  m.df <- melt(t(df), measure.vars=rownames(df))
  
  m.df$Var1 <- factor(m.df$Var1, levels = rev(levels(m.df$Var1)))
  m.df$Var3 <- sapply(m.df$Var2, function(x) {
    if(x=="0_+3kb"| x=="-3_0kb"){
      la <- "<3kb"
    } else if(x == "+3_+10kb" | x =="-10_-3kb") {
      la <- "3-10kb"
    } else if(x == "+10_+50kb" | x == "-50_-10kb") {
      la <- "10-50kb"
    } else if(x=="+50_+100kb" | x == "-100_-50kb") {
      la <- "50-100kb"
    } else if(x==">+100kb" | x == "<-100kb"){
      la <- ">100kb"
    }
    
    return(la)
  })
  
  m.df$Var3 <- factor(m.df$Var3, levels=c(">100kb", "50-100kb", "10-50kb", "3-10kb", "<3kb"))
  m.df$Var2 <- factor(m.df$Var2, levels =cat_level)
  print(head(m.df))
  
  
  
  num <- length(cat_level)
  #num <- length(levels(m.df$Var2))
  
  p <- ggplot(m.df)+
    geom_bar(aes(x=Var1, y=value, fill=Var2),stat="identity")+
    scale_fill_manual(values = choose_color(num))+
    guides(fill=guide_legend(reverse = T))+coord_flip()+theme_bw()+
    ylab("percentage")+xlab(NULL)+labs(fill="Feature")
  
  pdf(fname, w, h, useDingbats = F)
  print(p)
  dev.off()
  
  
  return(m.df)
  
  
}

########## function 6: conduct fisher exact test of genomic annotation ##################
### given two df (output from ls2df), two prefix (all, close, open or shared)
twodf2fisher <- function(df1, df2, pre1, pre2, adj_method="bonferroni") {
  name1 <- deparse(substitute(df1)); name1 <- unlist(strsplit(name1, "[.]"))[1]; name1 <- paste0(name1,"." ,pre1)
  name2 <- deparse(substitute(df2)); name2 <- unlist(strsplit(name2, "[.]"))[1]; name2 <- paste0(name2, ".",pre2)
  #percent1 <- paste0(pre1, ".percent")
  #percent2 <- paste0(pre2, ".percent")
  #df1 <- subset(df1, select=c(pre1, percent1))
  #df2 <- subset(df2, select=c(pre2, percent2))
  
  collapse_category <- function(df, pre) {
    percent <- paste0(pre, ".percent")
    df <- subset(df, select=c(pre, percent))
    
    all_exon <- sum(df[rownames(df) %in% c("3' UTR", "5' UTR", "exon", "non-coding"), pre])
    all_exon.per <- sum(df[rownames(df) %in% c("3' UTR", "5' UTR", "exon", "non-coding"), percent])
    
    df["all_exon", ] <- c(all_exon, all_exon.per)
    df <- df[c("all_exon","Intergenic", "intron", "promoter", "TTS"), ]
    return(df)
    
  }
  
  df1 <- collapse_category(df1, pre1)
  df2 <- collapse_category(df2, pre2)
  
  library(rcompanion)
  
  #print(name1); print(name2)
  
  stat.ls <- lapply(rownames(df1), function(anno) {
    df <- data.frame(V1=rep(0, 2),
                     V2=rep(0, 2))
    rownames(df) <- c(anno, "not")
    colnames(df) <- c(name1, name2)
    
    
    #sum1 <- df1[anno, "num"]*100/df1[anno, "percent"]
    #sum2 <- df2[anno, "num"]*100/df2[anno, "percent"]
    
    df[anno, name1] <- df1[[anno, pre1]]; df["not", name1] <- sum(df1[!(rownames(df1) %in% anno), pre1])
    df[anno, name2] <- df2[[anno, pre2]]; df["not", name2] <- sum(df2[!(rownames(df2) %in% anno), pre2])
    
    print(df)
    
    stat_res <- pairwiseNominalIndependence(as.matrix(df), compare = "column")
    
    return(stat_res)
  })
  names(stat.ls) <- rownames(df1)
  
  stats_df <- do.call(rbind, stat.ls)
  stats_df$p.adj.Fisher <- p.adjust(stats_df$p.Fisher, method = adj_method)
  stats_df$p.adj.Chisq <- p.adjust(stats_df$p.Chisq, method = adj_method)
  stats_df$p.adj.Gtest <- p.adjust(stats_df$p.Gtest, method = adj_method)
  
  
  
  stats_df$fish_pass <- sapply(stats_df$p.adj.Fisher, 
                               function(x) {
                                 if (x < 0.05 ) return("Y")
                                 else{return("N")}
                               })
  stats_df$chiseq_pass <- sapply(stats_df$p.adj.Chisq, 
                                 function(x) {
                                   if (x < 0.05 ) return("Y")
                                   else{return("N")}
                                 })
  stats_df$Gtest_pass <- sapply(stats_df$p.adj.Gtest, 
                                function(x) {
                                  if (x < 0.05 ) return("Y")
                                  else{return("N")}
                                })
  stats_df <- stats_df[order(stats_df$p.adj.Fisher, decreasing = F),]
  
  
  print(df1); print(df2)
  print(stats_df)
  return(stats_df)
  
}

######### function 7: fisher exact test on distance to TSS bins #######################
twodf2fisher_select <- function(df1, df2, pre1, pre2, adj_method="bonferroni") {
  
  
  name1 <- deparse(substitute(df1)); name1 <- unlist(strsplit(name1, "[.]"))[1]; name1 <- paste0(name1,"." ,pre1)
  name2 <- deparse(substitute(df2)); name2 <- unlist(strsplit(name2, "[.]"))[1]; name2 <- paste0(name2, ".",pre2)
  
  #name1 <- "g56.close"
  #name2 <- "g56.open"
  
  collapse_category <- function(df, pre) {
    
    #pre <- "close"
    #df <- g56.disTSS.df
    
    percent <- paste0(pre, ".percent")
    df <- subset(df, select=c(pre, percent))
    
    pro_enh <- sum(df[rownames(df) %in% c("-50_-10kb", "+10_+50kb", "-10_-3kb", "+3_+10kb"), pre])
    pro_enh.per <- sum(df[rownames(df) %in% c("-50_-10kb", "+10_+50kb", "-10_-3kb", "+3_+10kb"), percent])
    
    dis_enh <- sum(df[rownames(df) %in% c("<-100kb", ">+100kb"), pre])
    dis_enh.per <- sum(df[rownames(df) %in% c("<-100kb", ">+100kb"), percent])
    
    b0_3 <- sum(df[rownames(df) %in% c("-3_0kb", "0_+3kb"), pre])
    b0_3.per <- sum(df[rownames(df) %in% c("-3_0kb", "0_+3kb"), percent])
    
    b50_100 <- sum(df[rownames(df) %in% c("-100_-50kb","+50_+100kb"), pre])
    b50_100.per <- sum(df[rownames(df) %in% c("-100_-50kb","+50_+100kb"), percent])
    
    df["pro_enh", ] <- c(pro_enh, pro_enh.per)
    df["dis_enh", ] <- c(dis_enh, dis_enh.per)
    df["b50_100", ] <- c(b50_100, b50_100.per)
    df["b0_3", ] <- c(b0_3, b0_3.per)
    
    df <- df[c("pro_enh", "dis_enh", "b50_100", "b0_3"), ]
    
    return(df)
    
  }
  
  df1 <- collapse_category(df1, pre1)
  df2 <- collapse_category(df2, pre2)
  
  library(rcompanion)
  
  #print(name1); print(name2)
  
  stat.ls <- lapply(rownames(df1), function(anno) {
    df <- data.frame(V1=rep(0, 2),
                     V2=rep(0, 2))
    rownames(df) <- c(anno, "not")
    colnames(df) <- c(name1, name2)
    
    
    #sum1 <- df1[anno, "num"]*100/df1[anno, "percent"]
    #sum2 <- df2[anno, "num"]*100/df2[anno, "percent"]
    
    df[anno, name1] <- df1[[anno, pre1]]; df["not", name1] <- sum(df1[!(rownames(df1) %in% anno), pre1])
    df[anno, name2] <- df2[[anno, pre2]]; df["not", name2] <- sum(df2[!(rownames(df2) %in% anno), pre2])
    
    print(df)
    
    stat_res <- pairwiseNominalIndependence(as.matrix(df), compare = "column")
    
    return(stat_res)
  })
  names(stat.ls) <- rownames(df1)
  
  stats_df <- do.call(rbind, stat.ls)
  stats_df$p.adj.Fisher <- p.adjust(stats_df$p.Fisher, method = adj_method)
  stats_df$p.adj.Chisq <- p.adjust(stats_df$p.Chisq, method = adj_method)
  stats_df$p.adj.Gtest <- p.adjust(stats_df$p.Gtest, method = adj_method)
  
  
  
  stats_df$fish_pass <- sapply(stats_df$p.adj.Fisher, 
                               function(x) {
                                 if (x < 0.05 ) return("Y")
                                 else{return("N")}
                               })
  stats_df$chiseq_pass <- sapply(stats_df$p.adj.Chisq, 
                                 function(x) {
                                   if (x < 0.05 ) return("Y")
                                   else{return("N")}
                                 })
  stats_df$Gtest_pass <- sapply(stats_df$p.adj.Gtest, 
                                function(x) {
                                  if (x < 0.05 ) return("Y")
                                  else{return("N")}
                                })
  stats_df <- stats_df[order(stats_df$p.adj.Fisher, decreasing = F),]
  
  
  print(df1); print(df2)
  print(stats_df)
  return(stats_df)
  
}

########## depcreated function ########################
changeAnno <- function(name, pref) {
  
  
  anno <- paste0(pref, name)
  
  print(anno)
  
  anno.df <- read.table(file = anno, sep = "\t", quote = "", header = T )
  #prom.df <- read.table(file = prom, sep = "\t", quote = "", header = F)
  
  colnames(anno.df)[1] <- "ID"
  anno.df <- subset(anno.df, select=c(ID, Chr, Start, End, Annotation, Distance.to.TSS))
  #anno.df$ID <- sapply(rownames(anno.df), function(x) {
  #  
  #  start <- as.numeric(anno.df[x, "Start"])-1
  #  id <- paste0(anno.df[x, "Chr"], ":", start, "-", anno.df[x, "End"])
  #  return(id)
  #  
  #  #paste0(anno.df$Chr,":", anno.df$Start,"-", anno.df$End)
  #})
  rownames(anno.df) <- anno.df$ID
  
  anno.df$New_anno <- sapply(anno.df$Annotation, function(x) {
    category <- unlist(strsplit(as.character(x), "[()]"))[1]
    category <- trimws(category)
    
    
    return(category)
  })
  
  anno.df$New_anno <- sapply(rownames(anno.df), function(i){
    if(anno.df[i, "Distance.to.TSS"] < 1000 & anno.df[i, "Distance.to.TSS"] > -3000) {
      category <- "promoter"
    } else {
      category <- anno.df[i, "New_anno"]
    }
    
    return(category)
  })
  
  
  
  
  #print(head(anno.df))
  #print(head(prom.df))
  print(as.data.frame(table(anno.df$New_anno)))
  #print(dim(anno.df))
  df <- as.data.frame(table(anno.df$New_anno))
  return(df)
}