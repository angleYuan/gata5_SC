######### function 1: read in great output.tsv #################
readInGreatResult <- function(path, filter_by=list("BinomFdrQ"=0.05,
                                                "RegionFoldEnrich"=2,
                                                "HyperFdrQ"=0.05)) {
  #path <- ls[[1]]
  #filter_by <- list("BinomFdrQ"=0.05,
  #                  "RegionFoldEnrich"=2,
  #                  "HyperFdrQ"=0.05)
  
  df <- read.table(file = path, quote = "", sep = "\t")
  colnames(df) <- c("Ontology", "ID", "Desc", "BinomRank",
                    "BinomP", "BinomBonfP", "BinomFdrQ", "RegionFoldEnrich",
                    "ExpRegions", "ObsRegions", "GenomeFrac", "SetCov", 
                    "HyperRank", "HyperP", "HyperBonfP", "HyperFdrQ", 
                    "GeneFoldEnrich", "ExpGenes", "ObsGenes", "TotalGenes",
                    "GeneSetCov", "TermCov", "Regions", "Genes")
  

  df <- df[df$BinomFdrQ < filter_by$BinomFdrQ & df$RegionFoldEnrich > filter_by$RegionFoldEnrich & df$HyperFdrQ < filter_by$HyperFdrQ, ]
  print(path)
  
  df[, c("BinomBonfP", "ExpRegions", "ObsRegions", "GenomeFrac", 
         "SetCov","HyperBonfP", "ExpGenes", "ObsGenes",
         "Regions", "Genes")] <- NULL
  
  return(df)
}

######### function 2: retrieve the top N terms from certain ontology to plot ########
dotplotTopOntology <- function(df, topN=10, ontology=c("GO Biological Process", "GO Molecular Function",
                                                "InterPro", "Wiki Pathways"),
                               sort_by="BinomFdrQ", w=6,h=NULL, pre, size_break=NULL) {
  library(RColorBrewer)
  library(ggplot2)
  #df <- df.ls[[2]]
  #ontology <- c("GO Biological Process", "GO Molecular Function","InterPro", "Wiki Pathways")
  #sort_by <- "BinomFdrQ"
  #topN <- 10
  
  topN_df.ls <- lapply(ontology, function(x) {
    #x <- ontology[2]
    topN.df <- df[df$Ontology == x, ]
    topN.df <- topN.df[order(topN.df[, sort_by], decreasing = F), ]
    topN.df <- topN.df[1:min(topN, dim(topN.df)[1]), ]
    return(topN.df)
  })
  names(topN_df.ls) <- ontology
  
  topN_df <- Reduce(rbind, topN_df.ls)
  topN_df <- na.omit(topN_df)
  
  topN_df$Desc <- factor(topN_df$Desc, levels = rev(topN_df$Desc))
  
  color_breaks <- c("GO Biological Process"="#66C2A5",
                    "GO Molecular Function"="#FC8D62",
                    "InterPro"="#8DA0CB",
                    "Wiki Pathways"="#E78AC3")
  color_breaks <- color_breaks[as.character(unique(topN_df$Ontology))]
  
  if(is.null(size_break)){
    p <- ggplot(topN_df, aes(Desc, -log10(BinomFdrQ)))+
      geom_bar(aes(fill=Ontology), stat = "identity")+
      geom_point(aes(size=RegionFoldEnrich), shape=8)+
      coord_flip()+theme_bw()+
      #scale_fill_brewer(palette = "Set2")
      scale_fill_manual(values = color_breaks)
    
    p.name <- paste0(pre, ".topN_", topN, ".sby_", sort_by,
                     ".sizeDefault.pdf")
    
  } else {
    p <- ggplot(topN_df, aes(Desc, -log10(BinomFdrQ)))+
      geom_bar(aes(fill=Ontology), stat = "identity")+
      geom_point(aes(size=RegionFoldEnrich), shape=8)+
      coord_flip()+theme_bw()+
      #scale_fill_brewer(palette = "Set2")
      scale_fill_manual(values = color_breaks)+
      scale_size_continuous(breaks = size_break)
    
    p.name <- paste0(pre, ".topN_", topN, ".sby_", sort_by,
                     ".sizeAdjust.pdf")
  }
  
  if(is.null(h)) {
    h <- round(dim(topN_df)[1]/4)+1
  }
  
  pdf(p.name, w, h)
  print(p)
  dev.off()
  
}
  