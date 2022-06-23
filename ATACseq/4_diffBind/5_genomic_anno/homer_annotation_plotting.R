library(ggplot2)
library(rcompanion)
library(reshape2)



ls <- list("all"="all.anno", "close"="close.anno",
           "open"="open.anno", "shared"="shared.anno")

anno.df <- read.table("homer_results/gfpPosi_wtVSg56KD.edgeRGLM.pV_0.05.fold_1.open.anno", sep = "\t", quote = "", header = T)

### use the customized promoter criteria for annotation


###################################### Prepare the homer output ########################################
g56.anno.ls <- lapply(ls, function(x) extractAnnoDist(x, pref = "homer_results/gfpPosi_wtVSg56KD.edgeRGLM.pV_0.05.fold_1."))
gfp.anno.ls <- lapply(ls, function(x) extractAnnoDist(x, pref = "homer_results/wt_gfpPosiVSNega.edgeRGLM.pV_0.05.fold_1."))

#g56.ls <- lapply(ls, function(x) changeAnno(x, pref = "homer_results/gfpPosi_wtVSg56KD.edgeRGLM.pV_0.05.fold_1."))
#gfp.ls <- lapply(ls, function(x) changeAnno(x, pref = "homer_results/wt_gfpPosiVSNega.edgeRGLM.pV_0.05.fold_1."))

##################################### to plot genomic annotation  ########################################
g56.genomicA.ls <- lapply(g56.anno.ls, function(x) as.data.frame(table(x$New_anno)))
gfp.genomicA.ls <- lapply(gfp.anno.ls, function(x) as.data.frame(table(x$New_anno)))

## see some non-coding annotation, which are the regions mapped to the exon of non-coding genes
g56.genomicA.df <- ls2df(g56.genomicA.ls)
gfp.genomicA.df <- ls2df(gfp.genomicA.ls)

g56.m.genomicA.df <- df2plot(df = g56.genomicA.df, fname = "g56_homer_anno.pdf",
                     cat_level = rev(c("promoter", "TTS", "exon", "5' UTR", "3' UTR","non-coding","intron", "Intergenic")))
gfp.m.genomicA.df <- df2plot(df = gfp.genomicA.df, fname = "gfp_homer_anno.pdf",
                     cat_level = rev(c("promoter", "TTS", "exon", "5' UTR", "3' UTR","non-coding","intron", "Intergenic")))

###################################### to plot the distance to TSS ##########################################
g56.disTSS.ls <- lapply(g56.anno.ls, function(x) as.data.frame(table(x$distance_type_1)))
gfp.disTSS.ls <- lapply(gfp.anno.ls, function(x) as.data.frame(table(x$distance_type_1)))

g56.disTSS.df <- ls2df(g56.disTSS.ls)
gfp.disTSS.df <- ls2df(gfp.disTSS.ls)

#g56.m.disTSS.df <- df2plotTSS(df = g56.disTSS.df, fname = "g56_disTSS_anno.pdf",
#                    cat_level = rev(c("<3kb", "3-10kb", "10-50kb", "50-100kb", ">100kb")))
g56.m.disTSS.df <- df2plotTSS(df = g56.disTSS.df, fname = "g56_disTSS_anno.pdf",
                              cat_level = rev(c("<-100kb", "-100_-50kb", "-50_-10kb", 
                                                "-10_-3kb", "-3_0kb",  
                                                "0_+3kb", "+3_+10kb",
                                                "+10_+50kb", "+50_+100kb", ">+100kb")))
gfp.m.disTSS.df <- df2plotTSS(df = gfp.disTSS.df, fname = "gfp_disTSS_anno.pdf",
                              cat_level = rev(c("<-100kb", "-100_-50kb", "-50_-10kb", 
                                                "-10_-3kb", "-3_0kb",  
                                                "0_+3kb", "+3_+10kb",
                                                "+10_+50kb", "+50_+100kb", ">+100kb")))
################################## to conduct statistic on genomic annotation ############################
g56.shared_VS_g56.close <- twodf2fisher(g56.genomicA.df, g56.genomicA.df, "shared", "close")
g56.shared_VS_g56.open <- twodf2fisher(g56.genomicA.df, g56.genomicA.df, "shared", "open")
g56.close_VS_g56.open <- twodf2fisher(g56.genomicA.df, g56.genomicA.df, "close", "open")

gfp.shared_VS_gfp.close <- twodf2fisher(gfp.genomicA.df, gfp.genomicA.df, "shared", "close")
gfp.shared_VS_gfp.open <- twodf2fisher(gfp.genomicA.df, gfp.genomicA.df, "shared", "open")
gfp.close_VS_gfp.open <- twodf2fisher(gfp.genomicA.df, gfp.genomicA.df, "close", "open")


################################ conduct statistic on distances to TSS ##################################
#wilcox.test(x=g56.anno.ls[[2]]$Distance.to.TSS, y=g56.anno.ls[[3]]$Distance.to.TSS)
### KS would be a fairly standard to use for a two-sample distribution comparison
ks.test(x=g56.anno.ls[[2]]$Distance.to.TSS, y=g56.anno.ls[[3]]$Distance.to.TSS)

#Two-sample Kolmogorov-Smirnov test

#data:  g56.anno.ls[[2]]$Distance.to.TSS and g56.anno.ls[[3]]$Distance.to.TSS
#D = 0.12036, p-value = 0.0001952
#alternative hypothesis: two-sided

ks.test(x=g56.anno.ls[[2]]$Distance.to.TSS, y=g56.anno.ls[[4]]$Distance.to.TSS)

#Two-sample Kolmogorov-Smirnov test

#data:  g56.anno.ls[[2]]$Distance.to.TSS and g56.anno.ls[[4]]$Distance.to.TSS
#D = 0.11657, p-value < 2.2e-16
#alternative hypothesis: two-sided

ks.test(x=g56.anno.ls[[3]]$Distance.to.TSS, y=g56.anno.ls[[4]]$Distance.to.TSS)
#data:  g56.anno.ls[[3]]$Distance.to.TSS and g56.anno.ls[[4]]$Distance.to.TSS
#D = 0.078466, p-value = 0.02081
#alternative hypothesis: two-sided

ks.test(x=gfp.anno.ls[[2]]$Distance.to.TSS, y=gfp.anno.ls[[3]]$Distance.to.TSS)

#data:  gfp.anno.ls[[2]]$Distance.to.TSS and gfp.anno.ls[[3]]$Distance.to.TSS
#D = 0.013061, p-value = 0.5126
#alternative hypothesis: two-sided

ks.test(x=gfp.anno.ls[[2]]$Distance.to.TSS, y=gfp.anno.ls[[4]]$Distance.to.TSS)

#data:  gfp.anno.ls[[2]]$Distance.to.TSS and gfp.anno.ls[[4]]$Distance.to.TSS
#D = 0.089782, p-value < 2.2e-16
#alternative hypothesis: two-sided

ks.test(x=gfp.anno.ls[[3]]$Distance.to.TSS, y=gfp.anno.ls[[4]]$Distance.to.TSS)
#data:  gfp.anno.ls[[3]]$Distance.to.TSS and gfp.anno.ls[[4]]$Distance.to.TSS
#D = 0.084767, p-value < 2.2e-16
#alternative hypothesis: two-sided

### we can also try to conduct fisher exact test to 
g56TSS.close_VS_g56TSS.open <- twodf2fisher_select(g56.disTSS.df, g56.disTSS.df, "close", "open")
g56TSS.close_VS_gfpTSS.close <- twodf2fisher_select(g56.disTSS.df, gfp.disTSS.df, "close", "close")
g56TSS.close_VS_gfpTSS.open <- twodf2fisher_select(g56.disTSS.df, gfp.disTSS.df, "close", "open")

g56TSS.close_VS_g56TSS.shared <- twodf2fisher_select(g56.disTSS.df, g56.disTSS.df, "close", "shared")


gfpTSS.close_VS_gfpTSS.open <- twodf2fisher_select(gfp.disTSS.df, gfp.disTSS.df, "close", "open")


############################ specifically plot the proximal and distal enhancer percentage of DARs ################
proDis_df <- data.frame("closed_DARs"=c(0.32988764, 0.42876404),
                        "open_DARs"=c(0.16935484, 0.59677419),
                        "GFPposi"=c(0.24918899, 0.49309143),
                        "GFPnega"=c(0.25010036, 0.49270708))
rownames(proDis_df) <- c("distal", "proximal")

proDis_df.m <- melt(t(proDis_df), measure.vars=rownames(proDis_df))

pdf("distal_proximal.pdf", 5,5)
ggplot(proDis_df.m)+
  geom_bar(aes(x=Var2, y=value, fill=Var1),position="dodge", stat="identity")+
  scale_fill_manual(values = c("#F7941D", "#27AAE1", "#2BB673", "#231F20"))+
  #scale_fill_manual(values = choose_color(num))+
  #guides(fill=guide_legend(reverse = T))+
  theme_bw()+
  ylab("percentage")+xlab(NULL)+labs(fill="Feature")
dev.off()

############################## make violin plots regarding TSS distance ##########################
g56.disTSS_v.df <- lapply(names(g56.anno.ls)[2:4], function(x){
  df <- as.data.frame(g56.anno.ls[[x]]$Distance.to.TSS)
  df$type <- x
  colnames(df) <- c("distanceTSS", "type")
  df$distanceTSS <- abs(df$distanceTSS)
  
  return(df)
})
g56.disTSS_v.df <- Reduce(rbind, g56.disTSS_v.df)

pdf("vplot_g56_disTSS.pdf",2.5,3)
ggplot(g56.disTSS_v.df)+theme_bw()+
  geom_violin(aes(x=type, y=distanceTSS))+
  ylab("distance to TSS")+xlab("peak type")+
  ylim(0, 2e6)
dev.off()

gfp.disTSS_v.df <- lapply(names(gfp.anno.ls)[2:4], function(x){
  df <- as.data.frame(gfp.anno.ls[[x]]$Distance.to.TSS)
  df$type <- x
  colnames(df) <- c("distanceTSS", "type")
  df$distanceTSS <- abs(df$distanceTSS)
  
  return(df)
})
gfp.disTSS_v.df <- Reduce(rbind, gfp.disTSS_v.df)

pdf("vplot_gfp_disTSS.pdf",2.5,3)
ggplot(gfp.disTSS_v.df)+
  geom_violin(aes(x=type, y=distanceTSS))+theme_bw()+
  ylab("distance to TSS")+xlab("peak type")+
  ylim(0, 2e6)
dev.off()
