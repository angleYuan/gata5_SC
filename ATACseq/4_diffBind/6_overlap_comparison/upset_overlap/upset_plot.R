library(UpSetR)


makeUpsetPlots(path = "../all_interested/", f.n="all.upset.pdf")

### not sure why looking at a subset doesn't work here
makeUpsetPlots(path = "../all_interested/", f.n="gfp_g56.upset.pdf", 
               sets = c("gfpPosi_wtVSg56KD.close", "gfpPosi_wtVSg56KD.open",
                        "wt_gfpPosiVSNega.close", "wt_gfpPosiVSNega.open"))

############ analysis on gfp and gata56KD only ##############################
in.ls2 <- list.files(path = "../gfp_g56_intersect/", pattern = "u.bed$" , full.names = T)
n.ls2 <- list.files(path = "../gfp_g56_intersect/", pattern = "u.bed$" , full.names = F)
names(in.ls2) <- lapply(n.ls2, function(x) {
  v <- unlist(strsplit(x, '[.]'))
  n <- paste0(v[1],".", v[6])
  return(n)
})

in_bed.ls2 <- lapply(in.ls2, bedToVennlist)
names(in_bed.ls2) <- names(in.ls2)
in_bed.ls2 <- Filter(Negate(is.null),in_bed.ls2)

pdf("./upset_gfp_g56KD_6sample.pdf", 7,4.5, useDingbats = F)
upset(fromList(in_bed.ls2), nsets = length(in_bed.ls2), order.by = "freq", number.angles = 45, 
      sets.bar.color = "#A9A9A9",
      #empty.intersections = "on",
      #group.by = "sets",
      mb.ratio = c(0.6, 0.4),
      text.scale = 1.4, keep.order = T)
dev.off()


############# analysis on previous intersect that also include aplnrKD ####################
in.ls <- list.files(path = "../all_interested/", pattern = "u.bed$" , full.names = T)
n.ls <- list.files(path = "../all_interested/", pattern = "u.bed$" , full.names = F)
names(in.ls) <- lapply(n.ls, function(x) {
  v <- unlist(strsplit(x, '[.]'))
  n <- paste0(v[1],".", v[6])
  return(n)
})

in_bed.ls <- lapply(in.ls, bedToVennlist)
names(in_bed.ls) <- names(in.ls)
in_bed.ls <- Filter(Negate(is.null),in_bed.ls)

pdf("./upset_4sample.pdf", 7,5, useDingbats = F)
upset(fromList(in_bed.ls), nsets = length(in_bed.ls), order.by = "freq", number.angles = 45, 
      sets.bar.color = "#A9A9A9",
      #empty.intersections = "on",
      #group.by = "sets",
      mb.ratio = c(0.65, 0.35),
      text.scale = 1.4, keep.order = T, sets = c("gfpPosi_wtVSg56KD.close", "gfpPosi_wtVSg56KD.open",
                                                 "wt_gfpPosiVSNega.close", "wt_gfpPosiVSNega.open"))
dev.off()

library(VennDiagram)
venn.diagram(list("gfpPosi_wtVSg56KD.close"=in_bed.ls$gfpPosi_wtVSg56KD.close, 
                  "gfpPosi_wtVSg56KD.open"=in_bed.ls$gfpPosi_wtVSg56KD.open,
                  "wt_gfpPosiVSNega.close"=in_bed.ls$wt_gfpPosiVSNega.close,
                  "wt_gfpPosiVSNega.open"=in_bed.ls$wt_gfpPosiVSNega.open),
             filename = "./venn_4sample.png", 
             fill = c("red", "green", "blue", "yellow"),
             #alpha = c(0.5, 0.5), 
             #cex = 2,cat.fontface = 4,
             #lty =2, fontfamily =3,
             imagetype = "png")
