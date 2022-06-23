library(ggplot2)
source("function_greatPlot.r")


ls <- list.files(path="../", pattern = ".tsv", full.names = T)

names(ls) <- lapply(list.files(path="../", pattern = ".tsv", full.names = F), function(x) {
  name <- unlist(strsplit(x, "[.]"))[c(1,5,6)]
  name <- paste(name, collapse = ".")
  return(name)
})

df.ls <- lapply(ls, function(x) readInGreatResult(path = x))

lapply(names(df.ls), function(x) dotplotTopOntology(df=df.ls[[x]], w=8,
                                                    pre = paste0(x, ".default")))

lapply(names(df.ls), function(x) dotplotTopOntology(df=df.ls[[x]], w=8,topN = 15,
                                                    pre = paste0(x, ".default")))

dotplotTopOntology(df=df.ls[["GREAT_gfpPosi_wtVSg56KD.fold_1.close"]], w=8, topN = 15,
                   pre = paste0("GREAT_gfpPosi_wtVSg56KD.fold_1.close", ".default"),
                   size_break = c("2.1"=2.1,"2.5"=2.5, "3.0"=3.0,"3.5"=3.5))
dotplotTopOntology(df=df.ls[["GREAT_gfpPosi_wtVSg56KD.fold_1.close"]], w=7.5, topN = 10,
                   pre = paste0("GREAT_gfpPosi_wtVSg56KD.fold_1.close", ".default"),
                   size_break = c("2.04"=2.04,"2.5"=2.5, "3.0"=3.0,"3.5"=3.5))

dotplotTopOntology(df=df.ls[["GREAT_gfpPosi_wtVSg56KD.fold_1.open"]], w=7.3, topN = 10,
                   pre = paste0("GREAT_gfpPosi_wtVSg56KD.fold_1.open", ".default"),
                   size_break = c("2.01"=2.01,"3"=3, "4"=4,"5"=5))

dotplotTopOntology(df=df.ls[["GREAT_wt_gfpPosiVSNega.fold_1.close"]], w=7.2, topN = 10,
                   pre = paste0("GREAT_wt_gfpPosiVSNega.fold_1.close", ".default"),
                   size_break = c("2.01"=2.01,"2.5"=2.5, "3.0"=3.0,"3.5"=3.5))

dotplotTopOntology(df=df.ls[["GREAT_wt_gfpPosiVSNega.fold_1.open"]], w=6.8, topN = 10,
                   pre = paste0("GREAT_wt_gfpPosiVSNega.fold_1.open", ".default"),
                   size_break = c("2.04"=2.04,"2.5"=2.5, "3.0"=3.0,"3.5"=3.5))


dotplotTopOntology(df=df.ls[["GREAT_gfpPosi_wtVSg56KD.fold_1.close"]], w=5.5, topN = 10,
                   ontology = "GO Biological Process", 
                   pre = paste0("GREAT_gfpPosi_wtVSg56KD.fold_1.close", ".bioPro"))

dotplotTopOntology(df=df.ls[["GREAT_gfpPosi_wtVSg56KD.fold_1.close"]], w=5.7, topN = 10,
                   ontology = "GO Biological Process", 
                   pre = paste0("GREAT_gfpPosi_wtVSg56KD.fold_1.close", ".bioPro"),
                   size_break = c("2.1"=2.1,"2.4"=2.4, "2.7"=2.7,"3.0"=3.0, "3.3"=3.3))
