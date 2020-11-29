############### 1. Subset a seurat object, find variable genes, conduct PCA  and prepare datasets for Dynverse
## return a dataset ready for dynverse package
subset2Dynverse <- function(seurat_o = e, subset_id, project_n, n_vari_genes = 3000,
                            assay_type, marker_gene) {
  
  print(project_n)
  print(subset_id)
  ss <- subset(x = seurat_o, idents = subset_id)
  #ss <- subset(x=e.SCT, idents = id.ls[[1]])
  ss@project.name <- project_n
  table(ss@active.ident)
  
  #DefaultAssay(ss) <- "RNA"
  
  dir.create(paste0("dynverse_plot/", project_n, ".plot"), showWarnings = F)
  dir.create(paste0("dynverse_plot/", project_n, ".plot", "/analysis_check"), showWarnings = F)
  directory <- paste0("dynverse_plot/", project_n, ".plot", "/analysis_check/")
  
  ### another QC to check if the Counts of each cluster is comparable
  feature.p <- FeatureScatter(object = ss, feature1 = "nCount_RNA", feature2 = "percent.mt"
                              #cells = colnames(ss)[grep("6hpf", colnames(ss))]
  )
  pdf(paste0(directory, "nCount_mt.scatter.pdf"), 5, 5)
  print(feature.p)
  dev.off()
  
  
  ss <- FindVariableFeatures(ss, selection.method = "vst", nfeatures = n_vari_genes)
  print("variable genes finished")
  ss <- RunPCA(object = ss, verbose = T, npcs = 60)
  
  
  ######### format dataset to fit the dynverse analysis
  ### only keep the highly variable genes and gata5, gfp, marker_gene (in case not included)
  #assay_type <- "SCT"; marker_gene <- "nkx2.5"
  genes_k <- union(ss[[assay_type]]@var.features, c("gata5", "gfp", marker_gene))
  t.data <- t(as.matrix(ss[[assay_type]]@data))
  t.data <- t.data[, genes_k]
  t.data <- as(t.data, "dgCMatrix")
  
  t.counts <- t(as.matrix(ss[[assay_type]]@counts))
  t.counts <- t.counts[, genes_k]
  t.counts <- as(t.data, "dgCMatrix")
  
  dataset <- wrap_expression(
    expression = t.data,
    counts = t.counts)
  
  dataset <- add_grouping(
    dataset,
    ss@active.ident#,
    #group_ids = "orig_ident"
  ) 
  
  ### add start_id
  start_ID <- dataset[["cell_ids"]][grep("6hpf", dataset[["cell_ids"]])]
  ### add time points
  tp <- sapply(ss$orig.ident, function(x) as.numeric(unlist(strsplit(x, "\\_|\\h"))[2]))
  
  
  dataset <- dataset %>% add_prior_information(start_id = start_ID, 
                                               timecourse_discrete = tp,
                                               groups_id = ss$orig.ident)
  print("dataset finished")
  
  ## take the longest time to run, so put this step at the end
  ss <- JackStraw(object = ss, num.replicate = 100, dims = 30)
  print("JackStraw finished")
  ss <- ScoreJackStraw(object = ss, dims = 1:30)
  plot_PCnumbers(gbm = ss)
  
  return(dataset)
}

################ 2. Function used in subset2Dynverse for save JackStraw Plot and Elbow Plot ##################
plot_PCnumbers <- function(gbm=gbm, w1=8, h1=6, w2=5, h2=5, j.dim=30,  e.dim=60) {
  project_n <- gbm@project.name
  dir.create(paste0("dynverse_plot/", project_n, ".plot"), showWarnings = F)
  dir.create(paste0("dynverse_plot/", project_n, ".plot", "/analysis_check"), showWarnings = F)
  directory <- paste0("dynverse_plot/", project_n, ".plot", "/analysis_check/")
  
  j.name <- paste0(directory, "jackStrawPlot.dim", j.dim, ".pdf")
  e.name <- paste0(directory, "elbowPlot.dim", e.dim, ".pdf")
  
  p1 <- JackStrawPlot(object = gbm, dims = 1:j.dim)
  
  pdf(j.name, w1, h1, useDingbats = F)
  print(p1)
  dev.off()
  
  p2 <- ElbowPlot(object = gbm, ndims = e.dim)
  pdf(e.name, w2, h2, useDingbats = F)
  print(p2)
  #ElbowPlot(object = gbm, ndims = e.dim)
  dev.off()
  
  
}

################ 3. Function that run slingshot to infer trajectory on dynverse datasets and make plots ########
### return a trajectory model
dyn_slingshot_plot <- function(dataset, n_dataset, marker_gene,
                               ndim = 20L, shrink = 1L,
                               thresh = 0.001, maxit = 10, stretch = 2,
                               smoother = "smooth.spline", 
                               shrink.method = "cosine", cluster_m = "pam",
                               give_priors = c("start_id", "timecourse_discrete"),
                               verbose = T, seed = 1, 
                               w = 12, h=12) {
  
  dir.create(paste0("dynverse_plot/",n_dataset, ".plot"), showWarnings = F)
  directory <- paste0("dynverse_plot/", n_dataset, ".plot/")
  paras <- paste0("slingshot.ndim_",ndim, ".cluM_", cluster_m,".sh_", shrink, ".th_", thresh,
                  ".max_", maxit, ".str_", stretch,
                  ".sm_", smoother, ".shMeth_", shrink.method, ".pdf")
  
  f.name <- paste0(directory, paras)
  print(f.name)
  
  
  model <- infer_trajectory(dataset, seed = seed,
                            ti_slingshot(cluster_method = cluster_m,
                                         ndim = ndim,
                                         shrink = shrink,
                                         thresh = thresh,
                                         maxit = maxit,
                                         stretch = stretch,
                                         smoother = smoother,
                                         shrink.method = shrink.method),
                            give_priors = give_priors, verbose = verbose)
  print("model finished")
  
  
  p <- patchwork::wrap_plots(
    plot_dimred(model, grouping = dataset$prior_information$timecourse_discrete, 
                size_cells = 1) + ggtitle("timepoint"),
    plot_dimred(model, grouping = dataset$grouping, 
                size_cells = 1) + ggtitle("cluster"),
    plot_dimred(model, grouping = group_onto_nearest_milestones(model),
                size_cells = 1) + ggtitle("milestone_grouping"),
    plot_dimred(model, grouping = dataset$prior_information$groups_id, 
                size_cells = 1) + ggtitle("samples"),
    plot_dimred(model, feature_oi = "gata5", size_cells = 1,
                expression_source = dataset) + ggtitle("gata5 expression"),
    plot_dimred(model, feature_oi = "gfp", size_cells = 1,
                expression_source = dataset) + ggtitle("gfp expression"),
    plot_dimred(model, feature_oi = marker_gene, size_cells = 1,
                expression_source = dataset) + ggtitle("marker expression"),
    plot_dimred(model, "pseudotime", size_cells = 1,
                pseudotime = calculate_pseudotime(model)) + ggtitle("Pseudotime")
  )
  
  pdf(f.name, w, h)
  print(p)
  dev.off()
  
  print("plot finished")
  
  return(model)
  
}

############### 4. function to make geom_point plot of gene expression across a trajectory ###################
## return an expression-pseudotime dataframe
geneAcrossPseudoT <- function(model, dataset, n_dataset,ylims=NULL,
                              gene = "gata5",name, w=6.5, h=4.5) {
  dir.create(paste0("dynverse_plot/",n_dataset, ".plot"), showWarnings = F)
  directory <- paste0("dynverse_plot/",n_dataset, ".plot/")
  f_name <- paste0(directory,name, ".exprPseudoT.",gene,".pdf")
  print(f_name)
  
  #dataframe for ggplot2
  df <- data.frame("pseudoT" = calculate_pseudotime(trajectory = model),
                   "gene_exp" = dataset$expression[, gene],
                   "timepoint" = factor(dataset$prior_information$timecourse_discrete,
                                        levels = c("6", "8", "10", "13"))
  )
  
  if(is.null(ylims)) {
    p <- ggplot(df, aes(x = pseudoT, y = gene_exp))+
      geom_point(aes(color=timepoint, alpha=0.5))+
      geom_smooth()+
      #scale_color_brewer(palette = "Set1")+
      scale_color_manual(values = c("6"="#001b8f", "8" = "#06c8f9",
                                    "10" = "#FFBC01", "13" = "#7c000c"))+
      #ylim(...)+
      ylab(gene)+
      ggtitle(n_dataset)+
      theme_cowplot()
  } else {
    
    p <- ggplot(df, aes(x = pseudoT, y = gene_exp))+
      geom_point(aes(color=timepoint, alpha=0.5))+
      geom_smooth()+
      #scale_color_brewer(palette = "Set1")+
      scale_color_manual(values = c("6"="#001b8f", "8" = "#06c8f9",
                                    "10" = "#FFBC01", "13" = "#7c000c"))+
      ylim(ylims)+
      ylab(gene)+
      ggtitle(n_dataset)+
      theme_cowplot()
  }
  

  pdf(f_name, w, h)
  print(p)
  dev.off()
  
  print(p) ## for demonstrating in Rmd
  return(df)
}

############### 5. function to plot the trajectory for each lineage using customized color ##################
plotCustomize <- function(model, dataset, n_dataset, w=4.5, h=4.5, name,
                          colors_C = c("#FFBC01", "#7c000c","#001b8f", "#06c8f9")) {
  p <- plot_dimred(model, grouping = dataset$prior_information$timecourse_discrete,size_cells = 1,
                   size_trajectory = 0.5, 
                   groups = tibble(group_id=unique(dataset$prior_information$timecourse_discrete),
                                   color= colors_C)) + ggtitle(paste0(n_dataset, " color by timepoint"))
  
  dir.create(paste0("dynverse_plot/", n_dataset, ".plot"), showWarnings = F)
  directory <- paste0("dynverse_plot/", n_dataset, ".plot/")
  f.name <- paste0(directory, name, ".pdf")
  
  pdf(f.name, w, h)
  print(p)
  dev.off()
  return(p)
}

############### 6. function to bin the pseudotime at each timepoint and caculate mean gene expression for each bin #####
## return mean gene expression for each quantile bin
quantileBin_byTime <- function(df) {
  #df <- g5.tp.ls[[1]]
  ### split into different timepoint, bin into 10 quantile, calculate the mean g5 expression for each bin
  df.t.ls <- lapply(c(6,8,10,13), function(x) {
    
    df.t <- df[df$timepoint == x, ]
    df.t <- df.t[order(df.t$pseudoT, decreasing = F), ]
    df.t$rank <- c(1:dim(df.t)[1])
    #df.t$sc.pseudoT <- (df.t$pseudoT-min(df.t$pseudoT))*9.9999/(max(df.t$pseudoT)-min(df.t$pseudoT))
    df.t$quantile <- cut(df.t$rank,breaks=10,include.lowest = T, right = F)
    quantile_n <- levels(df.t$quantile)
    ### calculate the mean gata5 expression for each bin
    quantile_ave <- sapply(quantile_n, function(x) {
      ave_g5 <- mean(df.t[df.t$quantile == x, ]$gene_exp)
      return(ave_g5)
    })
    names(quantile_ave) <- sapply(c(1:10), function(q) paste0(x, "h.q",q))
    
    #c("q1", "q2", "q3", "q4", "q5", "q6", "q7", "q8", "q9","q10")
    return(quantile_ave)
  })
  g5.quantile <- do.call(c, df.t.ls)
  return(g5.quantile)
  
}



