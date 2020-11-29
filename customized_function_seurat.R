######################## 1. QC  function ##################################################### 
# ## This function plot the original or filtered number of genes (nGene), number of UMI counts (nUMI) and percent of mitochondria reads (percent.mito) before or after filtering 
# and save the plots in a folder called "analysis_check" (will automatically generate if the folder is not existing). 
# This function makes three plots, a combined violin plot (name1), a combined histogram plot (name2) and a combined scatter plot (name3).
# ## For the sake of demonstration in Rmd, it also returns the combined histogram plot 

plot_nGene_nUMI_mito <- function(gbm=gbm, break_gene=200, break_umi=200, break_mito=200,
                                 w1=15,h1=10, w2=6, h2=8, w3=10, h3=6, raw_data=T) {
  project_n <- gbm@project.name
  dir.create(paste0("seurat_plots/",project_n, ".plot"), showWarnings = F)
  dir.create(paste0("seurat_plots/",project_n, ".plot", "/analysis_check"), showWarnings = F)
  directory <- paste0("seurat_plots/",project_n, ".plot", "/analysis_check/")
  
  if(raw_data==T) {
    name1 <- paste0(directory, "vln.raw_nGene_nUMI_mito.pdf")
    name2 <- paste0(directory, "hist.raw_nGene_nUMI_mito.pdf")
    name3 <- paste0(directory, "scatter.raw_nGene_nUMI_mito.pdf")
  } else{
    name1 <- paste0(directory, "vln.filterd_nGene_nUMI_mito.pdf")
    name2 <- paste0(directory, "hist.filtered_nGene_nUMI_mito.pdf")
    name3 <- paste0(directory, "scatter.filtered_nGene_nUMI_mito.pdf")
  }
  
  
  
  ## To make a combined violin plot 
  p1 <- VlnPlot(object = gbm, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  pdf(name1, w1, h1, useDingbats = F)
  print(p1)
  dev.off()
  
  
  df <- data.frame(nCount_RNA=gbm$nCount_RNA,
                   nFeature_RNA=gbm$nFeature_RNA,
                   percent_mt=gbm$percent.mt)
  
  
  
  ##### To make a combined histogram
  #par(mfrow=c(1,3))
  p2.1 <- ggplot(df, aes(x=nCount_RNA))+
    geom_histogram(bins = break_umi,  color="#FF6666", fill="#FF6666",alpha=0.5)+
    theme_bw()
  p2.2 <- ggplot(df, aes(x=nFeature_RNA))+
    geom_histogram(bins=break_gene,  color="#FF6666", fill="#FF6666",alpha=0.5)+
    theme_bw()
  p2.3 <- ggplot(df, aes(x=percent_mt))+
    geom_histogram(bins = break_mito,  color="#FF6666", fill="#FF6666",alpha=0.5)+
    theme_bw()
  
  pdf(name2, w2, h2, useDingbats = F)
  plot2 <- CombinePlots(plots = list(p2.1, p2.2, p2.3), ncol = 1)
  print(plot2)
  dev.off()
  
  #### To make combined feature plot
  p3.1 <- FeatureScatter(object = gbm, feature1 = "nCount_RNA", feature2 = "percent.mt", group.by = "orig.ident")
  p3.2 <- FeatureScatter(object = gbm, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = "orig.ident")
  plot3 <- CombinePlots(plots = list(p3.1, p3.2))
  
  pdf(name3, w3, h3, useDingbats = F)
  print(plot3)
  dev.off()
  
  return(plot2)
  
}

####################### 2. plot heatmap and save the heatmap #################################
plotHeatmap <- function(gbm=gbm, test.use="MAST", diff.ls = diff.glob.ls,
                        topN=40, cells.use=NULL, w=18, h=topN*2) {
  library(dplyr)
  project_n <- gbm@project.name
  dir.create(paste0("seurat_plots/",project_n, ".plot"), showWarnings = F)
  directory <- paste0("seurat_plots/",project_n, ".plot/")
  
  diff.markers <- diff.ls[[test.use]]
  
  if(is.null(cells.use)){
    p.name <- paste0(directory,test.use, ".hm.allCells.top", topN,".pdf")
    gene.df <- diff.markers %>% group_by(cluster) %>% top_n(topN, avg_logFC)
    #print(head(gene.df))
    
    p <- DoHeatmap(object = gbm, features = gene.df$gene)
    pdf(p.name, w, h, useDingbats = F)
    print(p)
    dev.off()
  } else{
    gbm.use <- SubsetData(object = gbm, ident.use = cells.use)
    cell.name <- paste(cells.use, collapse = "-")
    p.name <- paste0(directory, "hm.Cells", cell.name,".top", topN,".pdf")
    
  }
  print(paste0("plot ", test.use, " heatmap finished"))
}

####################### 3. modified DotPlot that allows split by samples #####################
### dot plot that can split cells by samples while color the dots by expression levels
DotPlot.m <- function (object, assay = NULL, features, cols = c("lightgrey", "blue"), 
                       col.min = -2.5, col.max = 2.5, dot.min = 0, dot.scale = 6, 
                       group.by = NULL, split.by = NULL, scale.by = "radius", scale.min = NA, 
                       scale.max = NA) 
{
  #assay <- assay %||% DefaultAssay(object = object)
  #DefaultAssay(object = object) <- assay
  library(ggplot2)
  library(scales)
  library(cowplot)
  scale.func <- function (name = waiver(), breaks = waiver(), labels = waiver(), 
                          limits = NULL, range = c(1, 6), trans = "identity", guide = "legend") 
  {
    continuous_scale("size", "radius", rescale_pal(range), name = name, 
                     breaks = breaks, labels = labels, limits = limits, trans = trans, 
                     guide = guide)
  }
  PercentAbove <- function(x, threshold) {
    return(length(x = x[x > threshold]) / length(x = x))
  }
  data.features <- FetchData(object = object, vars = features)
  data.features$id <- if (is.null(x = group.by)) {
    Idents(object = object)
  }
  else {
    object[[group.by, drop = TRUE]]
  }
  if (!is.factor(x = data.features$id)) {
    data.features$id <- factor(x = data.features$id)
  }
  id.levels <- levels(x = data.features$id)
  data.features$id <- as.vector(x = data.features$id)
  if (!is.null(x = split.by)) {
    splits <- object[[split.by, drop = TRUE]]
    if (length(x = unique(x = splits)) > length(x = cols)) {
      stop("Not enought colors for the number of groups")
    }
    cols <- cols[1:length(x = unique(x = splits))]
    names(x = cols) <- unique(x = splits)
    data.features$id <- paste(data.features$id, splits, sep = "_")
    unique.splits <- rev(unique(x = splits))
    id.levels <- paste0(rep(x = id.levels, each = length(x = unique.splits)), 
                        "_", rep(x = unique.splits, times = length(x = id.levels)))
  }
  data.plot <- lapply(X = unique(x = data.features$id), FUN = function(ident) {
    data.use <- data.features[data.features$id == ident, 
                              1:(ncol(x = data.features) - 1), drop = FALSE]
    avg.exp <- apply(X = data.use, MARGIN = 2, FUN = function(x) {
      return(mean(x = expm1(x = x)))
    })
    pct.exp <- apply(X = data.use, MARGIN = 2, FUN = PercentAbove, 
                     threshold = 0)
    return(list(avg.exp = avg.exp, pct.exp = pct.exp))
  })
  names(x = data.plot) <- unique(x = data.features$id)
  data.plot <- lapply(X = names(x = data.plot), FUN = function(x) {
    data.use <- as.data.frame(x = data.plot[[x]])
    data.use$features.plot <- rownames(x = data.use)
    data.use$id <- x
    return(data.use)
  })
  data.plot <- do.call(what = "rbind", args = data.plot)
  if (!is.null(x = id.levels)) {
    data.plot$id <- factor(x = data.plot$id, levels = id.levels)
  }
  avg.exp.scaled <- sapply(X = unique(x = data.plot$features.plot), 
                           FUN = function(x) {
                             data.use <- data.plot[data.plot$features.plot == 
                                                     x, "avg.exp"]
                             data.use <- scale(x = data.use)
                             data.use <- MinMax(data = data.use, min = col.min, 
                                                max = col.max)
                             return(data.use)
                           })
  avg.exp.scaled <- as.vector(x = t(x = avg.exp.scaled))
  #if (!is.null(x = split.by)) {
  #  avg.exp.scaled <- as.numeric(x = cut(x = avg.exp.scaled, 
  #                                       breaks = 20))
  #}
  data.plot$avg.exp.scaled <- avg.exp.scaled
  data.plot$features.plot <- factor(x = data.plot$features.plot, 
                                    levels = rev(x = features))
  data.plot$pct.exp[data.plot$pct.exp < dot.min] <- NA
  data.plot$pct.exp <- data.plot$pct.exp * 100
  #if (!is.null(x = split.by)) {
  #  splits.use <- vapply(X = strsplit(x = as.character(x = data.plot$id), 
  #                                    split = "_"), FUN = "[[", FUN.VALUE = character(length = 1L), 
  #                       2)
  #  data.plot$colors <- mapply(FUN = function(color, value) {
  #    return(colorRampPalette(colors = c("grey", color))(20)[value])
  #  }, color = cols[splits.use], value = avg.exp.scaled)
  #}
  #color.by <- ifelse(test = is.null(x = split.by), yes = "avg.exp.scaled", 
  #                   no = "colors")
  color.by <- "avg.exp.scaled"
  if (!is.na(x = scale.min)) {
    data.plot[data.plot$pct.exp < scale.min, "pct.exp"] <- scale.min
  }
  if (!is.na(x = scale.max)) {
    data.plot[data.plot$pct.exp > scale.max, "pct.exp"] <- scale.max
  }
  plot <- ggplot(data = data.plot, mapping = aes_string(x = 'features.plot', y = 'id')) +
    geom_point(mapping = aes_string(size = 'pct.exp', color = color.by)) +
    scale.func(range = c(0, dot.scale), limits = c(scale.min, scale.max)) +
    theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
    guides(size = guide_legend(title = 'Percent Expressed')) +
    labs(
      x = 'Features',
      y = ifelse(test = is.null(x = split.by), yes = 'Identity', no = 'Split Identity')
    ) + 
    theme_cowplot()+ 
    coord_flip()+
    scale_colour_viridis_c(direction = -1)
  #if (!is.null(x = split.by)) {
  #  plot <- plot + scale_color_identity()
  #}
  #else if (length(x = cols) == 1) {
  #  plot <- plot + scale_color_distiller(palette = cols)
  #}
  #else {
  #  plot <- plot + scale_color_gradient(low = cols[1], high = cols[2])
  #}
  #if (is.null(x = split.by)) {
  #  plot <- plot + guides(color = guide_colorbar(title = "Average Expression"))
  #}
  return(plot)
}

###################### 4. modified DoHeatmap that allows split by samples ####################
### heatmap that allows split cells by both clusters and samples/conditions
DoHeatmap.m <- function (object, features = NULL, cells = NULL, group.by = "ident", split.by = "orig.ident",
                                 group.bar = TRUE, group.colors = NULL, disp.min = -2.5, 
                                 disp.max = NULL, slot = "scale.data", assay = NULL, label = TRUE, 
                                 size = 5.5, hjust = 0, angle = 45, raster = TRUE, draw.lines = TRUE, 
                                 lines.width = NULL, group.bar.height = 0.02, combine = TRUE) 
{
  library(scales)
  library(ggplot2)
  library(patchwork)
  
  CheckDots <- function(..., fxns = NULL) {
    args.names <- names(x = list(...))
    if (length(x = list(...)) == 0) {
      return(invisible(x = NULL))
    }
    if (is.null(x = args.names)) {
      stop("No named arguments passed")
    }
    if (length(x = fxns) == 1) {
      fxns <- list(fxns)
    }
    for (f in fxns) {
      if (!(is.character(x = f) || is.function(x = f))) {
        stop("CheckDots only works on characters or functions, not ", class(x = f))
      }
    }
    fxn.args <- suppressWarnings(expr = sapply(
      X = fxns,
      FUN = function(x) {
        x <- tryCatch(
          expr = if (isS3stdGeneric(f = x)) {
            as.character(x = methods(generic.function = x))
          } else {
            x
          },
          error = function(...) {
            return(x)
          }
        )
        x <- if (is.character(x = x)) {
          sapply(X = x, FUN = argsAnywhere, simplify = FALSE, USE.NAMES = TRUE)
        } else if (length(x = x) <= 1) {
          list(x)
        }
        return(sapply(
          X = x,
          FUN = function(f) {
            return(names(x = formals(fun = f)))
          },
          simplify = FALSE,
          USE.NAMES = TRUE
        ))
      },
      simplify = FALSE,
      USE.NAMES = TRUE
    ))
    fxn.args <- unlist(x = fxn.args, recursive = FALSE)
    fxn.null <- vapply(X = fxn.args, FUN = is.null, FUN.VALUE = logical(length = 1L))
    if (all(fxn.null) && !is.null(x = fxns)) {
      stop("None of the functions passed could be found")
    } else if (any(fxn.null)) {
      warning(
        "The following functions passed could not be found: ",
        paste(names(x = which(x = fxn.null)), collapse = ', '),
        call. = FALSE,
        immediate. = TRUE
      )
      fxn.args <- Filter(f = Negate(f = is.null), x = fxn.args)
    }
    dfxns <- vector(mode = 'logical', length = length(x = fxn.args))
    names(x = dfxns) <- names(x = fxn.args)
    for (i in 1:length(x = fxn.args)) {
      dfxns[i] <- any(grepl(pattern = '...', x = fxn.args[[i]], fixed = TRUE))
    }
    if (any(dfxns)) {
      dfxns <- names(x = which(x = dfxns))
      if (any(nchar(x = dfxns) > 0)) {
        fx <- vapply(
          X = Filter(f = nchar, x = dfxns),
          FUN = function(x) {
            if (isS3method(method = x)) {
              x <- unlist(x = strsplit(x = x, split = '\\.'))
              x <- x[length(x = x) - 1L]
            }
            return(x)
          },
          FUN.VALUE = character(length = 1L)
        )
        message(
          "The following functions and any applicable methods accept the dots: ",
          paste(unique(x = fx), collapse = ', ')
        )
        if (any(nchar(x = dfxns) < 1)) {
          message(
            "In addition, there is/are ",
            length(x = Filter(f = Negate(f = nchar), x = dfxns)),
            " other function(s) that accept(s) the dots"
          )
        }
      } else {
        message("There is/are ", length(x = dfxns), 'function(s) that accept(s) the dots')
      }
    } else {
      unused <- Filter(
        f = function(x) {
          return(!x %in% unlist(x = fxn.args))
        },
        x = args.names
      )
      if (length(x = unused) > 0) {
        msg <- paste0(
          "The following arguments are not used: ",
          paste(unused, collapse = ', ')
        )
        switch(
          EXPR = getOption(x = "Seurat.checkdots"),
          "warn" = warning(msg, call. = FALSE, immediate. = TRUE),
          "stop" = stop(msg),
          "silent" = NULL,
          stop("Invalid Seurat.checkdots option. Please choose one of warn, stop, silent")
        )
        unused.hints <- sapply(X = unused, FUN = OldParamHints)
        names(x = unused.hints) <- unused
        unused.hints <- na.omit(object = unused.hints)
        if (length(x = unused.hints) > 0) {
          message(
            "Suggested parameter: ",
            paste(unused.hints, "instead of", names(x = unused.hints), collapse = '; '),
            "\n"
          )
        }
      }
    }
  }
  SingleRasterMap <- function (data, raster = TRUE, cell.order = NULL, feature.order = NULL, 
                               colors = PurpleAndYellow(), disp.min = -2.5, disp.max = 2.5, 
                               limits = NULL, group.by = NULL) 
  {
    data <- MinMax(data = data, min = disp.min, max = disp.max)
    data <- Melt(x = t(x = data))
    colnames(x = data) <- c("Feature", "Cell", "Expression")
    if (!is.null(x = feature.order)) {
      data$Feature <- factor(x = data$Feature, levels = unique(x = feature.order))
    }
    if (!is.null(x = cell.order)) {
      data$Cell <- factor(x = data$Cell, levels = unique(x = cell.order))
    }
    if (!is.null(x = group.by)) {
      data$Identity <- group.by[data$Cell]
    }
    limits <- limits %||% c(min(data$Expression), max(data$Expression))
    if (length(x = limits) != 2 || !is.numeric(x = limits)) {
      stop("limits' must be a two-length numeric vector")
    }
    my_geom <- ifelse(test = raster, yes = geom_raster, no = geom_tile)
    plot <- ggplot(data = data) + my_geom(mapping = aes_string(x = "Cell", 
                                                               y = "Feature", fill = "Expression")) + theme(axis.text.x = element_blank(), 
                                                                                                            axis.ticks.x = element_blank()) + scale_fill_gradientn(limits = limits, 
                                                                                                                                                                   colors = colors, na.value = "white") + labs(x = NULL, 
                                                                                                                                                                                                               y = NULL, fill = group.by %iff% "Expression") + WhiteBackground() + 
      NoAxes(keep.text = TRUE)
    if (!is.null(x = group.by)) {
      plot <- plot + geom_point(mapping = aes_string(x = "Cell", 
                                                     y = "Feature", color = "Identity"), alpha = 0) + 
        guides(color = guide_legend(override.aes = list(alpha = 1)))
    }
    return(plot)
  }
  RandomName <- function(length = 5L, ...) {
    CheckDots(..., fxns = 'sample')
    return(paste(sample(x = letters, size = length, ...), collapse = ''))
  }
  Melt <- function(x) {
    if (!is.data.frame(x = x)) {
      x <- as.data.frame(x = x)
    }
    return(data.frame(
      rows = rep.int(x = rownames(x = x), times = ncol(x = x)),
      cols = unlist(x = lapply(X = colnames(x = x), FUN = rep.int, times = nrow(x = x))),
      vals = unlist(x = x, use.names = FALSE)
    ))
  }
  `%||%` <- function(lhs, rhs) {
    if (!is.null(x = lhs)) {
      return(lhs)
    } else {
      return(rhs)
    }
  }
  `%iff%` <- function(lhs, rhs) {
    if (!is.null(x = lhs)) {
      return(rhs)
    } else {
      return(lhs)
    }
  }
  
  cells <- cells %||% colnames(x = object)
  if (is.numeric(x = cells)) {
    cells <- colnames(x = object)[cells]
  }
  assay <- assay %||% DefaultAssay(object = object)
  DefaultAssay(object = object) <- assay
  features <- features %||% VariableFeatures(object = object)
  features <- rev(x = unique(x = features))
  disp.max <- disp.max %||% ifelse(test = slot == "scale.data", 
                                   yes = 2.5, no = 6)
  possible.features <- rownames(x = GetAssayData(object = object, 
                                                 slot = slot))
  if (any(!features %in% possible.features)) {
    bad.features <- features[!features %in% possible.features]
    features <- features[features %in% possible.features]
    if (length(x = features) == 0) {
      stop("No requested features found in the ", slot, 
           " slot for the ", assay, " assay.")
    }
    warning("The following features were omitted as they were not found in the ", 
            slot, " slot for the ", assay, " assay: ", paste(bad.features, 
                                                             collapse = ", "))
  }
  data <- as.data.frame(x = as.matrix(x = t(x = GetAssayData(object = object, 
                                                             slot = slot)[features, cells, drop = FALSE])))
  object <- suppressMessages(expr = StashIdent(object = object, 
                                               save.name = "ident"))
  group.by <- group.by %||% "ident"
  #groups.use <- object[[group.by]][cells, , drop = FALSE]
  
  old.group <- object[[group.by, drop = TRUE]]
  if (!is.factor(x = old.group)) {
    old.group <- factor(x = old.group)
  }
  id.levels <- levels(x = old.group)
  old.group <- as.vector(x = old.group)
  if (!is.null(x = split.by)) {
    splits <- object[[split.by, drop = TRUE]]
    groups.use <- paste(old.group, splits, sep = "_")
    groups.use <- data.frame(ident=groups.use)
    rownames(groups.use) <- rownames(object[[group.by]][cells, , drop = FALSE])
    unique.splits <- rev(unique(x = splits))
    id.levels <- paste0(rep(x = id.levels, each = length(x = unique.splits)), 
                        "_", rep(x = unique.splits, times = length(x = id.levels)))
    groups.use$ident <- factor(x=groups.use$ident, levels = id.levels)
  }
  
  plots <- vector(mode = "list", length = ncol(x = groups.use))
  for (i in 1:ncol(x = groups.use)) {
    data.group <- data
    group.use <- groups.use[, i, drop = TRUE]
    if (!is.factor(x = group.use)) {
      group.use <- factor(x = group.use)
    }
    names(x = group.use) <- cells
    if (draw.lines) {
      lines.width <- lines.width %||% ceiling(x = nrow(x = data.group) * 
                                                0.0025)
      placeholder.cells <- sapply(X = 1:(length(x = levels(x = group.use)) * 
                                           lines.width), FUN = function(x) {
                                             return(RandomName(length = 20))
                                           })
      placeholder.groups <- rep(x = levels(x = group.use), 
                                times = lines.width)
      group.levels <- levels(x = group.use)
      names(x = placeholder.groups) <- placeholder.cells
      group.use <- as.vector(x = group.use)
      names(x = group.use) <- cells
      group.use <- factor(x = c(group.use, placeholder.groups), 
                          levels = group.levels)
      na.data.group <- matrix(data = NA, nrow = length(x = placeholder.cells), 
                              ncol = ncol(x = data.group), dimnames = list(placeholder.cells, 
                                                                           colnames(x = data.group)))
      data.group <- rbind(data.group, na.data.group)
    }
    lgroup <- length(levels(group.use))
    plot <- SingleRasterMap(data = data.group, raster = raster, 
                            disp.min = disp.min, disp.max = disp.max, feature.order = features, 
                            cell.order = names(x = sort(x = group.use)), group.by = group.use)
    if (group.bar) {
      default.colors <- c(hue_pal()(length(x = levels(x = group.use))))
      cols <- group.colors[1:length(x = levels(x = group.use))] %||% 
        default.colors
      if (any(is.na(x = cols))) {
        cols[is.na(x = cols)] <- default.colors[is.na(x = cols)]
        cols <- Col2Hex(cols)
        col.dups <- sort(x = unique(x = which(x = duplicated(x = substr(x = cols, 
                                                                        start = 1, stop = 7)))))
        through <- length(x = default.colors)
        while (length(x = col.dups) > 0) {
          pal.max <- length(x = col.dups) + through
          cols.extra <- hue_pal()(pal.max)[(through + 
                                              1):pal.max]
          cols[col.dups] <- cols.extra
          col.dups <- sort(x = unique(x = which(x = duplicated(x = substr(x = cols, 
                                                                          start = 1, stop = 7)))))
        }
      }
      group.use2 <- sort(x = group.use)
      if (draw.lines) {
        na.group <- RandomName(length = 20)
        levels(x = group.use2) <- c(levels(x = group.use2), 
                                    na.group)
        group.use2[placeholder.cells] <- na.group
        cols <- c(cols, "#FFFFFF")
      }
      pbuild <- ggplot_build(plot = plot)
      names(x = cols) <- levels(x = group.use2)
      y.range <- diff(x = pbuild$layout$panel_params[[1]]$y.range)
      y.pos <- max(pbuild$layout$panel_params[[1]]$y.range) + 
        y.range * 0.015
      y.max <- y.pos + group.bar.height * y.range
      plot <- plot + annotation_raster(raster = t(x = cols[group.use2]), 
                                       xmin = -Inf, xmax = Inf, ymin = y.pos, ymax = y.max) + 
        coord_cartesian(ylim = c(0, y.max), clip = "off") + 
        scale_color_manual(values = cols)
      if (label) {
        x.max <- max(pbuild$layout$panel_params[[1]]$x.range)
        x.divs <- pbuild$layout$panel_params[[1]]$x.major %||% 
          pbuild$layout$panel_params[[1]]$x$break_positions()
        x <- data.frame(group = sort(x = group.use), 
                        x = x.divs)
        label.x.pos <- tapply(X = x$x, INDEX = x$group, 
                              FUN = median) * x.max
        label.x.pos <- data.frame(group = names(x = label.x.pos), 
                                  label.x.pos)
        plot <- plot + geom_text(stat = "identity", 
                                 data = label.x.pos, aes_string(label = "group", 
                                                                x = "label.x.pos"), y = y.max + y.max * 
                                   0.03 * 0.5, angle = angle, hjust = hjust, 
                                 size = size)
        plot <- suppressMessages(plot + coord_cartesian(ylim = c(0, 
                                                                 y.max + y.max * 0.002 * max(nchar(x = levels(x = group.use))) * 
                                                                   size), clip = "off"))
      }
    }
    plot <- plot + theme(line = element_blank())
    plots[[i]] <- plot
  }
  if (combine) {
    plots <- wrap_plots(plots)
  }
  return(plots)
}

##################### 5. caculate cell composition change and determine signficance ##########
### function to run stats to see if the change is significant
getStats_FC <- function(gbm=gbm, cutoff=0.05, allcells=NULL, toRemove=NULL) {
  library(rcompanion) 
  prop_num <- table(gbm@active.ident, gbm@meta.data[, "orig.ident"])
  prop_num.df <- data.frame("wt.13hpf"=prop_num[, 2],
                            "gata56MO.13hpf"=prop_num[, 1])
  prop_num.df$clu_id <- rownames(prop_num.df)
  prop_num.df <- prop_num.df[!(rownames(prop_num.df) %in% toRemove), ]
  
  prop_test_each <- function(cluster_id, df=prop_num.df) {
    if (is.null(allcells)) {
      wt_sum <- sum(df$wt.13hpf)#; print(wt_sum)
      mo_sum <- sum(df$gata56MO.13hpf)#; print(mo_sum)
    } else { ## for the case of subsetting
      total_sum <- as.data.frame(table(allcells$orig.ident))
      wt_sum <- total_sum[2, 2]#; print(wt_sum)
      mo_sum <- total_sum[1, 2]
    }
    
    
    new_df <- data.frame("cluster"=c(df[cluster_id,"wt.13hpf"], df[cluster_id,"gata56MO.13hpf"]),
                         "no" =c(wt_sum-df[cluster_id,"wt.13hpf"], mo_sum-df[cluster_id,"gata56MO.13hpf"])
    )
    rownames(new_df) <- c("wt", "mo")
    #print(new_df)
    stat_res <- pairwiseNominalIndependence(as.matrix(new_df))
    return(stat_res)
  }
  
  stats_ls <- lapply(prop_num.df$clu_id, prop_test_each)
  names(stats_ls) <- prop_num.df$clu_id
  stats_df <- do.call(rbind, stats_ls)
  stats_df$p.adj.Fisher <- p.adjust(stats_df$p.Fisher, method = "bonferroni")
  stats_df$p.adj.Chisq <- p.adjust(stats_df$p.Chisq, method = "bonferroni")
  stats_df$p.adj.Gtest <- p.adjust(stats_df$p.Gtest, method = "bonferroni")
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
  return(stats_df)
}