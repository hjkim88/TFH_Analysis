###
#   File name : TFH_Analyses_Both_Donors.R
#   Author    : Hyunjin Kim
#   Date      : Aug 16, 2020
#   Email     : hyunjin.kim@stjude.org
#   Purpose   : Now we have a new dataset from the other donor.
#               Perform the same analyses and find shared clones between the two donors.
#   
#   Instruction
#               1. Source("TFH_Analyses_Both_Donors.R")
#               2. Run the function "tfh_analyses_both_donors" - specify the input file paths and the output directory
#               3. The results will be generated under the output directory
#
#   Example
#               > source("The_directory_of_TFH_Analyses_Both_Donors.R/TFH_Analyses_Both_Donors.R")
#               > tfh_analyses_both_donors(Seurat_RObj_path="./data/SS_Tfh_BothDonors/BothDonors_Tfh_matched_clones_v2.rds",
#                                          outputDir="./results/v2/")
###

tfh_analyses_both_donors <- function(Seurat_RObj_path="./data/SS_Tfh_BothDonors/BothDonors_Tfh_matched_clones_v2.rds",
                                     outputDir="./results/v2/") {
  
  ### load libraries
  if(!require(Seurat, quietly = TRUE)) {
    install.packages("Seurat")
    require(Seurat, quietly = TRUE)
  }
  options(java.parameters = "-Xmx10240m")
  if(!require(xlsx, quietly = TRUE)) {
    install.packages("xlsx")
    require(xlsx, quietly = TRUE)
  }
  if(!require(ggplot2, quietly = TRUE)) {
    install.packages("ggplot2")
    require(ggplot2, quietly = TRUE)
  }
  if(!require(gridExtra, quietly = TRUE)) {
    install.packages("gridExtra")
    require(gridExtra, quietly = TRUE)
  }
  if(!require(ggalluvial, quietly = TRUE)) {
    install.packages("ggalluvial")
    require(ggalluvial, quietly = TRUE)
  }
  if(!require(ggpubr, quietly = TRUE)) {
    install.packages("ggpubr")
    require(ggpubr, quietly = TRUE)
  }
  if(!require(viridis, quietly = TRUE)) {
    install.packages("viridis")
    require(viridis, quietly = TRUE)
  }
  if(!require(ggsci, quietly = TRUE)) {
    install.packages("ggsci")
    require(ggsci, quietly = TRUE)
  }
  if(!require(org.Hs.eg.db, quietly = TRUE)) {
    if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
    BiocManager::install("org.Hs.eg.db")
    require(org.Hs.eg.db, quietly = TRUE)
  }
  if(!require(gplots, quietly = TRUE)) {
    install.packages("gplots")
    library(gplots, quietly = TRUE)
  }
  if(!require(slingshot, quietly = TRUE)) {
    if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
    BiocManager::install("slingshot")
    require(slingshot, quietly = TRUE)
  }
  if(!require(mclust, quietly = TRUE)) {
    install.packages("mclust")
    library(mclust, quietly = TRUE)
  }
  if(!require(scales, quietly = TRUE)) {
    install.packages("scales")
    library(scales, quietly = TRUE)
  }
  if(!require(tidymodels, quietly = TRUE)) {
    install.packages("tidymodels")
    library(tidymodels, quietly = TRUE)
  }
  if(!require(ranger, quietly = TRUE)) {
    install.packages("ranger")
    library(ranger, quietly = TRUE)
  }
  if(!require(viridis, quietly = TRUE)) {
    install.packages("viridis")
    require(viridis, quietly = TRUE)
  }
  if(!require(ggbeeswarm, quietly = TRUE)) {
    install.packages("ggbeeswarm")
    require(ggbeeswarm, quietly = TRUE)
  }
  if(!require(RColorBrewer, quietly = TRUE)) {
    install.packages("RColorBrewer")
    require(RColorBrewer, quietly = TRUE)
  }
  
  ### load Stefan's t-cell object v2 - both donors
  Seurat_Obj <- readRDS(Seurat_RObj_path)
  
  ### rownames in the meta.data should be in the same order as colnames in the counts
  Seurat_Obj@meta.data <- Seurat_Obj@meta.data[colnames(Seurat_Obj@assays$RNA@counts),]
  
  ### see how cells are distributed based on Donor info
  DimPlot(Seurat_Obj, reduction = "umap", group.by = "Donor", pt.size = 1)
  
  ### change the Day column of the meta.data
  Seurat_Obj@meta.data$Day <- factor(paste0("d", as.character(Seurat_Obj@meta.data$Day)),
                                     levels = c("d0", "d5", "d12", "d28", "d60", "d90", "d120", "d180"))
  
  ### shared clone between two donors
  match_cdr_04 <- Seurat_Obj@meta.data$match.cdr[which(Seurat_Obj@meta.data$Donor == "321-04")]
  match_cdr_05 <- Seurat_Obj@meta.data$match.cdr[which(Seurat_Obj@meta.data$Donor == "321-05")]
  print(length(intersect(match_cdr_04, match_cdr_05)))
  
  ### theme that draws dotted lines for each y-axis ticks
  ### this function is from "immunarch" package
  theme_cleveland2 <- function(rotate = TRUE) {
    if (rotate) {
      theme(
        panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_line(
          colour = "grey70",
          linetype = "dashed"
        )
      )
    }
    else {
      theme(
        panel.grid.major.x = element_line(
          colour = "grey70",
          linetype = "dashed"
        ), panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank()
      )
    }
  }
  
  # ******************************************************************************************
  # Pathway Analysis with clusterProfiler package
  # Input: geneList     = a vector of gene Entrez IDs for pathway analysis [numeric or character]
  #        org          = organism that will be used in the analysis ["human" or "mouse"]
  #                       should be either "human" or "mouse"
  #        database     = pathway analysis database (KEGG or GO) ["KEGG" or "GO"]
  #        title        = title of the pathway figure [character]
  #        pv_threshold = pathway analysis p-value threshold (not DE analysis threshold) [numeric]
  #        displayNum   = the number of pathways that will be displayed [numeric]
  #                       (If there are many significant pathways show the few top pathways)
  #        imgPrint     = print a plot of pathway analysis [TRUE/FALSE]
  #        dir          = file directory path of the output pathway figure [character]
  #
  # Output: Pathway analysis results in figure - using KEGG and GO pathways
  #         The x-axis represents the number of DE genes in the pathway
  #         The y-axis represents pathway names
  #         The color of a bar indicates adjusted p-value from the pathway analysis
  #         For Pathview Result, all colored genes are found DE genes in the pathway,
  #         and the color indicates log2(fold change) of the DE gene from DE analysis
  # ******************************************************************************************
  pathwayAnalysis_CP <- function(geneList,
                                 org,
                                 database,
                                 title="Pathway_Results",
                                 pv_threshold=0.05,
                                 displayNum=Inf,
                                 imgPrint=TRUE,
                                 dir="./") {
    
    ### load library
    if(!require(clusterProfiler, quietly = TRUE)) {
      if (!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager")
      BiocManager::install("clusterProfiler")
      require(clusterProfiler, quietly = TRUE)
    }
    if(!require(ggplot2)) {
      install.packages("ggplot2")
      library(ggplot2)
    }
    
    
    ### collect gene list (Entrez IDs)
    geneList <- geneList[which(!is.na(geneList))]
    
    if(!is.null(geneList)) {
      ### make an empty list
      p <- list()
      
      if(database == "KEGG") {
        ### KEGG Pathway
        kegg_enrich <- enrichKEGG(gene = geneList, organism = org, pvalueCutoff = pv_threshold)
        
        if(is.null(kegg_enrich)) {
          writeLines("KEGG Result does not exist")
          return(NULL)
        } else {
          kegg_enrich@result <- kegg_enrich@result[which(kegg_enrich@result$p.adjust < pv_threshold),]
          
          if(imgPrint == TRUE) {
            if((displayNum == Inf) || (nrow(kegg_enrich@result) <= displayNum)) {
              result <- kegg_enrich@result
              description <- kegg_enrich@result$Description
            } else {
              result <- kegg_enrich@result[1:displayNum,]
              description <- kegg_enrich@result$Description[1:displayNum]
            }
            
            if(nrow(kegg_enrich) > 0) {
              p[[1]] <- ggplot(result, aes(x=Description, y=Count)) + labs(x="", y="Gene Counts") + 
                theme_classic(base_size = 16) + geom_bar(aes(fill = p.adjust), stat="identity") + coord_flip() +
                scale_x_discrete(limits = rev(description)) +
                guides(fill = guide_colorbar(ticks=FALSE, title="P.Val", barheight=10)) +
                ggtitle(paste0("KEGG ", title))
              
              png(paste0(dir, "kegg_", title, "_CB.png"), width = 2000, height = 1000)
              print(p[[1]])
              dev.off()
            } else {
              writeLines("KEGG Result does not exist")
            }
          }
          
          return(kegg_enrich@result)
        }
      } else if(database == "GO") {
        ### GO Pathway
        if(org == "human") {
          go_enrich <- enrichGO(gene = geneList, OrgDb = 'org.Hs.eg.db', readable = T, ont = "BP", pvalueCutoff = pv_threshold)
        } else if(org == "mouse") {
          go_enrich <- enrichGO(gene = geneList, OrgDb = 'org.Mm.eg.db', readable = T, ont = "BP", pvalueCutoff = pv_threshold)
        } else {
          go_enrich <- NULL
          writeLines(paste("Unknown org variable:", org))
        }
        
        if(is.null(go_enrich)) {
          writeLines("GO Result does not exist")
          return(NULL)
        } else {
          go_enrich@result <- go_enrich@result[which(go_enrich@result$p.adjust < pv_threshold),]
          
          if(imgPrint == TRUE) {
            if((displayNum == Inf) || (nrow(go_enrich@result) <= displayNum)) {
              result <- go_enrich@result
              description <- go_enrich@result$Description
            } else {
              result <- go_enrich@result[1:displayNum,]
              description <- go_enrich@result$Description[1:displayNum]
            }
            
            if(nrow(go_enrich) > 0) {
              p[[2]] <- ggplot(result, aes(x=Description, y=Count)) + labs(x="", y="Gene Counts") + 
                theme_classic(base_size = 16) + geom_bar(aes(fill = p.adjust), stat="identity") + coord_flip() +
                scale_x_discrete(limits = rev(description)) +
                guides(fill = guide_colorbar(ticks=FALSE, title="P.Val", barheight=10)) +
                ggtitle(paste0("GO ", title))
              
              png(paste0(dir, "go_", title, "_CB.png"), width = 2000, height = 1000)
              print(p[[2]])
              dev.off()
            } else {
              writeLines("GO Result does not exist")
            }
          }
          
          return(go_enrich@result)
        }
      } else {
        stop("database prameter should be \"GO\" or \"KEGG\"")
      }
    } else {
      writeLines("geneList = NULL")
    }
  }
  
  ### a function that returns multi-figures from PCA, TSNE, AND UMAP
  ### One plot with all the groups + each group
  ### x: a numeric vector of first component
  ### y: a numeric vector of second component
  ### group_color: a factor vector of group1 info
  ### group_shape: a factor vector of group2 info
  ### type: type of the dimensionality reduction method
  multiReducPlot <- function(x, y, group_color, group_shape=NULL, type=c("PCA", "TSNE", "UMAP"), fName="Whole Group & Each", isConvex=FALSE, isPrint=FALSE) {
    
    ### load library
    if(!require(ggplot2, quietly = TRUE)) {
      install.packages("ggplot2")
      require(ggplot2, quietly = TRUE)
    }
    if(!require(gridExtra, quietly = TRUE)) {
      install.packages("gridExtra")
      require(gridExtra, quietly = TRUE)
    }
    if(!require(scales, quietly = TRUE)) {
      install.packages("scales")
      require(scales, quietly = TRUE)
    }
    if(!require(ggConvexHull, quietly = TRUE)) {
      if(!require(remotes, quietly = TRUE)) {
        install.packages("remotes")
      }
      remotes::install_github("cmartin/ggConvexHull")
      require(ggConvexHull, quietly = TRUE)
    }
    
    ### create a data frame for ggplot
    group_color <- as.factor(group_color)
    group_color <- factor(as.character(group_color),
                          levels = intersect(levels(group_color), unique(group_color)))
    if(is.null(group_shape)) {
      plot_df <- data.frame(X=x, Y=y, Group=group_color, Group2="")
    } else {
      group_shape <- as.factor(group_shape)
      group_shape <- factor(as.character(group_shape),
                            levels = intersect(levels(group_shape), unique(group_shape)))
      plot_df <- data.frame(X=x, Y=y, Group=group_color, Group2=group_shape)
    }
    
    ### set x & y axes labels
    if(type[1] == "PCA") {
      x_label <- "PC1"
      y_label <- "PC2"
    } else if(type[1] == "TSNE") {
      x_label <- "TSNE1"
      y_label <- "TSNE2"
    } else if(type[1] == "UMAP") {
      x_label <- "UMAP1"
      y_label <- "UMAP2"
    } else {
      stop("ERROR: type parameter should be either \"PCA\", \"TSNE\", or \"UMAP\".")
    }
    
    ### set colors for each group
    col_palette <- hue_pal()(length(levels(group_color)))
    names(col_palette) <- levels(group_color)
    
    ### 1. One plot with all the groups + each group
    p <- vector("list", length(levels(group_color))+1)
    x_range <- c(min(plot_df$X), max(plot_df$X))
    y_range <- c(min(plot_df$Y), max(plot_df$Y))
    if(isConvex) {
      p[[1]] <- ggplot(plot_df, aes_string(x="X", y="Y")) +
        geom_point(aes_string(col="Group", shape="Group2"), size=2, alpha=0.6) +
        xlab(x_label) + ylab(y_label) +
        xlim(x_range[1], x_range[2]) + ylim(y_range[1], y_range[2]) +
        ggtitle("All") +
        theme_classic(base_size = 16) +
        theme(legend.position = "none", legend.title = element_blank(),
              plot.title = element_text(hjust = 0.5),
              legend.spacing.y = unit(0, "mm"),
              legend.background = element_blank(),
              legend.box.background = element_rect(colour = "gray50")) +
        geom_convexhull(aes_string(col="Group", fill="Group"), size=1, alpha=0)
    } else {
      p[[1]] <- ggplot(plot_df, aes_string(x="X", y="Y")) +
        geom_point(aes_string(col="Group", shape="Group2"), size=2, alpha=0.6) +
        xlab(x_label) + ylab(y_label) +
        xlim(x_range[1], x_range[2]) + ylim(y_range[1], y_range[2]) +
        ggtitle("All") +
        theme_classic(base_size = 16) +
        theme(legend.position = "none", legend.title = element_blank(),
              plot.title = element_text(hjust = 0.5),
              legend.background = element_blank(),
              legend.spacing.y = unit(0, "mm"),
              legend.box.background = element_rect(colour = "gray50"))
    }
    for(i in 2:(length(levels(group_color))+1)) {
      p[[i]] <- ggplot(plot_df[which(plot_df$Group == levels(group_color)[i-1]),], aes_string(x="X", y="Y")) +
        geom_point(aes_string(shape="Group2"), col=col_palette[levels(group_color)[i-1]], size=2, alpha=0.8) +
        xlab(x_label) + ylab(y_label) +
        xlim(x_range[1], x_range[2]) + ylim(y_range[1], y_range[2]) +
        ggtitle(levels(group_color)[i-1]) +
        theme_classic(base_size = 16) +
        theme(legend.title = element_blank(),
              plot.title = element_text(hjust = 0.5, 
                                        color = col_palette[levels(group_color)[i-1]]),
              legend.background = element_blank(),
              legend.spacing.y = unit(0, "mm"),
              legend.box.background = element_rect(colour = "gray50"))
    }
    ### arrange the plots
    g <- arrangeGrob(grobs = p,
                     nrow = ceiling(sqrt(length(levels(group_color))+1)),
                     ncol = ceiling(sqrt(length(levels(group_color))+1)),
                     top = fName)
    
    if(isPrint) {
      plot(g)
    }
    
    return(g)
    
  }
  
  ###
  #   This function was downloaded from: https://raw.githubusercontent.com/obigriffith/biostar-tutorials/master/Heatmaps/heatmap.3.R
  ###
  heatmap.3 <- function(x,
                        Rowv = TRUE, Colv = if (symm) "Rowv" else TRUE,
                        distfun = dist,
                        hclustfun = hclust,
                        dendrogram = c("both","row", "column", "none"),
                        symm = FALSE,
                        scale = c("none","row", "column"),
                        na.rm = TRUE,
                        revC = identical(Colv,"Rowv"),
                        add.expr,
                        breaks,
                        symbreaks = max(x < 0, na.rm = TRUE) || scale != "none",
                        col = "heat.colors",
                        colsep,
                        rowsep,
                        sepcolor = "white",
                        sepwidth = c(0.05, 0.05),
                        cellnote,
                        notecex = 1,
                        notecol = "cyan",
                        na.color = par("bg"),
                        trace = c("none", "column","row", "both"),
                        tracecol = "cyan",
                        hline = median(breaks),
                        vline = median(breaks),
                        linecol = tracecol,
                        margins = c(5,5),
                        ColSideColors,
                        RowSideColors,
                        side.height.fraction=0.3,
                        cexRow = 0.2 + 1/log10(nr),
                        cexCol = 0.2 + 1/log10(nc),
                        labRow = NULL,
                        labCol = NULL,
                        key = TRUE,
                        keysize = 1.5,
                        density.info = c("none", "histogram", "density"),
                        denscol = tracecol,
                        symkey = max(x < 0, na.rm = TRUE) || symbreaks,
                        densadj = 0.25,
                        main = NULL,
                        xlab = NULL,
                        ylab = NULL,
                        lmat = NULL,
                        lhei = NULL,
                        lwid = NULL,
                        ColSideColorsSize = 1,
                        RowSideColorsSize = 1,
                        KeyValueName="Value",...){
    
    invalid <- function (x) {
      if (missing(x) || is.null(x) || length(x) == 0)
        return(TRUE)
      if (is.list(x))
        return(all(sapply(x, invalid)))
      else if (is.vector(x))
        return(all(is.na(x)))
      else return(FALSE)
    }
    
    x <- as.matrix(x)
    scale01 <- function(x, low = min(x), high = max(x)) {
      x <- (x - low)/(high - low)
      x
    }
    retval <- list()
    scale <- if (symm && missing(scale))
      "none"
    else match.arg(scale)
    dendrogram <- match.arg(dendrogram)
    trace <- match.arg(trace)
    density.info <- match.arg(density.info)
    if (length(col) == 1 && is.character(col))
      col <- get(col, mode = "function")
    if (!missing(breaks) && (scale != "none"))
      warning("Using scale=\"row\" or scale=\"column\" when breaks are",
              "specified can produce unpredictable results.", "Please consider using only one or the other.")
    if (is.null(Rowv) || is.na(Rowv))
      Rowv <- FALSE
    if (is.null(Colv) || is.na(Colv))
      Colv <- FALSE
    else if (Colv == "Rowv" && !isTRUE(Rowv))
      Colv <- FALSE
    if (length(di <- dim(x)) != 2 || !is.numeric(x))
      stop("`x' must be a numeric matrix")
    nr <- di[1]
    nc <- di[2]
    if (nr <= 1 || nc <= 1)
      stop("`x' must have at least 2 rows and 2 columns")
    if (!is.numeric(margins) || length(margins) != 2)
      stop("`margins' must be a numeric vector of length 2")
    if (missing(cellnote))
      cellnote <- matrix("", ncol = ncol(x), nrow = nrow(x))
    if (!inherits(Rowv, "dendrogram")) {
      if (((!isTRUE(Rowv)) || (is.null(Rowv))) && (dendrogram %in%
                                                   c("both", "row"))) {
        if (is.logical(Colv) && (Colv))
          dendrogram <- "column"
        else dedrogram <- "none"
        warning("Discrepancy: Rowv is FALSE, while dendrogram is `",
                dendrogram, "'. Omitting row dendogram.")
      }
    }
    if (!inherits(Colv, "dendrogram")) {
      if (((!isTRUE(Colv)) || (is.null(Colv))) && (dendrogram %in%
                                                   c("both", "column"))) {
        if (is.logical(Rowv) && (Rowv))
          dendrogram <- "row"
        else dendrogram <- "none"
        warning("Discrepancy: Colv is FALSE, while dendrogram is `",
                dendrogram, "'. Omitting column dendogram.")
      }
    }
    if (inherits(Rowv, "dendrogram")) {
      ddr <- Rowv
      rowInd <- order.dendrogram(ddr)
    }
    else if (is.integer(Rowv)) {
      hcr <- hclustfun(distfun(x))
      ddr <- as.dendrogram(hcr)
      ddr <- reorder(ddr, Rowv)
      rowInd <- order.dendrogram(ddr)
      if (nr != length(rowInd))
        stop("row dendrogram ordering gave index of wrong length")
    }
    else if (isTRUE(Rowv)) {
      Rowv <- rowMeans(x, na.rm = na.rm)
      hcr <- hclustfun(distfun(x))
      ddr <- as.dendrogram(hcr)
      ddr <- reorder(ddr, Rowv)
      rowInd <- order.dendrogram(ddr)
      if (nr != length(rowInd))
        stop("row dendrogram ordering gave index of wrong length")
    }
    else {
      rowInd <- nr:1
    }
    if (inherits(Colv, "dendrogram")) {
      ddc <- Colv
      colInd <- order.dendrogram(ddc)
    }
    else if (identical(Colv, "Rowv")) {
      if (nr != nc)
        stop("Colv = \"Rowv\" but nrow(x) != ncol(x)")
      if (exists("ddr")) {
        ddc <- ddr
        colInd <- order.dendrogram(ddc)
      }
      else colInd <- rowInd
    }
    else if (is.integer(Colv)) {
      hcc <- hclustfun(distfun(if (symm)
        x
        else t(x)))
      ddc <- as.dendrogram(hcc)
      ddc <- reorder(ddc, Colv)
      colInd <- order.dendrogram(ddc)
      if (nc != length(colInd))
        stop("column dendrogram ordering gave index of wrong length")
    }
    else if (isTRUE(Colv)) {
      Colv <- colMeans(x, na.rm = na.rm)
      hcc <- hclustfun(distfun(if (symm)
        x
        else t(x)))
      ddc <- as.dendrogram(hcc)
      ddc <- reorder(ddc, Colv)
      colInd <- order.dendrogram(ddc)
      if (nc != length(colInd))
        stop("column dendrogram ordering gave index of wrong length")
    }
    else {
      colInd <- 1:nc
    }
    retval$rowInd <- rowInd
    retval$colInd <- colInd
    retval$call <- match.call()
    x <- x[rowInd, colInd]
    x.unscaled <- x
    cellnote <- cellnote[rowInd, colInd]
    if (is.null(labRow))
      labRow <- if (is.null(rownames(x)))
        (1:nr)[rowInd]
    else rownames(x)
    else labRow <- labRow[rowInd]
    if (is.null(labCol))
      labCol <- if (is.null(colnames(x)))
        (1:nc)[colInd]
    else colnames(x)
    else labCol <- labCol[colInd]
    if (scale == "row") {
      retval$rowMeans <- rm <- rowMeans(x, na.rm = na.rm)
      x <- sweep(x, 1, rm)
      retval$rowSDs <- sx <- apply(x, 1, sd, na.rm = na.rm)
      x <- sweep(x, 1, sx, "/")
    }
    else if (scale == "column") {
      retval$colMeans <- rm <- colMeans(x, na.rm = na.rm)
      x <- sweep(x, 2, rm)
      retval$colSDs <- sx <- apply(x, 2, sd, na.rm = na.rm)
      x <- sweep(x, 2, sx, "/")
    }
    if (missing(breaks) || is.null(breaks) || length(breaks) < 1) {
      if (missing(col) || is.function(col))
        breaks <- 16
      else breaks <- length(col) + 1
    }
    if (length(breaks) == 1) {
      if (!symbreaks)
        breaks <- seq(min(x, na.rm = na.rm), max(x, na.rm = na.rm),
                      length = breaks)
      else {
        extreme <- max(abs(x), na.rm = TRUE)
        breaks <- seq(-extreme, extreme, length = breaks)
      }
    }
    nbr <- length(breaks)
    ncol <- length(breaks) - 1
    if (class(col) == "function")
      col <- col(ncol)
    min.breaks <- min(breaks)
    max.breaks <- max(breaks)
    x[x < min.breaks] <- min.breaks
    x[x > max.breaks] <- max.breaks
    if (missing(lhei) || is.null(lhei))
      lhei <- c(keysize, 4)
    if (missing(lwid) || is.null(lwid))
      lwid <- c(keysize, 4)
    if (missing(lmat) || is.null(lmat)) {
      lmat <- rbind(4:3, 2:1)
      
      if (!missing(ColSideColors)) {
        #if (!is.matrix(ColSideColors))
        #stop("'ColSideColors' must be a matrix")
        if (!is.character(ColSideColors) || nrow(ColSideColors) != nc)
          stop("'ColSideColors' must be a matrix of nrow(x) rows")
        lmat <- rbind(lmat[1, ] + 1, c(NA, 1), lmat[2, ] + 1)
        #lhei <- c(lhei[1], 0.2, lhei[2])
        lhei=c(lhei[1], side.height.fraction*ColSideColorsSize/2, lhei[2])
      }
      
      if (!missing(RowSideColors)) {
        #if (!is.matrix(RowSideColors))
        #stop("'RowSideColors' must be a matrix")
        if (!is.character(RowSideColors) || ncol(RowSideColors) != nr)
          stop("'RowSideColors' must be a matrix of ncol(x) columns")
        lmat <- cbind(lmat[, 1] + 1, c(rep(NA, nrow(lmat) - 1), 1), lmat[,2] + 1)
        #lwid <- c(lwid[1], 0.2, lwid[2])
        lwid <- c(lwid[1], side.height.fraction*RowSideColorsSize/2, lwid[2])
      }
      lmat[is.na(lmat)] <- 0
    }
    
    if (length(lhei) != nrow(lmat))
      stop("lhei must have length = nrow(lmat) = ", nrow(lmat))
    if (length(lwid) != ncol(lmat))
      stop("lwid must have length = ncol(lmat) =", ncol(lmat))
    op <- par(no.readonly = TRUE)
    on.exit(par(op))
    
    layout(lmat, widths = lwid, heights = lhei, respect = FALSE)
    
    if (!missing(RowSideColors)) {
      if (!is.matrix(RowSideColors)){
        par(mar = c(margins[1], 0, 0, 0.5))
        image(rbind(1:nr), col = RowSideColors[rowInd], axes = FALSE)
      } else {
        par(mar = c(margins[1], 0, 0, 0.5))
        rsc = t(RowSideColors[,rowInd, drop=F])
        rsc.colors = matrix()
        rsc.names = names(table(rsc))
        rsc.i = 1
        for (rsc.name in rsc.names) {
          rsc.colors[rsc.i] = rsc.name
          rsc[rsc == rsc.name] = rsc.i
          rsc.i = rsc.i + 1
        }
        rsc = matrix(as.numeric(rsc), nrow = dim(rsc)[1])
        image(t(rsc), col = as.vector(rsc.colors), axes = FALSE)
        if (length(rownames(RowSideColors)) > 0) {
          axis(1, 0:(dim(rsc)[2] - 1)/max(1,(dim(rsc)[2] - 1)), rownames(RowSideColors), las = 2, tick = FALSE)
        }
      }
    }
    
    if (!missing(ColSideColors)) {
      
      if (!is.matrix(ColSideColors)){
        par(mar = c(0.5, 0, 0, margins[2]))
        image(cbind(1:nc), col = ColSideColors[colInd], axes = FALSE)
      } else {
        par(mar = c(0.5, 0, 0, margins[2]))
        csc = ColSideColors[colInd, , drop=F]
        csc.colors = matrix()
        csc.names = names(table(csc))
        csc.i = 1
        for (csc.name in csc.names) {
          csc.colors[csc.i] = csc.name
          csc[csc == csc.name] = csc.i
          csc.i = csc.i + 1
        }
        csc = matrix(as.numeric(csc), nrow = dim(csc)[1])
        image(csc, col = as.vector(csc.colors), axes = FALSE)
        if (length(colnames(ColSideColors)) > 0) {
          axis(2, 0:(dim(csc)[2] - 1)/max(1,(dim(csc)[2] - 1)), colnames(ColSideColors), las = 2, tick = FALSE)
        }
      }
    }
    
    par(mar = c(margins[1], 0, 0, margins[2]))
    x <- t(x)
    cellnote <- t(cellnote)
    if (revC) {
      iy <- nr:1
      if (exists("ddr"))
        ddr <- rev(ddr)
      x <- x[, iy]
      cellnote <- cellnote[, iy]
    }
    else iy <- 1:nr
    image(1:nc, 1:nr, x, xlim = 0.5 + c(0, nc), ylim = 0.5 + c(0, nr), axes = FALSE, xlab = "", ylab = "", col = col, breaks = breaks, ...)
    retval$carpet <- x
    if (exists("ddr"))
      retval$rowDendrogram <- ddr
    if (exists("ddc"))
      retval$colDendrogram <- ddc
    retval$breaks <- breaks
    retval$col <- col
    if (!invalid(na.color) & any(is.na(x))) { # load library(gplots)
      mmat <- ifelse(is.na(x), 1, NA)
      image(1:nc, 1:nr, mmat, axes = FALSE, xlab = "", ylab = "",
            col = na.color, add = TRUE)
    }
    axis(1, 1:nc, labels = labCol, las = 2, line = -0.5, tick = 0,
         cex.axis = cexCol)
    if (!is.null(xlab))
      mtext(xlab, side = 1, line = margins[1] - 1.25)
    axis(4, iy, labels = labRow, las = 2, line = -0.5, tick = 0,
         cex.axis = cexRow)
    if (!is.null(ylab))
      mtext(ylab, side = 4, line = margins[2] - 1.25)
    if (!missing(add.expr))
      eval(substitute(add.expr))
    if (!missing(colsep))
      for (csep in colsep) rect(xleft = csep + 0.5, ybottom = rep(0, length(csep)), xright = csep + 0.5 + sepwidth[1], ytop = rep(ncol(x) + 1, csep), lty = 1, lwd = 1, col = sepcolor, border = sepcolor)
    if (!missing(rowsep))
      for (rsep in rowsep) rect(xleft = 0, ybottom = (ncol(x) + 1 - rsep) - 0.5, xright = nrow(x) + 1, ytop = (ncol(x) + 1 - rsep) - 0.5 - sepwidth[2], lty = 1, lwd = 1, col = sepcolor, border = sepcolor)
    min.scale <- min(breaks)
    max.scale <- max(breaks)
    x.scaled <- scale01(t(x), min.scale, max.scale)
    if (trace %in% c("both", "column")) {
      retval$vline <- vline
      vline.vals <- scale01(vline, min.scale, max.scale)
      for (i in colInd) {
        if (!is.null(vline)) {
          abline(v = i - 0.5 + vline.vals, col = linecol,
                 lty = 2)
        }
        xv <- rep(i, nrow(x.scaled)) + x.scaled[, i] - 0.5
        xv <- c(xv[1], xv)
        yv <- 1:length(xv) - 0.5
        lines(x = xv, y = yv, lwd = 1, col = tracecol, type = "s")
      }
    }
    if (trace %in% c("both", "row")) {
      retval$hline <- hline
      hline.vals <- scale01(hline, min.scale, max.scale)
      for (i in rowInd) {
        if (!is.null(hline)) {
          abline(h = i + hline, col = linecol, lty = 2)
        }
        yv <- rep(i, ncol(x.scaled)) + x.scaled[i, ] - 0.5
        yv <- rev(c(yv[1], yv))
        xv <- length(yv):1 - 0.5
        lines(x = xv, y = yv, lwd = 1, col = tracecol, type = "s")
      }
    }
    if (!missing(cellnote))
      text(x = c(row(cellnote)), y = c(col(cellnote)), labels = c(cellnote),
           col = notecol, cex = notecex)
    par(mar = c(margins[1], 0, 0, 0))
    if (dendrogram %in% c("both", "row")) {
      plot(ddr, horiz = TRUE, axes = FALSE, yaxs = "i", leaflab = "none")
    }
    else plot.new()
    par(mar = c(0, 0, if (!is.null(main)) 5 else 0, margins[2]))
    if (dendrogram %in% c("both", "column")) {
      plot(ddc, axes = FALSE, xaxs = "i", leaflab = "none")
    }
    else plot.new()
    if (!is.null(main))
      title(main, cex.main = 1.5 * op[["cex.main"]])
    if (key) {
      par(mar = c(5, 4, 2, 1), cex = 0.75)
      tmpbreaks <- breaks
      if (symkey) {
        max.raw <- max(abs(c(x, breaks)), na.rm = TRUE)
        min.raw <- -max.raw
        tmpbreaks[1] <- -max(abs(x), na.rm = TRUE)
        tmpbreaks[length(tmpbreaks)] <- max(abs(x), na.rm = TRUE)
      }
      else {
        min.raw <- min(x, na.rm = TRUE)
        max.raw <- max(x, na.rm = TRUE)
      }
      
      z <- seq(min.raw, max.raw, length = length(col))
      image(z = matrix(z, ncol = 1), col = col, breaks = tmpbreaks,
            xaxt = "n", yaxt = "n")
      par(usr = c(0, 1, 0, 1))
      lv <- pretty(breaks)
      xv <- scale01(as.numeric(lv), min.raw, max.raw)
      axis(1, at = xv, labels = lv)
      if (scale == "row")
        mtext(side = 1, "Row Z-Score", line = 2)
      else if (scale == "column")
        mtext(side = 1, "Column Z-Score", line = 2)
      else mtext(side = 1, KeyValueName, line = 2)
      if (density.info == "density") {
        dens <- density(x, adjust = densadj, na.rm = TRUE)
        omit <- dens$x < min(breaks) | dens$x > max(breaks)
        dens$x <- dens$x[-omit]
        dens$y <- dens$y[-omit]
        dens$x <- scale01(dens$x, min.raw, max.raw)
        lines(dens$x, dens$y/max(dens$y) * 0.95, col = denscol,
              lwd = 1)
        axis(2, at = pretty(dens$y)/max(dens$y) * 0.95, pretty(dens$y))
        title("Color Key\nand Density Plot")
        par(cex = 0.5)
        mtext(side = 2, "Density", line = 2)
      }
      else if (density.info == "histogram") {
        h <- hist(x, plot = FALSE, breaks = breaks)
        hx <- scale01(breaks, min.raw, max.raw)
        hy <- c(h$counts, h$counts[length(h$counts)])
        lines(hx, hy/max(hy) * 0.95, lwd = 1, type = "s",
              col = denscol)
        axis(2, at = pretty(hy)/max(hy) * 0.95, pretty(hy))
        title("Color Key\nand Histogram")
        par(cex = 0.5)
        mtext(side = 2, "Count", line = 2)
      }
      else title("Color Key")
    }
    else plot.new()
    retval$colorTable <- data.frame(low = retval$breaks[-length(retval$breaks)],
                                    high = retval$breaks[-1], color = retval$col)
    invisible(retval)
  }
  
  ### A function for scaling for heatmap
  scale_h <- function(data, type, na.rm=TRUE) {
    
    if(type == "row") {
      scaled <- t(scale(t(data)))
    } else if(type == "col") {
      scaled <- scale(data)
    } else {
      stop("Type is required: row or col")
    }
    
    if(na.rm == TRUE && (length(which(is.na(scaled))) > 0))  {
      scaled <- scaled[-unique(which(is.na(scaled), arr.ind = TRUE)[,1]),]
    }
    
    return(scaled)
  }
  
  ### a function for color brewer
  cell_pal <- function(cell_vars, pal_fun) {
    if (is.numeric(cell_vars)) {
      pal <- pal_fun(100)
      return(pal[cut(cell_vars, breaks = 100)])
    } else {
      categories <- sort(unique(cell_vars))
      pal <- setNames(pal_fun(length(categories)), categories)
      return(pal[cell_vars])
    }
  }
  
  #' @title Plot Slingshot output
  #' @name plot-SlingshotDataSet
  #' @aliases plot-SlingshotDataSet plot,SlingshotDataSet,ANY-method
  #'
  #' @description Tools for visualizing lineages inferred by \code{slingshot}.
  #'
  #' @param x a \code{SlingshotDataSet} with results to be plotted.
  #' @param type character, the type of output to be plotted, can be one of
  #'   \code{"lineages"}, \code{"curves"}, or \code{"both"} (by partial matching),
  #'   see Details for more.
  #' @param linInd integer, an index indicating which lineages should be plotted
  #'   (default is to plot all lineages). If \code{col} is a vector, it will be
  #'   subsetted by \code{linInd}.
  #' @param show.constraints logical, whether or not the user-specified initial
  #'   and terminal clusters should be specially denoted by green and red dots,
  #'   respectively.
  #' @param add logical, indicates whether the output should be added to an
  #'   existing plot.
  #' @param dims numeric, which dimensions to plot (default is \code{1:2}).
  #' @param asp numeric, the y/x aspect ratio, see \code{\link{plot.window}}.
  #' @param cex numeric, amount by which points should be magnified, see
  #'   \code{\link{par}}.
  #' @param lwd numeric, the line width, see \code{\link{par}}.
  #' @param col character or numeric, color(s) for lines, see \code{\link{par}}.
  #' @param ... additional parameters to be passed to \code{\link{lines}}.
  #'
  #' @details If \code{type == 'lineages'}, straight line connectors between
  #'   cluster centers will be plotted. If \code{type == 'curves'}, simultaneous
  #'   principal curves will be plotted.
  #'
  #' @details When \code{type} is not specified, the function will first check the
  #'   \code{curves} slot and plot the curves, if present. Otherwise,
  #'   \code{lineages} will be plotted, if present.
  #'
  #' @return returns \code{NULL}.
  #'
  #' @examples
  #' data("slingshotExample")
  #' rd <- slingshotExample$rd
  #' cl <- slingshotExample$cl
  #' sds <- slingshot(rd, cl, start.clus = "1")
  #' plot(sds, type = 'b')
  #'
  #' # add to existing plot
  #' plot(rd, col = 'grey50')
  #' lines(sds, lwd = 3)
  #'
  #' @import graphics
  #' @import grDevices
  #' @export
  setMethod(
    f = "plot",
    signature = signature(x = "SlingshotDataSet"),
    definition = function(x, type = NULL,
                          linInd = NULL,
                          show.constraints = FALSE,
                          constraints.col = NULL,
                          add = FALSE,
                          dims = seq_len(2),
                          asp = 1,
                          cex = 2,
                          lwd = 2,
                          col = 1,
                          ...) {
      col <- rep(col, length(slingLineages(x)))
      curves <- FALSE
      lineages <- FALSE
      if(is.null(type)){
        if(length(slingCurves(x)) > 0){
          type <- 'curves'
        }else if(length(slingLineages(x)) > 0){
          type <- 'lineages'
        }else{
          stop('No lineages or curves detected.')
        }
      }else{
        type <- c('curves','lineages','both')[pmatch(type,
                                                     c('curves','lineages','both'))]
        if(is.na(type)){
          stop('Unrecognized type argument.')
        }
      }
      
      if(type %in% c('lineages','both')){
        lineages <- TRUE
      }
      if(type %in% c('curves','both')){
        curves <- TRUE
      }
      
      if(lineages & (length(slingLineages(x))==0)){
        stop('No lineages detected.')
      }
      if(curves & (length(slingCurves(x))==0)){
        stop('No curves detected.')
      }
      
      if(is.null(linInd)){
        linInd <- seq_along(slingLineages(x))
      }else{
        linInd <- as.integer(linInd)
        if(!all(linInd %in% seq_along(slingLineages(x)))){
          if(any(linInd %in% seq_along(slingLineages(x)))){
            linInd.removed <-
              linInd[! linInd %in% seq_along(slingLineages(x))]
            linInd <-
              linInd[linInd %in% seq_along(slingLineages(x))]
            message('Unrecognized lineage indices (linInd): ',
                    paste(linInd.removed, collapse = ", "))
          }else{
            stop('None of the provided lineage indices',
                 ' (linInd) were found.')
          }
        }
      }
      
      if(lineages){
        X <- reducedDim(x)
        clusterLabels <- slingClusterLabels(x)
        connectivity <- slingAdjacency(x)
        clusters <- rownames(connectivity)
        nclus <- nrow(connectivity)
        centers <- t(vapply(clusters,function(clID){
          w <- clusterLabels[,clID]
          return(apply(X, 2, weighted.mean, w = w))
        }, rep(0,ncol(X))))
        rownames(centers) <- clusters
        X <- X[rowSums(clusterLabels) > 0, , drop = FALSE]
        clusterLabels <- clusterLabels[rowSums(clusterLabels) > 0, ,
                                       drop = FALSE]
        linC <- slingParams(x)
        clus2include <- unique(unlist(slingLineages(x)[linInd]))
      }
      
      if(!add){
        xs <- NULL
        ys <- NULL
        if(lineages){
          xs <- c(xs, centers[,dims[1]])
          ys <- c(ys, centers[,dims[2]])
        }
        if(curves){
          npoints <- nrow(slingCurves(x)[[1]]$s)
          xs <- c(xs, as.numeric(vapply(slingCurves(x),
                                        function(c){ c$s[,dims[1]] }, rep(0,npoints))))
          ys <- c(ys, as.numeric(vapply(slingCurves(x),
                                        function(c){ c$s[,dims[2]] }, rep(0,npoints))))
        }
        plot(x = NULL, y = NULL, asp = asp,
             xlim = range(xs), ylim = range(ys),
             xlab = colnames(reducedDim(x))[dims[1]],
             ylab = colnames(reducedDim(x))[dims[2]])
      }
      
      if(lineages){
        for(i in seq_len(nclus-1)){
          for(j in seq(i+1,nclus)){
            if(connectivity[i,j]==1 &
               all(clusters[c(i,j)] %in% clus2include)){
              lines(centers[c(i,j), dims],
                    lwd = lwd, col = col[1], ...)
            }
          }
        }
        points(centers[clusters %in% clus2include, dims],
               cex = cex, pch = 16, col = col[1])
        if(show.constraints && !is.null(constraints.col)){
          for(const in names(constraints.col)) {
            points(centers[clusters %in% const, dims,
                           drop=FALSE], cex = cex / 2,
                   col = constraints.col[const], pch = 16)
          }
        }
      }
      if(curves){
        for(ii in seq_along(slingCurves(x))[linInd]){
          c <- slingCurves(x)[[ii]]
          lines(c$s[c$ord, dims], lwd = lwd, col = col[ii], ...)
        }
      }
      invisible(NULL)
    }
  )
  
  #' @title Pairs plot of Slingshot output
  #' @name pairs-SlingshotDataSet
  #'
  #' @description A tool for quickly visualizing lineages inferred by
  #'   \code{slingshot}.
  #'
  #' @param x a \code{SlingshotDataSet} with results to be plotted.
  #' @param type character, the type of output to be plotted, can be one of
  #'   \code{"lineages"}, \code{curves}, or \code{both} (by partial matching), see
  #'   Details for more.
  #' @param show.constraints logical, whether or not the user-specified initial
  #'   and terminal clusters should be specially denoted by green and red dots,
  #'   respectively.
  #' @param col character, color vector for points.
  #' @param pch integer or character specifying the plotting symbol, see
  #'   \code{\link{par}}.
  #' @param cex numeric, amount by which points should be magnified, see
  #'   \code{\link{par}}.
  #' @param lwd numeric, the line width, see \code{\link{par}}.
  #' @param ... additional parameters for \code{plot} or \code{axis}, see
  #'   \code{\link{pairs}}.
  #' @param labels character, the names of the variables, see \code{\link{pairs}}.
  #' @param horInd see \code{\link{pairs}}.
  #' @param verInd see \code{\link{pairs}}.
  #' @param lower.panel see \code{\link{pairs}}.
  #' @param upper.panel see \code{\link{pairs}}.
  #' @param diag.panel see \code{\link{pairs}}.
  #' @param text.panel see \code{\link{pairs}}.
  #' @param label.pos see \code{\link{pairs}}.
  #' @param line.main see \code{\link{pairs}}.
  #' @param cex.labels see \code{\link{pairs}}.
  #' @param font.labels see \code{\link{pairs}}.
  #' @param row1attop see \code{\link{pairs}}.
  #' @param gap see \code{\link{pairs}}.
  #'
  #' @details If \code{type == 'lineages'}, straight line connectors between
  #'   cluster centers will be plotted. If \code{type == 'curves'}, simultaneous
  #'   principal curves will be plotted.
  #'
  #' @details When \code{type} is not specified, the function will first check the
  #'   \code{curves} slot and plot the curves, if present. Otherwise,
  #'   \code{lineages} will be plotted, if present.
  #'
  #' @return returns \code{NULL}.
  #'
  #' @examples
  #' data("slingshotExample")
  #' rd <- slingshotExample$rd
  #' cl <- slingshotExample$cl
  #' sds <- slingshot(rd, cl, start.clus = "1")
  #' pairs(sds, type = 'curves')
  #'
  #' @export
  pairs.SlingshotDataSet <-
    function (x, type = NULL, show.constraints = FALSE, col = NULL,
              constraints.col = NULL,
              pch = 16, cex=1, lwd=2, ...,
              labels, horInd = seq_len(nc), verInd = seq_len(nc),
              lower.panel = FALSE, upper.panel = TRUE,
              diag.panel = NULL, text.panel = textPanel,
              label.pos = 0.5 + has.diag/3, line.main = 3,
              cex.labels = NULL, font.labels = 1,
              row1attop = TRUE, gap = 1) {
      #####
      lp.sling <- lower.panel
      up.sling <- upper.panel
      panel <- points
      if(!up.sling){
        upper.panel <- NULL
      }else{
        upper.panel <- panel
      }
      if(!lower.panel){
        lower.panel <- NULL
      }else{
        lower.panel <- panel
      }
      log = ""
      sds <- x
      x <- reducedDim(sds)
      curves <- FALSE
      lineages <- FALSE
      if(is.null(type)){
        if(length(slingCurves(sds)) > 0){
          type <- 'curves'
        }else if(length(slingLineages(sds)) > 0){
          type <- 'lineages'
        }else{
          stop('No lineages or curves detected.')
        }
      }else{
        type <- c('curves','lineages','both')[pmatch(type,
                                                     c('curves','lineages',
                                                       'both'))]
        if(is.na(type)){
          stop('Unrecognized type argument.')
        }
      }
      if(type %in% c('lineages','both')){
        lineages <- TRUE
      }
      if(type %in% c('curves','both')){
        curves <- TRUE
      }
      if(lineages & (length(slingLineages(sds))==0)){
        stop('No lineages detected.')
      }
      if(curves & (length(slingCurves(sds))==0)){
        stop('No curves detected.')
      }
      if(lineages){
        forest <- slingAdjacency(sds)
        clusters <- rownames(forest)
        nclus <- nrow(forest)
        centers <- t(vapply(clusters,function(clID){
          w <- slingClusterLabels(sds)[,clID]
          return(apply(x, 2, weighted.mean, w = w))
        }, rep(0,ncol(reducedDim(sds)))))
        rownames(centers) <- clusters
        linC <- slingParams(sds)
      }
      range.max <- max(apply(x,2,function(xi){
        r <- range(xi, na.rm = TRUE)
        return(abs(r[2] - r[1]))
      }))
      plot.ranges <- apply(x,2,function(xi){
        mid <- (max(xi,na.rm = TRUE) + min(xi,na.rm = TRUE))/2
        return(c(mid - range.max/2, mid + range.max/2))
      })
      if(is.null(col)){
        if(requireNamespace("RColorBrewer", quietly = TRUE)) {
          cc <- c(RColorBrewer::brewer.pal(9, "Set1")[-c(1,3,6)],
                  RColorBrewer::brewer.pal(7, "Set2")[-2],
                  RColorBrewer::brewer.pal(6, "Dark2")[-5],
                  RColorBrewer::brewer.pal(8, "Set3")[-c(1,2)])
        } else {
          cc <- seq_len(100)
        }
        col <- cc[apply(slingClusterLabels(sds),1,which.max)]
      }
      #####
      if(doText <- missing(text.panel) || is.function(text.panel))
        textPanel <-
        function(x = 0.5, y = 0.5, txt, cex, font)
          text(x, y, txt, cex = cex, font = font)
      
      localAxis <- function(side, x, y, xpd, bg, col=NULL, lwd=NULL, main,
                            oma, ...) {
        ## Explicitly ignore any color argument passed in as
        ## it was most likely meant for the data points and
        ## not for the axis.
        xpd <- NA
        if(side %% 2L == 1L && xl[j]) xpd <- FALSE
        if(side %% 2L == 0L && yl[i]) xpd <- FALSE
        if(side %% 2L == 1L) Axis(x, side = side, xpd = xpd, ...)
        else Axis(y, side = side, xpd = xpd, ...)
      }
      
      localPlot <- function(..., main, oma, font.main, cex.main) plot(...)
      localLowerPanel <- function(..., main, oma, font.main, cex.main)
        lower.panel(...)
      localUpperPanel <- function(..., main, oma, font.main, cex.main)
        upper.panel(...)
      localDiagPanel <- function(..., main, oma, font.main, cex.main)
        diag.panel(...)
      
      dots <- list(...); nmdots <- names(dots)
      if (!is.matrix(x)) {
        x <- as.data.frame(x)
        for(i in seq_along(names(x))) {
          if(is.factor(x[[i]]) || is.logical(x[[i]]))
            x[[i]] <- as.numeric(x[[i]])
          if(!is.numeric(unclass(x[[i]])))
            stop("non-numeric argument to 'pairs'")
        }
      } else if (!is.numeric(x)) stop("non-numeric argument to 'pairs'")
      panel <- match.fun(panel)
      if((has.lower <- !is.null(lower.panel)) && !missing(lower.panel))
        lower.panel <- match.fun(lower.panel)
      if((has.upper <- !is.null(upper.panel)) && !missing(upper.panel))
        upper.panel <- match.fun(upper.panel)
      if((has.diag  <- !is.null( diag.panel)) && !missing( diag.panel))
        diag.panel <- match.fun( diag.panel)
      
      if(row1attop) {
        tmp <- lower.panel; lower.panel <- upper.panel; upper.panel <- tmp
        tmp <- has.lower; has.lower <- has.upper; has.upper <- tmp
      }
      
      nc <- ncol(x)
      if (nc < 2L) stop("only one column in the argument to 'pairs'")
      if(!all(horInd >= 1L & horInd <= nc))
        stop("invalid argument 'horInd'")
      if(!all(verInd >= 1L & verInd <= nc))
        stop("invalid argument 'verInd'")
      if(doText) {
        if (missing(labels)) {
          labels <- colnames(x)
          if (is.null(labels)) labels <- paste("var", 1L:nc)
        }
        else if(is.null(labels)) doText <- FALSE
      }
      oma <- if("oma" %in% nmdots) dots$oma
      main <- if("main" %in% nmdots) dots$main
      if (is.null(oma))
        oma <- c(4, 4, if(!is.null(main)) 6 else 4, 4)
      opar <- par(mfrow = c(length(horInd), length(verInd)),
                  mar = rep.int(gap/2, 4), oma = oma)
      on.exit(par(opar))
      dev.hold(); on.exit(dev.flush(), add = TRUE)
      
      xl <- yl <- logical(nc)
      if (is.numeric(log)) xl[log] <- yl[log] <- TRUE
      else {xl[] <- grepl("x", log); yl[] <- grepl("y", log)}
      for (i in if(row1attop) verInd else rev(verInd))
        for (j in horInd) {
          l <- paste0(ifelse(xl[j], "x", ""), ifelse(yl[i], "y", ""))
          localPlot(x[, j], x[, i], xlab = "", ylab = "",
                    axes = FALSE, type = "n", ..., log = l,
                    xlim = plot.ranges[,j], ylim = plot.ranges[,i])
          if(i == j || (i < j && has.lower) || (i > j && has.upper) ) {
            box()
            if(i == 1  && (!(j %% 2L) || !has.upper || !has.lower ))
              localAxis(1L + 2L*row1attop, x[, j], x[, i], ...)
            if(i == nc && (  j %% 2L  || !has.upper || !has.lower ))
              localAxis(3L - 2L*row1attop, x[, j], x[, i], ...)
            if(j == 1  && (!(i %% 2L) || !has.upper || !has.lower ))
              localAxis(2L, x[, j], x[, i], ...)
            if(j == nc && (  i %% 2L  || !has.upper || !has.lower ))
              localAxis(4L, x[, j], x[, i], ...)
            mfg <- par("mfg")
            if(i == j) {
              if (has.diag) localDiagPanel(as.vector(x[, i]), ...)
              if (doText) {
                par(usr = c(0, 1, 0, 1))
                if(is.null(cex.labels)) {
                  l.wid <- strwidth(labels, "user")
                  cex.labels <- max(0.8, min(2, .9 / max(l.wid)))
                }
                xlp <- if(xl[i]) 10^0.5 else 0.5
                ylp <- if(yl[j]) 10^label.pos else label.pos
                text.panel(xlp, ylp, labels[i],
                           cex = cex.labels, font = font.labels)
              }
            } else if(i < j){
              if(up.sling){
                points(as.vector(x[, j]), as.vector(x[, i]),
                       col = col, cex = cex, pch=pch, ...)
                if(lineages){
                  for(ii in seq_len(nclus-1)){
                    for(jj in seq(ii+1,nclus)){
                      if(forest[ii,jj]==1){
                        seg.col <- 1
                        lines(centers[c(ii,jj),j],
                              centers[c(ii,jj),i],
                              lwd = lwd, col = seg.col, ...)
                      }
                    }
                  }
                  points(centers[,j],centers[,i], pch = pch,
                         cex=2*cex)
                  if(show.constraints && is.null(constraints.col)){
                    if(any(linC$start.given)){
                      st.ind <- clusters %in%
                        linC$start.clus[linC$start.given]
                      points(centers[st.ind,j],
                             centers[st.ind,i], cex = cex,
                             col = 'green3',
                             pch = pch)
                    }
                    if(any(linC$end.given)){
                      en.ind <- clusters %in%
                        linC$end.clus[linC$end.given]
                      points(centers[en.ind,j],
                             centers[en.ind,i], cex = cex,
                             col = 'red2', pch = pch)
                    }
                  } else if(show.constraints && !is.null(constraints.col)){
                    for(const in names(constraints.col)) {
                      points(centers[clusters %in% const, j, drop=FALSE],
                             centers[clusters %in% const, i, drop=FALSE],
                             cex = cex, pch = 16,
                             col = constraints.col[const])
                    }
                  }
                }
                if(curves){
                  for(c in slingCurves(sds)){
                    lines(c$s[c$ord,c(j,i)], lwd = lwd,
                          col=1, ...)
                  }
                }
              }
            }
            else{
              if(lp.sling){
                points(as.vector(x[, j]), as.vector(x[, i]),
                       col = col, cex = cex, pch=pch, ...)
                if(lineages){
                  for(ii in seq_len(nclus-1)){
                    for(jj in seq(ii+1,nclus)){
                      if(forest[ii,jj]==1){
                        if(clusters[ii] %in%
                           linC$start.clus |
                           clusters[jj] %in%
                           linC$start.clus){
                          seg.col <- 'green3'
                        }else if(clusters[ii] %in%
                                 linC$end.clus[
                                   linC$end.given] |
                                 clusters[jj] %in%
                                 linC$end.clus[
                                   linC$end.given]){
                          seg.col <- 'red2'
                        }else{
                          seg.col <- 1
                        }
                        lines(centers[c(ii,jj),j],
                              centers[c(ii,jj),i],
                              lwd = lwd, col = seg.col,...)
                      }
                    }
                  }
                  points(centers[,j],centers[,i], pch = pch,
                         cex = 2*cex)
                }
                if(curves){
                  for(c in slingCurves(sds)){
                    lines(c$s[c$ord,c(j,i)],lwd = lwd,
                          col=1, ...)
                  }
                }
              }
            }
            if (any(par("mfg") != mfg))
              stop("the 'panel' function made a new plot")
          } else par(new = FALSE)
          
        }
      if (!is.null(main)) {
        font.main <- if("font.main" %in% nmdots){
          dots$font.main
        }else par("font.main")
        cex.main <- if("cex.main" %in%
                       nmdots) dots$cex.main else par("cex.main")
        mtext(main, 3, line.main, outer=TRUE, at = 0.5, cex = cex.main,
              font = font.main)
      }
      invisible(NULL)
    }
  
  #' @name plot3d-SlingshotDataSet
  #' @title Plot Slingshot output in 3D
  #'
  #' @description Tools for visualizing lineages inferred by \code{slingshot}.
  #'
  #' @param x a \code{SlingshotDataSet} with results to be plotted.
  #' @param type character, the type of output to be plotted, can be one of
  #'   \code{"lineages"}, \code{curves}, or \code{both} (by partial matching), see
  #'   Details for more.
  #' @param linInd integer, an index indicating which lineages should be plotted
  #'   (default is to plot all lineages). If \code{col} is a vector, it will be
  #'   subsetted by \code{linInd}.
  #' @param add logical, indicates whether the output should be added to an
  #'   existing plot.
  #' @param dims numeric, which dimensions to plot (default is \code{1:3}).
  #' @param aspect either a logical indicating whether to adjust the aspect ratio
  #'   or a new ratio, see \code{\link[rgl:plot3d]{plot3d}}.
  #' @param size numeric, size of points for MST (default is \code{10}), see
  #'   \code{\link[rgl:plot3d]{plot3d}}.
  #' @param col character or numeric, color(s) for lines, see \code{\link{par}}.
  #' @param ... additional parameters to be passed to \code{lines3d}.
  #'
  #' @details If \code{type == 'lineages'}, straight line connectors between
  #'   cluster centers will be plotted. If \code{type == 'curves'}, simultaneous
  #'   principal curves will be plotted.
  #'
  #' @details When \code{type} is not specified, the function will first check the
  #'   \code{curves} slot and plot the curves, if present. Otherwise,
  #'   \code{lineages} will be plotted, if present.
  #'
  #' @return returns \code{NULL}.
  #'
  #' @examples
  #' \dontrun{
  #' library(rgl)
  #' data("slingshotExample")
  #' rd <- slingshotExample$rd
  #' cl <- slingshotExample$cl
  #' rd <- cbind(rd, rnorm(nrow(rd)))
  #' sds <- slingshot(rd, cl, start.clus = "1")
  #' plot3d(sds, type = 'b')
  #'
  #' # add to existing plot
  #' plot3d(rd, col = 'grey50', aspect = 'iso')
  #' plot3d(sds, lwd = 3, add = TRUE)
  #' }
  # #' @importFrom rgl plot3d
  #' @export
  plot3d.SlingshotDataSet <- function(x,
                                      type = NULL,
                                      linInd = NULL,
                                      add = FALSE,
                                      dims = seq_len(3),
                                      aspect = 'iso',
                                      size = 10,
                                      col = 1,
                                      col2 = NULL,
                                      ...){
    if (!requireNamespace("rgl", quietly = TRUE)) {
      stop("Package 'rgl' is required for 3D plotting.",
           call. = FALSE)
    }
    col <- rep(col, length(slingLineages(x)))
    n <- nrow(reducedDim(x))
    curves <- FALSE
    lineages <- FALSE
    if(is.null(type)){
      if(length(slingCurves(x)) > 0){
        type <- 'curves'
      }else if(length(slingLineages(x)) > 0){
        type <- 'lineages'
      }else{
        stop('No lineages or curves detected.')
      }
    }else{
      type <- c('curves','lineages','both')[pmatch(type,c('curves','lineages',
                                                          'both'))]
      if(is.na(type)){
        stop('Unrecognized type argument.')
      }
    }
    
    if(type %in% c('lineages','both')){
      lineages <- TRUE
    }
    if(type %in% c('curves','both')){
      curves <- TRUE
    }
    
    if(lineages & (length(slingLineages(x))==0)){
      stop('No lineages detected.')
    }
    if(curves & (length(slingCurves(x))==0)){
      stop('No curves detected.')
    }
    
    if(is.null(linInd)){
      linInd <- seq_along(slingLineages(x))
    }else{
      linInd <- as.integer(linInd)
      if(!all(linInd %in% seq_along(slingLineages(x)))){
        if(any(linInd %in% seq_along(slingLineages(x)))){
          linInd.removed <-
            linInd[! linInd %in% seq_along(slingLineages(x))]
          linInd <-
            linInd[linInd %in% seq_along(slingLineages(x))]
          message('Unrecognized lineage indices (linInd): ',
                  paste(linInd.removed, collapse = ", "))
        }else{
          stop('None of the provided lineage indices',
               ' (linInd) were found.')
        }
      }
    }
    
    if(lineages){
      X <- reducedDim(x)
      clusterLabels <- slingClusterLabels(x)
      connectivity <- slingAdjacency(x)
      clusters <- rownames(connectivity)
      nclus <- nrow(connectivity)
      centers <- t(vapply(clusters,function(clID){
        w <- clusterLabels[,clID]
        return(apply(X, 2, weighted.mean, w = w))
      }, rep(0,ncol(X))))
      rownames(centers) <- clusters
      X <- X[rowSums(clusterLabels) > 0, , drop = FALSE]
      clusterLabels <- clusterLabels[rowSums(clusterLabels) > 0, ,
                                     drop = FALSE]
      clus2include <- unique(unlist(slingLineages(x)[linInd]))
    }
    
    if(!add){
      xs <- NULL
      ys <- NULL
      zs <- NULL
      if(lineages){
        xs <- c(xs, centers[,dims[1]])
        ys <- c(ys, centers[,dims[2]])
        zs <- c(zs, centers[,dims[3]])
      }
      if(curves){
        npoints <- nrow(slingCurves(x)[[1]]$s)
        xs <- c(xs, as.numeric(vapply(slingCurves(x), function(c){
          c$s[,dims[1]] }, rep(0,npoints))))
        ys <- c(ys, as.numeric(vapply(slingCurves(x), function(c){
          c$s[,dims[2]] }, rep(0,npoints))))
        zs <- c(zs, as.numeric(vapply(slingCurves(x), function(c){
          c$s[,dims[3]] }, rep(0,npoints))))
      }
      rgl::plot3d(x = NULL, y = NULL, z = NULL, aspect = aspect,
                  xlim = range(xs), ylim = range(ys), zlim = range(zs),
                  xlab = colnames(reducedDim(x))[dims[1]],
                  ylab = colnames(reducedDim(x))[dims[2]],
                  zlab = colnames(reducedDim(x))[dims[3]])
    }
    
    if(lineages){
      for(i in seq_len(nclus-1)){
        for(j in seq(i+1,nclus)){
          if(connectivity[i,j]==1 &
             all(clusters[c(i,j)] %in% clus2include)){
            rgl::lines3d(x = centers[c(i,j),dims[1]],
                         y = centers[c(i,j),dims[2]],
                         z = centers[c(i,j),dims[3]],
                         col = col[1], ...)
          }
        }
      }
      rgl::points3d(centers[clusters %in% clus2include, dims],
                    size = size/2, col = col2[clusters[clusters %in% clus2include]])
      rgl::points3d(centers[clusters %in% clus2include, dims],
                    size = size, col = col[1])
    }
    if(curves){
      for(ii in seq_along(slingCurves(x))[linInd]){
        c <- slingCurves(x)[[ii]]
        rgl::lines3d(c$s[c$ord,dims], col = col[ii], ...)
      }
    }
    invisible(NULL)
  }
  
  ### the plot3d.SlingshotDataSet of the slingshot package is incomplete and too simple,
  ### so, i'm implementing a 3d plot function myself
  slingshot_3d_lineages <- function(slingshot_obj, color, title,
                                    print=FALSE, outputDir=NULL,
                                    width=1200, height=800) {
    
    ### load libraries
    if(!require(Seurat, quietly = TRUE)) {
      install.packages("Seurat")
      require(Seurat, quietly = TRUE)
    }
    if(!require(slingshot, quietly = TRUE)) {
      if (!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager")
      BiocManager::install("slingshot")
      require(slingshot, quietly = TRUE)
    }
    if(!require(rgl, quietly = TRUE)) {
      install.packages("rgl")
      require(rgl, quietly = TRUE)
    }
    if(!require(rmarkdown, quietly = TRUE)) {
      install.packages("rmarkdown")
      require(rmarkdown, quietly = TRUE)
    }
    
    #
    ### 3D Slingshot
    #
    
    ### draw 3D PCA
    par3d(windowRect = c(50, 50, width+50, height+50))
    plot3d.SlingshotDataSet(slingshot_obj, dims = 1:3, col = "black", col2 = color, type = "lineages", add = TRUE)
    plot3d(slingshot_obj@reducedDim, col = apply(slingshot_obj@clusterLabels, 1, function(x) color[names(x)[which(x == 1)]]),
           size = 5, alpha = 0.5, aspect = FALSE, add = TRUE)
    axes3d(edges=c("x+-", "y+-", "z++"), lwd = 2,
           labels=TRUE, tick = FALSE, nticks = 3, box = TRUE, expand = 1.05)
    mtext3d(text = expression(bold("PC1")), edge="x+-", line = -2, at = min(slingshot_obj@reducedDim[,1]), pos = NA)
    mtext3d(text = expression(bold("PC2")), edge="y+-", line = -2, at = min(slingshot_obj@reducedDim[,2]), pos = NA)
    mtext3d(text = expression(bold("PC3")), edge="z++", line = -2, at = max(slingshot_obj@reducedDim[,3]), pos = NA)
    decorate3d(xlim = NULL, ylim = NULL, zlim = NULL, 
               xlab = "", ylab = "", zlab = "", 
               box = FALSE, axes = FALSE, main = title, sub = NULL,
               top = TRUE, aspect = FALSE, expand = 1.05)
    legend3d("topright", legend = names(color), title = "Clusters",
             col = color, pch = 19, cex=2)
    if(print) {
      writeWebGL(dir=outputDir, filename = paste0(outputDir, title, ".html"),
                 width=width, height = height)
    }
    
  }
  
  #'****************************************************************************************
  #' Compare different ranks of the same elements
  #' 
  #' It receives (1) two named numeric vectors that contain two different rank metrics
  #' OR (2) receives two character vectors of the ordered/ranked element names
  #' 
  #' In case of (1), if the two given vectors are A and B,
  #' names(A) and names(B) should have the same elements.
  #' length(intersect(names(A), names(B))) = length(A) = length(B).
  #' In case of (2), if the two given vectors are A and B,
  #' A and B should have the same elements.
  #' length(intersect(A, B)) = length(A) = length(B).
  #' 
  #' The function will generate 4 figures in one plot.
  #' If the two given vectors are A and B,
  #' 1. a barplot of ranked metric of A
  #' 2. a line plot that shows lines connect top same elements between ranked A and ranked B
  #' 3. a barplot of ranked metric of B
  #' 4. a scatter plot of A and B
  #'****************************************************************************************
  #' @title	compare_two_different_ranks
  #' 
  #' @param A           A vector of named numeric vector or a character vector
  #' @param B           A vector of named numeric vector or a character vector
  #' @param A_name      The title of the given rank A
  #'                    (Default = "A")      
  #' @param B_name      The title of the given rank B
  #'                    (Default = "B")
  #' @param ordering    Order of the rank metric ["decreasing" or "increasing"]
  #'                    (Default = "decreasing")
  #' @param top         The number of top hubs from each of A and B that will be
  #'                    presented in the line graph
  #'                    (Default = 100)
  #' @param alternative The alternative hypothesis for enrichment permutation ["greater" or "less"]
  #'                    (Default = "greater")
  #' @param permutation The number of permutation test
  #'                    (Default = 10000)
  #' @param fileName    The name of the plot file
  #'                    (Default = "Comparison_of_Ranks_Between_A_and_B")
  #' @param printPath   Print the plot in the designated path
  #'                    (Default = "./")
  #' @param width       The width of the plot file
  #'                    (Default = 24)
  #' @param height      The height of the plot file
  #'                    (Default = 12)
  #' 
  #' @return 	          It does not return any objects
  #' 
  compare_two_different_ranks <- function(A,
                                          B,
                                          A_name = "A",
                                          B_name = "B",
                                          ordering = "decreasing",
                                          top = 100,
                                          alternative = "greater",
                                          permutation = 10000,
                                          fileName = "Comparison_of_Ranks_Between_A_and_B",
                                          printPath = "./",
                                          width = 24,
                                          height = 10) {
    
    ### load required libraries
    if(!require(ggplot2, quietly = TRUE)) {
      install.packages("ggplot2")
      require(ggplot2, quietly = TRUE)
    }
    if(!require(gridExtra, quietly = TRUE)) {
      install.packages("gridExtra")
      require(gridExtra, quietly = TRUE)
    }
    if(!require("checkmate", quietly = TRUE)) {
      install.packages("checkmate")
      require("checkmate", quietly = TRUE)
    }
    
    ### argument checking
    assertVector(A)
    assertVector(B)
    assertString(A_name)
    assertString(B_name)
    assertChoice(ordering, c("decreasing", "increasing"))
    assertIntegerish(top)
    assertChoice(alternative, c("greater", "less"))
    assertIntegerish(permutation)
    assertString(fileName)
    assertString(printPath)
    assertIntegerish(width)
    assertIntegerish(height)
    if(length(A) != length(B)) {
      stop("ERROR: length(A) != length(B)")
    }
    if(class(A) == "numeric" && class(B) == "numeric") {
      if(length(intersect(names(A), names(B))) != length(A)) {
        stop("ERROR: A and B should have names, and names(sort(A)) and names(sort(B)) should be the same")
      }
      type <- "numeric"
    } else if(class(A) == "character" && class(B) == "character") {
      if(length(intersect(A, B)) != length(A)) {
        stop("ERROR: sort(A) and sort(B) should be the same")
      }
      type <- "character"
      names(A) <- A
      A <- length(A):1
      names(B) <- B
      B <- length(B):1
    } else {
      stop("ERROR: A and B should be either named numeric vectors or character vectors")
    }
    
    ### order the A and B
    if(ordering == "decreasing") {
      A <- A[order(-A)]
      B <- B[order(-B)]
    } else {
      A <- A[order(A)]
      B <- B[order(B)]
    }
    
    ### data frame of A and B
    df <- data.frame(idx=1:length(A), A=A[order(names(A))], B=B[order(names(B))],
                     stringsAsFactors = FALSE, check.names = FALSE)
    
    ### give colors to the df
    colors <- rep("lightgray", nrow(df))
    names(colors) <- rownames(df)
    colors[names(A)[1:top]] <- "skyblue"
    colors[names(B)[1:top]] <- "pink"
    colors[intersect(names(A)[1:top], names(B)[1:top])] <- "mediumpurple"
    df$colors <- factor(colors, levels = c("skyblue", "pink", "lightgray", "mediumpurple"))
    
    ### data frame for the line graph
    df2 <- data.frame(X=c(Reduce(function(x, y) {
      return(c(x, y, which(names(B) == names(A)[y])))
    }, 0:top)[-1],
    Reduce(function(x, y) {
      return(c(x, y, which(names(A) == names(B)[y])))
    }, 0:top)[-1]),
    Y=c(rep(c(1,0), top),
        rep(c(0,1), top)),
    Element=paste0("e", c(rbind(1:(top*2), 1:(top*2)))),
    Collection=c(rep("A", top*2),
                 rep("B", top*2)),
    stringsAsFactors = FALSE, check.names = FALSE)
    
    ### add top and bottom x-axis to the df2
    df2 <- rbind(df2, data.frame(X=c(0, length(A), 0, length(B)),
                                 Y=c(1, 1, 0, 0),
                                 Element=paste0("e", c(rbind((top*2+1):(top*2+2),(top*2+1):(top*2+2)))),
                                 Collection=rep("C", 4)))
    
    ### permutation test for getting a p-value (1-tail)
    set.seed(1234)
    permu_result <- sapply(1:(permutation-1), function(x) {
      random_hubs <- names(A)[sample(length(A), top)]
      return(sum(B[random_hubs]))
    })
    if(alternative == "greater") {
      pVal_A <- (length(which(permu_result > sum(B[names(A)[1:top]])))+1) / permutation
    } else {
      pVal_A <- (length(which(permu_result < sum(B[names(A)[1:top]])))+1) / permutation
    }
    permu_result <- sapply(1:(permutation-1), function(x) {
      random_hubs <- names(B)[sample(length(B), top)]
      return(sum(A[random_hubs]))
    })
    if(alternative == "greater") {
      pVal_B <- (length(which(permu_result > sum(A[names(B)[1:top]])))+1) / permutation
    } else {
      pVal_B <- (length(which(permu_result < sum(A[names(B)[1:top]])))+1) / permutation
    }
    
    if(ordering == "decreasing") {
      ### 1. a barplot of ranked metric of A
      p1 <- ggplot(data = df, aes(x = reorder(idx, -A), y = A, width = 1)) +
        ylab(paste(A_name, "(A)")) +
        theme(axis.title.x=element_blank(),
              axis.text.x=element_blank(),
              axis.ticks.x=element_blank(),
              plot.margin = margin(t = 0, r = 115, b = 0, l = 5)) +
        geom_bar(stat = "identity", color = "skyblue")
      
      ### 3. a barplot of ranked metric of B
      p3 <- ggplot(data = df, aes(x = reorder(idx, -B), y = B, width = 1)) +
        ylab(paste(B_name, "(B)")) +
        theme(axis.title.x=element_blank(),
              axis.text.x = element_blank(),
              axis.ticks.x=element_blank(),
              plot.margin = margin(t = 0, r = 115, b = 0, l = 5)) +
        geom_bar(stat = "identity", color = "pink")
    } else {
      ### 1. a barplot of ranked metric of A
      p1 <- ggplot(data = df, aes(x = reorder(idx, A), y = A, width = 1)) +
        ylab(paste(A_name, "(A)")) +
        theme(axis.title.x=element_blank(),
              axis.text.x=element_blank(),
              axis.ticks.x=element_blank(),
              plot.margin = margin(t = 0, r = 115, b = 0, l = 5)) +
        geom_bar(stat = "identity", color = "skyblue")
      
      ### 3. a barplot of ranked metric of B
      p3 <- ggplot(data = df, aes(x = reorder(idx, B), y = B, width = 1)) +
        ylab(paste(B_name, "(B)")) +
        theme(axis.title.x=element_blank(),
              axis.text.x = element_blank(),
              axis.ticks.x=element_blank(),
              plot.margin = margin(t = 0, r = 115, b = 0, l = 5)) +
        geom_bar(stat = "identity", color = "pink")
    }
    
    ### 2. a line plot that shows lines connect top same elements between ranked A and ranked B
    p2 <- ggplot(data = df2, aes(x = X, y = Y, group = Element, color = Collection)) +
      geom_line() +
      scale_color_manual(labels = c("A", "B", "C"), values = c("skyblue", "pink", "gray")) +
      theme_void() +
      theme(legend.position = "none", plot.margin = margin(t = 0, r = 5, b = 0, l = 5)) +
      annotate(geom = "text", x = length(A), y = 0.1,
               color = "skyblue", hjust = 0, vjust = 0, fontface = "bold",
               label = paste("Top", top, "(A)\nPermutation p-value\n=", pVal_A)) +
      annotate(geom = "text", x = length(B), y = 0.9,
               color = "pink", hjust = 0, vjust = 1, fontface = "bold",
               label = paste("Top", top, "(B)\nPermutation p-value\n=", pVal_B)) +
      expand_limits(x = length(A)*1.1)
    
    ### 4. a scatter plot of A and B
    p4 <- ggplot(data = df, aes(x=A, y=B)) +
      geom_point(aes(color=colors), size = 2) +
      labs(subtitle=paste0("Pearson Correlation = ", round(cor(df$A, df$B), 5),
                          " (P-value = ", round(cor.test(df$A, df$B)$p.value, 5), ")"),
           color=paste("Top", top)) +
      xlab(A_name) +
      ylab(B_name) +
      geom_smooth(method = lm, color="black", se=FALSE) +
      scale_color_manual(labels=c("A", "B", "Others", "Intersect(A, B)"), values = c("skyblue", "pink", "lightgray", "mediumpurple")) +
      theme_classic(base_size = 16)
    
    ### arrange the plots and print out
    g <- arrangeGrob(p1, p2, p3, p4, layout_matrix = rbind(c(1, 4),
                                                           c(2, 4),
                                                           c(3, 4)),
                     top = fileName)
    ggsave(file = paste0(printPath, fileName, ".png"), g, width = width, height = height)
    
  }
  
  ### for each donor run the same analyses
  for(donor in unique(Seurat_Obj@meta.data$Donor)) {
    
    ### new output directory
    outputDir2 <- paste0(outputDir, donor, "/")
    dir.create(outputDir2, showWarnings = FALSE, recursive = TRUE)
    
    ### check whether the orders are the same
    print(identical(names(Idents(object = Seurat_Obj)), rownames(Seurat_Obj@meta.data)))
    
    ### split the Seurat obj based on donor info
    Seurat_Obj <- SetIdent(object = Seurat_Obj,
                           cells = rownames(Seurat_Obj@meta.data),
                           value = Seurat_Obj@meta.data$Donor)
    subset_Seurat_Obj <- subset(Seurat_Obj, idents=donor)
    
    ### order the meta data by time
    subset_Seurat_Obj@meta.data <- subset_Seurat_Obj@meta.data[order(subset_Seurat_Obj@meta.data$Day),]
    
    ### rownames in the meta.data should be in the same order as colnames in the counts
    subset_Seurat_Obj@assays$RNA@counts <- subset_Seurat_Obj@assays$RNA@counts[,rownames(subset_Seurat_Obj@meta.data)]
    
    ### make a clone summary table
    unique_clone_idx <- which(!duplicated(subset_Seurat_Obj@meta.data$clone_id))
    clone_summary_table <- data.frame(clone_id=subset_Seurat_Obj@meta.data$clone_id[unique_clone_idx],
                                      cdr_ab=subset_Seurat_Obj@meta.data$match.cdr[unique_clone_idx])
    rownames(clone_summary_table) <- clone_summary_table$clone_id
    
    ### add time point counts and the total count of the clonotypes
    time_points <- c("d0", "d5", "d12", "d28", "d60", "d90", "d120", "d180")
    clone_summary_table[time_points] <- 0
    clone_summary_table$total_count <- 0
    
    ### fill out the counts
    tp_indicies <- lapply(time_points, function(x) which(subset_Seurat_Obj@meta.data$Day == x))
    names(tp_indicies) <- time_points
    for(clone in rownames(clone_summary_table)) {
      clone_idx <- which(subset_Seurat_Obj@meta.data$clone_id == clone)
      for(tp in time_points) {
        clone_summary_table[clone,tp] <- length(intersect(clone_idx, tp_indicies[[tp]]))
      }
    }
    clone_summary_table$total_count <- as.numeric(apply(clone_summary_table[,time_points], 1, sum))
    
    ### order by the total_count
    clone_summary_table <- clone_summary_table[order(-clone_summary_table$total_count),]
    
    ### save the table as Excel file
    write.xlsx2(clone_summary_table, file = paste0(outputDir2, donor, "_Clone_Count_Summary.xlsx"),
                sheetName = "CLONE_SUMMARY", row.names = FALSE)
    
    
    ### Alluvial plot - visualization of the lineage tracing
    
    ### only keep lineages from the clone_summary_table
    ### * lineage = a clone that appears in more than one (2, 3, ...) time points
    lineage_table <- clone_summary_table[apply(clone_summary_table[,time_points], 1, function(x) {
      return(length(which(x > 0)) > 1)
    }),]
    
    ### get an input data frame for the alluvial plot
    total_rows <- length(which(lineage_table[,time_points] > 0))
    plot_df <- data.frame(Time=rep("", total_rows),
                          Clone_Size=rep(0, total_rows),
                          Clone=rep("", total_rows),
                          CDR3=rep("", total_rows))
    cnt <- 1
    for(i in 1:nrow(lineage_table)) {
      for(tp in time_points) {
        if(lineage_table[i,tp] > 0) {
          plot_df[cnt,] <- c(tp,
                             lineage_table[i,tp],
                             lineage_table$clone_id[i],
                             lineage_table$cdr_ab[i])
          cnt <- cnt + 1
        }
      }
    }
    plot_df$Time <- factor(plot_df$Time, levels = time_points)
    
    ### numerize the clone_size column
    plot_df$Clone_Size <- as.numeric(plot_df$Clone_Size)
    
    ### draw the alluvial plot
    ggplot(plot_df,
           aes(x = Time, stratum = Clone, alluvium = Clone,
               y = Clone_Size,
               fill = Clone, label = CDR3)) +
      ggtitle(paste("Clonal Tracing of the TFH-related Cells", donor)) +
      geom_flow() +
      geom_stratum(alpha = 1) +
      geom_text(stat = "stratum", size = 0.8) +
      rotate_x_text(90) +
      theme_pubr(legend = "none") +
      theme(axis.title.x = element_blank()) +
      theme_cleveland2() +
      scale_fill_viridis(discrete = T) +
      scale_y_continuous(expand = c(0, 0), limits = c(0, NA))
    ggsave(file = paste0(outputDir2, donor, "_TFH_Clonal_Tracing.png"), width = 18, height = 9, dpi = 300)
    
    
    #
    ### PCA & UMAP with the top 9 clones from the TFH result
    #
    
    ### run PCA
    subset_Seurat_Obj <- RunPCA(subset_Seurat_Obj, npcs = 10)
    
    ### run UMAP
    subset_Seurat_Obj <- RunUMAP(subset_Seurat_Obj, dims = 1:5)
    
    ### if the donor is 321-05, flip the sign of the PC1 and PC2 for consistency with the previous results
    if(donor == "321-05") {
      subset_Seurat_Obj@reductions$pca@cell.embeddings <- -subset_Seurat_Obj@reductions$pca@cell.embeddings
      subset_Seurat_Obj@reductions$pca@feature.loadings <- -subset_Seurat_Obj@reductions$pca@feature.loadings
    }
    
    ### check whether the orders are the same
    print(identical(names(Idents(object = subset_Seurat_Obj)), rownames(subset_Seurat_Obj@meta.data)))
    
    ### subset preparation for the PCA & UMAP
    subset_Seurat_Obj <- SetIdent(object = subset_Seurat_Obj,
                                  cells = rownames(subset_Seurat_Obj@meta.data),
                                  value = subset_Seurat_Obj@meta.data$clone_id)
    subset_Seurat_Obj2 <- subset(subset_Seurat_Obj, idents=lineage_table$clone_id[1:9])
    subset_Seurat_Obj2@meta.data$Day <- factor(subset_Seurat_Obj2@meta.data$Day,
                                              levels = c("d0", "d5", "d12", "d28", "d60", "d90", "d120", "d180"))
    
    ### draw the PCA & UMAP
    DimPlot(subset_Seurat_Obj2, reduction = "pca", split.by = "clone_id", group.by = "Day", shape.by = "Tissue",
            ncol = 3, pt.size = 2) +
      labs(title = paste("PCA of the Top 9 Clones", donor))
    ggsave(file = paste0(outputDir2, donor, "_PCA_Top_9_Clones_TFH.png"), width = 20, height = 10, dpi = 300)
    DimPlot(subset_Seurat_Obj2, reduction = "umap", split.by = "clone_id", group.by = "Day", shape.by = "Tissue",
            ncol = 3, pt.size = 2) +
      labs(title = paste("UMAP of the Top 9 Clones", donor))
    ggsave(file = paste0(outputDir2, donor, "_UMAP_Top_9_Clones_TFH.png"), width = 20, height = 10, dpi = 300)
    
    
    ### set clone size threshold
    ### clones with clone size > clone size threshold will be used
    clone_size_thresh <- 5
    
    ### DE gene FDR threshold
    de_signif_thresh <- 0.05
    
    ### get clone names that will be used
    clone_names <- clone_summary_table$clone_id[which(clone_summary_table$total_count > clone_size_thresh)]
    
    ### set new output directory
    outputDir3 <- paste0(outputDir2, "Clone_Markers/")
    dir.create(path = outputDir3, showWarnings = FALSE, recursive = TRUE)
    
    ### for each interesting clone perform marker discovery
    for(clone in clone_names) {
      
      ### print progress
      writeLines(paste(clone))
      
      ### set new output directory
      outputDir4 <- paste0(outputDir3, clone, "/")
      dir.create(path = outputDir4, showWarnings = FALSE, recursive = TRUE)
      
      ### Ident configure
      new.ident <- rep(NA, nrow(subset_Seurat_Obj@meta.data))
      new.ident[which(subset_Seurat_Obj@meta.data$clone_id == clone)] <- "ident1"
      new.ident[which(subset_Seurat_Obj@meta.data$clone_id != clone)] <- "ident2"
      subset_Seurat_Obj <- SetIdent(object = subset_Seurat_Obj,
                                    cells = rownames(subset_Seurat_Obj@meta.data),
                                    value = new.ident)
      
      ### DE analysis with DESeq2
      de_result <- FindMarkers(subset_Seurat_Obj,
                               ident.1 = "ident1",
                               ident.2 = "ident2",
                               logfc.threshold = 0,
                               min.pct = 0.1,
                               test.use = "DESeq2")
      
      ### rearange the columns
      de_result <- data.frame(Gene_Symbol=rownames(de_result),
                              de_result[,-which(colnames(de_result) == "p_val_adj")],
                              FDR=p.adjust(de_result$p_val, method = "BH"),
                              stringsAsFactors = FALSE, check.names = FALSE)
      
      ### save in Excel
      write.xlsx2(de_result, file = paste0(outputDir4, donor, "_DESeq2_", clone, "_vs_Cluster17.xlsx"),
                  sheetName = paste0(clone, "_vs_Others"), row.names = FALSE)
      
      ### pathway analysis
      if(length(which(de_result$FDR < de_signif_thresh)) > 0) {
        pathway_result_GO <- pathwayAnalysis_CP(geneList = mapIds(org.Hs.eg.db,
                                                                  de_result[which(de_result$FDR < de_signif_thresh),
                                                                            "Gene_Symbol"],
                                                                  "ENTREZID", "SYMBOL"),
                                                org = "human", database = "GO",
                                                title = paste0(donor, "_Pathway_Results_", clone, "_vs_Others"),
                                                displayNum = 50, imgPrint = TRUE,
                                                dir = paste0(outputDir4))
        pathway_result_KEGG <- pathwayAnalysis_CP(geneList = mapIds(org.Hs.eg.db,
                                                                    de_result[which(de_result$FDR < de_signif_thresh),
                                                                              "Gene_Symbol"],
                                                                    "ENTREZID", "SYMBOL"),
                                                  org = "human", database = "KEGG",
                                                  title = paste0(donor, "_Pathway_Results_", clone, "_vs_Others"),
                                                  displayNum = 50, imgPrint = TRUE,
                                                  dir = paste0(outputDir4))
        if(!is.null(pathway_result_GO) && nrow(pathway_result_GO) > 0) {
          write.xlsx2(pathway_result_GO, file = paste0(outputDir4, donor, "_GO_pathway_results_", clone, "_vs_Others.xlsx"),
                      row.names = FALSE, sheetName = paste0("GO_Results_", clone, "_vs_Others"))
        }
        if(!is.null(pathway_result_KEGG) && nrow(pathway_result_KEGG) > 0) {
          write.xlsx2(pathway_result_KEGG, file = paste0(outputDir4, donor, "_KEGG_pathway_results_", clone, "_vs_Others.xlsx"),
                      row.names = FALSE, sheetName = paste0("KEGG_Results_", clone, "_vs_Others"))
        }
      }
      
    }
    
    
    ### PCA plot
    g <- multiReducPlot(x = subset_Seurat_Obj@reductions$pca@cell.embeddings[,"PC_1"],
                        y = subset_Seurat_Obj@reductions$pca@cell.embeddings[,"PC_2"],
                        group_color = subset_Seurat_Obj@meta.data$Day,
                        group_shape = subset_Seurat_Obj@meta.data$Tissue,
                        type = "PCA",
                        fName = paste0(donor, "_PCA_TFH"), isPrint = TRUE, isConvex = TRUE)
    ggsave(file = paste0(outputDir2, donor, "_PCA_TFH.png"), g, width = 20, height = 10, dpi = 300)
    
    ### UMAP plot
    g <- multiReducPlot(x = subset_Seurat_Obj@reductions$umap@cell.embeddings[,"UMAP_1"],
                        y = subset_Seurat_Obj@reductions$umap@cell.embeddings[,"UMAP_2"],
                        group_color = subset_Seurat_Obj@meta.data$Day,
                        group_shape = subset_Seurat_Obj@meta.data$Tissue,
                        type = "UMAP",
                        fName = paste0(donor, "_UMAP_TFH"), isPrint = TRUE, isConvex = TRUE)
    ggsave(file = paste0(outputDir2, donor, "_UMAP_TFH.png"), g, width = 20, height = 10, dpi = 300)
    
    
    ### find feature contributions of the PC1 
    pca_cos2 <- subset_Seurat_Obj@reductions$pca@feature.loadings * subset_Seurat_Obj@reductions$pca@feature.loadings
    pca_contb <- pca_cos2
    for(i in 1:ncol(pca_contb)) {
      s <- sum(pca_cos2[,i])
      for(j in 1:nrow(pca_contb)) {
        pca_contb[j,i] <- pca_cos2[j,i] * 100 / s
      }
    }
    pca_contb <- pca_contb[order(-pca_contb[,"PC_1"]),]
    write.xlsx2(data.frame(Gene=rownames(pca_contb), pca_contb,
                           stringsAsFactors = FALSE, check.names = FALSE),
                file = paste0(outputDir2, donor, "_PCA_Contributions.xlsx"),
                row.names = FALSE, sheetName = "PCA_Contributions")
    
    ### pathway analysis with the important genes of the PC1
    contb_threshold <- 0.1
    important_genes <- rownames(pca_contb)[which(pca_contb[,"PC_1"] > contb_threshold)]
    pathway_result_GO <- pathwayAnalysis_CP(geneList = mapIds(org.Hs.eg.db,
                                                              important_genes,
                                                              "ENTREZID", "SYMBOL"),
                                            org = "human", database = "GO",
                                            title = paste0(donor, "_Pathway_Results_", length(important_genes), "_PC1_Genes_", contb_threshold),
                                            displayNum = 50, imgPrint = TRUE,
                                            dir = paste0(outputDir2))
    pathway_result_KEGG <- pathwayAnalysis_CP(geneList = mapIds(org.Hs.eg.db,
                                                                important_genes,
                                                                "ENTREZID", "SYMBOL"),
                                              org = "human", database = "KEGG",
                                              title = paste0(donor, "_Pathway_Results_", length(important_genes), "_PC1_Genes_", contb_threshold),
                                              displayNum = 50, imgPrint = TRUE,
                                              dir = paste0(outputDir2))
    write.xlsx2(pathway_result_GO, file = paste0(outputDir2, donor, "_GO_pathway_results_", length(important_genes), "_PC1_Genes_", contb_threshold, ".xlsx"),
                row.names = FALSE, sheetName = paste0("GO_Results"))
    write.xlsx2(pathway_result_KEGG, file = paste0(outputDir2, donor, "_KEGG_pathway_results_", length(important_genes), "_PC1_Genes_", contb_threshold, ".xlsx"),
                row.names = FALSE, sheetName = paste0("KEGG_Results"))
    
    
    ### get genes that contributed to the PC1 the most
    contb_threshold <- 1
    important_genes <- rownames(pca_contb)[which(pca_contb[,"PC_1"] > contb_threshold)]
    if(length(important_genes) < 5) {
      contb_threshold <- 0.5
      important_genes <- rownames(pca_contb)[which(pca_contb[,"PC_1"] > contb_threshold)]
    }
    
    #
    ### a heatmap with the genes
    #
    
    ### get a matrix for the heatmap
    heatmap_mat <- data.frame(subset_Seurat_Obj@assays$RNA@counts[important_genes,], check.names = FALSE)
    
    ### scale the data
    heatmap_mat_scaled <- scale_h(heatmap_mat, type = "row")
    
    ### because there are some outliers in positive values
    ### we set the maximum as abs(minimum)
    heatmap_mat_scaled[which(heatmap_mat_scaled > abs(min(heatmap_mat_scaled)))] <- abs(min(heatmap_mat_scaled))
    
    ### set colside colors
    uniqueV <- levels(subset_Seurat_Obj@meta.data$Day)
    colors <- colorRampPalette(brewer.pal(9,"Blues"))(length(uniqueV))
    names(colors) <- uniqueV
    
    ### hierarchical clustering functions
    dist.spear <- function(x) as.dist(1-cor(t(x), method = "spearman"))
    hclust.ave <- function(x) hclust(x, method="average")
    
    ### heatmap
    png(paste0(outputDir2, donor, "_PC1_Genes_Heatmap.png"), width = 2000, height = 1000)
    par(oma=c(0,0,2,6))
    heatmap.3(as.matrix(heatmap_mat_scaled), main = paste0("PC1_Genes_Heatmap_(",
                                                           nrow(heatmap_mat_scaled), " Genes x ",
                                                           ncol(heatmap_mat_scaled), " Cells)"),
              xlab = "", ylab = "", col=greenred(100),
              scale="none", key=T, keysize=0.8, density.info="density",
              dendrogram = "none", trace = "none",
              labRow = rownames(heatmap_mat_scaled), labCol = FALSE,
              Rowv = TRUE, Colv = FALSE,
              distfun=dist.spear, hclustfun=hclust.ave,
              ColSideColors = cbind(colors[as.character(subset_Seurat_Obj@meta.data$Day)]),
              cexRow = 2, cexCol = 2, na.rm = TRUE)
    legend("left", inset = 0, xpd = TRUE, title = "Cell Collection Time", legend = names(colors), fill = colors, cex = 2, box.lty = 0)
    dev.off()
    
    
    #
    ### Diversity & Clonality
    #
    
    ### diversity
    
    ### data preparation
    seurat_diversity <- rep(0, length(levels(subset_Seurat_Obj@meta.data$Day)))
    seurat_cellNum <- rep(0, length(levels(subset_Seurat_Obj@meta.data$Day)))
    names(seurat_diversity) <- levels(subset_Seurat_Obj@meta.data$Day)
    names(seurat_cellNum) <- levels(subset_Seurat_Obj@meta.data$Day)
    for(tp in names(seurat_diversity)) {
      seurat_diversity[tp] <- length(unique(subset_Seurat_Obj@meta.data$clone_id[which(subset_Seurat_Obj@meta.data$Day == tp)]))
      seurat_cellNum[tp] <- length(which(subset_Seurat_Obj@meta.data$Day == tp))
    }
    
    ### data frame for barplot
    plot_df <- data.frame(Value=c(seurat_diversity, seurat_cellNum),
                          Time=c(names(seurat_diversity), names(seurat_cellNum)),
                          Type=c(rep("Clone #", length(seurat_diversity)), rep("Cell #", length(seurat_cellNum))),
                          Pct=c(paste0(signif((seurat_diversity / seurat_cellNum) * 100, digits = 4), "%"),
                                rep("", length(seurat_diversity))),
                          stringsAsFactors = FALSE, check.names = FALSE)
    
    ### factorize the time column
    plot_df$Time <- factor(plot_df$Time, levels = levels(subset_Seurat_Obj@meta.data$Day))
    
    ### draw the plot
    ggplot(data=plot_df, aes(x=Time, y=Value, fill=Type)) +
      geom_bar(stat="identity", position=position_dodge()) +
      geom_text(aes(label=Value), vjust=-1, color="black",
                position = position_dodge(0.9), size=3.5) +
      geom_text(aes(label=Pct, y=75), vjust=0, color="blue", size=3.5) +
      labs(title = paste(donor, "TFH Cluster Clonal Diversity Over Time"),
           subtitle = "(Clone # / Cell #) x 100") +
      scale_fill_brewer(palette="Paired") +
      ylim(0, max(plot_df$Value) * 1.1) +
      theme_classic(base_size = 16) +
      theme(axis.title.y = element_blank(),
            plot.subtitle=element_text(size=12, color="blue"))
    ggsave(file = paste0(outputDir2, donor, "_TFH_Cluster_Clonal_Diversity.png"), width = 12, height = 8, dpi = 300)
    
    
    ### pie plots for the clonality
    
    ### select top n clones
    top_clone_num <- 12
    
    ### data frame for the pie plot
    plot_df <- data.frame(freq=as.numeric(as.vector(t(clone_summary_table[1:top_clone_num,time_points]))),
                          clone_id=c(sapply(clone_summary_table$clone_id[1:top_clone_num], function(x) rep(x, length(time_points)))),
                          cdr_ab=c(sapply(clone_summary_table$cdr_ab[1:top_clone_num], function(x) rep(x, length(time_points)))),
                          time=c(rep(time_points, top_clone_num)),
                          stringsAsFactors = FALSE, check.names = FALSE)
    
    ### remove zero frequencies
    plot_df <- plot_df[which(plot_df$freq != 0),]
    
    ### create empty pie chart results
    p <- vector("list", length = length(unique(plot_df$clone_id)))
    names(p) <- unique(plot_df$clone_id)
    
    ### draw pie plots for the top clones
    for(clone in unique(plot_df$clone_id)) {
      ### get the clone indicies
      clone_idx <- which(plot_df$clone_id == clone)
      
      ### calculate percentages
      total_sample_num <- sum(plot_df$freq[clone_idx])
      pct <- sapply(plot_df$freq[clone_idx], function(x) signif(x*100/total_sample_num, digits = 3))
      pct <- paste0(plot_df$freq[clone_idx], "(", pct, "%)")
      
      ### ggplot drawing
      p[[clone]] <- ggplot(data = plot_df[clone_idx,],
                           aes(x = "", y = freq, fill = time)) +
        geom_bar(stat = "identity", width = 1) +
        coord_polar(theta="y") +
        geom_text(label = pct,
                  position = position_stack(vjust = 0.5), size = 3) +
        labs(x = NULL, y = NULL, title = clone) +
        theme_classic(base_size = 16) +
        theme(plot.title = element_text(hjust = 0.5, color = "black", size = 12),
              axis.line = element_blank(),
              axis.ticks = element_blank(),
              axis.text = element_blank())
    }
    
    ### arrange the plots and print out
    fName <- paste0(donor, "_TFH_Cluster_Clonality_Pie_Charts")
    g <- arrangeGrob(grobs = p,
                     nrow = 3,
                     ncol = 4,
                     top = fName)
    ggsave(file = paste0(outputDir2, fName, ".png"), g, width = 15, height = 9, dpi = 300)
    
    
    ### pseudotime analysis - Slingshot
    
    ### clustering on the TFH cluster
    subset_Seurat_Obj <- FindNeighbors(subset_Seurat_Obj, dims = 1:5, k.param = 5)
    subset_Seurat_Obj <- FindClusters(subset_Seurat_Obj, resolution = 0.4)
    
    ### run clustering on the PCA
    set.seed(1234)
    pca_map <- Embeddings(subset_Seurat_Obj, reduction = "pca")[rownames(subset_Seurat_Obj@meta.data),1:10]
    subset_Seurat_Obj@meta.data$mclust_clusters <- Mclust(pca_map)$classification[rownames(subset_Seurat_Obj@meta.data)]
    
    ### add "cluster" to the mclust_clusters column
    subset_Seurat_Obj@meta.data$mclust_clusters <- paste0("cluster", subset_Seurat_Obj@meta.data$mclust_clusters)
    
    ### factorize the mclust_clusters column
    ordered_clusters <- unique(subset_Seurat_Obj@meta.data$mclust_clusters)
    ordered_clusters <- ordered_clusters[order(ordered_clusters)]
    subset_Seurat_Obj@meta.data$mclust_clusters <- factor(subset_Seurat_Obj@meta.data$mclust_clusters,
                                                          levels = as.character(ordered_clusters))
    
    ### get slingshot object
    slingshot_obj <- slingshot(pca_map,
                               clusterLabels = subset_Seurat_Obj@meta.data$mclust_clusters, 
                               reducedDim = "PCA")
    
    ### get colors for the clustering result
    cell_colors_clust <- cell_pal(levels(subset_Seurat_Obj@meta.data$mclust_clusters), hue_pal())
    
    ### Trajectory inference
    png(paste0(outputDir2, donor, "_Trajectory_Inference_Mclust_PCA.png"), width = 2500, height = 1500, res = 200)
    plot(pca_map,
         main=paste(donor, "Trajectory Inference Based On Mclust Clusters (PCA)"),
         col = cell_colors_clust[as.character(subset_Seurat_Obj@meta.data$mclust_clusters)],
         pch = 19, cex = 1)
    lines(slingshot_obj, lwd = 2, type = "lineages", col = "black",
          show.constraints = TRUE, constraints.col = cell_colors_clust)
    legend("bottomleft", legend = names(cell_colors_clust), col = cell_colors_clust,
           pch = 19)
    dev.off()
    
    ### Trajectory inference on multi dimentional PCA
    png(paste0(outputDir2, donor, "_Trajectory_Inference_Mclust_Multi-PCA.png"), width = 2500, height = 1500, res = 200)
    pairs(slingshot_obj, type="lineages", col = apply(slingshot_obj@clusterLabels, 1, function(x) cell_colors_clust[names(x)[which(x == 1)]]),
          show.constraints = TRUE, constraints.col = cell_colors_clust, cex = 0.8,
          horInd = 1:5, verInd = 1:5, main = paste0(donor, "_Trajectory_Inference_Mclust"))
    par(xpd = TRUE)
    legend("bottomleft", legend = names(cell_colors_clust), col = cell_colors_clust,
           pch = 19, title = "Time")
    dev.off()
    
    ### 3D Slingshot
    slingshot_3d_lineages(slingshot_obj = slingshot_obj,
                          color = cell_colors_clust,
                          title = paste0(donor, "_Trajectory_Inference_3D_PCA_Mclust"),
                          print = TRUE,
                          outputDir = paste0(outputDir2, "Slingshot_3D/"),
                          width = 1200,
                          height = 800)
    rgl.close()
    
    
    #
    ### do the same thing with the time points
    #
    
    ### get slingshot object
    slingshot_obj <- slingshot(pca_map,
                               clusterLabels = subset_Seurat_Obj@meta.data$Day, 
                               reducedDim = "PCA")
    
    ### get colors for the clustering result
    cell_colors_clust <- cell_pal(levels(subset_Seurat_Obj@meta.data$Day), hue_pal())
    
    ### Trajectory inference
    png(paste0(outputDir2, donor, "_Trajectory_Inference_Time_PCA.png"), width = 2500, height = 1500, res = 200)
    plot(reducedDim(slingshot_obj),
         main=paste(donor, "Trajectory Inference Based On Time (PCA)"),
         col = cell_colors_clust[subset_Seurat_Obj@meta.data$Day],
         pch = 19, cex = 1)
    lines(slingshot_obj, lwd = 2, type = "lineages", col = "black",
          show.constraints = TRUE, constraints.col = cell_colors_clust)
    legend("bottomleft", legend = names(cell_colors_clust), col = cell_colors_clust,
           pch = 19)
    dev.off()
    
    ### Trajectory inference on multi dimentional PCA
    png(paste0(outputDir2, donor, "_Trajectory_Inference_Time_Multi-PCA.png"), width = 2500, height = 1500, res = 200)
    pairs(slingshot_obj, type="lineages", col = apply(slingshot_obj@clusterLabels, 1, function(x) cell_colors_clust[names(x)[which(x == 1)]]),
          show.constraints = TRUE, constraints.col = cell_colors_clust, cex = 0.8,
          horInd = 1:5, verInd = 1:5, main = paste0(donor, "_Trajectory_Inference_Time"))
    par(xpd = TRUE)
    legend("bottomleft", legend = names(cell_colors_clust), col = cell_colors_clust,
           pch = 19, title = "Time")
    dev.off()
    
    ### 3D Slingshot
    slingshot_3d_lineages(slingshot_obj = slingshot_obj,
                          color = cell_colors_clust,
                          title = paste0(donor, "_Trajectory_Inference_3D_PCA_Time"),
                          print = TRUE,
                          outputDir = paste0(outputDir2, "Slingshot_3D/"),
                          width = 1200,
                          height = 800)
    rgl.close()
    
  }
  
  ### get shared genes between the two donors that contributed to the PC1 a lot in each
  gene_list1 <- read.xlsx2(file = "./results/v2/321-04/321-04_PCA_Contributions.xlsx",
                           sheetIndex = 1)
  gene_list2 <- read.xlsx2(file = "./results/v2/321-05/321-05_PCA_Contributions.xlsx",
                           sheetIndex = 1)
  shared_genes <- intersect(gene_list1$Gene[1:100], gene_list2$Gene[1:100])
  write.table(data.frame(shared_genes), file = "./results/v2/62_Shared_PC1_Genes.txt",
              sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
  
  ### pathway analysis on the shared genes
  pathway_result_GO <- pathwayAnalysis_CP(geneList = mapIds(org.Hs.eg.db,
                                                            shared_genes,
                                                            "ENTREZID", "SYMBOL"),
                                          org = "human", database = "GO",
                                          title = "Pathway_Results_62_Shared_PC1_Genes",
                                          displayNum = 50, imgPrint = TRUE,
                                          dir = "./results/v2/")
  pathway_result_KEGG <- pathwayAnalysis_CP(geneList = mapIds(org.Hs.eg.db,
                                                              shared_genes,
                                                              "ENTREZID", "SYMBOL"),
                                            org = "human", database = "KEGG",
                                            title = "Pathway_Results_62_Shared_PC1_Genes",
                                            displayNum = 50, imgPrint = TRUE,
                                            dir = "./results/v2/")
  write.xlsx2(pathway_result_GO, file = "./results/v2/GO_Pathway_Results_62_Shared_PC1_Genes.xlsx",
              row.names = FALSE, sheetName = paste0("GO_Results"))
  write.xlsx2(pathway_result_KEGG, file = "./results/v2/KEGG_Pathway_Results_62_Shared_PC1_Genes.xlsx",
              row.names = FALSE, sheetName = paste0("KEGG_Results"))
  
  
  ### load the gene set from the previous dataset
  ### that contributed the PC1 the most
  ### then see the gene expression pattern in the new dataset
  
  ### load the gene set
  gene_set <- read.xlsx2(file = "./results/Cluster17-TFH/Interesting_Clones/PCA_Contributions.xlsx",
                         sheetIndex = 1, stringsAsFactors = FALSE, check.names = FALSE)
  gene_set[,-1] <- sapply(gene_set[,-1], as.numeric)
  rownames(gene_set) <- gene_set$Gene
  gene_set <- gene_set[,-1]
  
  ### get important genes from the previous result
  old_important_genes <- rownames(gene_set)[which(gene_set[,"PC_1"] > 1)]
  
  ### hierarchical clustering functions
  dist.spear <- function(x) as.dist(1-cor(t(x), method = "spearman"))
  hclust.ave <- function(x) hclust(x, method="average")
  
  ### for each donor run the same analyses
  for(donor in unique(Seurat_Obj@meta.data$Donor)) {
    
    ### check whether the orders are the same
    print(identical(names(Idents(object = Seurat_Obj)), rownames(Seurat_Obj@meta.data)))
    
    ### split the Seurat obj based on donor info
    Seurat_Obj <- SetIdent(object = Seurat_Obj,
                           cells = rownames(Seurat_Obj@meta.data),
                           value = Seurat_Obj@meta.data$Donor)
    subset_Seurat_Obj <- subset(Seurat_Obj, idents=donor)
    
    ### order the meta data by time
    subset_Seurat_Obj@meta.data <- subset_Seurat_Obj@meta.data[order(subset_Seurat_Obj@meta.data$Day),]
    
    ### rownames in the meta.data should be in the same order as colnames in the counts
    subset_Seurat_Obj@assays$RNA@counts <- subset_Seurat_Obj@assays$RNA@counts[,rownames(subset_Seurat_Obj@meta.data)]
    
    ### run PCA
    subset_Seurat_Obj <- RunPCA(subset_Seurat_Obj, npcs = 10)
    
    ### find feature contributions of the PC1 
    pca_cos2 <- subset_Seurat_Obj@reductions$pca@feature.loadings * subset_Seurat_Obj@reductions$pca@feature.loadings
    pca_contb <- pca_cos2
    for(i in 1:ncol(pca_contb)) {
      s <- sum(pca_cos2[,i])
      for(j in 1:nrow(pca_contb)) {
        pca_contb[j,i] <- pca_cos2[j,i] * 100 / s
      }
    }
    pca_contb <- pca_contb[order(-pca_contb[,"PC_1"]),]
    
    ### prepare the gene lists that will be compared
    old_gene_list <- rownames(gene_set)
    new_gene_list <- rownames(pca_contb)
    shared_genes <- intersect(old_gene_list, new_gene_list)
    old_gene_list <- old_gene_list[which(old_gene_list %in% shared_genes)]
    new_gene_list <- new_gene_list[which(new_gene_list %in% shared_genes)]
    
    ### prepare the signature for the comparison
    old_gene_signature <- gene_set[old_gene_list, "PC_1"]
    new_gene_signature <- pca_contb[new_gene_list, "PC_1"]
    names(old_gene_signature) <- old_gene_list
    names(new_gene_signature) <- new_gene_list
    
    ### compare the ranks
    compare_two_different_ranks(A = old_gene_signature,
                                B = new_gene_signature,
                                A_name = "Old_PC1_Genes_321-05",
                                B_name = paste0("New_PC1_Genes_", donor),
                                ordering = "decreasing",
                                top = 100,
                                alternative = "greater",
                                permutation = 10000,
                                fileName = paste0("Rank_Comparison_PC1_Genes_(Old_321-05)_vs_(New_", donor, ")"),
                                printPath = "./results/v2/",
                                width = 24,
                                height = 12)
    
    #
    ### a heatmap with the genes
    #
    
    ### get a matrix for the heatmap
    heatmap_mat <- data.frame(subset_Seurat_Obj@assays$RNA@counts[old_important_genes,], check.names = FALSE)
    
    ### scale the data
    heatmap_mat_scaled <- scale_h(heatmap_mat, type = "row")
    
    ### because there are some outliers in positive values
    ### we set the maximum as abs(minimum)
    heatmap_mat_scaled[which(heatmap_mat_scaled > abs(min(heatmap_mat_scaled)))] <- abs(min(heatmap_mat_scaled))
    
    ### set colside colors
    uniqueV <- levels(subset_Seurat_Obj@meta.data$Day)
    colors <- colorRampPalette(brewer.pal(9,"Blues"))(length(uniqueV))
    names(colors) <- uniqueV
    
    ### heatmap
    png(paste0("./results/v2/Old_Genes_in_New_", donor, "_Heatmap.png"), width = 2000, height = 1000)
    par(oma=c(0,0,2,6))
    heatmap.3(as.matrix(heatmap_mat_scaled), main = paste0("Old_PC1_Genes_Heatmap_(",
                                                           nrow(heatmap_mat_scaled), " Genes x ",
                                                           ncol(heatmap_mat_scaled), " Cells)"),
              xlab = "", ylab = "", col=greenred(100),
              scale="none", key=T, keysize=0.8, density.info="density",
              dendrogram = "none", trace = "none",
              labRow = rownames(heatmap_mat_scaled), labCol = FALSE,
              Rowv = TRUE, Colv = FALSE,
              distfun=dist.spear, hclustfun=hclust.ave,
              ColSideColors = cbind(colors[as.character(subset_Seurat_Obj@meta.data$Day)]),
              cexRow = 2, cexCol = 2, na.rm = TRUE)
    legend("left", inset = 0, xpd = TRUE, title = "Cell Collection Time", legend = names(colors), fill = colors, cex = 2, box.lty = 0)
    dev.off()
    
  }
  
  
  
}
