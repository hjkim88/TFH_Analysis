###
#   File name : TFH_Additional_Analyses.R
#   Author    : Hyunjin Kim
#   Date      : Jun 8, 2020
#   Email     : hyunjin.kim@stjude.org
#   Purpose   : Here are additional things to do:
#               1. PCA & UMAPs labeld with time points
#               2. Evidence of naive and recall, and their differences?
#                  * Naive as being size 1 and also located within the naive clusters.
#                    Resting would be clone size >= 1 and outside the naive clusters.
#                    Recalled would be size > 1 and located in a cluster with other highly activated cells.
#
#   Instruction
#               1. Source("TFH_Additional_Analyses.R")
#               2. Run the function "tfh_additional_analyses" - specify the input file path and the output directory
#               3. The results will be generated under the output directory
#
#   Example
#               > source("The_directory_of_TFH_Additional_Analyses.R/TFH_Additional_Analyses.R")
#               > tfh_additional_analyses(Seurat_RObj_path="./data/Ali_Tcell_combined.RDATA",
#                                         outputDir="./results/")
###

tfh_additional_analyses <- function(Seurat_RObj_path="./data/Ali_Tcell_combined.RDATA",
                                    outputDir="./results/") {
  
  ### load libraries
  if(!require(Seurat, quietly = TRUE)) {
    install.packages("Seurat")
    require(Seurat, quietly = TRUE)
  }
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
  
  ### load the Seurat object and save the object name
  tmp_env <- new.env()
  load(Seurat_RObj_path, tmp_env)
  obj_name <- ls(tmp_env)
  assign("Seurat_Obj", get(obj_name, envir = tmp_env))
  rm(tmp_env)
  gc()
  
  ### rownames in the meta.data should be in the same order as colnames in the counts
  Seurat_Obj@meta.data <- Seurat_Obj@meta.data[colnames(Seurat_Obj@assays$RNA@counts),]
  
  ### get unique clone ids in the cluster 17
  cluster_17_clone_ids <- unique(Seurat_Obj@meta.data$clone_id[which(Seurat_Obj@meta.data$seurat_clusters == "17")])
  cluster_17_clone_ids <- cluster_17_clone_ids[intersect(intersect(which(!is.na(cluster_17_clone_ids)),
                                                                   which(cluster_17_clone_ids != "NA")),
                                                         which(cluster_17_clone_ids != ""))]
  
  ### get meta.data of the cells that have the unique clone ids from the cluster 17
  cluster_17_clones_meta.data <- Seurat_Obj@meta.data[which(Seurat_Obj@meta.data$clone_id %in% cluster_17_clone_ids),]
  
  ### get a subset for the cluster 17
  new.ident <- rep(NA, nrow(Seurat_Obj@meta.data))
  new.ident[which(Seurat_Obj@meta.data$clone_id %in% cluster_17_clone_ids)] <- "Cluster17"
  Idents(object = Seurat_Obj) <- new.ident
  subset_Seurat_Obj <- subset(Seurat_Obj, idents=c("Cluster17"))
  
  ### factorize the Day column
  subset_Seurat_Obj@meta.data$Day <- factor(subset_Seurat_Obj@meta.data$Day,
                                            levels = c("d0", "d5", "d12", "d28", "d60", "d90", "d120", "d180"))
  
  ### a function that returns multi-figures from PCA, TSNE, AND UMAP
  ### One plot with all the groups + each group
  ### x: a numeric vector of first component
  ### y: a numeric vector of second component
  ### group: a factor vector of group info
  ### type: type of the dimensionality reduction method
  multiReducPlot <- function(x, y, group, type=c("PCA", "TSNE", "UMAP"), fName="Whole Group & Each", isConvex=FALSE, isPrint=FALSE) {
    
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
    group <- factor(as.character(group),
                    levels = intersect(levels(group), unique(group)))
    plot_df <- data.frame(X=x, Y=y, Group=group)
    
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
    col_palette <- hue_pal()(length(levels(group)))
    names(col_palette) <- levels(group)
    
    ### 1. One plot with all the groups + each group
    p <- vector("list", length(levels(group))+1)
    x_range <- c(min(plot_df$X), max(plot_df$X))
    y_range <- c(min(plot_df$Y), max(plot_df$Y))
    if(isConvex) {
      p[[1]] <- ggplot(plot_df, aes_string(x="X", y="Y")) +
        geom_point(aes_string(col="Group"), size=2, alpha=0.6) +
        xlab(x_label) + ylab(y_label) +
        xlim(x_range[1], x_range[2]) + ylim(y_range[1], y_range[2]) +
        ggtitle("All") +
        theme_classic(base_size = 16) +
        theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) +
        geom_convexhull(aes_string(col="Group", fill="Group"), size=1, alpha=0)
    } else {
      p[[1]] <- ggplot(plot_df, aes_string(x="X", y="Y")) +
        geom_point(aes_string(col="Group"), size=2, alpha=0.6) +
        xlab(x_label) + ylab(y_label) +
        xlim(x_range[1], x_range[2]) + ylim(y_range[1], y_range[2]) +
        ggtitle("All") +
        theme_classic(base_size = 16) +
        theme(legend.position = "none", plot.title = element_text(hjust = 0.5))
    }
    for(i in 2:(length(levels(group))+1)) {
      p[[i]] <- ggplot(plot_df[which(plot_df$Group == levels(group)[i-1]),], aes_string(x="X", y="Y")) +
        geom_point(col=col_palette[levels(group)[i-1]], size=2, alpha=0.8) +
        xlab(x_label) + ylab(y_label) +
        xlim(x_range[1], x_range[2]) + ylim(y_range[1], y_range[2]) +
        ggtitle(levels(group)[i-1]) +
        theme_classic(base_size = 16) +
        theme(plot.title = element_text(hjust = 0.5, 
                                        color = col_palette[levels(group)[i-1]]))
    }
    ### arrange the plots
    g <- arrangeGrob(grobs = p,
                     nrow = ceiling(sqrt(length(levels(group))+1)),
                     ncol = ceiling(sqrt(length(levels(group))+1)),
                     top = fName)
    
    if(isPrint) {
      plot(g)
    }
    
    return(g)
    
  }
  
  ### run PCA
  subset_Seurat_Obj <- RunPCA(subset_Seurat_Obj, npcs = 10)
  
  ### PCA plot
  g <- multiReducPlot(x = subset_Seurat_Obj@reductions$pca@cell.embeddings[,"PC_1"],
                      y = subset_Seurat_Obj@reductions$pca@cell.embeddings[,"PC_2"],
                      group = subset_Seurat_Obj@meta.data$Day, type = "PCA",
                      fName = "PCA_TFH_Cluster_17", isPrint = TRUE, isConvex = TRUE)
  ggsave(file = paste0(outputDir, "PCA_TFH_Cluster_17.png"), g, width = 20, height = 10, dpi = 300)
  
  # DimPlot(subset_Seurat_Obj, reduction = "pca", group.by = "Day", pt.size = 2) +
  #   labs(title = "PCA TFH Cluster17")
  # ggsave(file = paste0(outputDir, "PCA_TFH_Cluster_17_All.png"), width = 20, height = 10, dpi = 300)
  # DimPlot(subset_Seurat_Obj, reduction = "pca", split.by = "Day", group.by = "Day",
  #         pt.size = 2, ncol = 4) +
  #   labs(title = "PCA TFH Cluster17 by Time")
  
  ### run UMAP
  subset_Seurat_Obj <- RunUMAP(subset_Seurat_Obj, dims = 1:5)
  
  ### UMAP plot
  g <- multiReducPlot(x = subset_Seurat_Obj@reductions$umap@cell.embeddings[,"UMAP_1"],
                      y = subset_Seurat_Obj@reductions$umap@cell.embeddings[,"UMAP_2"],
                      group = subset_Seurat_Obj@meta.data$Day, type = "UMAP",
                      fName = "UMAP_TFH_Cluster_17", isPrint = TRUE, isConvex = TRUE)
  ggsave(file = paste0(outputDir, "UMAP_TFH_Cluster_17.png"), g, width = 20, height = 10, dpi = 300)
  
  # DimPlot(subset_Seurat_Obj, reduction = "umap", group.by = "Day", pt.size = 2) +
  #   labs(title = "UMAP TFH Cluster17")
  # DimPlot(subset_Seurat_Obj, reduction = "umap", split.by = "Day", group.by = "Day",
  #         pt.size = 2, ncol = 4) +
  #   labs(title = "UMAP TFH Cluster17 by Time")
  
  
  #
  ### In PCA, there are no differences among the time points,
  ### but in UMAP, most of the plots have 2 clusters
  ### and I want to know what makes the separation
  #
  
  ### test if it is explained by existing info
  p <- list()
  p[[1]] <- DimPlot(subset_Seurat_Obj, reduction = "umap", group.by = "Library", pt.size = 2) +
    labs(title = "Library")
  p[[2]] <- DimPlot(subset_Seurat_Obj, reduction = "umap", group.by = "Type", pt.size = 2) +
    labs(title = "Type")
  p[[3]] <- DimPlot(subset_Seurat_Obj, reduction = "umap", group.by = "Chip", pt.size = 2) +
    labs(title = "Chip")
  p[[4]] <- DimPlot(subset_Seurat_Obj, reduction = "umap", group.by = "LibNum", pt.size = 2) +
    labs(title = "LibNum")
  p[[5]] <- DimPlot(subset_Seurat_Obj, reduction = "umap", group.by = "Phase", pt.size = 2) +
    labs(title = "Phase")
  p[[6]] <- DimPlot(subset_Seurat_Obj, reduction = "umap", group.by = "seurat_clusters", pt.size = 2) +
    labs(title = "seurat_clusters")
  p[[7]] <- DimPlot(subset_Seurat_Obj, reduction = "umap", group.by = "ident", pt.size = 2) +
    labs(title = "ident")
  p[[8]] <- DimPlot(subset_Seurat_Obj, reduction = "umap", group.by = "clone_call", pt.size = 2) +
    labs(title = "clone_call")
  p[[9]] <- DimPlot(subset_Seurat_Obj, reduction = "umap", group.by = "Day", pt.size = 2) +
    labs(title = "Day")
  g <- arrangeGrob(grobs = p,
                   nrow = 3,
                   ncol = 3,
                   top = "UMAP coloring with every possible column")
  ggsave(file = paste0(outputDir, "UMAP_TFH_Cluster_17_Various.png"), g, width = 20, height = 12, dpi = 300)
  
  ### see differentially expressed genes between two groups
  group1 <- which(subset_Seurat_Obj@reductions$umap@cell.embeddings[,"UMAP_1"] < 0)
  group2 <- which(subset_Seurat_Obj@reductions$umap@cell.embeddings[,"UMAP_1"] >= 0)
  new.ident <- rep(NA, nrow(subset_Seurat_Obj@meta.data))
  new.ident[group1] <- "ident1"
  new.ident[group2] <- "ident2"
  Idents(object = subset_Seurat_Obj) <- new.ident
  de_result <- FindMarkers(subset_Seurat_Obj,
                           ident.1 = "ident1",
                           ident.2 = "ident2",
                           logfc.threshold = 0,
                           min.pct = 0.1,
                           test.use = "DESeq2")
  de_result <- data.frame(Gene_Symbol=rownames(de_result),
                          de_result[,-which(colnames(de_result) == "p_val_adj")],
                          FDR=p.adjust(de_result$p_val, method = "BH"),
                          stringsAsFactors = FALSE, check.names = FALSE)
  
  ### write out the DE result
  write.xlsx2(de_result, file = paste0(outputDir, "DESeq2_Two_Interesting_Clusters.xlsx"),
              sheetName = "DESeq2", row.names = FALSE)
  
  ### DE Gene Expressions in UMAP
  FeaturePlot(object = subset_Seurat_Obj,
              features = c("COTL1", "H3F3B", "TMSB10", "IL32", "RPS12", "TPT1", "HLA-A", "PPP1CC", "PTPRCAP"),
              cols = c("grey", "red"), reduction = "umap", pt.size = 2) +
    labs(plot.title = "DE Gene Expressions in UMAP")
  ggsave(file = paste0(outputDir, "Two_Cluster_DE_Gene_Expressions_in_UMAP.png"), width = 20, height = 12, dpi = 300)
  
  
  #
  ### compare naive and recall
  #
  
  ### Adding cell types to the meta.data
  Seurat_Obj@meta.data$CD4_CD8 <- NA
  Seurat_Obj@meta.data$CD4_CD8[which(Seurat_Obj@meta.data$seurat_clusters %in% c(0,1,3,4,5,8,9,10,13,14,15,17))] <- "CD4"
  Seurat_Obj@meta.data$CD4_CD8[which(Seurat_Obj@meta.data$seurat_clusters %in% c(2,6,7,11,12,16,18))] <- "CD8"
  DimPlot(Seurat_Obj, reduction = "umap", group.by = "CD4_CD8", pt.size = 1, label = TRUE)
  Seurat_Obj@meta.data$Cell_Type <- NA
  Seurat_Obj@meta.data$Cell_Type[which(Seurat_Obj@meta.data$seurat_clusters %in% c(0,1,2,3,9,10,11,14,16))] <- "Naive"
  Seurat_Obj@meta.data$Cell_Type[which(Seurat_Obj@meta.data$seurat_clusters %in% c(4,5,7,12,13,15))] <- "Eff-Mem"
  Seurat_Obj@meta.data$Cell_Type[which(Seurat_Obj@meta.data$seurat_clusters %in% c(6))] <- "MAIT-NKT"
  Seurat_Obj@meta.data$Cell_Type[which(Seurat_Obj@meta.data$seurat_clusters %in% c(18))] <- "Hobits"
  Seurat_Obj@meta.data$Cell_Type[which(Seurat_Obj@meta.data$seurat_clusters %in% c(8))] <- "Treg"
  Seurat_Obj@meta.data$Cell_Type[which(Seurat_Obj@meta.data$seurat_clusters %in% c(17))] <- "TFH"
  DimPlot(Seurat_Obj, reduction = "umap", group.by = "Cell_Type", pt.size = 1, label = TRUE)
  
  length(intersect(which(Seurat_Obj@meta.data$Cell_Type == "Naive"),
                   which(Seurat_Obj@meta.data$Day == "d0")))
  
  ### Are there clones that have the naive phenotype at d0 and expanded at later?
  ### 1. True naive at d0 vs expanded clones
  ### 2. True naive at d0 that will be expanded vs will not be expanded
  
  
  ### "Recall" vs "Resting"
  ### Expanded at Eff-Mem vs the other Eff-Mem 
  
  
  
}
