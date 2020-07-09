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
#               3. Classifier to classify "recall" and "resting"
#   
#   Instruction
#               1. Source("TFH_Additional_Analyses.R")
#               2. Run the function "tfh_additional_analyses" - specify the input file path and the output directory
#               3. The results will be generated under the output directory
#
#   Example
#               > source("The_directory_of_TFH_Additional_Analyses.R/TFH_Additional_Analyses.R")
#               > tfh_additional_analyses(Seurat_RObj_path="./data/Ali_Tcell_combined_NEW.RDATA",
#                                         outputDir="./results/")
###

tfh_additional_analyses <- function(Seurat_RObj_path="./data/Ali_Tcell_combined_NEW.RDATA",
                                    outputDir="./results/") {
  
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
  if(!require(caret, quietly = TRUE)) {
    install.packages("caret")
    require(caret, quietly = TRUE)
  }
  if(!require(pROC, quietly = TRUE)) {
    install.packages("pROC")
    require(pROC, quietly = TRUE)
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
  
  # DimPlot(Seurat_Obj, reduction = "umap", group.by = "seurat_clusters", pt.size = 1, label = TRUE)
  
  ### Adding more informative columns to the meta.data
  Seurat_Obj@meta.data$CD4_CD8 <- NA
  Seurat_Obj@meta.data$CD4_CD8[which(Seurat_Obj@meta.data$seurat_clusters %in% c(0,1,3,4,7,8,9,11,13,14,16,17))] <- "CD4"
  Seurat_Obj@meta.data$CD4_CD8[which(Seurat_Obj@meta.data$seurat_clusters %in% c(2,5,6,10,12,15,18))] <- "CD8"
  DimPlot(Seurat_Obj, reduction = "umap", group.by = "CD4_CD8", pt.size = 1, label = TRUE)
  Seurat_Obj@meta.data$Cell_Type <- NA
  Seurat_Obj@meta.data$Cell_Type[which(Seurat_Obj@meta.data$seurat_clusters %in% c(0,1,2,4,8,10,11,15,16))] <- "Naive"
  Seurat_Obj@meta.data$Cell_Type[which(Seurat_Obj@meta.data$seurat_clusters %in% c(3,6,7,12,13,14))] <- "Eff-Mem"
  Seurat_Obj@meta.data$Cell_Type[which(Seurat_Obj@meta.data$seurat_clusters %in% c(5))] <- "MAIT-NKT"
  Seurat_Obj@meta.data$Cell_Type[which(Seurat_Obj@meta.data$seurat_clusters %in% c(18))] <- "Hobits"
  Seurat_Obj@meta.data$Cell_Type[which(Seurat_Obj@meta.data$seurat_clusters %in% c(9))] <- "Treg"
  Seurat_Obj@meta.data$Cell_Type[which(Seurat_Obj@meta.data$seurat_clusters %in% c(17))] <- "TFH"
  DimPlot(Seurat_Obj, reduction = "umap", group.by = "Cell_Type", pt.size = 1, label = TRUE)
  Seurat_Obj@meta.data$Clone_Size_At_The_Time <- NA
  for(i in 1:nrow(Seurat_Obj@meta.data)) {
    Seurat_Obj@meta.data$Clone_Size_At_The_Time[i] <- length(intersect(which(Seurat_Obj@meta.data$clone_id == Seurat_Obj@meta.data$clone_id[i]),
                                                                       which(Seurat_Obj@meta.data$Day == Seurat_Obj@meta.data$Day[i])))
  }
  
  ### Are there clones that have the naive phenotype at d0 and expanded at later?
  ### 1. True naive at d0 vs expanded clones
  ### 2. True naive at d0 that will be expanded vs will not be expanded
  
  ### new output directory for the results
  outputDir2 <- paste0(outputDir, "Naive_at_d0_also_appear_later/")
  dir.create(outputDir2, showWarnings = FALSE, recursive = TRUE)
  
  ### get indicies of the true naive
  true_naive_idx <- intersect(intersect(which(Seurat_Obj@meta.data$Cell_Type == "Naive"),
                                        which(Seurat_Obj@meta.data$Day == "d0")),
                              which(Seurat_Obj@meta.data$Clone_Size_At_The_Time == 1))
  
  ### unique clones of the true naive
  true_naive_unique_clones <- unique(Seurat_Obj@meta.data$clone_id[true_naive_idx])
  
  ### load the clone summary table
  clone_summary_table <- read.xlsx2(file = paste0(outputDir, "All_Clones_Count_Summary.xlsx"), sheetIndex = 1,
                                    stringsAsFactors=FALSE, check.names=FALSE)
  clone_summary_table[3:11] <- sapply(clone_summary_table[3:11], as.numeric)
  rownames(clone_summary_table) <- clone_summary_table$clone_id
  
  ### clone summary table subset with the unique clones of the true naive
  true_naive_at_d0_clone_summary_table <- clone_summary_table[true_naive_unique_clones,]
  
  ### clones that also appear in the later time points
  result_table <- true_naive_at_d0_clone_summary_table[which(true_naive_at_d0_clone_summary_table$total_count > 1),]
  
  ### save the result
  ### * there are only 2 samples in Group1, so it's impossible to run DE analysis
  write.xlsx2(result_table, file = paste0(outputDir2, "Clones_Naive_at_d0_and_appear_later.xlsx"),
              sheetName = "Naive at d0 and also appear later", row.names = FALSE)
  for(clone in result_table$clone_id) {
    write.xlsx2(Seurat_Obj@meta.data[which(Seurat_Obj@meta.data$clone_id == clone),],
                file = paste0(outputDir2, "Clones_Naive_at_d0_and_appear_later.xlsx"),
                sheetName = clone, row.names = FALSE, append = TRUE)
  }
  
  ### UMAPs with Clusters and with Cells Type
  p <- list()
  p[[1]] <- DimPlot(Seurat_Obj, reduction = "umap", group.by = "seurat_clusters", pt.size = 1, label = TRUE) +
              labs(title = "UMAP with Clusters")
  p[[2]] <- DimPlot(Seurat_Obj, reduction = "umap", group.by = "Cell_Type", pt.size = 1, label = TRUE) +
              labs(title = "UMAP with Cell Types")
  ggsave(file = paste0(outputDir2, "UMAP_Clones_Naive_at_d0_and_appear_later.png"),
         arrangeGrob(p[[1]], p[[2]], nrow = 2, ncol = 1),
         width = 16, height = 12, dpi = 300)
  
  ### "Recall" vs "Resting"
  ### Expanded at Eff-Mem vs the other Eff-Mem 
  
  ### new output directory for the results
  outputDir2 <- paste0(outputDir, "Recall_vs_Resting/")
  dir.create(outputDir2, showWarnings = FALSE, recursive = TRUE)
  
  ### expansion threshold
  ### if clone size > expansion threshold, the clone is regarded as "clonaly expanded"
  expansion_threshold <- 2
  
  ### UMAP of "recall" vs "resting" in each time point
  Seurat_Obj@meta.data$Recall_Resting <- NA
  Seurat_Obj@meta.data$Recall_Resting[intersect(which(Seurat_Obj@meta.data$Cell_Type != "Naive"),
                                                which(Seurat_Obj@meta.data$Clone_Size_At_The_Time > expansion_threshold))] <- "Recall"
  Seurat_Obj@meta.data$Recall_Resting[intersect(which(Seurat_Obj@meta.data$Cell_Type != "Naive"),
                                                which(Seurat_Obj@meta.data$Clone_Size_At_The_Time <= expansion_threshold))] <- "Resting"
  p <- list()
  p[[1]] <- DimPlot(Seurat_Obj, reduction = "umap", group.by = "Recall_Resting", pt.size = 2, label = TRUE) +
              labs(title = paste0("UMAP (Recall - Clone Size > ", expansion_threshold, " or Resting)"))
  p[[2]] <- DimPlot(Seurat_Obj, reduction = "umap", group.by = "Cell_Type", pt.size = 2, label = TRUE) +
              labs(title = "UMAP with Cell Types")
  alpha_v <- rep(0.7, nrow(Seurat_Obj@meta.data))
  alpha_v[intersect(which(Seurat_Obj@meta.data$Cell_Type != "Naive"),
                 which(Seurat_Obj@meta.data$Clone_Size_At_The_Time > expansion_threshold))] <- 1
  p[[1]]$layers[[1]]$aes_params$alpha <- alpha_v
  size_v <- rep(1, nrow(Seurat_Obj@meta.data))
  size_v[intersect(which(Seurat_Obj@meta.data$Cell_Type != "Naive"),
                   which(Seurat_Obj@meta.data$Clone_Size_At_The_Time > expansion_threshold))] <- 3
  p[[1]]$layers[[1]]$aes_params$size <- size_v
  ggsave(file = paste0(outputDir2, "UMAP_Recall_vs_Resting.png"),
         arrangeGrob(p[[1]], p[[2]], nrow = 2, ncol = 1),
         width = 16, height = 12, dpi = 300)
  
  
  ### DE analysis of the "recall" and the "resting"
  
  ### get indicies of the "recall" and the "resting"
  recall_idx <- intersect(which(Seurat_Obj@meta.data$Clone_Size_At_The_Time > expansion_threshold),
                          which(Seurat_Obj@meta.data$Cell_Type != "Naive"))
  resting_idx <- intersect(which(Seurat_Obj@meta.data$Clone_Size_At_The_Time <= expansion_threshold),
                           which(Seurat_Obj@meta.data$Cell_Type != "Naive"))
  
  ### Ident configure
  new.ident <- rep(NA, nrow(Seurat_Obj@meta.data))
  new.ident[recall_idx] <- "ident1"
  new.ident[resting_idx] <- "ident2"
  Idents(object = Seurat_Obj) <- new.ident
  
  ### DE analysis
  de_result <- FindMarkers(Seurat_Obj,
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
  write.xlsx2(de_result, file = paste0(outputDir2, "Recall_vs_Resting_DE_Genes_DESeq2.xlsx"),
              sheetName = "DESeq2", row.names = FALSE)
  
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
  
  ### pathway analysis
  target_genes <- de_result$Gene_Symbol[which(de_result$FDR < 0.0001)]
  pathway_result_GO <- pathwayAnalysis_CP(geneList = mapIds(org.Hs.eg.db, target_genes, "ENTREZID", "SYMBOL"),
                                          org = "human", database = "GO",
                                          title = paste0("Pathway_Results_Recall_vs_Resting_DESeq2"),
                                          displayNum = 50, imgPrint = TRUE,
                                          dir = paste0(outputDir2))
  write.xlsx2(pathway_result_GO, file = paste0(outputDir2, "GO_pathway_results_recall_vs_resting_DESeq2.xlsx"),
              row.names = FALSE, sheetName = "GO")
  pathway_result_KEGG <- pathwayAnalysis_CP(geneList = mapIds(org.Hs.eg.db, target_genes, "ENTREZID", "SYMBOL"),
                                            org = "human", database = "KEGG",
                                            title = paste0("Pathway_Results_Recall_vs_Resting_DESeq2"),
                                            displayNum = 50, imgPrint = TRUE,
                                            dir = paste0(outputDir2))
  write.xlsx2(pathway_result_KEGG, file = paste0(outputDir2, "KEGG_pathway_results_recall_vs_resting_DESeq2.xlsx"),
              row.names = FALSE, sheetName = "KEGG")
  
  
  #
  ### Classifier to classify "recall" and "resting"
  #
  
  ### new output directory for the results
  outputDir2 <- paste0(outputDir2, "Classifier/")
  dir.create(outputDir2, showWarnings = FALSE, recursive = TRUE)
  
  ### set parameters
  set.seed(1234)
  featureSelectionNum <- 100
  random_sampleNum <- 500
  methodTypes <- c("svmLinear", "svmRadial", "gbm", "rf", "LogitBoost", "knn")
  methodNames <- c("SVMLinear", "SVMRadial", "GBM", "RandomForest", "LogitBoost", "K-NN")
  train_control <- trainControl(method="LOOCV", classProbs = TRUE, savePredictions = TRUE, verboseIter = FALSE)
  
  ### bar plot to show 'Day' distribution of the recall and random resting samples
  time_points <- c("d0", "d5", "d12", "d28", "d60", "d90", "d120", "d180")
  recall_day_dist <- sapply(time_points, function(x) length(which(Seurat_Obj@meta.data$Day[recall_idx] == x)))
  resting_day_dist <- sapply(time_points, function(x) length(which(Seurat_Obj@meta.data$Day[resting_idx] == x)))
  
  ### draw the bar plot
  png(paste0(outputDir2, "Barplot_Sample_Numbers_Distribution.png"), width = 2000, height = 1000, res = 110)
  par(mfrow = c(1,2))
  recall_bp <- barplot(recall_day_dist, main = paste0("The Number of Recall Samples (", length(recall_idx), " in total)"))
  text(recall_bp, 0, recall_day_dist, cex=1, pos=3)
  resting_bp <- barplot(resting_day_dist, main = paste0("The Number of Resting Samples (", length(resting_idx), " in total)"))
  text(resting_bp, 0, resting_day_dist, cex=1, pos=3)
  dev.off()
  
  ### because there are so many recall & resting samples, we randomly choose samples
  random_recall_idx <- sample(x = recall_idx, size = random_sampleNum)
  random_resting_idx <- sample(x = resting_idx, size = random_sampleNum)
  random_recall_day_dist <- sapply(time_points, function(x) length(which(Seurat_Obj@meta.data$Day[random_recall_idx] == x)))
  random_resting_day_dist <- sapply(time_points, function(x) length(which(Seurat_Obj@meta.data$Day[random_resting_idx] == x)))
  
  ### draw the bar plot with random samples that will be used for the classifier
  png(paste0(outputDir2, "Barplot_Random_Sample_Numbers_Distribution.png"), width = 2000, height = 1000, res = 110)
  par(mfrow = c(1,2))
  recall_bp <- barplot(random_recall_day_dist, main = paste0("The Number of Randomly Chosen Recall Samples (", length(random_recall_idx), " in total)"))
  text(recall_bp, 0, random_recall_day_dist, cex=1, pos=3)
  resting_bp <- barplot(random_resting_day_dist, main = paste0("The Number of Randomly Chosen Resting Samples (", length(random_resting_idx), " in total)"))
  text(resting_bp, 0, random_resting_day_dist, cex=1, pos=3)
  dev.off()
  
  
  #'******************************************************************************
  #' A function to transform RNA-Seq data with VST in DESeq2 package
  #' readCount: RNA-Seq rawcounts in a matrix or in a data frame form
  #'            Rows are genes and columns are samples
  #' filter_thresh: The function filters out genes that have at least one sample
  #'                with counts larger than the 'filter_thresh' value
  #'                e.g., if the 'filter_thresh' = 1, then it removes genes
  #'                that have counts <= 1 across all the samples
  #'                if 0, then there will be no filtering
  #'******************************************************************************
  normalizeRNASEQwithVST <- function(readCount, filter_thresh=1) {
    
    ### load library
    if(!require(DESeq2, quietly = TRUE)) {
      if(!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager")
      BiocManager::install("DESeq2")
      require(DESeq2, quietly = TRUE)
    }
    
    ### make a design matrix for DESeq2 data
    condition <- data.frame(factor(rep("OneClass", ncol(readCount))))
    
    ### Data preparation for DESeq2 format
    deSeqData <- DESeqDataSetFromMatrix(countData=readCount, colData=condition, design= ~0)
    
    if(filter_thresh > 0) {
      ### Remove rubbish rows - this will decrease the number of rows
      keep = apply(counts(deSeqData), 1, function(r){
        return(sum(r > filter_thresh) > 0)
      })
      deSeqData <- deSeqData[keep,]
    }
    
    ### VST
    vsd <- varianceStabilizingTransformation(deSeqData)
    transCnt <- data.frame(assay(vsd), check.names = FALSE)
    
    return (transCnt)
    
  }
  
  ### a function to select genes based on variance
  selectTopV <- function(x, selectNum) {
    v <- apply(x, 1, var)
    x <- x[order(-v),]
    x <- x[1:selectNum,]
    
    return (x)
  }
  
  ### normalize the read counts
  ### before the normalization, only keep the recall and random resting samples
  input_data <- normalizeRNASEQwithVST(readCount = data.frame(Seurat_Obj@assays$RNA@counts[,c(random_recall_idx,
                                                                                              random_resting_idx)],
                                                              stringsAsFactors = FALSE, check.names = FALSE))
  
  ### reduce the gene size based on variance
  ### only select high variance genes
  input_data <- selectTopV(input_data, featureSelectionNum)
  
  ### annotate class for the input data
  input_data <- data.frame(t(input_data), stringsAsFactors = FALSE, check.names = FALSE)
  input_data$Type <- factor(c(rep("Recall", length(random_recall_idx)),
                              rep("Resting", length(random_recall_idx))), levels = c("Resting", "Recall"))
  
  ### build classifier and test
  ### LOOCV
  p <- list()
  acc <- NULL
  for(i in 1:length(methodTypes)) {
    writeLines(paste(methodTypes[i]))
    model <- train(Type~., data=input_data, trControl=train_control, method=methodTypes[i])
    roc <- roc(model$pred$obs, model$pred$Recall)
    acc <- c(acc, round(mean(model$results$Accuracy), 3))
    p[[i]] <- plot.roc(roc, main = paste(methodNames[i], "Using Gene Expressions\n",
                                         "Accuracy =", acc[i]),
                       legacy.axes = TRUE, print.auc = TRUE, auc.polygon = TRUE,
                       xlim = c(1,0), ylim = c(0,1), grid = TRUE, cex.main = 1)
    gc()
  }
  
  ### draw ROC curves
  png(paste0(outputDir2, "Classifier_Recall_vs_Resting_AUCs_", featureSelectionNum, ".png"),
      width = 2000, height = 2000, res = 350)
  par(mfrow=c(3, 2))
  for(i in 1:length(methodTypes)) {
    plot.roc(p[[i]], main = paste(methodNames[i], "Using Gene Expressions\n",
                               "Accuracy =", acc[i]),
             legacy.axes = TRUE, print.auc = TRUE, auc.polygon = TRUE,
             xlim = c(1,0), ylim = c(0,1), grid = TRUE, cex.main = 1)
  }
  dev.off()
  
  
  ### in the "PCA_TFH_Cluster_17.png", pick cells with PC1 > 15 in d180
  ### and see how they changed over time
  ### also, which genes contribute a lot in the PC1? + pathway analysis
  
  ### new output directory for the results
  outputDir2 <- paste0(outputDir, "Cluster17-TFH/Interesting_Clones/")
  dir.create(outputDir2, showWarnings = FALSE, recursive = TRUE)
  
  ### get names of the cells of interest
  interesting_cells <- rownames(subset_Seurat_Obj@reductions$pca@cell.embeddings)[intersect(which(subset_Seurat_Obj@meta.data$Day == "d180"),
                                                                                            which(subset_Seurat_Obj@reductions$pca@cell.embeddings[,"PC_1"] > 15))]
  
  ### get clones of the cells
  interesting_clones <- unique(subset_Seurat_Obj@meta.data[interesting_cells,"clone_id"])
  
  ### get all the indicies of the clones
  all_interesting_clone_idx <- which(subset_Seurat_Obj@meta.data$clone_id %in% interesting_clones)
  
  ### make a table for lineage tracing
  interesting_tp <- unique(subset_Seurat_Obj@meta.data$Day[all_interesting_clone_idx])
  interesting_tp <- as.character(interesting_tp[order(interesting_tp)])
  interesting_clone_table <- matrix(0, length(interesting_clones), length(interesting_tp))
  rownames(interesting_clone_table) <- interesting_clones
  colnames(interesting_clone_table) <- interesting_tp
  
  ### fill out the table
  for(idx in all_interesting_clone_idx) {
    interesting_clone_table[subset_Seurat_Obj@meta.data$clone_id[idx],
                            as.character(subset_Seurat_Obj@meta.data$Day[idx])] <- interesting_clone_table[subset_Seurat_Obj@meta.data$clone_id[idx],
                                                                                                           as.character(subset_Seurat_Obj@meta.data$Day[idx])] + 1
  }
  
  ### save the table
  write.xlsx2(data.frame(Clone_ID=rownames(interesting_clone_table), interesting_clone_table,
                         stringsAsFactors = FALSE, check.names = FALSE),
              file = paste0(outputDir2, "Interesting_Clones_Lineage_Tracing.xlsx"),
              row.names = FALSE, sheetName = "Interesting_Clones")
  
  ### choose clones that were appeared in more than one time points
  target_clones <- interesting_clones[apply(interesting_clone_table, 1, function(x) {
    if(length(which(x > 0)) > 1) {
      return(TRUE)
    } else {
      return(FALSE)
    }
  })]
  
  ### draw PCA plot to describe changes of the clones over time
  p <- vector("list", length = length(target_clones))
  names(p) <- target_clones
  for(clone in target_clones) {
    plot_df <- data.frame(PC1=subset_Seurat_Obj@reductions$pca@cell.embeddings[which(subset_Seurat_Obj@meta.data$clone_id == clone),"PC_1"],
                          PC2=subset_Seurat_Obj@reductions$pca@cell.embeddings[which(subset_Seurat_Obj@meta.data$clone_id == clone),"PC_2"],
                          Day=subset_Seurat_Obj@meta.data$Day[which(subset_Seurat_Obj@meta.data$clone_id == clone)],
                          stringsAsFactors = FALSE, check.names = FALSE)
    p[[clone]] <- ggplot(plot_df, aes_string(x="PC1", y="PC2")) +
      geom_point(aes_string(col="Day"), size=2, alpha=0.8) +
      ggtitle(clone) +
      theme_classic(base_size = 16) +
      theme(plot.title = element_text(hjust = 0.5))
  }
  g <- arrangeGrob(grobs = p,
                   nrow = 2,
                   ncol = 2,
                   top = "Interesting-Clones Lineages")
  ggsave(file = paste0(outputDir2, "Interesting_Clones_Lineage_Tracing.png"), g, width = 10, height = 6, dpi = 500)
  
  ### interesting clone PCA ALL
  ### one PCA - all clones coloring with days
  plot_df <- data.frame(PC1=subset_Seurat_Obj@reductions$pca@cell.embeddings[all_interesting_clone_idx,"PC_1"],
                        PC2=subset_Seurat_Obj@reductions$pca@cell.embeddings[all_interesting_clone_idx,"PC_2"],
                        Day=subset_Seurat_Obj@meta.data$Day[all_interesting_clone_idx],
                        Clone=subset_Seurat_Obj@meta.data$clone_id[all_interesting_clone_idx],
                        stringsAsFactors = FALSE, check.names = FALSE)
  plot_df$Clone[which(!(plot_df$Clone %in% target_clones))] <- "unique_clones"
  ggplot(plot_df, aes_string(x="PC1", y="PC2")) +
    geom_point(aes_string(col="Day", shape="Clone"), size=2, alpha=1) +
    ggtitle("All_Interesting_Clones") +
    theme_classic(base_size = 16) +
    theme(plot.title = element_text(hjust = 0.5))
  ggsave(file = paste0(outputDir2, "All_Interesting_Clones_Over_Time.png"), width = 7, height = 4, dpi = 500)
  
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
              file = paste0(outputDir2, "PCA_Contributions.xlsx"),
              row.names = FALSE, sheetName = "PCA_Contributions")
  
  ### pathway analysis with the important genes of the PC1
  contb_threshold <- 0.1
  important_genes <- rownames(pca_contb)[which(pca_contb[,"PC_1"] > contb_threshold)]
  pathway_result_GO <- pathwayAnalysis_CP(geneList = mapIds(org.Hs.eg.db,
                                                            important_genes,
                                                            "ENTREZID", "SYMBOL"),
                                          org = "human", database = "GO",
                                          title = paste0("Pathway_Results_", length(important_genes), "_PC1_Genes_", contb_threshold),
                                          displayNum = 50, imgPrint = TRUE,
                                          dir = paste0(outputDir2))
  pathway_result_KEGG <- pathwayAnalysis_CP(geneList = mapIds(org.Hs.eg.db,
                                                              important_genes,
                                                              "ENTREZID", "SYMBOL"),
                                            org = "human", database = "KEGG",
                                            title = paste0("Pathway_Results_", length(important_genes), "_PC1_Genes_", contb_threshold),
                                            displayNum = 50, imgPrint = TRUE,
                                            dir = paste0(outputDir2))
  write.xlsx2(pathway_result_GO, file = paste0(outputDir2, "GO_pathway_results_", length(important_genes), "_PC1_Genes_", contb_threshold, ".xlsx"),
              row.names = FALSE, sheetName = paste0("GO_Results"))
  write.xlsx2(pathway_result_KEGG, file = paste0(outputDir2, "KEGG_pathway_results_", length(important_genes), "_PC1_Genes_", contb_threshold, ".xlsx"),
              row.names = FALSE, sheetName = paste0("KEGG_Results"))
  
  #
  ### same analysis with matched_SC and matched_bulk only with "yes"
  #
  
  ### new output directory for the results
  outputDir2 <- paste0(outputDir, "bulk_and_scTCR(PCR)/")
  dir.create(outputDir2, showWarnings = FALSE, recursive = TRUE)
  
  ### get "yes" indicies
  matched_indicies <- union(union(which(Seurat_Obj@meta.data$matched_SC_alpha_only == "yes"),
                                  which(Seurat_Obj@meta.data$matched_SC_beta_only == "yes")),
                            union(which(Seurat_Obj@meta.data$matched_bulk_alpha == "yes"),
                                  which(Seurat_Obj@meta.data$matched_bulk_beta == "yes")))
  
  ### see UMAP of the "yes" cells
  plot_df <- data.frame(UMAP1=Seurat_Obj@reductions$umap@cell.embeddings[,"UMAP_1"],
                        UMAP2=Seurat_Obj@reductions$umap@cell.embeddings[,"UMAP_2"],
                        Day=Seurat_Obj@meta.data$Day,
                        Clone=Seurat_Obj@meta.data$clone_id,
                        Cluster17=ifelse(Seurat_Obj@meta.data$clone_id %in% cluster_17_clone_ids, "yes", "no"),
                        Matched=ifelse(1:nrow(Seurat_Obj@meta.data) %in% matched_indicies, "yes", "no"),
                        stringsAsFactors = FALSE, check.names = FALSE)
  p <- list()
  p[[1]] <- ggplot(plot_df, aes_string(x="UMAP1", y="UMAP2")) +
    geom_point(aes_string(col="Cluster17"), size=1, alpha=1) +
    ggtitle("Cluster17 Clones") +
    theme_classic(base_size = 16) +
    theme(plot.title = element_text(hjust = 0.5))
  p[[2]] <- ggplot(plot_df, aes_string(x="UMAP1", y="UMAP2")) +
    geom_point(aes_string(col="Matched"), size=1, alpha=1) +
    ggtitle("Bulk TCR and scTCR(PCR) Clones") +
    theme_classic(base_size = 16) +
    theme(plot.title = element_text(hjust = 0.5))
  g <- arrangeGrob(grobs = p,
                   nrow = 1,
                   ncol = 2,
                   top = "Cluster17 vs bulk_and_scTCR(PCR)")
  ggsave(file = paste0(outputDir2, "UMAP_comparison.png"), g, width = 14, height = 6, dpi = 200)
  
  ### start to perform the same analysis as for the cluster 17
  
  ### get a subset for the bulk_and_scTCR(PCR)
  new.ident <- rep(NA, nrow(Seurat_Obj@meta.data))
  new.ident[matched_indicies] <- "matched"
  Idents(object = Seurat_Obj) <- new.ident
  subset_Seurat_Obj <- subset(Seurat_Obj, idents=c("matched"))
  
  ### factorize the Day column
  subset_Seurat_Obj@meta.data$Day <- factor(subset_Seurat_Obj@meta.data$Day,
                                            levels = c("d0", "d5", "d12", "d28", "d60", "d90", "d120", "d180"))
  
  ### run PCA
  subset_Seurat_Obj <- RunPCA(subset_Seurat_Obj, npcs = 10)
  
  ### PCA plot
  g <- multiReducPlot(x = subset_Seurat_Obj@reductions$pca@cell.embeddings[,"PC_1"],
                      y = subset_Seurat_Obj@reductions$pca@cell.embeddings[,"PC_2"],
                      group = subset_Seurat_Obj@meta.data$Day, type = "PCA",
                      fName = "PCA_bulk_and_scTCR(PCR)", isPrint = TRUE, isConvex = TRUE)
  ggsave(file = paste0(outputDir2, "PCA_bulk_and_scTCR(PCR).png"), g, width = 20, height = 10, dpi = 300)
  
  ### run UMAP
  subset_Seurat_Obj <- RunUMAP(subset_Seurat_Obj, dims = 1:5)
  
  ### UMAP plot
  g <- multiReducPlot(x = subset_Seurat_Obj@reductions$umap@cell.embeddings[,"UMAP_1"],
                      y = subset_Seurat_Obj@reductions$umap@cell.embeddings[,"UMAP_2"],
                      group = subset_Seurat_Obj@meta.data$Day, type = "UMAP",
                      fName = "UMAP_bulk_and_scTCR(PCR)", isPrint = TRUE, isConvex = TRUE)
  ggsave(file = paste0(outputDir2, "UMAP_bulk_and_scTCR(PCR).png"), g, width = 20, height = 10, dpi = 300)
  
  
  ### make a clone summary table
  ### here I checked same match.cdr always has the same clone_id
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
  write.xlsx2(clone_summary_table, file = paste0(outputDir2, "Clone_Count_Summary.xlsx"),
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
  
  ### draw the alluvial plot
  ggplot(plot_df,
         aes(x = Time, stratum = Clone, alluvium = Clone,
             y = Clone_Size,
             fill = Clone, label = CDR3)) +
    ggtitle("Clonal Tracing of the bulk_and_scTCR(PCR) Matched Cells") +
    geom_flow() +
    geom_stratum(alpha = 1) +
    geom_text(stat = "stratum", size = 0.8) +
    rotate_x_text(90) +
    theme_pubr(legend = "none") +
    theme(axis.title.x = element_blank()) +
    theme_cleveland2() +
    scale_fill_viridis(discrete = T) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, NA))
  ggsave(file = paste0(outputDir2, "bulk_and_scTCR(PCR)_Clonal_Tracing.png"), width = 18, height = 9, dpi = 300)
  
  
  ### separate the LN and PB cells
  clone_summary_table_LNPB <- data.frame(sapply(clone_summary_table, function(x) c(rbind(x, x, x))),
                                         stringsAsFactors = FALSE, check.names = FALSE)
  clone_summary_table_LNPB <- data.frame(clone_summary_table_LNPB[,c("clone_id", "cdr_ab")],
                                         cell_type=rep(c("LN", "PB", "ALL"), nrow(clone_summary_table)),
                                         sapply(clone_summary_table_LNPB[,c(time_points, "total_count")],
                                                as.numeric),
                                         stringsAsFactors = FALSE, check.names = FALSE)
  rownames(clone_summary_table_LNPB) <- paste(c(rbind(rownames(clone_summary_table),
                                                      rownames(clone_summary_table),
                                                      rownames(clone_summary_table))),
                                              c("LN", "PB", "ALL"), sep = "_")
  
  ### fill out the new table
  for(i in 1:nrow(clone_summary_table_LNPB)) {
    if(clone_summary_table_LNPB$cell_type[i] != "ALL") {
      clone_idx <- intersect(which(subset_Seurat_Obj@meta.data$clone_id == clone_summary_table_LNPB$clone_id[i]),
                             which(subset_Seurat_Obj@meta.data$Tissue == clone_summary_table_LNPB$cell_type[i]))
      for(tp in time_points) {
        clone_summary_table_LNPB[i,tp] <- length(intersect(clone_idx, tp_indicies[[tp]]))
      }
      clone_summary_table_LNPB$total_count[i] <- sum(clone_summary_table_LNPB[i,time_points])
    }
  }
  
  ### pb-associated table
  ### there are only 12 PB cells in this subset
  pb_idx <- intersect(which(clone_summary_table_LNPB$cell_type == "PB"),
                      which(clone_summary_table_LNPB$total_count > 0))
  clone_summary_table_PB <- clone_summary_table_LNPB[c(rbind(pb_idx-1,
                                                             pb_idx,
                                                             pb_idx+1)),]
  
  ### only keep lineages from the pb-associated table
  line_idx <- intersect(which(clone_summary_table_PB$total_count > 1),
                        which(clone_summary_table_PB$cell_type == "ALL"))
  lineage_table_PB <- clone_summary_table_PB[c(rbind(line_idx-2,
                                                     line_idx-1,
                                                     line_idx)),]
  
  ### save it in Excel file
  write.xlsx2(lineage_table_PB, file = paste0(outputDir2, "PB_Associated_Lineages.xlsx"),
              sheetName = "PB_Lineages", row.names = FALSE)
  
  ### get the number of cells for each group
  ln_cellNum_subset <- sapply(time_points, function(x) {
    return(length(intersect(which(subset_Seurat_Obj@meta.data$Tissue == "LN"),
                            which(subset_Seurat_Obj@meta.data$Day == x))))
  })
  pb_cellNum_subset <- sapply(time_points, function(x) {
    return(length(intersect(which(subset_Seurat_Obj@meta.data$Tissue == "PB"),
                            which(subset_Seurat_Obj@meta.data$Day == x))))
  })
  ln_cellNum_all <- sapply(time_points, function(x) {
    return(length(intersect(intersect(which(Seurat_Obj@meta.data$Tissue == "LN"),
                                      which(Seurat_Obj@meta.data$Day == x)),
                            which(!is.na(Seurat_Obj@meta.data$match.cdr)))))
  })
  pb_cellNum_all <- sapply(time_points, function(x) {
    return(length(intersect(intersect(which(Seurat_Obj@meta.data$Tissue == "PB"),
                                      which(Seurat_Obj@meta.data$Day == x)),
                            which(!is.na(Seurat_Obj@meta.data$match.cdr)))))
  })
  
  ### Alluvial plot - visualization of the lineage tracing (PB-associated lineages only)
  
  ### get an input data frame for the alluvial plot
  plot_df <- plot_df[which(plot_df$Clone %in% unique(lineage_table_PB$clone_id)),]
  
  ### add the number of PB cells
  plot_df$PB_Num <- ""
  for(i in 1:nrow(plot_df)) {
    num <- lineage_table_PB[paste0(plot_df$Clone[i], "_PB"),as.character(plot_df$Time[i])]
    if(num > 0) {
      plot_df$PB_Num[i] <- num
    }
  }
  
  ### draw the alluvial plot
  ggplot(plot_df,
         aes(x = Time, stratum = Clone, alluvium = Clone,
             y = Clone_Size,
             fill = CDR3, label = PB_Num)) +
    ggtitle("Clonal Tracing of the bulk_and_scTCR(PCR) Matched Cells (PB-Associated Lineages)") +
    geom_stratum(alpha = 1) +
    geom_text(stat = "stratum", size = 3, col = "black") +
    geom_flow() +
    rotate_x_text(90) +
    theme_pubr(legend = "none") +
    theme(axis.title.x = element_blank()) +
    theme_cleveland2() +
    scale_fill_viridis(discrete = T) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, NA))
  ggsave(file = paste0(outputDir2, "bulk_and_scTCR(PCR)_Clonal_Tracing_PB.png"), width = 18, height = 9, dpi = 300)
  
  
  #
  ### PCA & UMAP with the top 9 clones from the bulk_and_scTCR(PCR) result
  #
  
  ### subset preparation for the PCA & UMAP
  Idents(object = Seurat_Obj) <- Seurat_Obj@meta.data$clone_id
  subset_Seurat_Obj2 <- subset(Seurat_Obj, idents=lineage_table$clone_id[1:9])
  subset_Seurat_Obj2@meta.data$Day <- factor(subset_Seurat_Obj2@meta.data$Day,
                                            levels = c("d0", "d5", "d12", "d28", "d60", "d90", "d120", "d180"))
  
  ### draw the PCA & UMAP
  DimPlot(subset_Seurat_Obj2, reduction = "pca", split.by = "clone_id", group.by = "Day",
          ncol = 3, pt.size = 2) +
    labs(title = "PCA of the Top 9 Clones")
  ggsave(file = paste0(outputDir2, "PCA_Top_9_Clones_bulk_and_scTCR(PCR).png"), width = 20, height = 10, dpi = 300)
  DimPlot(subset_Seurat_Obj2, reduction = "umap", split.by = "clone_id", group.by = "Day",
          ncol = 3, pt.size = 2) +
    labs(title = "UMAP of the Top 9 Clones")
  ggsave(file = paste0(outputDir2, "UMAP_Top_9_Clones_bulk_and_scTCR(PCR).png"), width = 20, height = 10, dpi = 300)
  
  
  #
  ### Markers associated with the interesting clones in the bulk_and_scTCR(PCR) matched cells
  #
  
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
    new.ident <- rep(NA, nrow(Seurat_Obj@meta.data))
    new.ident[intersect(matched_indicies,
                        which(Seurat_Obj@meta.data$clone_id == clone))] <- "ident1"
    new.ident[setdiff(matched_indicies,
                      which(Seurat_Obj@meta.data$clone_id == clone))] <- "ident2"
    Idents(object = Seurat_Obj) <- new.ident
    
    ### DE analysis with DESeq2
    de_result <- FindMarkers(Seurat_Obj,
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
    write.xlsx2(de_result, file = paste0(outputDir4, "DESeq2_", clone, "_vs_others.xlsx"),
                sheetName = paste0(clone, "_vs_others"), row.names = FALSE)
    
    ### pathway analysis
    if(length(which(de_result$FDR < de_signif_thresh)) > 0) {
      pathway_result_GO <- pathwayAnalysis_CP(geneList = mapIds(org.Hs.eg.db,
                                                                de_result[which(de_result$FDR < de_signif_thresh),
                                                                          "Gene_Symbol"],
                                                                "ENTREZID", "SYMBOL"),
                                              org = "human", database = "GO",
                                              title = paste0("Pathway_Results_", clone, "_vs_others"),
                                              displayNum = 50, imgPrint = TRUE,
                                              dir = paste0(outputDir4))
      pathway_result_KEGG <- pathwayAnalysis_CP(geneList = mapIds(org.Hs.eg.db,
                                                                  de_result[which(de_result$FDR < de_signif_thresh),
                                                                            "Gene_Symbol"],
                                                                  "ENTREZID", "SYMBOL"),
                                                org = "human", database = "KEGG",
                                                title = paste0("Pathway_Results_", clone, "_vs_others"),
                                                displayNum = 50, imgPrint = TRUE,
                                                dir = paste0(outputDir4))
      if(!is.null(pathway_result_GO) && nrow(pathway_result_GO) > 0) {
        write.xlsx2(pathway_result_GO, file = paste0(outputDir4, "GO_pathway_results_", clone, "_vs_others.xlsx"),
                    row.names = FALSE, sheetName = paste0("GO_Results_", clone, "_vs_others"))
      }
      if(!is.null(pathway_result_KEGG) && nrow(pathway_result_KEGG) > 0) {
        write.xlsx2(pathway_result_KEGG, file = paste0(outputDir4, "KEGG_pathway_results_", clone, "_vs_others.xlsx"),
                    row.names = FALSE, sheetName = paste0("KEGG_Results_", clone, "_vs_others"))
      }
    }
    
  }
  
  
  ### in the "PCA_bulk_and_scTCR(PCR).png", pick cells with PC1 > 5 in d180
  ### and see how they changed over time
  ### also, which genes contribute a lot in the PC1? + pathway analysis
  
  ### new output directory for the results
  outputDir3 <- paste0(outputDir2, "Interesting_Clones/")
  dir.create(outputDir3, showWarnings = FALSE, recursive = TRUE)
  
  ### get names of the cells of interest
  interesting_cells <- rownames(subset_Seurat_Obj@reductions$pca@cell.embeddings)[intersect(which(subset_Seurat_Obj@meta.data$Day == "d180"),
                                                                                            which(subset_Seurat_Obj@reductions$pca@cell.embeddings[,"PC_1"] > 5))]
  
  ### get clones of the cells
  interesting_clones <- unique(subset_Seurat_Obj@meta.data[interesting_cells,"clone_id"])
  
  ### get all the indicies of the clones
  all_interesting_clone_idx <- which(subset_Seurat_Obj@meta.data$clone_id %in% interesting_clones)
  
  ### make a table for lineage tracing
  interesting_tp <- unique(subset_Seurat_Obj@meta.data$Day[all_interesting_clone_idx])
  interesting_tp <- as.character(interesting_tp[order(interesting_tp)])
  interesting_clone_table <- matrix(0, length(interesting_clones), length(interesting_tp))
  rownames(interesting_clone_table) <- interesting_clones
  colnames(interesting_clone_table) <- interesting_tp
  
  ### fill out the table
  for(idx in all_interesting_clone_idx) {
    interesting_clone_table[subset_Seurat_Obj@meta.data$clone_id[idx],
                            as.character(subset_Seurat_Obj@meta.data$Day[idx])] <- interesting_clone_table[subset_Seurat_Obj@meta.data$clone_id[idx],
                                                                                                           as.character(subset_Seurat_Obj@meta.data$Day[idx])] + 1
  }
  
  ### save the table
  write.xlsx2(data.frame(Clone_ID=rownames(interesting_clone_table), interesting_clone_table,
                         stringsAsFactors = FALSE, check.names = FALSE),
              file = paste0(outputDir3, "Interesting_Clones_Lineage_Tracing.xlsx"),
              row.names = FALSE, sheetName = "Interesting_Clones")
  
  ### choose clones that were appeared in more than one time points
  target_clones <- interesting_clones[apply(interesting_clone_table, 1, function(x) {
    if(length(which(x > 0)) > 1) {
      return(TRUE)
    } else {
      return(FALSE)
    }
  })]
  
  ### draw PCA plot to describe changes of the clones over time
  p <- vector("list", length = length(target_clones))
  names(p) <- target_clones
  for(clone in target_clones) {
    plot_df <- data.frame(PC1=subset_Seurat_Obj@reductions$pca@cell.embeddings[which(subset_Seurat_Obj@meta.data$clone_id == clone),"PC_1"],
                          PC2=subset_Seurat_Obj@reductions$pca@cell.embeddings[which(subset_Seurat_Obj@meta.data$clone_id == clone),"PC_2"],
                          Day=subset_Seurat_Obj@meta.data$Day[which(subset_Seurat_Obj@meta.data$clone_id == clone)],
                          stringsAsFactors = FALSE, check.names = FALSE)
    p[[clone]] <- ggplot(plot_df, aes_string(x="PC1", y="PC2")) +
      geom_point(aes_string(col="Day"), size=2, alpha=0.8) +
      ggtitle(clone) +
      theme_classic(base_size = 16) +
      theme(plot.title = element_text(hjust = 0.5))
  }
  g <- arrangeGrob(grobs = p,
                   nrow = 3,
                   ncol = 3,
                   top = "Interesting-Clones Lineages")
  ggsave(file = paste0(outputDir3, "Interesting_Clones_Lineage_Tracing.png"), g, width = 10, height = 6, dpi = 500)
  
  ### interesting clone PCA ALL
  ### one PCA - all clones coloring with days
  plot_df <- data.frame(PC1=subset_Seurat_Obj@reductions$pca@cell.embeddings[all_interesting_clone_idx,"PC_1"],
                        PC2=subset_Seurat_Obj@reductions$pca@cell.embeddings[all_interesting_clone_idx,"PC_2"],
                        Day=subset_Seurat_Obj@meta.data$Day[all_interesting_clone_idx],
                        Clone=subset_Seurat_Obj@meta.data$clone_id[all_interesting_clone_idx],
                        stringsAsFactors = FALSE, check.names = FALSE)
  plot_df$Clone[which(!(plot_df$Clone %in% target_clones))] <- "unique_clones"
  plot_df$Clone <- factor(plot_df$Clone)
  ggplot(plot_df, aes_string(x="PC1", y="PC2")) +
    geom_point(aes_string(col="Day", shape="Clone"), size=2, alpha=1) +
    scale_shape_manual(values=1:nlevels(plot_df$Clone)) +
    ggtitle("All_Interesting_Clones") +
    theme_classic(base_size = 16) +
    theme(plot.title = element_text(hjust = 0.5))
  ggsave(file = paste0(outputDir3, "All_Interesting_Clones_Over_Time.png"), width = 10, height = 6, dpi = 500)
  
  
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
              file = paste0(outputDir3, "PCA_Contributions.xlsx"),
              row.names = FALSE, sheetName = "PCA_Contributions")
  
  ### pathway analysis with the important genes of the PC1
  contb_threshold <- 0.1
  important_genes <- rownames(pca_contb)[which(pca_contb[,"PC_1"] > contb_threshold)]
  pathway_result_GO <- pathwayAnalysis_CP(geneList = mapIds(org.Hs.eg.db,
                                                            important_genes,
                                                            "ENTREZID", "SYMBOL"),
                                          org = "human", database = "GO",
                                          title = paste0("Pathway_Results_", length(important_genes), "_PC1_Genes_", contb_threshold),
                                          displayNum = 50, imgPrint = TRUE,
                                          dir = paste0(outputDir3))
  pathway_result_KEGG <- pathwayAnalysis_CP(geneList = mapIds(org.Hs.eg.db,
                                                              important_genes,
                                                              "ENTREZID", "SYMBOL"),
                                            org = "human", database = "KEGG",
                                            title = paste0("Pathway_Results_", length(important_genes), "_PC1_Genes_", contb_threshold),
                                            displayNum = 50, imgPrint = TRUE,
                                            dir = paste0(outputDir3))
  write.xlsx2(pathway_result_GO, file = paste0(outputDir3, "GO_pathway_results_", length(important_genes), "_PC1_Genes_", contb_threshold, ".xlsx"),
              row.names = FALSE, sheetName = paste0("GO_Results"))
  write.xlsx2(pathway_result_KEGG, file = paste0(outputDir3, "KEGG_pathway_results_", length(important_genes), "_PC1_Genes_", contb_threshold, ".xlsx"),
              row.names = FALSE, sheetName = paste0("KEGG_Results"))
  
}
