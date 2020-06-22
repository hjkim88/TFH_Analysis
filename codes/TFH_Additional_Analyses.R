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
  
  ### Adding more informative columns to the meta.data
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
  expansion_threshold <- 5
  
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
  methodTypes <- c("svmLinear", "svmRadial", "gbm", "rf", "LogitBoost", "knn")
  methodNames <- c("SVMLinear", "SVMRadial", "GBM", "RandomForest", "LogitBoost", "K-NN")
  train_control <- trainControl(method="LOOCV", classProbs = TRUE, savePredictions = TRUE, verboseIter = FALSE)
  
  ### because there are so many resting samples, we randomly choose the resting samples
  ### as the same number of the recall samples
  random_resting_idx <- sample(x = resting_idx, size = length(recall_idx))
  
  ### bar plot to show 'Day' distribution of the recall and random resting samples
  time_points <- c("d0", "d5", "d12", "d28", "d60", "d90", "d120", "d180")
  recall_day_dist <- sapply(time_points, function(x) length(which(Seurat_Obj@meta.data$Day[recall_idx] == x)))
  resting_day_dist <- sapply(time_points, function(x) length(which(Seurat_Obj@meta.data$Day[random_resting_idx] == x)))
  
  ### draw the bar plot
  png(paste0(outputDir2, "Barplot_Sample_Numbers_Classifier.png"), width = 2000, height = 1000, res = 110)
  par(mfrow = c(1,2))
  recall_bp <- barplot(recall_day_dist, main = "The Number of Recall Samples (1092 in total) for Classifier")
  text(recall_bp, 0, recall_day_dist, cex=1, pos=3)
  resting_bp <- barplot(resting_day_dist, main = "The Number of Resting Samples (1092 in total) for Classifier")
  text(resting_bp, 0, resting_day_dist, cex=1, pos=3)
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
  input_data <- normalizeRNASEQwithVST(readCount = data.frame(Seurat_Obj@assays$RNA@counts[,c(recall_idx,
                                                                                              random_resting_idx)],
                                                              stringsAsFactors = FALSE, check.names = FALSE))
  
  ### reduce the gene size based on variance
  ### only select high variance genes
  input_data <- selectTopV(input_data, featureSelectionNum)
  
  ### annotate class for the input data
  input_data <- data.frame(t(input_data), stringsAsFactors = FALSE, check.names = FALSE)
  input_data$Type <- factor(c(rep("Recall", length(recall_idx)),
                              rep("Resting", length(recall_idx))), levels = c("Resting", "Recall"))
  
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
  
}
