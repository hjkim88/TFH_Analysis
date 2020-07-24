###
#   File name : TFH_Additional_Analyses2.R
#   Author    : Hyunjin Kim
#   Date      : July 24, 2020
#   Email     : hyunjin.kim@stjude.org
#   Purpose   : Here are additional things to do:
#               1. A heatmap of gene expressions of the genes that were contributed to the PC1 the most in the PCA
#                  to see if there is a change over time - ordered by time or sub-side color for time
#               2. Diversity & Clonality - The diversity means how many clones in each time point (time point based)
#                  and the clonality means time point distribution in each clone (clone based)
#               3. Apply trajectory inference method (Pseudotime analysis) (i.e., Slingshot) on
#                  the gene expression data
#   
#   Instruction
#               1. Source("TFH_Additional_Analyses.R")
#               2. Run the function "tfh_additional_analyses" - specify the input file path and the output directory
#               3. The results will be generated under the output directory
#
#   Example
#               > source("The_directory_of_TFH_Additional_Analyses.R/TFH_Additional_Analyses.R")
#               > tfh_additional_analyses2(Seurat_RObj_path="./data/Ali_Tcell_combined_NEW.RDATA",
#                                          outputDir="./results/")
###

tfh_additional_analyses2 <- function(Seurat_RObj_path="./data/Ali_Tcell_combined_NEW.RDATA",
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
  
  
  
  
  
}
