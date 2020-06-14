###
#   File name : Clonal_Tracking.R
#   Author    : Hyunjin Kim
#   Date      : Jun 1, 2020
#   Email     : hyunjin.kim@stjude.org
#   Purpose   : Here's the question:
#               How do the Tfh cells of the same clone changing over time with regard to their phenotype and numbers?
#               Can we detect known Tfh clones from the LN in the blood, does this have some sort of pattern
#               with maybe a 1 cell in the blood early on turns into many in the LN, or the other way,
#               where we have expanded clones in the LN become detectable in blood later, and is there an alteration
#               in these cases?
#
#   Instruction
#               1. Source("Clonal_Tracking.R")
#               2. Run the function "clonal_tracking" - specify the input file path and the output directory
#               3. The results will be generated under the output directory
#
#   Example
#               > source("The_directory_of_Clonal_Tracking.R/Clonal_Tracking.R")
#               > clonal_tracking(Seurat_RObj_path="./data/Ali_Tcell_combined.RDATA",
#                                 outputDir="./results/")
###

clonal_tracking <- function(Seurat_RObj_path="./data/Ali_Tcell_combined.RDATA",
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
  
  ### make a clone summary table
  ### here I checked same match.cdr always has the same clone_id
  unique_clone_idx <- which(!duplicated(cluster_17_clones_meta.data$clone_id))
  clone_summary_table <- data.frame(clone_id=cluster_17_clones_meta.data$clone_id[unique_clone_idx],
                                    cdr_ab=cluster_17_clones_meta.data$match.cdr[unique_clone_idx])
  rownames(clone_summary_table) <- clone_summary_table$clone_id
  
  ### add time point counts and the total count of the clonotypes
  ### first you should check unique(cluster_17_clones_meta.data$Day)
  time_points <- c("d0", "d5", "d12", "d28", "d60", "d90", "d120", "d180")
  clone_summary_table[time_points] <- 0
  clone_summary_table$total_count <- 0
  
  ### fill out the counts
  tp_indicies <- lapply(time_points, function(x) which(cluster_17_clones_meta.data$Day == x))
  names(tp_indicies) <- time_points
  for(clone in rownames(clone_summary_table)) {
    clone_idx <- which(cluster_17_clones_meta.data$clone_id == clone)
    for(tp in time_points) {
      clone_summary_table[clone,tp] <- length(intersect(clone_idx, tp_indicies[[tp]]))
    }
  }
  clone_summary_table$total_count <- as.numeric(apply(clone_summary_table[,time_points], 1, sum))
  
  ### order by the total_count
  clone_summary_table <- clone_summary_table[order(-clone_summary_table$total_count),]
  
  ### save the table as Excel file
  write.xlsx2(clone_summary_table, file = paste0(outputDir, "Clone_Count_Summary.xlsx"),
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
    ggtitle("Clonal Tracing of the TFH-related Cells (The Cluster 17)") +
    geom_flow() +
    geom_stratum(alpha = 1) +
    geom_text(stat = "stratum", size = 0.8) +
    rotate_x_text(90) +
    theme_pubr(legend = "none") +
    theme(axis.title.x = element_blank()) +
    theme_cleveland2() +
    scale_fill_viridis(discrete = T) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, NA))
  ggsave(file = paste0(outputDir, "TFH_Cluster_17_Clonal_Tracing.png"), width = 18, height = 9, dpi = 300)
  
  
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
      clone_idx <- intersect(which(cluster_17_clones_meta.data$clone_id == clone_summary_table_LNPB$clone_id[i]),
                             which(cluster_17_clones_meta.data$Tissue == clone_summary_table_LNPB$cell_type[i]))
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
  write.xlsx2(lineage_table_PB, file = paste0(outputDir, "PB_Associated_Lineages.xlsx"),
              sheetName = "PB_Lineages", row.names = FALSE)
  
  ### get the number of cells for each group
  ln_cellNum_subset <- sapply(time_points, function(x) {
    return(length(intersect(which(cluster_17_clones_meta.data$Tissue == "LN"),
                            which(cluster_17_clones_meta.data$Day == x))))
  })
  pb_cellNum_subset <- sapply(time_points, function(x) {
    return(length(intersect(which(cluster_17_clones_meta.data$Tissue == "PB"),
                            which(cluster_17_clones_meta.data$Day == x))))
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
  lnpb_only_idx <- which(lineage_table_PB$cell_type != "ALL")
  total_rows <- length(which(lineage_table_PB[lnpb_only_idx,time_points] > 0))
  plot_df <- data.frame(Tissue=rep("", total_rows),
                        Time=rep("", total_rows),
                        Clone_Size=rep(0, total_rows),
                        Clone=rep("", total_rows),
                        CDR3=rep("", total_rows))
  cnt <- 1
  for(i in lnpb_only_idx) {
    for(tp in time_points) {
      if(lineage_table_PB[i,tp] > 0) {
        plot_df[cnt,] <- c(lineage_table_PB$cell_type[i],
                           tp,
                           lineage_table_PB[i,tp],
                           lineage_table_PB$clone_id[i],
                           lineage_table_PB$cdr_ab[i])
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
             fill = CDR3, label = CDR3)) +
    ggtitle("Clonal Tracing of the TFH-related Cells (PB-Associated Lineages)") +
    geom_stratum(alpha = 1) +
    geom_flow() +
    rotate_x_text(90) +
    theme_pubr(legend = "right") +
    theme(axis.title.x = element_blank()) +
    theme_cleveland2() +
    scale_fill_jco(name="CDR3 (TCRa:TCRb)") +
    scale_y_continuous(expand = c(0, 0), limits = c(0, NA))
  ggsave(file = paste0(outputDir, "TFH_Cluster_17_Clonal_Tracing.png"), width = 18, height = 9, dpi = 300)
  
  
  #
  ### PCA & UMAP with the top 9 clones from the TFH result
  #
  
  ### subset preparation for the PCA & UMAP
  Idents(object = Seurat_Obj) <- Seurat_Obj@meta.data$clone_id
  subset_Seurat_Obj <- subset(Seurat_Obj, idents=lineage_table$clone_id[1:9])
  subset_Seurat_Obj@meta.data$Day <- factor(subset_Seurat_Obj@meta.data$Day,
                                            levels = c("d0", "d5", "d12", "d28", "d60", "d90", "d120", "d180"))
  
  ### draw the PCA & UMAP
  DimPlot(subset_Seurat_Obj, reduction = "pca", split.by = "clone_id", group.by = "Day",
          ncol = 3, pt.size = 2) +
    labs(title = "PCA of the Top 9 Clones")
  ggsave(file = paste0(outputDir, "PCA_Top_9_Clones_TFH.png"), width = 20, height = 10, dpi = 300)
  DimPlot(subset_Seurat_Obj, reduction = "umap", split.by = "clone_id", group.by = "Day",
          ncol = 3, pt.size = 2) +
    labs(title = "UMAP of the Top 9 Clones")
  ggsave(file = paste0(outputDir, "UMAP_Top_9_Clones_TFH.png"), width = 20, height = 10, dpi = 300)
  
  
  #
  ### Markers associated with the interesting clones in the TFH cluster
  #
  
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
  
  ### set clone size threshold
  ### clones with clone size > clone size threshold will be used
  clone_size_thresh <- 5
  
  ### DE gene FDR threshold
  de_signif_thresh <- 0.05
  
  ### get clone names that will be used
  clone_names <- clone_summary_table$clone_id[which(clone_summary_table$total_count > clone_size_thresh)]
  
  ### set new output directory
  outputDir2 <- paste0(outputDir, "Clone_Markers/")
  dir.create(path = outputDir2, showWarnings = FALSE, recursive = TRUE)
  
  ### for each interesting clone perform marker discovery
  for(clone in clone_names) {
    
    ### print progress
    writeLines(paste(clone))
    
    ### set new output directory
    outputDir3 <- paste0(outputDir2, clone, "/")
    dir.create(path = outputDir3, showWarnings = FALSE, recursive = TRUE)
    
    ### Ident configure
    new.ident <- rep(NA, nrow(Seurat_Obj@meta.data))
    new.ident[which(Seurat_Obj@meta.data$clone_id == clone)] <- "ident1"
    new.ident[setdiff(which(Seurat_Obj@meta.data$clone_id %in% cluster_17_clone_ids),
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
    write.xlsx2(de_result, file = paste0(outputDir3, "DESeq2_", clone, "_vs_Cluster17.xlsx"),
                sheetName = paste0(clone, "_vs_Cluster17"), row.names = FALSE)
    
    ### pathway analysis
    if(length(which(de_result$FDR < de_signif_thresh)) > 0) {
      pathway_result_GO <- pathwayAnalysis_CP(geneList = mapIds(org.Hs.eg.db,
                                                                de_result[which(de_result$FDR < de_signif_thresh),
                                                                          "Gene_Symbol"],
                                                                "ENTREZID", "SYMBOL"),
                                              org = "human", database = "GO",
                                              title = paste0("Pathway_Results_", clone, "_vs_Cluster17"),
                                              displayNum = 50, imgPrint = TRUE,
                                              dir = paste0(outputDir3))
      pathway_result_KEGG <- pathwayAnalysis_CP(geneList = mapIds(org.Hs.eg.db,
                                                                  de_result[which(de_result$FDR < de_signif_thresh),
                                                                            "Gene_Symbol"],
                                                                  "ENTREZID", "SYMBOL"),
                                                org = "human", database = "KEGG",
                                                title = paste0("Pathway_Results_", clone, "_vs_Cluster17"),
                                                displayNum = 50, imgPrint = TRUE,
                                                dir = paste0(outputDir3))
      if(!is.null(pathway_result_GO) && nrow(pathway_result_GO) > 0) {
        write.xlsx2(pathway_result_GO, file = paste0(outputDir3, "GO_pathway_results_", clone, "_vs_Cluster17.xlsx"),
                    row.names = FALSE, sheetName = paste0("GO_Results_", clone, "_vs_Cluster17"))
      }
      if(!is.null(pathway_result_KEGG) && nrow(pathway_result_KEGG) > 0) {
        write.xlsx2(pathway_result_KEGG, file = paste0(outputDir3, "KEGG_pathway_results_", clone, "_vs_Cluster17.xlsx"),
                    row.names = FALSE, sheetName = paste0("KEGG_Results_", clone, "_vs_Cluster17"))
      }
    }
    
  }
  
  
  #
  ### Now this is not a lineage tracing only for the cluster 17 but for all the cells
  #
  
  ### make a clone summary table
  ### here I checked same match.cdr always has the same clone_id
  unique_clone_idx <- intersect(which(!duplicated(Seurat_Obj@meta.data$clone_id)),
                                which(!is.na(Seurat_Obj@meta.data$clone_id)))
  clone_summary_table <- data.frame(clone_id=Seurat_Obj@meta.data$clone_id[unique_clone_idx],
                                    cdr_ab=Seurat_Obj@meta.data$match.cdr[unique_clone_idx])
  rownames(clone_summary_table) <- clone_summary_table$clone_id
  
  ### add time point counts and the total count of the clonotypes
  ### first you should check unique(cluster_17_clones_meta.data$Day)
  time_points <- c("d0", "d5", "d12", "d28", "d60", "d90", "d120", "d180")
  clone_summary_table[time_points] <- 0
  clone_summary_table$total_count <- 0
  
  ### start time
  start_time <- Sys.time()
  
  ### set progress bar
  cnt <- 0
  pb <- txtProgressBar(min = 0, max = nrow(clone_summary_table), style = 3)
  
  ### fill out the counts
  tp_indicies <- lapply(time_points, function(x) which(Seurat_Obj@meta.data$Day == x))
  names(tp_indicies) <- time_points
  for(clone in rownames(clone_summary_table)) {
    clone_idx <- which(Seurat_Obj@meta.data$clone_id == clone)
    for(tp in time_points) {
      clone_summary_table[clone,tp] <- length(intersect(clone_idx, tp_indicies[[tp]]))
    }
    cnt <- cnt + 1
    setTxtProgressBar(pb, cnt)
  }
  close(pb)
  clone_summary_table$total_count <- as.numeric(apply(clone_summary_table[,time_points], 1, sum))
  
  ### end time
  end_time <- Sys.time()
  
  ### print out the running time
  cat(paste("Running Time:",
            signif(as.numeric(difftime(end_time, start_time, units = "mins")), digits = 3),
            "mins"))
  
  ### order by the total_count
  clone_summary_table <- clone_summary_table[order(-clone_summary_table$total_count),]
  
  ### save the table as Excel file
  write.xlsx2(clone_summary_table, file = paste0(outputDir, "All_Clones_Count_Summary.xlsx"),
              sheetName = "ALL_CLONES_SUMMARY", row.names = FALSE)
  
  
  ### Alluvial plot with the top 30 clones
  
  ### select the top 30 clones
  lineage_table <- clone_summary_table[1:30,]
  
  ### get an input data frame for the alluvial plot
  time_points <- c("d0", "d5", "d12", "d28", "d60", "d90", "d120", "d180")
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
    ggtitle("Clonal Tracing of the Top 30 Clones from All the Cells") +
    geom_flow() +
    geom_stratum(alpha = 1) +
    # geom_text(stat = "stratum", size = 0.8) +
    rotate_x_text(90) +
    theme_pubr(legend = "right") +
    theme(axis.title.x = element_blank()) +
    theme_cleveland2() +
    scale_fill_viridis(discrete = T) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, NA))
  ggsave(file = paste0(outputDir, "Top_30_All_Cell_Clonal_Tracing.png"), width = 20, height = 10, dpi = 300)
  
  
  #
  ### PCA & UMAP with the top 9 clones from the "All cell" result
  #
  
  ### subset preparation for the PCA & UMAP
  Idents(object = Seurat_Obj) <- Seurat_Obj@meta.data$clone_id
  subset_Seurat_Obj <- subset(Seurat_Obj, idents=lineage_table$clone_id[1:9])
  subset_Seurat_Obj@meta.data$Day <- factor(subset_Seurat_Obj@meta.data$Day,
                                            levels = c("d0", "d5", "d12", "d28", "d60", "d90", "d120", "d180"))

  ### draw the PCA & UMAP
  DimPlot(subset_Seurat_Obj, reduction = "pca", split.by = "clone_id", group.by = "Day",
          ncol = 3, pt.size = 2) +
    labs(title = "PCA of the Top 9 Clones")
  ggsave(file = paste0(outputDir, "PCA_Top_9_Clones_All_Cell.png"), width = 20, height = 10, dpi = 300)
  DimPlot(subset_Seurat_Obj, reduction = "umap", split.by = "clone_id", group.by = "Day",
          ncol = 3, pt.size = 2) +
    labs(title = "UMAP of the Top 9 Clones")
  ggsave(file = paste0(outputDir, "UMAP_Top_9_Clones_All_Cell.png"), width = 20, height = 10, dpi = 300)
  
  
  
  # ### Adding cell types to the meta.data
  # Seurat_Obj@meta.data$CD4_CD8 <- NA
  # Seurat_Obj@meta.data$CD4_CD8[which(Seurat_Obj@meta.data$seurat_clusters %in% c(0,1,3,4,5,8,9,10,13,14,15,17))] <- "CD4"
  # Seurat_Obj@meta.data$CD4_CD8[which(Seurat_Obj@meta.data$seurat_clusters %in% c(2,6,7,11,12,16,18))] <- "CD8"
  # Seurat_Obj@meta.data$Cell_Type <- NA
  # Seurat_Obj@meta.data$Cell_Type[which(Seurat_Obj@meta.data$seurat_clusters %in% c(0,1,2,3,9,10,11,14,16))] <- "Naive"
  # Seurat_Obj@meta.data$Cell_Type[which(Seurat_Obj@meta.data$seurat_clusters %in% c(4,5,7,12,13,15))] <- "Eff-Mem"
  # Seurat_Obj@meta.data$Cell_Type[which(Seurat_Obj@meta.data$seurat_clusters %in% c(6))] <- "MAIT-NKT"
  # Seurat_Obj@meta.data$Cell_Type[which(Seurat_Obj@meta.data$seurat_clusters %in% c(18))] <- "Hobits"
  # Seurat_Obj@meta.data$Cell_Type[which(Seurat_Obj@meta.data$seurat_clusters %in% c(8))] <- "Treg"
  # Seurat_Obj@meta.data$Cell_Type[which(Seurat_Obj@meta.data$seurat_clusters %in% c(17))] <- "TFH"
  # 
  # ### draw the UMAP
  # p <- vector("list", length = 9)
  # names(p) <- lineage_table$clone_id[1:length(p)]
  # for(i in 1:length(p)) {
  #   clones <- Seurat_Obj@meta.data$clone_id
  #   clones[which(Seurat_Obj@meta.data$clone_id != lineage_table$clone_id[i])] <- NA
  #   days <- Seurat_Obj@meta.data$Day
  #   days[which(!Seurat_Obj@meta.data$clone_id != lineage_table$clone_id[i])] <- NA
  #   plot_df <- data.frame(X=Seurat_Obj@reductions$umap@cell.embeddings[,"UMAP_1"],
  #                         Y=Seurat_Obj@reductions$umap@cell.embeddings[,"UMAP_2"],
  #                         Clone=clones,
  #                         Day=days,
  #                         Cell_Type=Seurat_Obj@meta.data$Cell_Type,
  #                         stringsAsFactors = FALSE, check.names = FALSE)
  #   p[[i]] <- ggplot(plot_df, aes_string(x="X", y="Y")) +
  #     geom_point(aes_string(col="Day", shape="Cell_Type"), size=2, alpha=0.8) +
  #     xlab("UMAP_1") + ylab("UMAP_2") +
  #     ggtitle(paste0(lineage_table$clone_id[i])) +
  #     theme_classic(base_size = 16) +
  #     theme(legend.position = "right", plot.title = element_text(hjust = 0.5))
  # }
  
}
