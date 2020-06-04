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
  
  ###
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
  
  
  
  
  
}
