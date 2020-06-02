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
  
  
  
  ### ggalluvial example
  data(vaccinations)
  levels(vaccinations$response) <- rev(levels(vaccinations$response))
  ggplot(vaccinations,
         aes(x = survey, stratum = response, alluvium = subject,
             y = freq,
             fill = response, label = response)) +
    scale_x_discrete(expand = c(.1, .1)) +
    geom_flow() +
    geom_stratum(alpha = .5) +
    geom_text(stat = "stratum", size = 3) +
    theme(legend.position = "none") +
    ggtitle("vaccination survey responses at three points in time")
  
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
  
  ### test
  temp <- data.frame(height=rep(c(1:5), 2),
                     label=c("A", "B", "C")[sample(3, 10, replace = TRUE)],
                     Time=c("1", "5", "10")[sample(3, 10, replace = TRUE)],
                     connection=sample(5, 10, replace = TRUE))
  ggplot(vaccinations,
         aes(x = survey, stratum = response, alluvium = subject,
             y = freq,
             fill = response, label = response)) +
    scale_x_discrete(expand = c(.1, .1)) +
    geom_flow() +
    geom_stratum(alpha = 1) +
    geom_text(stat = "stratum", size = 3) +
    theme_pubr(legend = "right") + rotate_x_text(90) + theme_cleveland2() +
    scale_fill_viridis(discrete = T)
  
  
  
}
