###
#   File name : Combine_GEX_And_TCR.R
#   Author    : Hyunjin Kim
#   Date      : May 27, 2020
#   Email     : hyunjin.kim@stjude.org
#   Purpose   : The 10X dataset contains lymph node and PBMC T cells from one donor over 8 different time points
#               after receiving a Flu vaccine on day 0. This makes it a great dataset to map these
#               clonal trajectories in response to a timed stimulus. We also have the bulk and scTCR sequencing
#               for this donor and two others from the study.
#               Here, the GEX and the TCR data are combined.
#
#   Instruction
#               1. Source("Combine_GEX_And_TCR.R")
#               2. Run the function "combine_gex_tcr" - specify the input file paths and the output directory
#               3. The results will be generated under the output directory
#
#   Example
#               > source("The_directory_of_Combine_GEX_And_TCR.R/Combine_GEX_And_TCR.R")
#               > combine_gex_tcr(Seurat_RObj_path="./data/Ali_Tcell_agg.rds",
#                                 TCR_data_dirs="./data/JCC280_VDJoutputs/",
#                                 outputDir="./results/")
###

combine_gex_tcr <- function(Seurat_RObj_path="./data/JCC243_JCC280_Aggregregress.Robj",
                            TCR_data_dirs="./data/JCC280_VDJoutputs/",
                            outputDir="./data/") {
  
  ### load library
  if(!require(Seurat, quietly = TRUE)) {
    install.packages("Seurat")
    require(Seurat, quietly = TRUE)
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
  
  ### get TCR info file paths
  TCR_data_dirs <- list.files(path = TCR_data_dirs, pattern = "filtered_contig_annotations.csv$",
                              full.names = TRUE, recursive = TRUE)
  
  ### get Library and libNum info from the TCR files
  tcr_libr_lib <- data.frame(matrix("", length(TCR_data_dirs), 2),
                             stringsAsFactors = FALSE, check.names = FALSE)
  rownames(tcr_libr_lib) <- basename(TCR_data_dirs)
  colnames(tcr_libr_lib) <- c("Library", "LibNum")
  for(i in 1:length(TCR_data_dirs)) {
    extracted <- basename(TCR_data_dirs[i])
    extracted <- substr(extracted, 15, nchar(extracted)-36)
    extracted <- strsplit(extracted, split = "-lib", fixed = TRUE)[[1]]
    
    ### "B_cell" to "Bcell"
    ### "PBMC_2" TO "PBMC-2"
    if(grepl("B_cell", extracted[1])) {
      extracted[1] <- paste0(strsplit(extracted[1], split = "B_cell", fixed = TRUE)[[1]][1], "Bcell")
    } else if(grepl("PBMC_2", extracted[1])) {
      extracted[1] <- paste0(strsplit(extracted[1], split = "PBMC_2", fixed = TRUE)[[1]][1], "PBMC-2")
    }
    
    tcr_libr_lib[i,"Library"] <- extracted[1]
    tcr_libr_lib[i,"LibNum"] <- extracted[2]
  }
  
  ### check every TCR info has matched Library & LibNum in the Seurat object
  for(i in 1:nrow(tcr_libr_lib)) {
    if(length(which(Seurat_Obj@meta.data$Library == tcr_libr_lib[i,"Library"])) == 0) {
      writeLines(paste(i, ":", tcr_libr_lib[i,"Library"]))
    }
  }
  
  
  
  
  
  
  
  
  
  ### load the Seurat object and save the object name
  Seurat_Obj <- readRDS(Seurat_RObj_path)
  
  ### load filtered contig annotation data
  tcr_data <- read.table(file = TCR_data_path, sep = "\t", header = TRUE,
                         stringsAsFactors = FALSE, check.names = FALSE)
  rownames(tcr_data) <- tcr_data$barcode
  
  ### check if NT sequences and clone ids are the same between the GEX and TCR
  a <- Seurat_Obj@meta.data[rownames(tcr_data),"match.cdr"]
  b <- tcr_data$match.cdr
  print(identical(a, b))
  a <- Seurat_Obj@meta.data[rownames(tcr_data),"clone_id"]
  b <- tcr_data$clone_id
  print(identical(a, b))
  
  ### attach some columns in the TCR to the GEX
  added_cols <- c("MHCi_score_cdr3a",
                  "MHCi_score_cdr3b",
                  "knn_MHCi_cdr3a",       
                  "knn_MHCi_cdr3b",
                  "matched_SC",
                  "matched_SC_alpha_only",
                  "matched_SC_beta_only", 
                  "matched_bulk_alpha",
                  "matched_bulk_beta")
  Seurat_Obj@meta.data[added_cols] <- NA
  Seurat_Obj@meta.data[rownames(tcr_data),added_cols] <- tcr_data[,added_cols]
  
  ### save the new combined RDATA file
  save(list = c("Seurat_Obj"), file = paste0(outputDir, "Ali_Tcell_combined.RDATA"))
  
}
