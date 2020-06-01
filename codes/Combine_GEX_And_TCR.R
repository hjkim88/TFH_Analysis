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
#                                 TCR_data_dirs="./data/Ali_clones_mapped_TCR_newMay11.txt",
#                                 outputDir="./results/")
###

combine_gex_tcr <- function(Seurat_RObj_path="./data/Ali_Tcell_agg.rds",
                            TCR_data_path="./data/Ali_clones_mapped_TCR_newMay11.txt",
                            outputDir="./data/") {
  
  ### load library
  if(!require(Seurat, quietly = TRUE)) {
    install.packages("Seurat")
    require(Seurat, quietly = TRUE)
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
