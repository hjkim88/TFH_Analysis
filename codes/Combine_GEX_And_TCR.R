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
  tcr_libs <- rep("", length(TCR_data_dirs))
  names(tcr_libs) <- basename(TCR_data_dirs)
  for(i in 1:length(TCR_data_dirs)) {
    extracted <- basename(TCR_data_dirs[i])
    extracted <- substr(extracted, 15, nchar(extracted)-41)
    
    ### "B_cell" to "Bcell"
    ### "PBMC_2" TO "PBMC-2"
    if(grepl("B_cell", extracted)) {
      extracted <- paste0(strsplit(extracted, split = "B_cell", fixed = TRUE)[[1]][1], "Bcell")
    } else if(grepl("PBMC_2", extracted)) {
      extracted <- paste0(strsplit(extracted, split = "PBMC_2", fixed = TRUE)[[1]][1], "PBMC-2")
    }
    
    tcr_libs[i] <- extracted
  }
  
  ### check if every TCR info has matched Library in the Seurat object
  keep_idx <- NULL
  for(i in 1:length(tcr_libs)) {
    if(length(which(Seurat_Obj@meta.data$Library == tcr_libs[i])) == 0) {
      writeLines(paste(i, ":", tcr_libs[i]))
    } else {
      keep_idx <- c(keep_idx, i)
    }
  }
  
  ### retain the TCR files that are already in the Seurat object 
  TCR_data_dirs <- TCR_data_dirs[keep_idx]
  
  ### remove BCR files - only include TCR files
  TCR_data_dirs <- TCR_data_dirs[grep("TCR", TCR_data_dirs)]
  
  ### combine the TCR info
  tcr <- NULL
  for(i in 1:length(TCR_data_dirs)) {
    
    ### write progress
    writeLines(paste(TCR_data_dirs[i]))
    
    ### load filtered contig annotation data
    tcr_data <- read.csv(file = TCR_data_dirs[i],
                         stringsAsFactors = FALSE, check.names = FALSE)
    
    ### remove the -* at the end of each barcode.
    tcr_data$barcode <- sapply(strsplit(tcr_data$barcode, split = "-", fixed = TRUE), function(x) x[1])
    
    ### remove the rows that do not have CDR3 sequence
    tcr_data <- tcr_data[which(tcr_data$cdr3 != "None"),]
    
    ### remove redundant rows
    tcr_data <- tcr_data[!duplicated(tcr_data[c("barcode", "cdr3")]),]
    
    ### order by "chain" so that TRA rows come first than TRB rows
    ### and secondly, order by "CDR3 Nucleotide" sequence
    tcr_data <- tcr_data[order(as.character(tcr_data$cdr3_nt)),]
    
    ### merge different TRA & TRB info to one row
    dups <- which(duplicated(tcr_data$barcode))
    if(length(dups) > 0) {
      temp <- tcr_data[dups,]
      tcr_data <- tcr_data[-dups,]
      rownames(tcr_data) <- tcr_data$barcode
      for(barcode in tcr_data$barcode) {
        idx <- which(temp$barcode == barcode)
        tcr_data[barcode,"cdr3"] <- paste(c(tcr_data[barcode,"cdr3"], temp[idx, "cdr3"]), collapse = ";")
        tcr_data[barcode,"cdr3_nt"] <- paste(c(tcr_data[barcode,"cdr3_nt"], temp[idx, "cdr3_nt"]), collapse = ";")
        tcr_data[barcode,"productive"] <- paste(c(tcr_data[barcode,"productive"], temp[idx, "productive"]), collapse = ";")
        tcr_data[barcode,"reads"] <- paste(c(tcr_data[barcode,"reads"], temp[idx, "reads"]), collapse = ";")
        tcr_data[barcode,"umis"] <- paste(c(tcr_data[barcode,"umis"], temp[idx, "umis"]), collapse = ";")
      }
    }
    
    ### only retain informative columns
    tcr_data <- tcr_data[,c("barcode", "raw_clonotype_id", "cdr3", "cdr3_nt", "reads", "umis", "productive")]
    colnames(tcr_data) <- c("barcode", "raw_clonotype_id", "cdr3_aa", "cdr3_nt", "tcr_reads", "tcr_umis", "tcr_productive")
    
    ### add library info
    tcr_data$library <- tcr_libs[basename(TCR_data_dirs[i])]
    
    ### combine the TCR data for all the data
    if(is.null(tcr)) {
      tcr <- tcr_data
    } else {
      tcr <- rbind(tcr, tcr_data)
    }
    
  }
  
  ### attach the TCR info to the metadata of the Seurat object
  obj_colnames <- colnames(Seurat_Obj@meta.data)
  tcr_colnames <- colnames(tcr)
  Seurat_Obj@meta.data$GexCellFull <- rownames(Seurat_Obj@meta.data)
  Seurat_Obj@meta.data$GexCellShort <-  sapply(strsplit(rownames(Seurat_Obj@meta.data), split = "-", fixed = TRUE), function(x) x[1])
  Seurat_Obj@meta.data <- merge(Seurat_Obj@meta.data,
                                tcr,
                                by.x = c("GexCellShort", "Library"),
                                by.y = c("barcode", "library"),
                                all.x = TRUE)
  Seurat_Obj@meta.data <- Seurat_Obj@meta.data[,c("GexCellFull", "GexCellShort", obj_colnames, setdiff(tcr_colnames, c("barcode", "library")))]
  rownames(Seurat_Obj@meta.data) <- Seurat_Obj@meta.data[,"GexCellFull"]
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
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
