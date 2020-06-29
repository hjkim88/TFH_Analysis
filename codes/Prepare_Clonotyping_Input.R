###
#   File name : Prepare_Clonotyping_Input.R
#   Author    : Hyunjin Kim
#   Date      : Jun 28, 2020
#   Email     : hyunjin.kim@stjude.org
#   Purpose   : We run TCRdist for clonotyping and it needs meta data info
#               First column: CellRanger VDJ file path
#               Second column: Barcode suffix
#
#   Instruction
#               1. Source("Prepare_Clonotyping_Input.R")
#               2. Run the function "prepare_input" - specify the input file paths and the output directory
#               3. The results will be generated under the output directory
#
#   Example
#               > source("The_directory_of_Prepare_Clonotyping_Input.R/Prepare_Clonotyping_Input.R")
#               > prepare_input(Seurat_RObj_path="./data/JCC243_JCC280_Aggregregress.Robj",
#                               TCR_data_dirs="./data/JCC280_VDJoutputs/",
#                               outputDir="./data/")
###

prepare_input <- function(Seurat_RObj_path="./data/JCC243_JCC280_Aggregregress.Robj",
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
  
  ### get Library info from the TCR files
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
  
  ### for each file extract the info
  result_mat <- matrix("", 1, 2)
  colnames(result_mat) <- c("file", "suffix")
  for(i in 1:length(TCR_data_dirs)) {
    
    ### write progress
    writeLines(paste(TCR_data_dirs[i]))
    
    ### load filtered contig annotation data
    tcr_data <- read.csv(file = TCR_data_dirs[i],
                         stringsAsFactors = FALSE, check.names = FALSE)
    
    ### check if this is TCR info
    if(length(which(unique(tcr_data$chain) %in% c("TRA", "TRB"))) > 0) {
      ### get file path
      filePath <- paste0(dirname(normalizePath(TCR_data_dirs[i])), "/", basename(TCR_data_dirs[i]))
      
      ### get suffix of the barcodes
      bcs <- rownames(Seurat_Obj@meta.data)[which(Seurat_Obj@meta.data$Library == tcr_libs[basename(TCR_data_dirs[i])])]
      sufs <- sapply(strsplit(bcs, split = "-", fixed = TRUE), function(x) x[2])
      sufs <- unique(sufs)
      
      ### get one suffix
      if(length(sufs) == 1) {
        suf <- sufs
      } else {
        stop("ERROR: suffix number not 1")
      }
      
      ### add to the result matrix
      result_mat <- rbind(result_mat, c(filePath, suf))
    }
    
  }
  
  ### remove the first null row
  result_mat <- result_mat[-1,]
  
  ### write out the input meta data file
  write.csv(result_mat, file = paste0(outputDir, "TFH_metadata_for_tcrdist.csv"), row.names = FALSE)
  
}
