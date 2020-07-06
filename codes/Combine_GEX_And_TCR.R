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
#               Since Stefan's object does not have some TCR info (single chain or multiple chains),
#               I attached those info additonally, but the basic idea of having only one alpha and
#               only one beta for clonal lineage make sense, so I used the clone_id for the lineage tracing.
#               single chain clonoal lineage does not make sense in two-chain data, and even if there are
#               two in-frame alpha chains, if one alpha chain's UMI count is significantly larger, then
#               that alpha and beta pair should be considered in the lineage tracing.
#
#   Instruction
#               1. Source("Combine_GEX_And_TCR.R")
#               2. Run the function "combine_gex_tcr" - specify the input file paths and the output directory
#               3. The results will be generated under the output directory
#
#   Example
#               > source("The_directory_of_Combine_GEX_And_TCR.R/Combine_GEX_And_TCR.R")
#               > combine_gex_tcr(Seurat_RObj_path="./data/JCC243_JCC280_Aggregregress.Robj",
#                                 TCR_data_dirs="./data/JCC280_VDJoutputs/",
#                                 stefan_obj_path="./data/Ali_all_agg_wTCR.rds",
#                                 scTCR_PCR_info_path="./data/scTCR(PCR)/Ali_paired_sc_clones_JCC283.tsv",
#                                 bulkTCR_info_path="./data/bulk/Ali_bulk_allDonor_clones.tsv",
#                                 outputDir="./data/")
###

combine_gex_tcr <- function(Seurat_RObj_path="./data/JCC243_JCC280_Aggregregress.Robj",
                            TCR_data_dirs="./data/JCC280_VDJoutputs/",
                            stefan_obj_path="./data/Ali_all_agg_wTCR.rds",
                            scTCR_PCR_info_path="./data/scTCR(PCR)/Ali_paired_sc_clones_JCC283.tsv",
                            bulkTCR_info_path="./data/bulk/Ali_bulk_allDonor_clones.tsv",
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
  
  ### combine the TCR info
  tcr <- NULL
  for(i in 1:length(TCR_data_dirs)) {
    
    ### write progress
    writeLines(paste(TCR_data_dirs[i]))
    
    ### load filtered contig annotation data
    tcr_data <- read.csv(file = TCR_data_dirs[i],
                         stringsAsFactors = FALSE, check.names = FALSE)
    
    ### check if this is TCR info
    if(length(which(unique(tcr_data$chain) %in% c("TRA", "TRB"))) > 0) {
      ### remove the -* at the end of each barcode.
      tcr_data$barcode <- sapply(strsplit(tcr_data$barcode, split = "-", fixed = TRUE), function(x) x[1])
      
      ### remove the rows that do not have CDR3 sequence
      tcr_data <- tcr_data[which(tcr_data$cdr3 != "None"),]
      
      ### remove the rows with productive == False
      tcr_data <- tcr_data[which(tcr_data$productive == "True"),]
      
      ### remove redundant rows
      tcr_data <- tcr_data[!duplicated(tcr_data[c("barcode", "cdr3")]),]
      
      ### order by "chain" so that TRA rows come first than TRB rows
      ### and secondly, order by "CDR3 Nucleotide" sequence
      tcr_data <- tcr_data[order(as.character(tcr_data$chain), as.character(tcr_data$cdr3_nt)),]
      
      ### annotate TRA & TRB info to the cdr3 sequences
      tcr_data$cdr3 <- paste0(tcr_data$chain, ":", tcr_data$cdr3)
      tcr_data$cdr3_nt <- paste0(tcr_data$chain, ":", tcr_data$cdr3_nt)
      
      ### now merge different TRA & TRB info to one row
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
    
  }
  
  # ### keep the TCR info only if there are both alpha & beta chains
  # is_keep <- sapply(tcr$cdr3_aa, function(x) {
  #   if(grepl("TRA:", x) && grepl("TRB:", x)) {
  #     return (TRUE)
  #   } else {
  #     return (FALSE)
  #   }
  # })
  # tcr <- tcr[which(is_keep),]
  
  
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
  
  ### remove the "productive" column since we only used the TRUE chains
  Seurat_Obj@meta.data <- Seurat_Obj@meta.data[,-which(colnames(Seurat_Obj@meta.data) == "tcr_productive")]
  
  
  ### load Stefan's Seurat object
  stefan_seurat_obj <- readRDS(stefan_obj_path)
  
  ### make the object order same as the main Seurat object
  stefan_seurat_obj@assays$RNA@counts <- stefan_seurat_obj@assays$RNA@counts[,rownames(Seurat_Obj@meta.data)]
  stefan_seurat_obj@meta.data <- stefan_seurat_obj@meta.data[colnames(stefan_seurat_obj@assays$RNA@counts),]
  
  ### get information from the Stefan's object
  new_cols <- setdiff(colnames(stefan_seurat_obj@meta.data), colnames(Seurat_Obj@meta.data))
  Seurat_Obj@meta.data[new_cols] <- stefan_seurat_obj@meta.data[,new_cols]
  
  
  ### load additional info (bulk TCR & scTCR[PCR])
  scTCR_PCR <- read.table(file = scTCR_PCR_info_path, header = TRUE, sep = "\t",
                          stringsAsFactors = FALSE, check.names = FALSE)
  bulk_TCR <- read.table(file = bulkTCR_info_path, header = TRUE, sep = "\t",
                         stringsAsFactors = FALSE, check.names = FALSE)
  
  ### attach scTCR[PCR] info to the Seurat object
  Seurat_Obj@meta.data$matched_SC <- ifelse(Seurat_Obj@meta.data$match.cdr %in% scTCR_PCR$match.cdr, "yes", "no" )
  Seurat_Obj@meta.data$matched_SC_alpha_only <- ifelse(Seurat_Obj@meta.data$cdr3a %in% scTCR_PCR$cdr3a, "yes", "no" )
  Seurat_Obj@meta.data$matched_SC_beta_only <- ifelse(Seurat_Obj@meta.data$cdr3b %in% scTCR_PCR$cdr3b, "yes", "no" )
  
  ### attach bulk TCR info to the Seurat object
  bulk_TCR <- bulk_TCR[which(bulk_TCR$donor == "321-05"),]
  Seurat_Obj@meta.data$matched_bulk_alpha <- ifelse(Seurat_Obj@meta.data$cdr3a %in% bulk_TCR$cdr3aa, "yes", "no" )
  Seurat_Obj@meta.data$matched_bulk_beta <- ifelse(Seurat_Obj@meta.data$cdr3b %in% bulk_TCR$cdr3aa, "yes", "no" )
  
  ### check same match.cdr regardless of library belong to the same clone_id
  dups <- unique(Seurat_Obj@meta.data$clone_id[which(duplicated(Seurat_Obj@meta.data$clone_id))])
  dups <- dups[which(!is.na(dups))]
  for(dup in dups) {
    target_idx <- which(Seurat_Obj@meta.data$clone_id == dup)
    if(length(unique(Seurat_Obj@meta.data$match.cdr[target_idx])) > 1) {
      writeLines(paste("clone_id =", dup))
    }
  }
  
  ### save the new combined RDATA file
  save(list = c("Seurat_Obj"), file = paste0(outputDir, "Ali_Tcell_combined_NEW.RDATA"))
  
}
