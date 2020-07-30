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
#                                          clone_count_path="./results/Cluster17-TFH/Clone_Count_Summary.xlsx",
#                                          outputDir="./results/")
###

tfh_additional_analyses2 <- function(Seurat_RObj_path="./data/Ali_Tcell_combined_NEW.RDATA",
                                     clone_count_path="./results/Cluster17-TFH/Clone_Count_Summary.xlsx",
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
  if(!require(RColorBrewer, quietly = TRUE)) {
    install.packages("RColorBrewer")
    require(RColorBrewer, quietly = TRUE)
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
  if(!require(gplots, quietly = TRUE)) {
    install.packages("gplots")
    library(gplots, quietly = TRUE)
  }
  if(!require(slingshot, quietly = TRUE)) {
    if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
    BiocManager::install("slingshot")
    require(slingshot, quietly = TRUE)
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
  
  ###
  #   This function was downloaded from: https://raw.githubusercontent.com/obigriffith/biostar-tutorials/master/Heatmaps/heatmap.3.R
  ###
  heatmap.3 <- function(x,
                        Rowv = TRUE, Colv = if (symm) "Rowv" else TRUE,
                        distfun = dist,
                        hclustfun = hclust,
                        dendrogram = c("both","row", "column", "none"),
                        symm = FALSE,
                        scale = c("none","row", "column"),
                        na.rm = TRUE,
                        revC = identical(Colv,"Rowv"),
                        add.expr,
                        breaks,
                        symbreaks = max(x < 0, na.rm = TRUE) || scale != "none",
                        col = "heat.colors",
                        colsep,
                        rowsep,
                        sepcolor = "white",
                        sepwidth = c(0.05, 0.05),
                        cellnote,
                        notecex = 1,
                        notecol = "cyan",
                        na.color = par("bg"),
                        trace = c("none", "column","row", "both"),
                        tracecol = "cyan",
                        hline = median(breaks),
                        vline = median(breaks),
                        linecol = tracecol,
                        margins = c(5,5),
                        ColSideColors,
                        RowSideColors,
                        side.height.fraction=0.3,
                        cexRow = 0.2 + 1/log10(nr),
                        cexCol = 0.2 + 1/log10(nc),
                        labRow = NULL,
                        labCol = NULL,
                        key = TRUE,
                        keysize = 1.5,
                        density.info = c("none", "histogram", "density"),
                        denscol = tracecol,
                        symkey = max(x < 0, na.rm = TRUE) || symbreaks,
                        densadj = 0.25,
                        main = NULL,
                        xlab = NULL,
                        ylab = NULL,
                        lmat = NULL,
                        lhei = NULL,
                        lwid = NULL,
                        ColSideColorsSize = 1,
                        RowSideColorsSize = 1,
                        KeyValueName="Value",...){
    
    invalid <- function (x) {
      if (missing(x) || is.null(x) || length(x) == 0)
        return(TRUE)
      if (is.list(x))
        return(all(sapply(x, invalid)))
      else if (is.vector(x))
        return(all(is.na(x)))
      else return(FALSE)
    }
    
    x <- as.matrix(x)
    scale01 <- function(x, low = min(x), high = max(x)) {
      x <- (x - low)/(high - low)
      x
    }
    retval <- list()
    scale <- if (symm && missing(scale))
      "none"
    else match.arg(scale)
    dendrogram <- match.arg(dendrogram)
    trace <- match.arg(trace)
    density.info <- match.arg(density.info)
    if (length(col) == 1 && is.character(col))
      col <- get(col, mode = "function")
    if (!missing(breaks) && (scale != "none"))
      warning("Using scale=\"row\" or scale=\"column\" when breaks are",
              "specified can produce unpredictable results.", "Please consider using only one or the other.")
    if (is.null(Rowv) || is.na(Rowv))
      Rowv <- FALSE
    if (is.null(Colv) || is.na(Colv))
      Colv <- FALSE
    else if (Colv == "Rowv" && !isTRUE(Rowv))
      Colv <- FALSE
    if (length(di <- dim(x)) != 2 || !is.numeric(x))
      stop("`x' must be a numeric matrix")
    nr <- di[1]
    nc <- di[2]
    if (nr <= 1 || nc <= 1)
      stop("`x' must have at least 2 rows and 2 columns")
    if (!is.numeric(margins) || length(margins) != 2)
      stop("`margins' must be a numeric vector of length 2")
    if (missing(cellnote))
      cellnote <- matrix("", ncol = ncol(x), nrow = nrow(x))
    if (!inherits(Rowv, "dendrogram")) {
      if (((!isTRUE(Rowv)) || (is.null(Rowv))) && (dendrogram %in%
                                                   c("both", "row"))) {
        if (is.logical(Colv) && (Colv))
          dendrogram <- "column"
        else dedrogram <- "none"
        warning("Discrepancy: Rowv is FALSE, while dendrogram is `",
                dendrogram, "'. Omitting row dendogram.")
      }
    }
    if (!inherits(Colv, "dendrogram")) {
      if (((!isTRUE(Colv)) || (is.null(Colv))) && (dendrogram %in%
                                                   c("both", "column"))) {
        if (is.logical(Rowv) && (Rowv))
          dendrogram <- "row"
        else dendrogram <- "none"
        warning("Discrepancy: Colv is FALSE, while dendrogram is `",
                dendrogram, "'. Omitting column dendogram.")
      }
    }
    if (inherits(Rowv, "dendrogram")) {
      ddr <- Rowv
      rowInd <- order.dendrogram(ddr)
    }
    else if (is.integer(Rowv)) {
      hcr <- hclustfun(distfun(x))
      ddr <- as.dendrogram(hcr)
      ddr <- reorder(ddr, Rowv)
      rowInd <- order.dendrogram(ddr)
      if (nr != length(rowInd))
        stop("row dendrogram ordering gave index of wrong length")
    }
    else if (isTRUE(Rowv)) {
      Rowv <- rowMeans(x, na.rm = na.rm)
      hcr <- hclustfun(distfun(x))
      ddr <- as.dendrogram(hcr)
      ddr <- reorder(ddr, Rowv)
      rowInd <- order.dendrogram(ddr)
      if (nr != length(rowInd))
        stop("row dendrogram ordering gave index of wrong length")
    }
    else {
      rowInd <- nr:1
    }
    if (inherits(Colv, "dendrogram")) {
      ddc <- Colv
      colInd <- order.dendrogram(ddc)
    }
    else if (identical(Colv, "Rowv")) {
      if (nr != nc)
        stop("Colv = \"Rowv\" but nrow(x) != ncol(x)")
      if (exists("ddr")) {
        ddc <- ddr
        colInd <- order.dendrogram(ddc)
      }
      else colInd <- rowInd
    }
    else if (is.integer(Colv)) {
      hcc <- hclustfun(distfun(if (symm)
        x
        else t(x)))
      ddc <- as.dendrogram(hcc)
      ddc <- reorder(ddc, Colv)
      colInd <- order.dendrogram(ddc)
      if (nc != length(colInd))
        stop("column dendrogram ordering gave index of wrong length")
    }
    else if (isTRUE(Colv)) {
      Colv <- colMeans(x, na.rm = na.rm)
      hcc <- hclustfun(distfun(if (symm)
        x
        else t(x)))
      ddc <- as.dendrogram(hcc)
      ddc <- reorder(ddc, Colv)
      colInd <- order.dendrogram(ddc)
      if (nc != length(colInd))
        stop("column dendrogram ordering gave index of wrong length")
    }
    else {
      colInd <- 1:nc
    }
    retval$rowInd <- rowInd
    retval$colInd <- colInd
    retval$call <- match.call()
    x <- x[rowInd, colInd]
    x.unscaled <- x
    cellnote <- cellnote[rowInd, colInd]
    if (is.null(labRow))
      labRow <- if (is.null(rownames(x)))
        (1:nr)[rowInd]
    else rownames(x)
    else labRow <- labRow[rowInd]
    if (is.null(labCol))
      labCol <- if (is.null(colnames(x)))
        (1:nc)[colInd]
    else colnames(x)
    else labCol <- labCol[colInd]
    if (scale == "row") {
      retval$rowMeans <- rm <- rowMeans(x, na.rm = na.rm)
      x <- sweep(x, 1, rm)
      retval$rowSDs <- sx <- apply(x, 1, sd, na.rm = na.rm)
      x <- sweep(x, 1, sx, "/")
    }
    else if (scale == "column") {
      retval$colMeans <- rm <- colMeans(x, na.rm = na.rm)
      x <- sweep(x, 2, rm)
      retval$colSDs <- sx <- apply(x, 2, sd, na.rm = na.rm)
      x <- sweep(x, 2, sx, "/")
    }
    if (missing(breaks) || is.null(breaks) || length(breaks) < 1) {
      if (missing(col) || is.function(col))
        breaks <- 16
      else breaks <- length(col) + 1
    }
    if (length(breaks) == 1) {
      if (!symbreaks)
        breaks <- seq(min(x, na.rm = na.rm), max(x, na.rm = na.rm),
                      length = breaks)
      else {
        extreme <- max(abs(x), na.rm = TRUE)
        breaks <- seq(-extreme, extreme, length = breaks)
      }
    }
    nbr <- length(breaks)
    ncol <- length(breaks) - 1
    if (class(col) == "function")
      col <- col(ncol)
    min.breaks <- min(breaks)
    max.breaks <- max(breaks)
    x[x < min.breaks] <- min.breaks
    x[x > max.breaks] <- max.breaks
    if (missing(lhei) || is.null(lhei))
      lhei <- c(keysize, 4)
    if (missing(lwid) || is.null(lwid))
      lwid <- c(keysize, 4)
    if (missing(lmat) || is.null(lmat)) {
      lmat <- rbind(4:3, 2:1)
      
      if (!missing(ColSideColors)) {
        #if (!is.matrix(ColSideColors))
        #stop("'ColSideColors' must be a matrix")
        if (!is.character(ColSideColors) || nrow(ColSideColors) != nc)
          stop("'ColSideColors' must be a matrix of nrow(x) rows")
        lmat <- rbind(lmat[1, ] + 1, c(NA, 1), lmat[2, ] + 1)
        #lhei <- c(lhei[1], 0.2, lhei[2])
        lhei=c(lhei[1], side.height.fraction*ColSideColorsSize/2, lhei[2])
      }
      
      if (!missing(RowSideColors)) {
        #if (!is.matrix(RowSideColors))
        #stop("'RowSideColors' must be a matrix")
        if (!is.character(RowSideColors) || ncol(RowSideColors) != nr)
          stop("'RowSideColors' must be a matrix of ncol(x) columns")
        lmat <- cbind(lmat[, 1] + 1, c(rep(NA, nrow(lmat) - 1), 1), lmat[,2] + 1)
        #lwid <- c(lwid[1], 0.2, lwid[2])
        lwid <- c(lwid[1], side.height.fraction*RowSideColorsSize/2, lwid[2])
      }
      lmat[is.na(lmat)] <- 0
    }
    
    if (length(lhei) != nrow(lmat))
      stop("lhei must have length = nrow(lmat) = ", nrow(lmat))
    if (length(lwid) != ncol(lmat))
      stop("lwid must have length = ncol(lmat) =", ncol(lmat))
    op <- par(no.readonly = TRUE)
    on.exit(par(op))
    
    layout(lmat, widths = lwid, heights = lhei, respect = FALSE)
    
    if (!missing(RowSideColors)) {
      if (!is.matrix(RowSideColors)){
        par(mar = c(margins[1], 0, 0, 0.5))
        image(rbind(1:nr), col = RowSideColors[rowInd], axes = FALSE)
      } else {
        par(mar = c(margins[1], 0, 0, 0.5))
        rsc = t(RowSideColors[,rowInd, drop=F])
        rsc.colors = matrix()
        rsc.names = names(table(rsc))
        rsc.i = 1
        for (rsc.name in rsc.names) {
          rsc.colors[rsc.i] = rsc.name
          rsc[rsc == rsc.name] = rsc.i
          rsc.i = rsc.i + 1
        }
        rsc = matrix(as.numeric(rsc), nrow = dim(rsc)[1])
        image(t(rsc), col = as.vector(rsc.colors), axes = FALSE)
        if (length(rownames(RowSideColors)) > 0) {
          axis(1, 0:(dim(rsc)[2] - 1)/max(1,(dim(rsc)[2] - 1)), rownames(RowSideColors), las = 2, tick = FALSE)
        }
      }
    }
    
    if (!missing(ColSideColors)) {
      
      if (!is.matrix(ColSideColors)){
        par(mar = c(0.5, 0, 0, margins[2]))
        image(cbind(1:nc), col = ColSideColors[colInd], axes = FALSE)
      } else {
        par(mar = c(0.5, 0, 0, margins[2]))
        csc = ColSideColors[colInd, , drop=F]
        csc.colors = matrix()
        csc.names = names(table(csc))
        csc.i = 1
        for (csc.name in csc.names) {
          csc.colors[csc.i] = csc.name
          csc[csc == csc.name] = csc.i
          csc.i = csc.i + 1
        }
        csc = matrix(as.numeric(csc), nrow = dim(csc)[1])
        image(csc, col = as.vector(csc.colors), axes = FALSE)
        if (length(colnames(ColSideColors)) > 0) {
          axis(2, 0:(dim(csc)[2] - 1)/max(1,(dim(csc)[2] - 1)), colnames(ColSideColors), las = 2, tick = FALSE)
        }
      }
    }
    
    par(mar = c(margins[1], 0, 0, margins[2]))
    x <- t(x)
    cellnote <- t(cellnote)
    if (revC) {
      iy <- nr:1
      if (exists("ddr"))
        ddr <- rev(ddr)
      x <- x[, iy]
      cellnote <- cellnote[, iy]
    }
    else iy <- 1:nr
    image(1:nc, 1:nr, x, xlim = 0.5 + c(0, nc), ylim = 0.5 + c(0, nr), axes = FALSE, xlab = "", ylab = "", col = col, breaks = breaks, ...)
    retval$carpet <- x
    if (exists("ddr"))
      retval$rowDendrogram <- ddr
    if (exists("ddc"))
      retval$colDendrogram <- ddc
    retval$breaks <- breaks
    retval$col <- col
    if (!invalid(na.color) & any(is.na(x))) { # load library(gplots)
      mmat <- ifelse(is.na(x), 1, NA)
      image(1:nc, 1:nr, mmat, axes = FALSE, xlab = "", ylab = "",
            col = na.color, add = TRUE)
    }
    axis(1, 1:nc, labels = labCol, las = 2, line = -0.5, tick = 0,
         cex.axis = cexCol)
    if (!is.null(xlab))
      mtext(xlab, side = 1, line = margins[1] - 1.25)
    axis(4, iy, labels = labRow, las = 2, line = -0.5, tick = 0,
         cex.axis = cexRow)
    if (!is.null(ylab))
      mtext(ylab, side = 4, line = margins[2] - 1.25)
    if (!missing(add.expr))
      eval(substitute(add.expr))
    if (!missing(colsep))
      for (csep in colsep) rect(xleft = csep + 0.5, ybottom = rep(0, length(csep)), xright = csep + 0.5 + sepwidth[1], ytop = rep(ncol(x) + 1, csep), lty = 1, lwd = 1, col = sepcolor, border = sepcolor)
    if (!missing(rowsep))
      for (rsep in rowsep) rect(xleft = 0, ybottom = (ncol(x) + 1 - rsep) - 0.5, xright = nrow(x) + 1, ytop = (ncol(x) + 1 - rsep) - 0.5 - sepwidth[2], lty = 1, lwd = 1, col = sepcolor, border = sepcolor)
    min.scale <- min(breaks)
    max.scale <- max(breaks)
    x.scaled <- scale01(t(x), min.scale, max.scale)
    if (trace %in% c("both", "column")) {
      retval$vline <- vline
      vline.vals <- scale01(vline, min.scale, max.scale)
      for (i in colInd) {
        if (!is.null(vline)) {
          abline(v = i - 0.5 + vline.vals, col = linecol,
                 lty = 2)
        }
        xv <- rep(i, nrow(x.scaled)) + x.scaled[, i] - 0.5
        xv <- c(xv[1], xv)
        yv <- 1:length(xv) - 0.5
        lines(x = xv, y = yv, lwd = 1, col = tracecol, type = "s")
      }
    }
    if (trace %in% c("both", "row")) {
      retval$hline <- hline
      hline.vals <- scale01(hline, min.scale, max.scale)
      for (i in rowInd) {
        if (!is.null(hline)) {
          abline(h = i + hline, col = linecol, lty = 2)
        }
        yv <- rep(i, ncol(x.scaled)) + x.scaled[i, ] - 0.5
        yv <- rev(c(yv[1], yv))
        xv <- length(yv):1 - 0.5
        lines(x = xv, y = yv, lwd = 1, col = tracecol, type = "s")
      }
    }
    if (!missing(cellnote))
      text(x = c(row(cellnote)), y = c(col(cellnote)), labels = c(cellnote),
           col = notecol, cex = notecex)
    par(mar = c(margins[1], 0, 0, 0))
    if (dendrogram %in% c("both", "row")) {
      plot(ddr, horiz = TRUE, axes = FALSE, yaxs = "i", leaflab = "none")
    }
    else plot.new()
    par(mar = c(0, 0, if (!is.null(main)) 5 else 0, margins[2]))
    if (dendrogram %in% c("both", "column")) {
      plot(ddc, axes = FALSE, xaxs = "i", leaflab = "none")
    }
    else plot.new()
    if (!is.null(main))
      title(main, cex.main = 1.5 * op[["cex.main"]])
    if (key) {
      par(mar = c(5, 4, 2, 1), cex = 0.75)
      tmpbreaks <- breaks
      if (symkey) {
        max.raw <- max(abs(c(x, breaks)), na.rm = TRUE)
        min.raw <- -max.raw
        tmpbreaks[1] <- -max(abs(x), na.rm = TRUE)
        tmpbreaks[length(tmpbreaks)] <- max(abs(x), na.rm = TRUE)
      }
      else {
        min.raw <- min(x, na.rm = TRUE)
        max.raw <- max(x, na.rm = TRUE)
      }
      
      z <- seq(min.raw, max.raw, length = length(col))
      image(z = matrix(z, ncol = 1), col = col, breaks = tmpbreaks,
            xaxt = "n", yaxt = "n")
      par(usr = c(0, 1, 0, 1))
      lv <- pretty(breaks)
      xv <- scale01(as.numeric(lv), min.raw, max.raw)
      axis(1, at = xv, labels = lv)
      if (scale == "row")
        mtext(side = 1, "Row Z-Score", line = 2)
      else if (scale == "column")
        mtext(side = 1, "Column Z-Score", line = 2)
      else mtext(side = 1, KeyValueName, line = 2)
      if (density.info == "density") {
        dens <- density(x, adjust = densadj, na.rm = TRUE)
        omit <- dens$x < min(breaks) | dens$x > max(breaks)
        dens$x <- dens$x[-omit]
        dens$y <- dens$y[-omit]
        dens$x <- scale01(dens$x, min.raw, max.raw)
        lines(dens$x, dens$y/max(dens$y) * 0.95, col = denscol,
              lwd = 1)
        axis(2, at = pretty(dens$y)/max(dens$y) * 0.95, pretty(dens$y))
        title("Color Key\nand Density Plot")
        par(cex = 0.5)
        mtext(side = 2, "Density", line = 2)
      }
      else if (density.info == "histogram") {
        h <- hist(x, plot = FALSE, breaks = breaks)
        hx <- scale01(breaks, min.raw, max.raw)
        hy <- c(h$counts, h$counts[length(h$counts)])
        lines(hx, hy/max(hy) * 0.95, lwd = 1, type = "s",
              col = denscol)
        axis(2, at = pretty(hy)/max(hy) * 0.95, pretty(hy))
        title("Color Key\nand Histogram")
        par(cex = 0.5)
        mtext(side = 2, "Count", line = 2)
      }
      else title("Color Key")
    }
    else plot.new()
    retval$colorTable <- data.frame(low = retval$breaks[-length(retval$breaks)],
                                    high = retval$breaks[-1], color = retval$col)
    invisible(retval)
  }
  
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
  
  ### order the meta data by time
  subset_Seurat_Obj@meta.data <- subset_Seurat_Obj@meta.data[order(subset_Seurat_Obj@meta.data$Day),]
  
  ### rownames in the meta.data should be in the same order as colnames in the counts
  Seurat_Obj@assays$RNA@counts <- Seurat_Obj@assays$RNA@counts[,rownames(subset_Seurat_Obj@meta.data)]
  
  ### run PCA
  subset_Seurat_Obj <- RunPCA(subset_Seurat_Obj, npcs = 10)
  
  ### find feature contributions of the PC1 
  pca_cos2 <- subset_Seurat_Obj@reductions$pca@feature.loadings * subset_Seurat_Obj@reductions$pca@feature.loadings
  pca_contb <- pca_cos2
  for(i in 1:ncol(pca_contb)) {
    s <- sum(pca_cos2[,i])
    for(j in 1:nrow(pca_contb)) {
      pca_contb[j,i] <- pca_cos2[j,i] * 100 / s
    }
  }
  pca_contb <- pca_contb[order(-pca_contb[,"PC_1"]),]
  
  ### get genes that contributed to the PC1 the most
  contb_threshold <- 1
  important_genes <- rownames(pca_contb)[which(pca_contb[,"PC_1"] > contb_threshold)]
  
  #
  ### a heatmap with the genes - 22 genes x 2127 cells
  #
  
  ### get a matrix for the heatmap
  heatmap_mat <- data.frame(subset_Seurat_Obj@assays$RNA@counts[important_genes,], check.names = FALSE)
  
  ### A function for scaling for heatmap
  scale_h <- function(data, type, na.rm=TRUE) {
    
    if(type == "row") {
      scaled <- t(scale(t(data)))
    } else if(type == "col") {
      scaled <- scale(data)
    } else {
      stop("Type is required: row or col")
    }
    
    if(na.rm == TRUE && (length(which(is.na(scaled))) > 0))  {
      scaled <- scaled[-unique(which(is.na(scaled), arr.ind = TRUE)[,1]),]
    }
    
    return(scaled)
  }
  
  ### scale the data
  heatmap_mat_scaled <- scale_h(heatmap_mat, type = "row")
  
  ### see the distribution of the gene expressions
  # plot(density(heatmap_mat_scaled))
  
  ### because there are some outliers in positive values
  ### we set the maximum as abs(minimum)
  heatmap_mat_scaled[which(heatmap_mat_scaled > abs(min(heatmap_mat_scaled)))] <- abs(min(heatmap_mat_scaled))
  
  ### set colside colors
  uniqueV <- levels(subset_Seurat_Obj@meta.data$Day)
  colors <- colorRampPalette(brewer.pal(9,"Blues"))(length(uniqueV))
  names(colors) <- uniqueV
  
  ### hierarchical clustering functions
  dist.spear <- function(x) as.dist(1-cor(t(x), method = "spearman"))
  hclust.ave <- function(x) hclust(x, method="average")
  
  ### heatmap
  png(paste0(outputDir, "PC1_Genes_Heatmap.png"), width = 2000, height = 1000)
  par(oma=c(0,0,2,6))
  heatmap.3(as.matrix(heatmap_mat_scaled), main = paste0("PC1_Genes_Heatmap_(",
                                                         nrow(heatmap_mat_scaled), " Genes x ",
                                                         ncol(heatmap_mat_scaled), " Cells)"),
            xlab = "", ylab = "", col=greenred(100),
            scale="none", key=T, keysize=0.8, density.info="density",
            dendrogram = "none", trace = "none",
            labRow = rownames(heatmap_mat_scaled), labCol = FALSE,
            Rowv = TRUE, Colv = FALSE,
            distfun=dist.spear, hclustfun=hclust.ave,
            ColSideColors = cbind(colors[as.character(subset_Seurat_Obj@meta.data$Day)]),
            cexRow = 2, cexCol = 2, na.rm = TRUE)
  legend("left", inset = 0, xpd = TRUE, title = "Cell Collection Time", legend = names(colors), fill = colors, cex = 2, box.lty = 0)
  dev.off()
  
  
  #
  ### Diversity & Clonality
  #
  
  ### diversity
  
  ### data preparation
  seurat_diversity <- rep(0, length(levels(subset_Seurat_Obj@meta.data$Day)))
  seurat_cellNum <- rep(0, length(levels(subset_Seurat_Obj@meta.data$Day)))
  names(seurat_diversity) <- levels(subset_Seurat_Obj@meta.data$Day)
  names(seurat_cellNum) <- levels(subset_Seurat_Obj@meta.data$Day)
  for(tp in names(seurat_diversity)) {
    seurat_diversity[tp] <- length(unique(subset_Seurat_Obj@meta.data$clone_id[which(subset_Seurat_Obj@meta.data$Day == tp)]))
    seurat_cellNum[tp] <- length(which(subset_Seurat_Obj@meta.data$Day == tp))
  }
  
  ### data frame for barplot
  plot_df <- data.frame(Value=c(seurat_diversity, seurat_cellNum),
                        Time=c(names(seurat_diversity), names(seurat_cellNum)),
                        Type=c(rep("Clone #", length(seurat_diversity)), rep("Cell #", length(seurat_cellNum))),
                        Pct=c(paste0(signif((seurat_diversity / seurat_cellNum) * 100, digits = 4), "%"),
                              rep("", length(seurat_diversity))),
                        stringsAsFactors = FALSE, check.names = FALSE)
  
  ### factorize the time column
  plot_df$Time <- factor(plot_df$Time, levels = levels(subset_Seurat_Obj@meta.data$Day))
  
  ### draw the plot
  ggplot(data=plot_df, aes(x=Time, y=Value, fill=Type)) +
    geom_bar(stat="identity", position=position_dodge()) +
    geom_text(aes(label=Value), vjust=-1, color="black",
              position = position_dodge(0.9), size=3.5) +
    geom_text(aes(label=Pct, y=75), vjust=0, color="blue", size=3.5) +
    labs(title = "TFH Cluster Clonal Diversity Over Time",
         subtitle = "(Clone # / Cell #) x 100") +
    scale_fill_brewer(palette="Paired") +
    ylim(0, max(plot_df$Value) * 1.1) +
    theme_classic(base_size = 16) +
    theme(axis.title.y = element_blank(),
          plot.subtitle=element_text(size=12, color="blue"))
  ggsave(file = paste0(outputDir, "TFH_Cluster_Clonal_Diversity.png"), width = 12, height = 8, dpi = 300)
  
  
  ### pie plots for the clonality
  
  ### select top n clones
  top_clone_num <- 12
  
  ### load the clone summary table
  clone_summary_table <- read.xlsx2(file = clone_count_path, sheetIndex = 1,
                                    stringsAsFactors=FALSE, check.names=FALSE)
  
  ### set time point values
  time_points <- c("d0", "d5", "d12", "d28", "d60", "d90", "d120", "d180")
  
  ### data frame for the pie plot
  plot_df <- data.frame(freq=as.numeric(as.vector(t(clone_summary_table[1:top_clone_num,time_points]))),
                        clone_id=c(sapply(clone_summary_table$clone_id[1:top_clone_num], function(x) rep(x, length(time_points)))),
                        cdr_ab=c(sapply(clone_summary_table$cdr_ab[1:top_clone_num], function(x) rep(x, length(time_points)))),
                        time=c(rep(time_points, top_clone_num)),
                        stringsAsFactors = FALSE, check.names = FALSE)
  
  ### remove zero frequencies
  plot_df <- plot_df[which(plot_df$freq != 0),]
  
  ### create empty pie chart results
  p <- vector("list", length = length(unique(plot_df$clone_id)))
  names(p) <- unique(plot_df$clone_id)
  
  ### draw pie plots for the top clones
  for(clone in unique(plot_df$clone_id)) {
    ### get the clone indicies
    clone_idx <- which(plot_df$clone_id == clone)
    
    ### calculate percentages
    total_sample_num <- sum(plot_df$freq[clone_idx])
    pct <- sapply(plot_df$freq[clone_idx], function(x) signif(x*100/total_sample_num, digits = 3))
    pct <- paste0(plot_df$freq[clone_idx], "(", pct, "%)")
    
    ### ggplot drawing
    p[[clone]] <- ggplot(data = plot_df[clone_idx,],
           aes(x = "", y = freq, fill = time)) +
      geom_bar(stat = "identity", width = 1) +
      coord_polar(theta="y") +
      geom_text(label = pct,
                position = position_stack(vjust = 0.5), size = 3) +
      labs(x = NULL, y = NULL, title = clone) +
      theme_classic(base_size = 16) +
      theme(plot.title = element_text(hjust = 0.5, color = "black", size = 12),
            axis.line = element_blank(),
            axis.ticks = element_blank(),
            axis.text = element_blank())
  }
  
  ### arrange the plots and print out
  fName <- paste0("TFH_Cluster_Clonality_Pie_Charts")
  g <- arrangeGrob(grobs = p,
                   nrow = 3,
                   ncol = 4,
                   top = fName)
  ggsave(file = paste0(outputDir, fName, ".png"), g, width = 15, height = 9, dpi = 300)
  
  
  ### pseudotime analysis - Slingshot
  
  # DimPlot(subset_Seurat_Obj, reduction = "pca", group.by = "Tissue", pt.size = 2)
  
  ###
  subset_Seurat_Obj <- FindNeighbors(subset_Seurat_Obj, dims = 1:5, k.param = 5)
  subset_Seurat_Obj <- FindClusters(subset_Seurat_Obj, resolution = 0.4)
  
  ###
  slingshot_obj <- slingshot(Embeddings(subset_Seurat_Obj, "pca"), clusterLabels = subset_Seurat_Obj@meta.data$Day, 
                             start.clus = "d0", reducedDim = "PCA")
  
  
  cell_pal <- function(cell_vars, pal_fun,...) {
if (is.numeric(cell_vars)) {
pal <- pal_fun(100, ...)
return(pal[cut(cell_vars, breaks = 100)])
} else {
categories <- sort(unique(cell_vars))
pal <- setNames(pal_fun(length(categories), ...), categories)
return(pal[cell_vars])
}
}

  
  library(scales)
cell_colors_clust <- cell_pal(subset_Seurat_Obj@meta.data$Day, hue_pal())
cell_colors_clust <- cell_pal(levels(subset_Seurat_Obj@meta.data$Day), hue_pal())
plot(reducedDim(slingshot_obj), col = cell_colors_clust, pch = 16, cex = 0.5)
lines(slingshot_obj, lwd = 2, type = 'lineages', col = 'black')

  
  
  
  
  
}
