#_ -------------------------------------------------------

# fun: AbsSummary CycifStack ----

#' A graphical summary of cycles and antibodies of a CycifStack class.
#'
#' Generates an absolute summary plot for a CycifStack object.
#'
#' @param x A CycifStack object.
#' @param show.cycles.in.row Logical, indicating whether to display cycle information in row labels.
#' FALSE by default. If TRUE, the output plot shows the number of cycles each sample was stained for.
#' @param ... Additional parameters to be passed to the \code{graphics::image} function.
#'
#' @return
#' The function generates a graphical absolute summary plot and returns \code{NULL}.
#'
#' @details
#' This function generates a plot that summarizes the antibodies used in a CycifStack object.
#' It visualizes the antibodies across different samples and cycles.
#'
#' @seealso
#' \code{\link{abs_list}}
#'
#' @importFrom graphics image box abline axis par
#' @importFrom grDevices grey
#'
#' @export
AbsSummary <- function(x,show.cycles.in.row=FALSE,...){
  if(!is(x,"CycifStack")){
    stop("input should be a CycifStack obj.")
  }
  uniq.abs <- as.character(x@abs_list$ab)
  n1 <- data.frame(do.call(rbind,lapply(x@samples,function(y){
    this.abs <- as.character(abs_list(y)$ab)
    tested <- uniq.abs %in% this.abs
    return(tested)
  })))

  idx <- rep(seq(ncol(n1)/3),each=3)
  n2 <- n1 * idx[col(n1)]
  n2[n2==0] <- NA
  pcols <- grDevices::grey(c(.8,.6))[(seq(max(n2,na.rm=T)) %% 2) + 1]

  if(show.cycles.in.row){
    smpl.labs <- paste0(x@names,"\n(",x@n_cycles, " cycles)")
  }else{
    smpl.labs <- x@names
  }
  col.labs <- paste0("Cycle\n",seq(x@max_cycles))
  graphics::image(seq(ncol(n2)),seq(nrow(n2)),t(n2[rev(seq(nrow(n2))),]),col=pcols,
                  xlab="",ylab="",axes=F,...)
  graphics::box()
  graphics::abline(h=c(0,seq(nrow(n1)))+.5,lwd=.5)
  graphics::abline(v=c(0,seq(ncol(n1)))+.5,lwd=.5)
  graphics::axis(2,at=rev(seq(nrow(n1))),labels=smpl.labs,las=1)
  graphics::axis(1,at=seq(ncol(n2)),labels=uniq.abs,las=2)
  graphics::par(tcl=0)
  graphics::axis(3,at=seq(x@max_cycles)*3-1,labels=col.labs)
  graphics::par(tcl=-.5)

  return(invisible(NULL))
}
# fun: h3/heatmap3 matrix, data.frame ----

#' 3D Heatmap Plot (extension of heatmap3)
#'
#' This function creates a 3D heatmap plot to visualize a matrix of data, extending \code{heatmap3::heatmap3}.
#'
#' @param x A numeric matrix containing the data to be plotted.
#' @param Rowv Specifies the row dendrogram or clustering (default is NULL).
#' @param Colv Specifies the column dendrogram or clustering (default is NULL).
#' @param distfun A function to calculate the distance matrix for rows (default is dist).
#' @param distfunC A function to calculate the distance matrix for columns (default is NULL).
#' @param distfunR A function to calculate the distance matrix for rows (default is NULL).
#' @param balanceColor Logical, indicating whether to balance colors (default is FALSE).
#' @param ColSideLabs Labels for the columns' side (default is NULL).
#' @param RowSideLabs Labels for the rows' side (default is NULL).
#' @param showColDendro Logical, indicating whether to show the column dendrogram (default is TRUE).
#' @param showRowDendro Logical, indicating whether to show the row dendrogram (default is TRUE).
#' @param col A color palette for the heatmap (default is a gradient from navy to firebrick3).
#' @param legendfun A custom function to create the legend (default is NULL).
#'
#' @details
#' - The `h3` function creates a 3D heatmap plot to visualize a matrix of data, extending and customizing \code{heatmap3::heatmap3} function.
#' - It allows customization of various aspects of the plot, including colors, labels, and dendrograms.
#' - Additional functionalities like distance matrix calculations and custom legends are available.
#'
#' @export
h3 <- function (x, Rowv = NULL, Colv = if (symm) "Rowv" else NULL,
                # distfun = function(x) as.dist(1 - cor(t(x), use = "pa")),
                distfun = dist, distfunC, distfunR, balanceColor = F, ColSideLabs, RowSideLabs,
                showColDendro = T, showRowDendro = T, col = colorRampPalette(c("navy","white", "firebrick3"))(1024),
                legendfun, method = "complete", ColAxisColors = 0, RowAxisColors = 0, hclustfun = hclust,
                reorderfun = function(d, w) reorder(d, w), add.expr, symm = FALSE,
                revC = identical(Colv, "Rowv"), scale = c("row", "column","none"),
                na.rm = TRUE, ColSideFun, ColSideAnn, ColSideWidth = 0.4,
                ColSideCut, colorCell, highlightCell, file = "heatmap3.pdf",
                topN = NA, filterFun = sd, returnDistMatrix = FALSE, margins = c(5,5),
                ColSideColors, RowSideColors, cexRow = 0.2 + 1/log10(nrow(x)),
                cexCol = 0.2 + 1/log10(ncol(x)), lasRow = 2, lasCol = 2,
                labRow = NULL, labCol = NULL, main = NULL, xlab = NULL, ylab = NULL, na.color = "grey",
                keep.dendro = FALSE, verbose = getOption("verbose"),
                useRaster = if (ncol(x) * nrow(x) >= 50000) TRUE else FALSE, ...)
{
  hcc <- NULL
  hcr <- NULL
  if (!all(is.na(topN))) {
    temp <- apply(x, 1, filterFun)
    pdf(file)
    for (n in topN) {
      xSub <- x[rev(order(temp))[1:n], , drop = F]
      if (!missing(RowSideColors)) {
        RowSideColorsBak <- RowSideColors
        RowSideColors <- RowSideColors[rev(order(temp))[1:n],
                                       , drop = F]
      }
      result[[paste0(n)]] <- heatmap3(xSub, Rowv = Rowv,
                                      Colv = Colv, distfun = distfun, balanceColor = balanceColor,
                                      ColSideLabs = ColSideLabs, RowSideLabs = RowSideLabs,
                                      showColDendro = showColDendro, showRowDendro = showRowDendro,
                                      col = col, legendfun = legendfun, method = "complete",
                                      ColAxisColors = 0, RowAxisColors = 0, hclustfun = hclust,
                                      reorderfun = reorderfun, add.expr = add.expr,
                                      symm = symm, revC = revC, scale = scale, na.rm = na.rm,
                                      ColSideFun = ColSideFun, ColSideAnn = ColSideAnn,
                                      ColSideWidth = ColSideWidth, ColSideCut = ColSideCut,
                                      margins = margins, ColSideColors = ColSideColors,
                                      RowSideColors = RowSideColors, cexRow = cexRow,
                                      cexCol = cexCol, labRow = labRow, labCol = labCol,
                                      main = paste0("top ", n), xlab = xlab, ylab = ylab,
                                      keep.dendro = keep.dendro, verbose = verbose,
                                      ...)
      if (!missing(RowSideColors)) {
        RowSideColors <- RowSideColorsBak
      }
    }
    temp <- dev.off()
    cat(paste0("The heatmaps were generated at ", file, "\n"))
    return(invisible(result))
  }
  scale <- if (symm && missing(scale))
    "none"
  else match.arg(scale)
  if (is.data.frame(x)) {
    x <- as.matrix(x)
  }
  ## a ---

  if (!missing(ColSideColors)) {
    if (is.vector(ColSideColors)) {
      ColSideColors <- cbind(ColSideColors)
    }
  }
  if (!missing(RowSideColors)) {
    if (is.vector(RowSideColors)) {
      RowSideColors <- cbind(RowSideColors)
    }
  }
  if (missing(distfunC)) {
    distfunC <- distfun
  }
  if (missing(distfunR)) {
    distfunR <- distfun
  }
  distMatrixC = NULL
  distMatrixR = NULL
  if (length(di <- dim(x)) != 2 || !is.numeric(x))
    stop("'x' must be a numeric matrix")
  nr <- di[1L]
  nc <- di[2L]
  if (nr <= 1 || nc <= 1)
    stop("'x' must have at least 2 rows and 2 columns")
  if (!is.numeric(margins) || length(margins) != 2L)
    stop("'margins' must be a numeric vector of length 2")
  doRdend <- !identical(Rowv, NA)
  doCdend <- !identical(Colv, NA)
  if (!doRdend && identical(Colv, "Rowv"))
    doCdend <- FALSE
  if (is.null(Rowv))
    Rowv <- rowMeans(x, na.rm = na.rm)
  if (is.null(Colv))
    Colv <- colMeans(x, na.rm = na.rm)
  if (doRdend) {
    if (inherits(Rowv, "dendrogram"))
      ddr <- Rowv
    else {
      distMatrixR = distfunR(x)
      hcr <- hclustfun(distMatrixR, method = method)
      ddr <- as.dendrogram(hcr)
      if (!is.logical(Rowv) || Rowv)
        ddr <- reorderfun(ddr, Rowv)
    }
    if (nr != length(rowInd <- order.dendrogram(ddr)))
      stop("row dendrogram ordering gave index of wrong length")
  }
  else rowInd <- 1L:nr
  if (doCdend) {
    if (inherits(Colv, "dendrogram"))
      ddc <- Colv
    else if (identical(Colv, "Rowv")) {
      if (nr != nc)
        stop("Colv = \"Rowv\" but nrow(x) != ncol(x)")
      ddc <- ddr
    }
    else {
      distMatrixC = distfunC(if (symm)
        x
        else t(x))
      hcc <- hclustfun(distMatrixC, method = method)
      ddc <- as.dendrogram(hcc)
      if (!is.logical(Colv) || Colv)
        ddc <- reorderfun(ddc, Colv)
    }
    if (nc != length(colInd <- order.dendrogram(ddc)))
      stop("column dendrogram ordering gave index of wrong length")
  }
  else colInd <- 1L:nc
  x <- x[rowInd, colInd]
  labRow <- if (is.null(labRow))
    if (is.null(rownames(x)))
      (1L:nr)[rowInd]
  else rownames(x)
  else labRow[rowInd]
  labCol <- if (is.null(labCol))
    if (is.null(colnames(x)))
      (1L:nc)[colInd]
  else colnames(x)
  else labCol[colInd]
  if (scale == "row") {
    x <- sweep(x, 1L, rowMeans(x, na.rm = na.rm), check.margin = FALSE)
    sx <- apply(x, 1L, sd, na.rm = na.rm)
    x <- sweep(x, 1L, sx, "/", check.margin = FALSE)
  }
  else if (scale == "column") {
    x <- sweep(x, 2L, colMeans(x, na.rm = na.rm), check.margin = FALSE)
    sx <- apply(x, 2L, sd, na.rm = na.rm)
    x <- sweep(x, 2L, sx, "/", check.margin = FALSE)
  }
  lmat <- rbind(c(NA, 3), 2:1)
  lwid <- c(1, 4)
  lhei <- c(1 + if (!is.null(main)) 0.2 else 0, 4)
  if (!missing(ColSideFun)) {
    lmat <- rbind(lmat[1, ] + 1, c(NA, 1), lmat[2, ] + 1)
    lhei <- c(lhei[1L], ColSideWidth, lhei[2L])
  }
  else if (!missing(ColSideColors)) {
    if (!is.character(ColSideColors) & nrow(ColSideColors) !=
        nc)
      stop("'ColSideColors' must be a character vector or matrix of length ncol(x)")
    lmat <- rbind(lmat[1, ] + 1, c(NA, 1), lmat[2, ] + 1)
    lhei <- c(lhei[1L], 0.2 * round(ncol(ColSideColors)/2 +
                                      0.1), lhei[2L])
  }
  if (!missing(RowSideColors)) {
    if (!is.character(RowSideColors) || nrow(RowSideColors) !=
        nr)
      stop("'RowSideColors' must be a character vector or matrix of length nrow(x)")
    lmat <- cbind(lmat[, 1] + 1, c(rep(NA, nrow(lmat) - 1),
                                   1), lmat[, 2] + 1)
    lwid <- c(lwid[1L], 0.2 * round(ncol(RowSideColors)/2 +
                                      0.1), lwid[2L])
  }
  lmat <- lmat + 1
  lmat[is.na(lmat)] <- 0
  lmat[1, 1] <- 1
  dev.hold()
  on.exit(dev.flush())
  op <- par(no.readonly = TRUE)
  on.exit(par(op), add = TRUE)
  if (balanceColor) {
    if (abs(max(x, na.rm = T)) >= abs(min(x, na.rm = T))) {
      cut.off <- round(quantile(1:length(col), probs = 1 -
                                  (abs(max(x, na.rm = T)) + abs(min(x, na.rm = T)))/(2 *
                                                                                       abs(max(x, na.rm = T)))))
      col <- col[cut.off:length(col)]
    }
    else {
      cut.off <- round(quantile(1:length(col), probs = (abs(max(x,
                                                                na.rm = T)) + abs(min(x, na.rm = T)))/(2 * abs(min(x,
                                                                                                                   na.rm = T)))))
      col <- col[1:cut.off]
    }
  }
  graphics::layout(lmat, widths = lwid, heights = lhei, respect = TRUE)
  # layout.show(l)
  # stop()
  if (!missing(legendfun)) {
    par(mar = c(0, 0, 0, 0))
    par(xpd = NA)
    legendfun()
  }
  else {
    par(mar = c(5, 1, 1, 0))
    dummy.x <- seq(min(x, na.rm = TRUE), max(x, na.rm = TRUE),
                   length = length(col))
    dummy.z <- matrix(dummy.x, ncol = 1)
    image(x = dummy.x, y = 1, z = dummy.z, yaxt = "n", col = col,
          cex.axis = cexCol, xlab = "")
  }
  if (!missing(RowSideColors)) {
    par(mar = c(margins[1L], 0, 0, 0.5))
    if (revC) {
      rsc = RowSideColors[rev(rowInd), , drop = F]
    }
    else {
      rsc = RowSideColors[rowInd, , drop = F]
    }
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
    if (missing(RowSideLabs)) {
      if (ncol(RowSideColors) == 1 & colnames(RowSideColors)[1] ==
          "") {
        RowSideLabs <- ""
      }
      else {
        RowSideLabs <- colnames(RowSideColors)
      }
    }
    if (dim(rsc)[2] == 1) {
      axis(1, 0, RowSideLabs, las = 2, tick = FALSE)
    }
    else {
      axis(1, 0:(dim(rsc)[2] - 1)/(dim(rsc)[2] - 1), RowSideLabs,
           las = 2, tick = FALSE)
    }
  }
  if (!missing(ColSideCut)) {
    ColSideCutResult <- cut(ddc, ColSideCut)$lower
    ColSideCutResultSubIndList <- list()
    for (i in 1:length(ColSideCutResult)) {
      ColSideCutResultSubInd <- order.dendrogram(ColSideCutResult[[i]])
      ColSideCutResultSubIndList[[i]] <- ColSideCutResultSubInd
    }
    cutTable <- NULL
    if (verbose) {
      cat(paste0("The samples could be cut into ", length(ColSideCutResult),
                 " parts with height ", ColSideCut))
      cat("\n")
      if (!missing(ColSideAnn)) {
        for (i in 1:ncol(ColSideAnn)) {
          if (is.factor(ColSideAnn[, i])) {
            cutTable[[i]] <- sapply(ColSideCutResultSubIndList,
                                    function(x) table(ColSideAnn[x, i]))
            colnames(cutTable[[i]]) <- paste0("Cluster ",
                                              1:length(ColSideCutResult))
            names(cutTable)[i] <- colnames(ColSideAnn)[i]
            pvalue <- chisq.test(cutTable[[i]])$p.value
            cat(paste0("Differential distribution for ",
                       colnames(ColSideAnn)[i], ", p value by chi-squared test: ",
                       round(pvalue, 3), "\n"))
            cutTable[[i]] <- rbind(cutTable[[i]], round(cutTable[[i]][1,
            ]/colSums(cutTable[[i]]), 2))
            row.names(cutTable[[i]])[nrow(cutTable[[i]])] <- paste0(row.names(cutTable[[i]])[1],
                                                                    "_Percent")
            cutTable[[i]] <- cbind(cutTable[[i]], pValue = c(pvalue,
                                                             rep(NA, nrow(cutTable[[i]]) - 1)))
          }
          else {
            cutTable[[i]] <- sapply(split(ColSideAnn[unlist(ColSideCutResultSubIndList),
                                                     i], rep(1:length(ColSideCutResultSubIndList),
                                                             sapply(ColSideCutResultSubIndList, length))),
                                    function(x) summary(na.omit(x)))
            colnames(cutTable[[i]]) <- paste0("Cluster ",
                                              1:length(ColSideCutResult))
            names(cutTable)[i] <- colnames(ColSideAnn)[i]
            temp <- aov(ColSideAnn[unlist(ColSideCutResultSubIndList),
                                   i] ~ as.factor(rep(1:length(ColSideCutResultSubIndList),
                                                      sapply(ColSideCutResultSubIndList, length))))
            pvalue <- summary(temp)[[1]]$"Pr(>F)"[1]
            cat(paste0("Differential distribution for ",
                       colnames(ColSideAnn)[i], ", p value by ANOVA: ",
                       round(pvalue, 3), "\n"))
            cutTable[[i]] <- cbind(cutTable[[i]], pValue = c(pvalue,
                                                             rep(NA, 5)))
          }
        }
      }
    }
    ColSideCutResultCol <- rainbow(length(ColSideCutResult),
                                   alpha = 0.2)
    ColNumber <- (ncol(x) - 1)
  }
  if (!missing(ColSideFun)) {
    par(mar = c(0.5, 0, 0, margins[2L]))
    ColSideAnn <- ColSideAnn[colInd, , drop = F]
    ColAnnHeight <- ColSideFun(ColSideAnn)
    if (!exists("ColAnnHeight")) {
      ColAnnHeight <- par("usr")[3:4]
    }
    if (!missing(ColSideCut)) {
      rect(c(0 - 1/ColNumber/2, (0 - 1/ColNumber/2) + 1/ColNumber *
               cumsum(sapply(ColSideCutResult, function(x) length(unlist(x))))[-length(ColSideCutResult)]),
           ColAnnHeight[1], c((0 - 1/ColNumber/2) + 1/ColNumber *
                                cumsum(sapply(ColSideCutResult, function(x) length(unlist(x))))),
           ColAnnHeight[2], col = ColSideCutResultCol)
    }
  }
  else if (!missing(ColSideColors)) {
    par(mar = c(0.5, 0, 0, margins[2L]))
    csc = ColSideColors[colInd, , drop = F]
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
    if (missing(ColSideLabs)) {
      if (ncol(ColSideColors) == 1 & colnames(ColSideColors)[1] ==
          "") {
        ColSideLabs <- ""
      }
      else {
        ColSideLabs <- colnames(ColSideColors)
      }
    }
    if (dim(csc)[2] == 1) {
      axis(4, 0, ColSideLabs, las = 2, tick = FALSE)
    }
    else {
      axis(4, 0:(dim(csc)[2] - 1)/(dim(csc)[2] - 1), ColSideLabs,
           las = 2, tick = FALSE)
    }
  }
  par(mar = c(margins[1L], 0, 0, margins[2L]))
  if (!symm || scale != "none")
    x <- t(x)
  if (revC) {
    iy <- nr:1
    if (doRdend)
      ddr <- rev(ddr)
    x <- x[, iy]
  }
  else iy <- 1L:nr
  # browser()
  image(1L:nc, 1L:nr, x, xlim = 0.5 + c(0, nc), ylim = 0.5 +
          c(0, nr), axes = FALSE, xlab = "", ylab = "", col = col,
        useRaster = useRaster, ...)
  if (any(is.na(x))) {
    mmat <- ifelse(is.na(x), 1, NA)
    image(1L:nc, 1L:nr, mmat, axes = FALSE, xlab = "", ylab = "",
          col = na.color, useRaster=useRaster, add = TRUE)
  }
  if (!missing(colorCell)) {
    colorCell[, 1] <- match(colorCell[, 1], rowInd)
    colorCell[, 2] <- match(colorCell[, 2], colInd)
    rect(colorCell[, 2] - 0.5, colorCell[, 1] - 0.5, colorCell[,
                                                               2] + 0.5, colorCell[, 1] + 0.5, col = as.character(colorCell[,
                                                                                                                            3]), border = NA)
  }
  if (!missing(highlightCell)) {
    if (ncol(highlightCell) == 3) {
      highlightCell$lwd <- 1
    }
    highlightCell[, 1] <- match(highlightCell[, 1], rowInd)
    highlightCell[, 2] <- match(highlightCell[, 2], colInd)
    rect(highlightCell[, 2] - 0.5, highlightCell[, 1] - 0.5,
         highlightCell[, 2] + 0.5, highlightCell[, 1] + 0.5,
         border = as.character(highlightCell[, 3]), lwd = as.integer(highlightCell[,
                                                                                   4]))
  }
  if (!missing(ColSideColors) & ColAxisColors != 0) {
    mtext(1, at = 1L:nc, text = labCol, las = lasCol, line = 0.5,
          cex = cexCol, col = ColSideColors[colInd, ColAxisColors])
  }
  else {
    axis(1, 1L:nc, labels = labCol, las = lasCol, line = -0.5,
         tick = 0, cex.axis = cexCol)
  }
  if (!is.null(xlab))
    mtext(xlab, side = 1, line = margins[1L] - 1.25)
  if (!missing(RowSideColors) & RowAxisColors != 0) {
    mtext(4, at = iy, text = labRow, las = lasRow, line = 0.5,
          cex = cexRow, col = RowSideColors[rowInd, RowAxisColors])
  }
  else {
    axis(4, iy, labels = labRow, las = lasRow, line = -0.5,
         tick = 0, cex.axis = cexRow)
  }
  if (!is.null(ylab))
    mtext(ylab, side = 4, line = margins[2L] - 1.25)
  if (!missing(add.expr))
    eval(substitute(add.expr))
  par(mar = c(margins[1L], 0, 0, 0))
  if (doRdend & showRowDendro)
    plot(ddr, horiz = TRUE, axes = FALSE, yaxs = "i", leaflab = "none")
  else frame()
  par(mar = c(0, 0, if (!is.null(main)) 1 else 0, margins[2L]))
  if (doCdend & showColDendro) {
    plot(ddc, axes = FALSE, xaxs = "i", leaflab = "none")
    if (!missing(ColSideCut)) {
      rect(c(0.5, 0.5 + cumsum(sapply(ColSideCutResult,
                                      function(x) length(unlist(x))))[-length(ColSideCutResult)]),
           0, cumsum(sapply(ColSideCutResult, function(x) length(unlist(x)))) +
             0.5, ColSideCut, col = ColSideCutResultCol)
    }
  }
  else if (!is.null(main))
    frame()
  if (!is.null(main)) {
    par(xpd = NA)
    title(main, cex.main = 1.5 * op[["cex.main"]])
  }
  invisible(list(rowInd = rowInd, colInd = colInd, Rowv = if (keep.dendro &&
                                                              doRdend) ddr, Colv = if (keep.dendro && doCdend) ddc,
                 cutTable = if (!missing(ColSideAnn) && !missing(ColSideCut)) cutTable,
                 cutColoumIndList = if (!missing(ColSideCut)) ColSideCutResultSubIndList,
                 DistMatrixC = if (returnDistMatrix) distMatrixC, DistMatrixR = if (returnDistMatrix) distMatrixR,
                 hcr = hcr, hcc = hcc))
}

#_ -------------------------------------------------------

# fun: heatmapSummAb (Cycif),CycifStack ----

#' Create Heatmap Summary of Protein Expression
#'
#' This function generates a heatmap summarizing protein expression data from a `CycifStack` object.
#'
#' @param x A `CycifStack` object containing protein expression data.
#' @param norm_type Normalization type for the data, one of "log" or "logTh" (default is "log").
#' @param ab The name of the protein to be used for summarization.
#' @param sum_type The type of summarization to be performed, one of "freq", "mean", "median", or "x percentile" (e.g., "95 percentile").
#' @param ct_name The name of the cell type column (default is "default").
#' @param uniq_cts Vector of unique cell types to include in the heatmap.
#' @param uniq.smpls Vector of unique samples to include in the heatmap.
#' @param strict Logical, indicating strict filtering of cell types (default is FALSE).
#' @param scale Scaling method for the heatmap, one of "none", "row", or "column" (default is "none").
#'
#' @details
#' - The `heatmapSummAb` function creates a heatmap summarizing protein expression data for the specified protein.
#' - Users can choose from different normalization types and summarization methods.
#' - Cell types and samples can be filtered and customized for the heatmap.
#'
#' @export
#' @rdname heatmapSummAb
setGeneric("heatmapSummAb", function(x,...) standardGeneric("heatmapSummAb"))

#' @rdname heatmapSummAb
#' @export
setMethod("heatmapSummAb", "CycifStack",
          function(x,norm_type=c("log","logTh"),
                   summBy=c("ab","heatmap"),ab,
                   sum_type=c("freq","mean","median","x percentile"),
                   ct_name="default",
                   uniq_cts,
                   uniq.smpls,
                   strict=FALSE,
                   p_thres=0.5,
                   scale="none",...){
    options(dplyr.summarise.inform = FALSE)

    if(missing(sum_type)){
      stop("'sum_type' should be specified")
    }
    if(missing(ab)){
      stop("ab should be always specified")
    }

    # norm_type
    if(sum_type=="freq"){
      norm_type="logTh"
    }else if(missing(norm_type)){
      norm_type <- "log"
    }

    # sum_fun
    if(grepl("percentile$",sum_type)){
      pct <- strsplit(sum_type," ")[[1]]
      p.ile <- as.numeric(pct[[1]])
      sum_fun=function(y,...)quantile(y,p.ile/100,na.rm=T)
    }else if(sum_type=="median"){
      sum_fun=function(y,...)quantile(y,.5,na.rm=T)
    }else if(sum_type=="mean"){
      sum_fun=function(y,...)mean(y,na.rm=T)
    }else if(sum_type=="freq"){
      sum_fun=function(y,th){
        if(length(y)==0){
          return(NA)
        }else{
          return(sum(y > th, na.rm=T)/sum(!is.na(y)))
        }
      }
    }

    # balanceColor
    if(sum_type=="freq" || norm_type=="log"){
      balanceColor <- FALSE
    }else if(norm_type=="logTh"){
      balanceColor <- TRUE
    }

    ## cell types
    cts <- cell_types(x,ct_name=ct_name,strict=strict)

    df <- exprs(x,type=norm_type) %>%
      cbind(cts) %>%
      filter(cell_types != "outOfROI") %>%
      rename(smpl = sample) %>%
      rename(cell_type = cell_types)

    if(missing(uniq.smpls)){
      uniq.smpls <- levels(df$smpl)
    }

    if(missing(uniq_cts)){
      uniq_cts <- levels(df$cell_type)
    }

    df <- df %>%
      mutate(smpl = factor(smpl,levels=uniq.smpls)) %>%
      filter(smpl %in% uniq.smpls) %>%
      mutate(cell_type=factor(cell_type,levels=uniq_cts)) %>%
      mutate(cell_type=factor(cell_type)) %>%
      filter(!is.na(cell_type)) %>%
      dplyr::rename(this_ab=!!sym(ab)) %>%
      group_by(smpl,cell_type) %>%
      dplyr::summarise(sum_ab=sum_fun(this_ab,th=p_thres))

    ttl <- paste0(ab,",",sum_type,",",norm_type,",(scale:",scale,")")

    m1 <- as.matrix(with(df, Matrix::sparseMatrix(as.integer(smpl), as.integer(cell_type), x=sum_ab)))

    m1[apply(m1,1,function(x)any(is.nan(x))),] <- NA
    if(sum_type!="freq"){
      m1[m1==0] <- NA
    }

    rownames(m1) <- levels(df$smpl)
    colnames(m1) <- levels(df$cell_type)
    m1 <- m1[rev(seq(nrow(m1))),]

    if(sum_type!="freq" && norm_type=="logTh"){
      m1 <- m1 - 0.5 ## hardcode!!!!
    }

    ## cols
    if(sum_type=="mean"){
      cols <- viridis::viridis(1024)
    }else if(sum_type == "freq"){
      cols <- colorRampPalette(RColorBrewer::brewer.pal(9,"YlOrRd"))(1024)
    }
    h3(m1,
       margins=c(10,5),
       main=ttl,
       cexRow=.7,
       na.rm = T,
       col = cols,
       balanceColor=balanceColor,
       scale=scale,...)

    invisible(df)
  })
