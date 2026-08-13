#_ -------------------------------------------------------

# fun: selectionSummary ----

#' Report how many eligible cells are available per sample, before selecting any
#'
#' A diagnostic step, meant to run before \code{computeSelection()}: reports how
#' many cells per sample match \code{used.cts}, and the percentile distribution of
#' those per-sample counts, so you can pick an \code{n_per_sample} value with the
#' actual numbers in front of you rather than guessing. Uses the exact same
#' eligibility filter \code{computeSelection()} will use, so the numbers here are
#' guaranteed to match what selection actually does.
#'
#' @param x A CycifStack object.
#' @param used.cts A character vector of eligible cell types.
#' @param ct_name Which cell-typing to use (default "default").
#' @param ... Additional arguments (unused).
#'
#' @return Invisibly, a list with \code{per_sample} (a data.frame of sample x
#' n_eligible) and \code{percentiles} (the 10th/20th/.../100th percentile of
#' those counts). Also prints both as a side effect.
#'
#' @seealso \code{\link{computeSelection}}
#'
#' @rdname selectionSummary
#' @export
setGeneric("selectionSummary", function(x,...) standardGeneric("selectionSummary"))

#' @rdname selectionSummary
#' @export
setMethod("selectionSummary", "CycifStack", function(x, used.cts, ct_name="default"){
  cts <- cell_types(x, ct_name=ct_name)
  eligible <- cts$cell_types %in% used.cts

  per_sample <- as.data.frame(table(sample = as.character(cts$sample)[eligible]))
  names(per_sample) <- c("sample","n_eligible")

  ## table() silently omits samples with zero eligible cells -- make sure
  ## every sample in x appears, with 0 rather than being missing entirely.
  all_smpls <- data.frame(sample = names(x), stringsAsFactors = FALSE)
  per_sample <- merge(all_smpls, per_sample, by = "sample", all.x = TRUE)
  per_sample$n_eligible[is.na(per_sample$n_eligible)] <- 0
  per_sample <- per_sample[match(names(x), per_sample$sample), ]
  rownames(per_sample) <- NULL

  pct <- seq(10, 100, by = 10)
  pct_vals <- stats::quantile(per_sample$n_eligible, probs = pct/100, type = 7)
  percentiles <- data.frame(percentile = paste0(pct, "%"), n_cells = round(pct_vals), row.names = NULL)

  cat("== selectionSummary:", nrow(per_sample), "samples, used.cts =",
      paste(used.cts, collapse=", "), "==\n")
  print(per_sample)
  cat("\n")
  print(percentiles)

  invisible(list(per_sample = per_sample, percentiles = percentiles))
})

#_ -------------------------------------------------------

# fun: computeSelection ----

#' Choose a reproducible subset of cells for a specific downstream analysis
#'
#' Builds a \code{CellSelection}: an explicit, named, reusable record of which
#' cells were chosen (which cell types, how many per sample, which seed), rather
#' than an inline \code{sample()} call repeated by hand in every analysis. Fresh
#' per analysis is the intended usage -- a Mac-cluster selection and a
#' Cancer-cluster selection are two different \code{CellSelection}s, not one
#' shared subsample.
#'
#' \code{n_per_sample} may be a single number (same cap applied to every sample,
#' capped further by however many eligible cells that sample actually has), or a
#' NAMED vector with one entry per sample -- matched by name against
#' \code{names(x)}, never by position, so a reordered or mismatched vector errors
#' loudly instead of silently mis-assigning caps to the wrong samples.
#'
#' If \code{n_per_sample} is not supplied, the 20th percentile of eligible cells
#' per sample (via \code{selectionSummary()}) is used, and this is printed --
#' never a silent default.
#'
#' @param x A CycifStack object.
#' @param used.cts A character vector of eligible cell types.
#' @param n_per_sample A single number, or a named vector (names = names(x)). If
#' missing, defaults to the 20th percentile of eligible-cells-per-sample (see
#' Details), and prints what value was used.
#' @param seed Random seed, for reproducibility.
#' @param ct_name Which cell-typing to use (default "default").
#' @param ... Additional arguments (unused).
#'
#' @return A CellSelection object. Pass it to \code{subsetCells()}.
#'
#' @seealso \code{\link{selectionSummary}}, \code{\link{subsetCells}}
#'
#' @rdname computeSelection
#' @export
setGeneric("computeSelection", function(x,...) standardGeneric("computeSelection"))

#' @rdname computeSelection
#' @export
setMethod("computeSelection", "CycifStack",
  function(x, used.cts, n_per_sample=NULL, seed=12345, ct_name="default"){
    call1 <- sys.call()
    uniq_smpls <- names(x)

    cts <- cell_types(x, ct_name=ct_name)
    eligible <- cts$cell_types %in% used.cts
    smpl_of_cell <- as.character(cts$sample)

    if (is.null(n_per_sample)) {
      summ <- selectionSummary(x, used.cts=used.cts, ct_name=ct_name)
      n_per_sample <- summ$percentiles$n_cells[summ$percentiles$percentile == "20%"]
      message("n_per_sample not specified -- using ", n_per_sample,
              " (20th percentile of ", length(uniq_smpls),
              " samples' eligible counts). Pass n_per_sample explicitly to override.")
    }

    if (length(n_per_sample) == 1) {
      n_per_sample_named <- stats::setNames(rep(n_per_sample, length(uniq_smpls)), uniq_smpls)
    } else {
      if (is.null(names(n_per_sample)) || !setequal(names(n_per_sample), uniq_smpls)) {
        stop("n_per_sample must be a single number, or a NAMED vector with one entry ",
             "per sample, names exactly matching names(x). Got names: ",
             paste(names(n_per_sample), collapse=", "))
      }
      n_per_sample_named <- n_per_sample[uniq_smpls] # reorder by name -- never trust position
    }

    set.seed(seed)
    is_used <- rep(FALSE, length(eligible))
    for (nm in uniq_smpls) {
      idx <- which(eligible & smpl_of_cell == nm)
      n_take <- min(n_per_sample_named[[nm]], length(idx))
      if (n_take > 0) {
        is_used[sample(idx, n_take)] <- TRUE
      }
    }

    new("CellSelection",
        used.cts = used.cts,
        n_per_sample = n_per_sample_named,
        seed = seed,
        is_used = is_used,
        source_samples = uniq_smpls,
        call = call1)
  }
)

#_ -------------------------------------------------------

# fun: subsetCells ----

#' Assemble a compact, analysis-ready table for a CellSelection
#'
#' Pulls together center-cell expression, neighborhood composition/expression
#' (if \code{computeCN()} has been run), and clinical metadata for exactly the
#' cells named by \code{selection}, into one \code{CellFeatures} object -- a
#' thin wrapper around a single \code{data.table} plus the selection's
#' provenance. This is the one place row-alignment across those different
#' sources is handled; nothing downstream should need to re-derive it.
#'
#' @param x A CycifStack object. If it has a \code{cell_neighborhood} computed
#' (via \code{computeCN()}), neighbor_* columns are included; cells in
#' \code{selection} that aren't covered by that CellNeighborhood (e.g. outside
#' its \code{used.cts}) get NA in those columns, with a warning.
#' @param selection A CellSelection, from \code{computeSelection()}.
#' @param norm_type Which expression slot to use for center-cell expression: "log" or "logTh".
#' @param ct_name Which cell-typing to use (default "default").
#' @param ... Additional arguments (unused).
#'
#' @return A CellFeatures object.
#'
#' @seealso \code{\link{computeSelection}}
#'
#' @rdname subsetCells
#' @export
setGeneric("subsetCells", function(x,...) standardGeneric("subsetCells"))

#' @rdname subsetCells
#' @export
setMethod("subsetCells", "CycifStack",
  function(x, selection, norm_type=c("log","logTh"), ct_name="default"){
    stopifnot(is(selection, "CellSelection"))
    norm_type <- match.arg(norm_type)

    sel_idx <- which(selection@is_used) # full-object indices of selected cells

    cts <- cell_types(x, ct_name=ct_name)
    stopifnot(length(selection@is_used) == nrow(cts))

    xy <- data.frame(data.table::rbindlist(cyApply(x, xys)))
    expr <- exprs(x, type=norm_type)
    stopifnot(nrow(xy) == nrow(cts), nrow(expr) == nrow(cts))

    dt <- data.table::data.table(
      global_id  = sel_idx,
      sample     = as.character(cts$sample[sel_idx]),
      cell_types = cts$cell_types[sel_idx]
    )
    dt <- cbind(dt, xy[sel_idx, , drop=FALSE])
    dt <- cbind(dt, expr[sel_idx, , drop=FALSE])
    stopifnot(nrow(dt) == length(sel_idx))

    ## neighbor composition/expression, if available -- joined by matching each
    ## selected cell's global_id against cell_neighborhood's OWN is_used mask
    ## (a different, larger population than `selection`), never by row position.
    cn <- tryCatch(cell_neighborhood(x), error=function(e) NULL)
    if (!is.null(cn) && length(cn@is_used) == length(selection@is_used)) {
      cn_row <- match(sel_idx, which(cn@is_used))
      has_cn <- !is.na(cn_row)
      if (any(!has_cn)) {
        warning(sum(!has_cn), " selected cell(s) have no cell_neighborhood result ",
                "(not eligible under computeCN()'s used.cts) -- their neighbor_* columns are NA")
      }

      nbhd_counts <- matrix(NA_real_, nrow=length(sel_idx), ncol=ncol(cn@neighbor_counts),
                             dimnames=list(NULL, paste0("nbhd_", colnames(cn@neighbor_counts))))
      nbhd_counts[has_cn, ] <- cn@neighbor_counts[cn_row[has_cn], , drop=FALSE]
      dt <- cbind(dt, as.data.frame(nbhd_counts))

      expr_cols <- setdiff(colnames(cn@mean_expr_per_nbhd), "cell_id")
      nbhd_expr <- matrix(NA_real_, nrow=length(sel_idx), ncol=length(expr_cols),
                           dimnames=list(NULL, paste0("nbhd_expr_", expr_cols)))
      nbhd_expr[has_cn, ] <- as.matrix(cn@mean_expr_per_nbhd[cn_row[has_cn], expr_cols, with=FALSE])
      dt <- cbind(dt, as.data.frame(nbhd_expr))
    }

    ## phenoData, joined by sample -- left_join preserves dt's row order
    ph <- pData(x)
    if (is.data.frame(ph) && nrow(ph) > 0 && "id" %in% names(ph)) {
      dt <- dplyr::left_join(dt, ph, by = c("sample" = "id"))
    }
    data.table::setDT(dt)
    stopifnot(nrow(dt) == length(sel_idx))

    new("CellFeatures", data = dt, selection = selection)
  }
)
