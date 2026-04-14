#' @import Matrix
#' @keywords internal
aggregate.Matrix <- function (x, groupings = NULL, form = NULL, fun = "sum", ...)
{
  if (!is(x, "Matrix"))
    x <- Matrix(as.matrix(x), sparse = TRUE)
  if (fun == "count")
    x <- x != 0
  groupings2 <- groupings
  if (!is(groupings2, "data.frame"))
    groupings2 <- as(groupings2, "data.frame")
  groupings2 <- data.frame(lapply(groupings2, as.factor))
  groupings2 <- data.frame(interaction(groupings2, sep = "_"))
  colnames(groupings2) <- "A"
  if (is.null(form))
    form <- as.formula("~0+.")
  form <- as.formula(form)
  mapping <- dMcast(groupings2, form)
  colnames(mapping) <- substring(colnames(mapping), 2)
  result <- t(mapping) %*% x
  if (fun == "mean")
    result@x <- result@x/(aggregate.Matrix(x, groupings2,
                                           fun = "count"))@x
  attr(result, "crosswalk") <- grr::extract(groupings, match(rownames(result),
                                                             groupings2$A))
  return(result)
}

#' @import Matrix
#' @keywords internal
dMcast <- function (data, formula, fun.aggregate = "sum", value.var = NULL, as.factors = FALSE, factor.nas = TRUE, drop.unused.levels = TRUE)
{
  values <- 1
  if (!is.null(value.var))
    values <- data[, value.var]
  alltms <- terms(formula, data = data)
  response <- rownames(attr(alltms, "factors"))[attr(alltms, "response")]
  tm <- attr(alltms, "term.labels")
  interactionsIndex <- grep(":", tm)
  interactions <- tm[interactionsIndex]
  simple <- setdiff(tm, interactions)
  i2 <- strsplit(interactions, ":")
  newterms <- unlist(lapply(i2, function(x) paste("paste(",paste(x, collapse = ","), ",", "sep='_'", ")")))
  newterms <- c(simple, newterms)
  newformula <- as.formula(paste("~0+", paste(newterms, collapse = "+")))
  allvars <- all.vars(alltms)
  data <- data[, c(allvars), drop = FALSE]
  if (as.factors)
    data <- data.frame(lapply(data, as.factor))
  characters <- unlist(lapply(data, is.character))
  data[, characters] <- lapply(data[, characters, drop = FALSE], as.factor)
  factors <- unlist(lapply(data, is.factor))
  data[, factors] <- lapply(data[, factors, drop = FALSE],
                            function(x) {
                              if (factor.nas)
                                if (any(is.na(x))) {
                                  levels(x) <- c(levels(x), "NA")
                                  x[is.na(x)] <- "NA"
                                }
                              if (drop.unused.levels)
                                if (nlevels(x) != length(na.omit(unique(x))))
                                  x <- factor(as.character(x))
                              y <- contrasts(x, contrasts = FALSE, sparse = TRUE)
                              attr(x, "contrasts") <- y
                              return(x)
                            })
  attr(data, "na.action") <- na.pass
  result <- Matrix::sparse.model.matrix(newformula, data, drop.unused.levels = FALSE, row.names = FALSE)
  brokenNames <- grep("paste(", colnames(result), fixed = TRUE)
  colnames(result)[brokenNames] <- lapply(colnames(result)[brokenNames],
                                          function(x) {
                                            x <- gsub("paste(", replacement = "", x = x, fixed = TRUE)
                                            x <- gsub(pattern = ", ", replacement = "_", x = x,
                                                      fixed = TRUE)
                                            x <- gsub(pattern = "_sep = \"_\")", replacement = "",
                                                      x = x, fixed = TRUE)
                                            return(x)
                                          })
  result <- result * values
  if (isTRUE(response > 0)) {
    responses = all.vars(terms(as.formula(paste(response, "~0"))))
    result <- aggregate.Matrix(result, data[, responses, drop = FALSE], fun = fun.aggregate)
  }
  return(result)
}

#' Compute pct_in and pct_out for arbitrary cell groups
#'
#' @param seurat Seurat object
#' @param cells_in Character vector of cell barcodes defining the "in" group
#' @param cells_out Character vector of cell barcodes defining the "out" group
#' @param assay Assay name to use (default: "RNA")
#' @param slot Slot with expression data (usually "data" or "counts")
#' @param features Optional vector of features to restrict to (default: all)
#'
#' @return A data.frame with columns:
#'   feature, pct_in, pct_out
#'   where:
#'     pct_in  = fraction of cells in cells_in with expr > 0
#'     pct_out = fraction of cells in cells_out with expr > 0
#'
#' @keywords internal
CalcPctInOutByCells <- function(seurat,
                                cells_in,
                                cells_out,
                                assay   = "RNA",
                                slot    = "counts",
                                features = NULL) {
  # Basic checks
  if (!inherits(seurat, "Seurat")) {
    stop("'seurat' must be a Seurat object")
  }
  if (!assay %in% names(seurat@assays)) {
    stop(paste0("Assay '", assay, "' not found in Seurat object"))
  }
  valid_slots <- c("counts", "data", "scale.data")
  if (!slot %in% valid_slots) {
    stop("Slot '", slot, "' is not one of: ",
         paste(valid_slots, collapse = ", "))
  }

  all_cells <- colnames(seurat)
  cells_in  <- intersect(cells_in, all_cells)
  cells_out <- intersect(cells_out, all_cells)

  if (length(cells_in) == 0) {
    stop("No valid cells found in 'cells_in'")
  }
  if (length(cells_out) == 0) {
    stop("No valid cells found in 'cells_out'")
  }

  # Pull expression matrix
  mat <- Seurat::GetAssayData(seurat, assay = assay, slot = slot)
  #if (!is.matrix(mat)) {
  #  mat <- as.matrix(mat)
  #}

  # Restrict to requested features if provided
  if (!is.null(features)) {
    features <- intersect(features, rownames(mat))
    if (length(features) == 0) {
      stop("None of the requested features found in assay '", assay, "'")
    }
    mat <- mat[features, , drop = FALSE]
  }

  # Subset matrix to in/out cells
  mat_in  <- mat[, cells_in,  drop = FALSE]
  mat_out <- mat[, cells_out, drop = FALSE]

  # Handle degenerate cases
  if (ncol(mat_in) == 0 || ncol(mat_out) == 0) {
    warning("One of the groups has no cells; setting pct_in/out to NA.")
    pct_in  <- rep(NA_real_, nrow(mat))
    pct_out <- rep(NA_real_, nrow(mat))
  } else {
    pct_in  <- Matrix::rowMeans(mat_in  > 0)
    pct_out <- Matrix::rowMeans(mat_out > 0)
  }

  data.frame(
    feature = rownames(mat),
    pct_in  = pct_in,
    pct_out = pct_out,
    row.names = NULL,
    stringsAsFactors = FALSE
  )
}
