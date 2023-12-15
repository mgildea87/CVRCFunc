#' @import Matrix
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
