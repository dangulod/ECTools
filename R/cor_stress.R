cor.stress = function(x, y, period = period, N_SIMUL = N_SIMUL, abs = abs, ...) {

  if (length(x) != length(y)) stop("Vectors do not have same length")
  if(!is.logical(abs)) stop("abs must be logical")

  if (!is.numeric(x) | !is.numeric(y)) return(NA)
  if (identical(x, y)) return(1)

  z = complete.cases(x, y)

  x1 = x[z]
  y1 = y[z]

  if ((length(x1) - period) < 1) return(cor(x, y, use = "pairwise.complete.obs"))

  n = length(x1) - period

  v_cor = rep(0, n)

  for (i in 1:n) {

    x2 = x1[i:(i + period)]
    y2 = y1[i:(i + period)]

    v_cor[i] = cor(x2, y2, ...)
  }

  cor.stress = if (abs) max(abs(v_cor)) else max(v_cor)

  return(cor.stress)
}

#' Stressed correlation
#'
#' @param x a numeric vector, matrix or data frame.
#' @param y NULL (default) or a vector, matrix or data frame with compatible
#' dimensions to x. The default is equivalent to y = x (but more efficient).
#' @param period numeric number of periods to use
#' @param abs logical. Should absolute values be used
#' @param use an optional character string giving a method for computing
#' covariances in the presence of missing values. This must be (an
#' abbreviation of) one of the strings "everything", "all.obs",
#' "complete.obs", "na.or.complete", or "pairwise.complete.obs"
#' @param method a character string indicating which correlation coefficient
#' (or covariance) is to be computed. One of "pearson" (default), "kendall",
#' or "spearman": can be abbreviated.
#' @param ...
#'
#' @return Return the stressed correlation of x and y in the selected period
#' if these are vectors. If x and y are matrices then the correlations
#' between the columns of x and the columns of y are returned
#'
#' @export
#'
cor_stress = function (x, y = NULL, period = period, abs = F,
                       use = "everything",
                       method = c("pearson", "kendall", "spearman"), ...) {

  if (is.null(y)) {
    ncy <- ncx <- ncol(x)
    if (ncx == 0)
      stop("'x' is empty")
    r <- matrix(0, nrow = ncx, ncol = ncy)
    for (i in seq_len(ncx)) {
      for (j in seq_len(ncy)) {
        r[i,j] = cor.stress(x[[i]], x[[j]], period = period, abs = abs, ...)
      }
    }
    rownames(r) <- colnames(x)
    colnames(r) <- colnames(x)
    return(r)
  } else {
    return(cor.stress(x, y, period = period, ...))
  }
}
