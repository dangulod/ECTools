#' pairwise.diss
#'
pairwise.diss <- function(series, dissfun, ...) {
  n <- length(series)
  distances <- matrix(0, n, n)
  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      tryCatch( {
        d <- dissfun( series, i, j, ...)
        distances[i,j] <- d
        distances[j,i] <- d
      }, error = function (e) {
        stop( paste("Applying diss, series (",i,",",j,") produced the following error: ", e) )
      })
    }
  }
  as.dist((distances))
}

.common.ts.sanity.check <- function(x) {
  if (missing(x)) {
    stop("At least one series is missing!")
  }
  if (!is.numeric(x)) {
    stop("Series must be numeric")
  }
  #check length
  if (length(x) < 2) {
    stop("Incorrect length of the series")
  }
  if (!is.null(dim(x))) {
    stop("Incorrect dimension of the series, please input univarate series")
  }
}

.ts.freq.check <- function(x, y) {
  if (is.ts(x) && is.ts(y)) { #check their frequencies
    cbind(x,y)
  }
}

.ts.sanity.check = function(x,y) {
  .common.ts.sanity.check(x)
  .common.ts.sanity.check(y)
  .ts.freq.check(x,y)
}

#' diss.COR
#'
#' Funcion modificada del paquete TSclust para que adminta series con NAs
#'
diss.COR = function(x, y, beta = NULL) {
  .ts.sanity.check(x, y)
  correl <- cor(x,y, use = "pairwise.complete.obs")
  if (is.na(correl) | correl == -1) correl = -0.99
  if (is.null(beta)) {
    return(sqrt(2*(1- correl)))
  } else {
    if (beta<0) {
      stop("beta must be greater than 0")
    }
    return(sqrt(((1-correl)/(1+correl ))**beta ))
  }
}

#' noindicesdiss
#'
noindicesdiss <- function( fun ) {
  function(series, i, j, ...) {
    fun(series[[i]], series[[j]], ...)
  }
}

#' Dissimilarity Index Combining Temporal Correlation and Raw Values
#'
#' @description Computes an adaptive dissimilarity index between two time series that covers both dissimilarity on
#' raw values and dissimilarity on temporal correlation behaviors.
#' @param SERIES data
#' @param ... beta = NULL
#'
#' @details ¡¡FILL!!
#'
#' @return The computed distance.
#'
#' @references Chouakria-Douzal, A. and Nagabhushan P. N. (2007) Adaptive dissimilarity index for measuring
#' time series proximity. Adv. Data Anal. Classif., 1(1), 5–21.
#'
#' Montero, P and Vilar, J.A. (2014) TSclust: An R Package for Time Series Clustering. Journal of
#' Statistical Software, 62(1), 1-43. http://www.jstatsoft.org/v62/i01/.
#'
#' @export
#'
diss_cor = function(SERIES, ...) {

  if (!is.matrix(SERIES) && !is.list(SERIES) && !is.mts(SERIES)) {
    stop("list, mts, matrix or data.frame object is required for SERIES ")
  }

  mat.ser = SERIES

  if (is.mts(SERIES)) {
    SERIES = t(as.matrix(SERIES))
  }

  if (!is.list(SERIES)) {
    tmpser = SERIES
    SERIES = list()
    for (i in 1:nrow(tmpser)) {
      SERIES[[i]] = tmpser[i,]
    }
    names(SERIES) = rownames(tmpser)
  }
  if (length(SERIES) < 2) {
    stop("Only one series provided")
  }

  list.to.matrix = function(series) {
    n = length(series)
    k = length(series[[n]])
    mat.ser = matrix(0, n, k)
    for (i in 1:n) {
      if ( length( series[[i]]) != k ) {
        stop("diss method requires same length series")
      }
      mat.ser[i,] = series[[i]]
    }
    rownames(mat.ser) = names(series)
    mat.ser
  }

  out.dist = NULL

  diss.fun = diss.COR

  out.dist = pairwise.diss(SERIES, noindicesdiss(diss.fun), ...)

  out.dist = as.matrix(out.dist)
  rownames(out.dist) = names(SERIES)
  colnames(out.dist) = names(SERIES)
  out.dist = as.dist(out.dist)

  out.dist
}
