#' Anderson-Darling Test of Goodness-of-Fit
#'
#' @description Performs the Anderson-Darling test of goodness-of-fit to a specified continuous univariate proba
#' bility distribution
#'
#' @usage ad.test(x, null = "punif", ..., nullname)
#'
#' @param x Numeric vector of data values.
#' @param null A function, or a character string giving the name of a function, to compute the cumulative
#' distribution function for the null distribution
#' @param ... Additional arguments for the cumulative distribution function.
#' @param nullname Optional character string describing the null distribution. The default is "uniform distribution".
#'
#' @details  This command performs the Anderson-Darling test of goodness-of-fit to the distribution specified
#' by the argument null. It is assumed that the values in x are independent and identically distributed
#' random values, with some cumulative distribution function F. The null hypothesis is that F is the
#' function specified by the argument null, while the alternative hypothesis is that F is some other
#' function.
#'
#' @return An object of class "htest" representing the result of the hypothesis test
#'
#' @author Original C code by George Marsaglia and John Marsaglia. R interface by Adrian Baddeley
#'
#' @references
#'
#' Anderson, T.W. and Darling, D.A. (1952) Asymptotic theory of certain ’goodness-of-fit’ criteria
#' based on stochastic processes. Annals of Mathematical Statistics 23, 193–212.
#'
#' Anderson, T.W. and Darling, D.A. (1954) A test of goodness of fit. Journal of the American Statistical Association 49, 765–769.
#'
#' Marsaglia, G. and Marsaglia, J. (2004) Evaluating the Anderson-Darling Distribution. Journal of Statistical Software 9 (2), 1–5.
#' February 2004. http://www.jstatsoft.org/v09/i02
#'
#' @examples
#'
#' x <- rnorm(10, mean=2, sd=1)
#' ad.test(x, "pnorm", mean=2, sd=1)
#'
#' @export
ad.test <- function(x, null="punif", ..., nullname) {
  xname <- deparse(substitute(x))
  nulltext <- deparse(substitute(null))
  if(is.character(null)) nulltext <- null
  if(missing(nullname) || is.null(nullname)) {
    reco <- recogniseCdf(nulltext)
    nullname <- if(!is.null(reco)) reco else
      paste("distribution", sQuote(nulltext))
  }
  stopifnot(is.numeric(x))
  x <- as.vector(x)
  n <- length(x)
  F0 <- if(is.function(null)) null else
    if(is.character(null)) get(null, mode="function") else
      stop("Argument 'null' should be a function, or the name of a function")
  U <- F0(x, ...)
  if(any(U < 0 | U > 1))
    stop("null distribution function returned values outside [0,1]")
  U <- sort(U)
  k <- seq_len(n)
  ## call Marsaglia C code
  z <- .Call("ECTools_ADtestR",
          x = as.double(U),
          n = as.integer(n),
          adstat = as.double(numeric(1)),
          pvalue = as.double(numeric(1)),
          PACKAGE="ECTools"
  )
  STATISTIC <- z$adstat
  names(STATISTIC) <- "An"
  PVAL <- z$pvalue
  METHOD <- c("Anderson-Darling test of goodness-of-fit",
              paste("Null hypothesis:", nullname))
  extras <- list(...)
  parnames <- intersect(names(extras), names(formals(F0)))
  if(length(parnames) > 0) {
    pars <- extras[parnames]
    pard <- character(0)
    for(i in seq_along(parnames))
      pard[i] <- paste(parnames[i], "=", paste(pars[[i]], collapse=" "))
    pard <- paste("with",
                  ngettext(length(pard), "parameter", "parameters"),
                  "  ",
                  paste(pard, collapse=", "))
    METHOD <- c(METHOD, pard)
  }
  out <- list(statistic = STATISTIC,
              p.value = PVAL,
              method = METHOD,
              data.name = xname)
  class(out) <- "htest"
  return(out)
}

pAD <- function(q, n=Inf, lower.tail=TRUE, fast=TRUE) {
  q <- as.numeric(q)
  p <- rep(NA_real_, length(q))
  if(any(ones <- is.infinite(q) & (q == Inf)))
    p[ones] <- 1
  if(any(zeroes <- (is.finite(q) & q <= 0) | (is.infinite(q) & (q == -Inf))))
    p[zeroes] <- 0
  ok <- is.finite(q) & (q > 0)
  nok <- sum(ok)
  if(nok > 0) {
    if(is.finite(n)) {
      z <- .Call("ECTools_ADprobN",
              a       = as.double(q[ok]),
              na      = as.integer(nok),
              nsample = as.integer(n),
              prob    = as.double(numeric(nok)),
              PACKAGE="ECTools")
      p[ok] <- z$prob
    } else if(fast) {
      ## fast version adinf()
      z <- .Call("ECTools_ADprobApproxInf",
              a    = as.double(q[ok]),
              na   = as.integer(nok),
              prob = as.double(numeric(nok)),
              PACKAGE="ECTools")
      p[ok] <- z$prob
    } else {
      ## slow, accurate version ADinf()
      z <- .Call("ECTools_ADprobExactInf",
              a    = as.double(q[ok]),
              na   = as.integer(nok),
              prob = as.double(numeric(nok)),
              PACKAGE="ECTools")
      p[ok] <- z$prob
    }

  }
  if(!lower.tail)
    p <- 1 - p
  return(p)
}

qAD <- local({

  f <- function(x, N, P, Fast) {
    pAD(x, N, fast=Fast) - P
  }

  qAD <- function(p, n=Inf, lower.tail=TRUE, fast=TRUE) {
    ## quantiles of null distribution of Anderson-Darling test statistic
    stopifnot(all(p >= 0))
    stopifnot(all(p <= 1))
    if(!lower.tail) p <- 1-p
    ans <- rep(NA_real_, length(p))
    for(i in which(p >= 0 & p < 1))
      ans[i] <- uniroot(f, c(0, 1), N=n, P=p[i], Fast=fast, extendInt="up")$root
    return(ans)
  }

  qAD
})


recogniseCdf <- function(s="punif") {
  if(!is.character(s) || length(s) != 1) return(NULL)
  if(nchar(s) <= 1 || substr(s,1,1) != "p") return(NULL)
  root <- substr(s, 2, nchar(s))
  a <- switch(root,
              beta     = "beta",
              binom    = "binomial",
              birthday = "birthday coincidence",
              cauchy   = "Cauchy",
              chisq    = "chi-squared",
              exp      = "exponential",
              f        = "F",
              gamma    = "Gamma",
              geom     = "geometric",
              hyper    = "hypergeometric",
              lnorm    = "log-normal",
              logis    = "logistic",
              nbinom   = "negative binomial",
              norm     = "Normal",
              pois     = "Poisson",
              t        = "Student's t",
              tukey    = "Tukey (Studentized range)",
              unif     = "uniform",
              weibull  = "Weibull",
              NULL)
  if(!is.null(a))
    return(paste(a, "distribution"))
  b <- switch(root,
              AD     = "Anderson-Darling",
              CvM    = "Cramer-von Mises",
              wilcox = "Wilcoxon Rank Sum",
              NULL)
  if(!is.null(b))
    return(paste("null distribution of", b, "Test Statistic"))
  return(NULL)
}
