# Classes ----

setClass(Class = "skewt",
         slots = c(
           call = "language",
           data = "numeric",
           n_days = "numeric",
           expDf = "numeric",
           value = "numeric",
           shape = "numeric",
           df = "numeric",
           escala = "numeric",
           location = "numeric",
           convergence = "numeric",
           quantiles = "data.frame"),
         package = "ECTools")

#' Moment estimation of the skew-t distribution
#'
#' @usage
#'
#' st = fit_skewt(x)
#' st
#' plot(st)
#' summary(st)
#' qst(0.99, df = get_param(st))
#' ks_test(st)
#' qqPlot(st)
#'
#' @param x Numeric vector to fit
#' @param n_days Number of days, to convolve the distribution
#'
#' @return This functions returns the estimated parameters of skew-t using the moment method and plot the observed and the fitted distribution
#'
#' * skewtpar functions return a list with the following values:
#'
#' * value, value of the function to minimize
#'
#' * shape shape parameter to skew-t distribution
#'
#' * df, nu parameter to skew-t distribution
#'
#' * escala, omega parameter to skew-t distribution
#'
#' * location, xi parameter to skew-t distribution
#'
#' * convergence, logical value indicating whether the optimization has converged
#' or not. 1 indicates successful convergence and 0 no convergence. A successful
#' convergence is considered when the value of objetive the funcion is less
#' than 0.0001
#'
#' @references Technical Note No. 2. Skew t Distribution of Economic Capital Models Department Methodology & Model Development Area
#'
#' @examples
#'
#' st = fit_skewt(x)
#' st
#' plot(st)
#' summary(st)
#' qst(0.99, df = get_param(st))
#' ks_test(st)
#' qqPlot(st)
#'
#' @export
#'
fit_skewt = function(x = x, n_days = 1) {

  if(!is.numeric(x)) stop("x must be numeric")

  call = match.call()

  q = c(0.0001, 0.0005, 0.001, 0.05, 0.1, 0.5, 0.9, 0.95, 0.999, 0.9995, 0.9999)

  # powel method in the powel package

  result = optim(c(0,1),
                 minfunc,
                 x = x,
                 n_days = n_days,
                 # method = "CG",
                 control = list(reltol = 1e-09))

  objet = list()

  object = new("skewt",
               call = call,
               n_days = n_days,
               data = x)

  object@expDf = result$par[2]
  object@value = result$value
  object@shape = result$par[1]
  object@df = exp(object@expDf) + 4

  st = st_cumulants(location = 0, escala = 1, shape = object@shape, df = object@df)

  object@escala = sqrt(k2(x) * n_days / st[2])
  object@location = (k1(x) * n_days)  - object@escala * st[1]

  object@convergence = ifelse(object@value > 0.0001, 1, 0)

  object@quantiles = data.frame(
    p = paste0(q * 100, "%"),
    empirical = unname(quantile(x, q)),
    theorical = qst(q, xi = object@location,
                    omega = object@escala,
                    nu = object@df,
                    alpha = object@shape),
    stringsAsFactors = F)

  if(object@convergence == 1) message("Warning the vector has not converged")

  return(object)
}

# methods ----

plot.skewt = function(x = x, n = 1e6) {

  object = x
  ggplot() +
    geom_density(aes(object@data,
                     color = "observed"),
                 size = 1) +
    geom_density(aes(rst(n = n,
                         dp = get_param(object)),
                     color = "estimated"),
                 size = 1) +
    xlab(NULL) +
    theme_bw()

}

qqPlot.skewt = function(x = x, ...) {

  qqPlot(x@data,
         dist="st",
         pch = 20,
         dp = get_param(x),
         ...)
}

show.skewt = function(object) {

  cat("An object of class \"skew t\"\n")
  cat("\nCall:\n", deparse(object@call), "\n\n",sep="")
  cat("Available slots:\n")
  print(slotNames(object))

}

#' Get parameters of skewt class
#'
#' @export
get_param = function(object, ...) {

  UseMethod("get_param", object)

}

get_param.skewt = function(object, ...) {

  param = c(object@location, object@escala, object@shape, object@df)
  names(param) = c("xi", "omega", "alpha", "nu")

  return(param)

}

#' Kolmogorov-Smirnov Tests
#'
#' @export
#'
ks_test = function(object, ...) {

  UseMethod("ks_test", object)

}

ks_test.skewt = function(object, ...) {

  ks.test(object@data,
          "pst",
          dp = get_param(object),
          ...)

}

# ad & rst

summary.skewt <- function(object, ...) {

  cat("+-----------------------------------+\n")
  cat("|            Skew t fitting         |\n")
  cat("+-----------------------------------+\n\n")
  cat("The vector", ifelse(object@convergence == 1, "has NOT", "has"), "converged\n\n")
  cat("Skew t settings: \n")
  cat(paste("  Number of days        = ", object@n_days, "\n"))
  cat(paste("  Convergence           = ", format(object@value, digits = 6, scientific = F), "\n\n"))
  cat("Estimation results:\n")
  cat(paste("  Location              = ", format(object@location, digits = 4, scientific = F), "\n"))
  cat(paste("  Escala                = ", format(object@escala, digits = 4, scientific = F), "\n"))
  cat(paste("  Shape                 = ", format(object@shape, digits = 4, scientific = F), "\n"))
  cat(paste("  df                    = ", format(object@df, digits = 4, scientific = F), "\n\n"))
  cat("Quantiles:\n\n")
  print(object@quantiles, digits = 4, row.names = F, scipen=999)

}

# setMethods & export ----

setMethod("show", "skewt", show.skewt)
setMethod("get_param", "skewt", get_param.skewt)
setMethod("print", "skewt", function(x, ...) str(x))
setMethod("ks_test", "skewt", ks_test.skewt)
setMethod("qqPlot", "skewt", qqPlot.skewt)

#' @export
setMethod("summary", signature(object = "skewt"), summary.skewt)

#' @export
setMethod("plot", signature(x = "skewt", y = "missing"), plot.skewt)
