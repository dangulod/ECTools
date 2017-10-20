# Class definition -------------------------------------------------------------------

setClass(Class = "sensitivities",
         slots = c(
           call = "language",
           hist = "matrix",
           map = "data.frame",
           lim = "numeric",
           maxiter = "numeric",
           sensitivities = "data.frame",
           error = "numeric",
           iterations = "numeric",
           CD = "matrix",
           suggestions = "numeric"),
         package = "ECTools")

# Optimization function ---------------------------------------------------

#' Global & Local factor optimization
#'
#' @param map Dataframe with the correspondence
#' @param hist Matrix with historical correlation
#' @param CD Matrix with credit drivers correlation
#' @param lim Numeric it acts as a limit of FG ^ 2 + FL ^ 2, by default is 0.6
#' @param maxiter Numeric the maximum number of iterations to run
#' @param parallel Logical, argument specifying if parallel computing should be used
#' (TRUE) or not (FALSE, default) for evaluating the fitness function. This argument
#' could also be used to specify the number of cores to employ; by default, this is
#' taken from detectCores. Finally, the functionality of parallelization depends on
#' system OS: on Windows only 'snow' type functionality is available, while on Unix/Linux/Mac
#' OSX both 'snow' and 'multicore' (default) functionalities are available.
#' @param ... Arguments to be passed to ga function
#'
#' @examples
#' x = cor_optim(map = map,
#'               hist = hist,
#'               CD = CD,
#'               lim = 1,
#'               maxiter = 500,
#'               run = 10)
#'
#' x = cor_optim(map = map,
#'               hist = hist,
#'               CD = CD,
#'               lim = 1,
#'               maxiter = 500,
#'               run = 10,
#'               suggestions = get_suggestions(x))
#'
#' @export
#'
cor_optim = function(map = map, hist = hist, CD = CD, lim = 1, maxiter = 1e4, parallel = F, ...) {

  # COMPROBACIONES INICIALES
    # esta instalada la libreria

  if (!"GA" %in% (installed.packages()[,"Package"])) install.packages("GA") else suppressPackageStartupMessages(require(GA))
  if (!"GA" %in% (installed.packages()[,"Package"])) stop("Library GA must be installed") else suppressPackageStartupMessages(require(GA))

    # matrices

  if (!is.matrix(hist)) stop("'hist' must be a matrix")
  if (!isSymmetric.matrix(hist)) stop ("'hist' is not a valid correlation matrix")

  if (!is.matrix(CD)) stop("'CD' must be a matrix")
  if (!isSymmetric.matrix(CD)) stop ("'CD' is not a valid correlation matrix")

  call = match.call()

  object = new("sensitivities",
               call = call,
               hist = hist,
               map = map,
               CD = CD,
               maxiter = maxiter,
               lim = lim)

  n = colnames(CD)
  nm = names(map)[-1]

  R_2 = function(FG = FG, FL = FL, RU = map, hist = hist, CD = CD, col = col) {

    matr = mat(FG = FG, FL = as.matrix(FL), RU = map, col = col)
    matr = matr %*% CD %*% t(matr)
    diag(matr) = 1
    if (max(matr >1)) return(1e10)

    return(sum((hist - matr) ^ 2))
  }

  fitness = function(x, lim = lim) {

    # seleccion

    y = fgyfl(x, lim = lim, n = nm)

    fit = -R_2(FG = y[,1], FL = y[,-1], RU = map, hist = hist, CD = CD, col = n)

    return(fit)

  }

  GA = ga(type = "real-valued",
          fitness = fitness,
          min = rep(0, nrow(map) * ncol(map)),
          max = rep(1, nrow(map) * ncol(map)),
          lim = lim,
          maxiter = maxiter,
          optim = T,
          parallel = parallel,
          ...)

  x = GA@solution[1,]

  x = fgyfl(x, lim = lim, n = nm)

  object@error = GA@fitnessValue
  object@iterations = GA@iter

  object@sensitivities = cbind(map[,1],
                               as.data.frame(x))

  object@suggestions = as.numeric(as.matrix(object@sensitivities[,-1]))

  return(object)

}

# Object method -----------------------------------------------------------

show.sensitivities = function(object) {

  cat("An object of class \"sensitities\"\n")
  cat("\nCall:\n", deparse(object@call), "\n\n",sep="")
  cat("Available slots:\n")
  print(slotNames(object))

}

#' Get R squared
#'
#' @export
r_squared = function(object, ...) {

  UseMethod("r_squared", object)

}


r_squared.sensitivities = function(object) {

  r2 = object@sensitivities
  # r2$r2 = r2$FG ^ 2 + r2$FL ^ 2

  r2$r2 = apply(object@sensitivities[,-1], 1, function(x) sum(x ^ 2))

  return(r2)

}

summary.sensitivities <- function(object, ...) {

  cat("+-----------------------------------+\n")
  cat("|Global & Local factor optimization |\n")
  cat("+-----------------------------------+\n\n")
  cat("Settings: \n")
  cat(paste("  lim                   = ", object@lim, "\n"))
  cat(paste("  maxiter               = ", object@maxiter, "\n\n"))
  cat("Estimation results:\n")
  cat(paste("  error                 = ", format(object@error, digits = 4, scientific = F), "\n"))
  cat(paste("  iterations            = ", object@iterations, "\n\n"))
  cat("Sensitivities:\n\n")
  print(r_squared(object), digits = 4, row.names = F, scipen=999)

}

#' Get Global & local factor sensitivities
#'
#' @export
get_sensitivities = function(object, ...) {

  UseMethod("get_sensitivities", object)

}

get_sensitivities.sensitivities = function(object) {

  return(object@sensitivities)

}

#' Get Global & local factor sensitivities
#'
#' @export
get_suggestions = function(object, ...) {

  UseMethod("get_sensitivities", object)

}

get_suggestions.sensitivities = function(object) {

  return(object@suggestions)

}

#' Fitted matrix
#'
#' @export
fitted_cor = function(object, ...) {

  UseMethod("fit_matrix", object)

}

fitted_cor.sensitivities = function(object) {

  n = colnames(object@CD)

  matr = mat(FG = object@sensitivities[,2], FL = as.matrix(object@sensitivities[,-(1:2)]), RU = object@map, col = n)
  matr = matr %*% object@CD %*% t(matr)
  diag(matr) = 1

  return(matr)
}


# Set methods -------------------------------------------------------------

setMethod("show", "sensitivities", show.sensitivities)
setMethod("print", "sensitivities", function(x, ...) str(x))
setMethod("r_squared", "sensitivities", r_squared.sensitivities)
setMethod("get_sensitivities", "sensitivities", get_sensitivities.sensitivities)
setMethod("get_suggestions", "sensitivities", get_suggestions.sensitivities)
setMethod("fitted_cor", "sensitivities", fitted_cor.sensitivities)

#' @export
setMethod("summary", signature(object = "sensitivities"), summary.sensitivities)


