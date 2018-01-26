# Class definition -------------------------------------------------------------------

setClass(Class = "sensitivities",
         slots = c(
           call = "language",
           hist = "matrix",
           map = "data.frame",
           lim = "numeric",
           minFG = "numeric",
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
#' @param lim Numeric, limit of the sum squared of the factors, by default 1
#' @param maxiter Numeric, the maximum number of iterations to run
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
#' summary(x)
#'
#' fitted_cor(x)
#'
#' get_sensitivities(x)
#'
#' get_suggestions(x)
#'
#' r_squared(x)
#'
#' @export
cor_optim = function(hist, ...) UseMethod("cor_optim", hist)

cor_optim.matrix = function(hist, ...) coroptim(hist, ...)

cor_optim.list = function(hist) {

  return(lapply(hist, function(x) cor_optim(x)))

}

cor_optim.sensitivities = function(hist, ...) {

  coroptim(hist = hist@hist,
           CD = hist@CD,
           map = hist@map,
           lim = hist@lim,
           maxiter = hist@maxiter,
           minFG = hist@minFG,
           suggestions = get_suggestions(hist),
           ...)
}

coroptim = function(hist = hist, CD = CD, map = map, lim = 1, maxiter = 1e4, parallel = F,  minFG = 0, seed = 1, ...) {

  # COMPROBACIONES INICIALES
    # esta instalada la libreria

  if (!"GA" %in% (installed.packages()[,"Package"])) install.packages("GA") else suppressPackageStartupMessages(require(GA))
  if (!"GA" %in% (installed.packages()[,"Package"])) stop("Library GA must be installed") else suppressPackageStartupMessages(require(GA))

    # matrices

  if (!is.matrix(hist)) stop("'hist' must be a matrix")
  if (any(hist != t(hist))) stop ("'hist' is not a valid correlation matrix")

  if (!is.matrix(CD)) stop("'CD' must be a matrix")
  if (any(hist != t(hist))) stop ("'CD' is not a valid correlation matrix")

  call = match.call()

  object = new("sensitivities",
               call = call,
               hist = hist,
               map = map,
               CD = CD,
               maxiter = maxiter,
               minFG = minFG,
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

  if (ncol(map) > 2 & minFG != 0) {stop("ERROR: this option is not allowed")} else {maxFL = sqrt(lim - (minFG ^ 2))}

  GA = ga(type = "real-valued",
          fitness = fitness,
          min = c(rep(minFG, nrow(map)), rep(0, nrow(map) * (ncol(map) - 1))), #  c(rep(minFG, nrow(map)), rep(minFG, nrow(map) * (ncol(map) - 1)))         rep(minFG, nrow(map) * ncol(map))
          max = c(rep(1, nrow(map)), rep(maxFL, nrow(map) * (ncol(map) - 1))),
          lim = lim,
          maxiter = maxiter,
          optim = T,
          parallel = parallel,
          seed = seed,
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

r_squared.list = function(object) {

  return(lapply(object, function(x) r_squared(x)))

}

summary.sensitivities <- function(object, ...) {

  cat("+-----------------------------------+\n")
  cat("|Global & Local factor optimization |\n")
  cat("+-----------------------------------+\n\n")
  cat("Settings: \n")
  cat(paste("  lim                   = ", object@lim, "\n"))
  cat(paste("  maxiter               = ", object@maxiter, "\n"))
  cat(paste("  minFG                 = ", object@minFG, "\n\n"))
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

get_sensitivities.list = function(object) {

  return(lapply(object, function(x) get_sensitivities(x)))

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

get_suggestions.list = function(object) {

  return(lapply(object, function(x) get_suggestions(x)))

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

fitted_cor.list = function(object) {

  return(lapply(object, function(x) fitted_cor(x)))

}

#' Get equation
#'
#' @export
get_equation = function(object, ...) {

  UseMethod("get_equation", object)

}

get_equation.sensitivities = function(object) {

  if (ncol(object@map) > 2) return("Not supported")

  rows = object@map$RU
  cols = colnames(object@CD)
  fl = as.data.frame(object@map)[,2]

  equ = matrix(0,
               nrow = length(rows),
               ncol = length(cols),
               dimnames = list(
                 rows,
                 cols
               ))

  equ[, 1] = object@sensitivities[,2]

  for (i in 1:length(rows)) {
    equ[i, fl[i]] = object@sensitivities[i,3]
  }

  equ = rownames_to_column(as.data.frame(equ), var = "Equation")

  return(equ)
}

get_equation.list = function(object) {

  return(bind_rows(lapply(object, function(x) get_equation(x))))

}

# Set methods -------------------------------------------------------------

setMethod("show", "sensitivities", show.sensitivities)
setMethod("print", "sensitivities", function(x, ...) str(x))
setMethod("r_squared", "sensitivities", r_squared.sensitivities)
setMethod("r_squared", "list", r_squared.list)
setMethod("get_sensitivities", "sensitivities", get_sensitivities.sensitivities)
setMethod("get_sensitivities", "list", get_sensitivities.list)
setMethod("get_suggestions", "sensitivities", get_suggestions.sensitivities)
setMethod("get_suggestions", "list", get_suggestions.list)
setMethod("fitted_cor", "sensitivities", fitted_cor.sensitivities)
setMethod("fitted_cor", "list", fitted_cor.list)

setMethod("get_equation", "sensitivities", get_equation.sensitivities)
setMethod("get_equation", "list", get_equation.list)


#' @export
setMethod("summary", signature(object = "sensitivities"), summary.sensitivities)

setMethod("cor_optim", "matrix", cor_optim.matrix)
setMethod("cor_optim", "list", cor_optim.list)
setMethod("cor_optim", "sensitivities", cor_optim.sensitivities)
