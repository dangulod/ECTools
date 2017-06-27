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
#'
#' @export
#'
cor_optim = function(map = map, hist = hist, CD = CD, lim = 1, maxiter = 1e4, parallel = F) {

  # COMPROBACIONES INICIALES
    # esta instalada la libreria

  if (!"GA" %in% (installed.packages()[,"Package"])) install.packages("GA") else suppressPackageStartupMessages(require(GA))
  if (!"GA" %in% (installed.packages()[,"Package"])) stop("Library GA must be installed") else suppressPackageStartupMessages(require(GA))

    # matrices

  if (!is.matrix(hist)) stop("'hist' must be a matrix")
  if (!isSymmetric.matrix(hist)) stop ("'hist' is not a valid correlation matrix")

  if (!is.matrix(CD)) stop("'CD' must be a matrix")
  if (!isSymmetric.matrix(CD)) stop ("'CD' is not a valid correlation matrix")

  n = colnames(CD)

  R_2 = function(FG = FG, FL = FL, RU = map$cor_cd, hist = hist, CD = CD) {

    mat = mat(FG = FG, FL = FL, RU = map$cor_cd, col = n)
    mat = mat %*% CD %*% t(mat)
    diag(mat) = 1

    sum((hist - mat) ^ 2)
  }

  fitness = function(x, lim = lim) {

    # seleccion

    y = fgyfl(x, lim = lim)

    fit = -R_2(FG = y$FG, FL = y$FL, RU = map$cor_cd, hist = hist, CD = CD)

    return(fit)

  }

  GA = ga(type = "real-valued",
          fitness = fitness,
          min = rep(0, length(map$cor_cd) * 2),
          max = rep(1, length(map$cor_cd) * 2),
          lim = lim,
          maxiter = maxiter,
          optim = T,
          parallel = parallel)

  x = GA@solution[1,]

  x = fgyfl(x, lim = lim)

  cor_optim = data.frame(
    RU = map$cor_his,
    FG = x$FG,
    FL = x$FL
  ) %>%
    unique()

  return(cor_optim)

}
