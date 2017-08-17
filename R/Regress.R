#' Get the residuals from a regression
#'
#' @param y dataframe or vector
#' @param x.var character variable, representing the dependent variable
#' @param scale logical, it the result should be standardie
#'
#' @return This function perform regressions for every variable inside 'y',
#' 'y[i] ~ x.var' and get the residuals 'y[i] - b0 * b1 x.var'
#'
#' In this way we get orthogonal vectors to 'x.var'
#'
#' @export
#'
flresi = function(y = y, x.var = x.var, scale = T) {

  if(!is.logical(scale)) stop("scale must be logical")

  if (!is.data.frame(y)) {

    flresi = flresid(Y = y, X = x.var, scale = scale)

  } else {

    FG = y %>% .[[x.var]]

    y = y %>%
      lapply(function(z)
        if(!is.numeric(z)) {
          z
        } else {
          flresid(Y = z, X = FG, scale = scale)})

    y[[x.var]] = FG

    flresi = bind_cols(y)
  }

  return(flresi)
}

flresid = function(Y = Y, X = X, scale = scale) {

  lm = lm(Y ~ X )
  flresid = Y - (lm$coefficients[1] + lm$coefficients[2] * X)
  if (scale) flresid = as.numeric(scale(flresid))

  return(flresid)
}


#' Optimization weights which maximize the correlation between two matrix
#'
#' @param x Vector of factors
#' @param y Vector of local factors
#'
#' @return weights which maximize the correlation between two matrix
#' @export
#'
#' @examples
#' library(ECTools)
#' X = data.frame(
#'  x1 = rnorm(100),
#'  x2 = rnorm(100),
#'  x3 = rnorm(100)
#'  )
#'
#'Y = data.frame(
#'  y1 = rnorm(100),
#'  y2 = rnorm(100)
#'  )
#'
# weopt(X, Y)
weopt = function(x = x, y = y) {

  NlcOptim = getFromNamespace("NlcOptim", "NlcOptim")

  mcor = function(b = rep(1, length(y))) {

    y = y %>% as.matrix()
    b = b %>% as.matrix()
    FF = y %*% b

    mcor = (apply(x, 2, function(x) cor(FF, x, use = "pairwise.complete.obs")) %>% mean)
    return(mcor)
  }

  confun=function(x) {
    f = NULL
    f = (x %>% abs %>% sum %>% abs- 1)
    return(list(ceq=f,c=NULL))
  }

  opt = NlcOptim(rep(1, length(y)), objfun=mcor ,confun=confun)
  weopt = opt$p %>% as.numeric
  return(weopt)
}


