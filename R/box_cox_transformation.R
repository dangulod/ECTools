box.cox <- function(s, lambda = 0) {

  if(!is.numeric(s)) {
    bc = s
  } else {
    if(lambda == 0) {
      bc = -log(s)
    } else {
      bc = -((s ^ lambda) - 1) / lambda}
  }
  return(bc)
}

loglik.box = function(x, lambda) {

  x_box = box.cox(x, lambda)
  x_box = (x_box - mean(x_box, na.rm = T)) / sd(x_box, na.rm = T)
  x_box = x_box[order(x_box)]

  max(abs(pnorm(x_box) - ((row_number(x_box) - 1) / length(x_box[!is.na(x_box)]))),
      abs(((row_number(x_box)) / length(x_box[!is.na(x_box)])) - pnorm(x_box)),
      na.rm =  T)

}

#' Fit the lambda which optimize the loglikehood function
#'
#' @param x vector or dataframe with data
#' @param interval a vector containing the end-points of the interval to be searched for the minimum., default c(-20,20)
#' @param tol the desired accuracy, default = 0.00001
#'
#' @return Return the lambda parameter which minimizes the loglikehood function
#' @export
#'
lambda.box = function(x = x, interval = c(-20,20), tol = 0.00001) {

  if (is.data.frame(x)) {
    x = lapply(x, function(x) x)
  }

  if (is.list(x)) {

    y = lapply(x, function(x)
      if (is.numeric(x)) {
        optimize(loglik.box,
                 interval,
                 tol = tol,
                 x = x)$minimum
      } else {NA})

    if (length(y) > 1) {

      z = data.frame(
        variable = names(y),
        lambda = as.vector(unlist(y))
      )
    }} else {
      if (is.numeric(x)) {
        z = optimize(loglik.box,
                     interval,
                     tol = tol,
                     x = x)$minimum
      } else {z = NA}
    }

  return(z)
}

#' boxcox transformation and evaluation
#'
#' @param x vector
#' @param lambda (Optional) lambda, if NULL, lambda is optimized
#' @param plot (Optional) Logical, plot qqplots
#' @param ks.test (Optional) Logical, kolmogorov-smirnov test
#'
#' @return Compute the box-cox transformation based of a vector. Plot and ks.test argument can be set to TRUE in order to asses the normality of the new vector
#' @export
#'
boxcox = function(x = x, lambda = NULL, plot = F, ks.test = F, ...) {

  if (is.data.frame(x)) {stop("x must be a vector")}

  if (missing("lambda")) {
    lambda = lambda.box(x)}

  boxcox = list()

  # x = abs(qnorm(x))

  boxcox$x.box = -box.cox(x,
                           lambda = lambda)

  if (plot == T) {

    if ("car" %in% rownames(installed.packages())) {
      qqPlot = getFromNamespace("qqPlot", "car")

      par(mfrow = c(1, 2))

      qqPlot(x,
             dist="norm",
             pch = 20, main = "Original", ylab = "")

      qqPlot(boxcox$x.box,
             dist="norm",
             pch = 20, main = paste("Transformed l = ", round(lambda, digits = 2)), ylab = "")
    } else {

      message("Warning: you must install 'car' package if 'plot = T'")

    }
  }

  if (ks.test == T) {

    kstest = list()

    kstest$Dtd = ks.test(x,
                         "pnorm",
                         mean = mean(x, na.rm = T),
                         sd = sd(x, na.rm = T), ...)

    kstest$transformed = ks.test(boxcox$x.box,
                                 "pnorm",
                                 mean = mean(boxcox$x.box, na.rm = T),
                                 sd = sd(boxcox$x.box, na.rm = T), ...)

    nombre = c("statistic", "p-value")
    Dtd = c(kstest$Dtd$statistic %>% as.numeric,
            kstest$Dtd$p.value %>% as.numeric)

    transformed = c(kstest$transformed$statistic %>% as.numeric,
                    kstest$transformed$p.value %>% as.numeric)

    boxcox$test = rbind(Dtd, transformed) %>% as.data.frame()

    names(boxcox$test) = nombre
  }

  if (plot == F & ks.test == F) {
    return(boxcox$x.box)
  } else if (ks.test == T) {
    return(boxcox[-1])}

}


#' Distance to default
#'
#' @param x Frequency Default Observation
#'
#' @return Distance to Default
#' @export
#'
DtD = function(x) {

  if (is.data.frame(x)) {
    DtD = lapply(x, function(x) if(is.numeric(x)) {-abs(qnorm(x))} else {x})
    DtD = bind_cols(DtD)
  } else {DtD = -abs(qnorm(x))
  }
  return(DtD)
}

#' scaledf
#'
#' @param x dataframe
#' @param center either a logical value or a numeric vector of length equal to the number of columns of x.
#' @param scale either a logical value or a numeric vector of length equal to the number of columns of x
#'
#' @return The value of center determines how column centering is performed. If center is a numeric vector with length equal to the number of columns of x, then each column of x has the corresponding value from center subtracted from it. If center is TRUE then centering is done by subtracting the column means (omitting NAs) of x from their corresponding columns, and if center is FALSE, no centering is done.
#'
#' The value of scale determines how column scaling is performed (after centering). If scale is a numeric vector with length equal to the number of columns of x, then each column of x is divided by the corresponding value from scale. If scale is TRUE then scaling is done by dividing the (centered) columns of x by their standard deviations if center is TRUE, and the root mean square otherwise. If scale is FALSE, no scaling is done.
#'
#' @export
#'
scaledf = function(x, center = TRUE, scale = TRUE) {

  x = lapply(x, function(x) x)
  x = lapply(x, function(x)
    if (is.numeric(x)) {
      as.numeric(scale(x, center = center, scale = T))
    } else {
      x
    })
  scaledf = bind_cols(x)
  return(scaledf)
}
