#' Modified lag function
#'
#' @param x vector
#' @param lag if possitive retard, negative delay
#'
#' @export
#'
lg = function(x = x , lag = lag) {

  if (lag < 0) {
    lg = x %>% rev %>% lag(abs(lag)) %>% rev
  } else if (lag == 0) {
    lg = x
  } else if (lag > 0) {
    lg = x %>% lag(lag)
  }

  return(lg)
}

#' Correlations bewteen lags
#'
#' @param x vector
#' @param y vector
#' @param max.lag max.lag
#' @param allow.negative (Optional) logical, if negative lags are allowed, by default TRUE
#'
#' @return Return the correlation between x and y lags
#'
#' @export
#'
llcor = function(x = x, y = y, max.lag = 6, allow.negative = T) {

  if(allow.negative) {
    z = matrix(c(-max.lag:max.lag,
                 rep(0,
                     length(-max.lag:max.lag))),
               ncol = 2)
  } else {
    z = matrix(c(0:max.lag,
                 rep(0, length(0:max.lag))),
               ncol = 2)
  }

  for(i in 1:nrow(z)) {

    j = z[i,1]

    if (j == 0) {
      z[i,2] = cor(x,
                   y,
                   use = "pairwise.complete.obs")
    } else {
      lg1 = y %>% lg(j)
      z[i,2] = cor(x,
                   lg1,
                   use = "pairwise.complete.obs")
  }}
  return(z)

}

#' maxccf
#'
#' @param x vector
#' @param y vector
#' @param lag.max (Optional) maximum lag
#' @param allow.negative (Optional) logical, if negative lags are allowed, by default TRUE
#' @param max (Optional) logical, Should maximum be returned, in other case minimum is returned
#' @param abs (Optional) logical, Should be the maximum in absolute value, by defaukt TRUE
#'
#' @return Return the lag with the maximum/minimum correlations bewteen two time series
#'
#' if x is a data frame or matrix, u must use maxccfdf fucntion
#'
#' @export
#'
maxccf = function(x = x, y = y, lag.max = 6, allow.negative = T, max = T, abs = T) {

  ccf = llcor(x,
              y,
              max.lag = lag.max,
              allow.negative = allow.negative) %>%
    as.data.frame()


  # ccf = ccf[order(abs(ccf$V2), decreasing = T),]
  # ccf = ccf[order(ccf$V2,
  #                 decreasing = max),]
  maxcorlag = if(abs) ccf[which.max(abs(ccf$V2)),"V1"] else ccf[which.max(ccf$V2),"V1"]

  return(maxcorlag)
}

### CONSIDERAR UNA FUNCIÓN PARA SACAR LA CORRELACIÓN MAXIMA TMB ADEMAS DEL LAG ###

#' maxccfdf
#'
#' @param x data.frame with data
#' @param lag.max maximum lag
#' @param allow.negative (Optional) logical, if negative lags are allowed, by default TRUE
#' @param max (Optional) logical, Should maximum be returned, in other case minimum is returned
#'
#' @return matrix with maxccf
#'
#' @export
#'
maxccfdf = function(x = x, lag.max = 6, allow.negative = T, max = T, abs = T) {

  x = apply(x, 2, function (x) as.numeric(x)) %>% as.data.frame()

  if (is.data.frame(x))
    x <- as.matrix(x)

  ncy <- ncx <- ncol(x)
  if (ncx == 0)
    stop("'x' is empty")
  r <- matrix(0,
              nrow = ncx,
              ncol = ncy)
  for (i in seq_len(ncx)) {
    for (j in seq_len(ncx)) {
      x2 <- x[, i]
      y2 <- x[, j]
      ok <- complete.cases(x2, y2)
      # x2 <- x2[ok]
      # y2 <- y2[ok]
      r[i, j] <- ifelse(any(ok), maxccf(x2,
                                        y2,
                                        lag.max = lag.max,
                                        allow.negative = allow.negative,
                                        max = max,
                                        abs = abs),
                        NA_real_)
    }
  }
  rownames(r) <- colnames(x)
  colnames(r) <- colnames(x)
  return(r)

}

#' corscale
#'
#' @param DATA data.frame with the data
#' @param ref a string with the reference variable
#' @param variables variables to change
#' @param lag.max maximum number of lags to asses
#' @param allow.negative (Optional) logical, if negative lags are allowed, by default TRUE
#' @param max (Optional) logical, Should maximum be returned, in other case minimum is returned
#'
#' @return Return a dataframe with the variables lagged or delayed with the max correlation to reference variable
#' @export
#'
corscale = function(DATA = DATA, ref = ref, variables = variables, lag.max = 6, allow.negative = T, max = T, abs = T) {

  if(!is.data.frame(DATA)) {
    stop("DATA must be a data.frame")
  }

  if(!ref %in% names(DATA)) {
    stop(paste(ref, "is not a valid variable"))
  }

  variables = variables[variables != ref]

  rf = DATA %>% .[[ref]]

  x = list()

  for(i in names(DATA)) {

    if(i %in% variables) {
      j = maxccf(rf,
                 DATA %>% .[[i]],
                 lag.max = lag.max,
                 allow.negative = allow.negative,
                 max = max,
                 abs = abs)
      x[[i]] =  DATA %>% .[[i]] %>% lg(j)

    } else {
      x[[i]]  = DATA %>% .[[i]]
    }
  }
  x = bind_cols(x)
  # x = do.call(cbind, x) %>% as.data.frame()
  # names(x) = names(DATA)
  return(x)
}

#' matrix_hclust
#'
#' @param groups vector of groups of the cluster
#'
#' @export
#'
matrix_hclust = function(groups = groups) {

  m = matrix(0,
             nrow = groups %>% length() ,
             ncol = groups %>% length())

  rownames(m) = groups %>% names
  colnames(m) = groups %>% names
  i =  1; j = 1;
  for (i in seq(nrow(m))) {
    for (j in seq(ncol(m))) {
      if (i != j) {
        x1 = rownames(m)[i]
        y1 = colnames(m)[j]
        x2 = groups[x1]
        y2 = groups[y1]
        if (x2 == 0 | y2 == 0) {
          m[i, j] = 1
        } else {
          m[i, j] = as.numeric(x2 == y2)
        }
      } else {
        m[i, j] = 0
      }
    }
  }

  return(m)

}

#' cutcor
#'
#' @param data dataframe
#' @param FG character
#'
#' @return return a dataframe with more correlation between the global factor and the distance to default vector
#' @export
#'
cutcor = function(data = data, FG = FG) {

  if (!(is.numeric(data) | is.data.frame(data))) stop("data should be a vector or a dataframe")

  if (!is.character(FG)) stop("ref should be a character")

  x = lapply(data, function(x) x)

  y = x[[FG]]

  cuti = function(xx = xx, FG = FG) {

    if (!is.numeric(xx)) return(xx)

    if (identical(xx, FG)) return(xx)

    l = length(na.omit(xx))

    if (l < 16) return(xx)# Four year 4 * 4

    z = data.frame(
      x = xx,
      y = FG
    )

    i = 0
    z = z %>% na.omit()
    while (cor(z$x[-1], z$y[-1], use =  "pairwise.complete.obs") < cor(z$x, z$y, use =  "pairwise.complete.obs")) {

      z$x[1] = NA
      z$y[1] = NA
      z = z %>% na.omit()
      i = i + 1
    }
    if (is.na(xx[1])) {

      s = min(which(!is.na(xx)))

    } else {

      s = 1

    }

    if ((l - i) < 16 & i > 0) i = l - 16

    if (i > 0) {

      xx[s:(s + i - 1)] = NA

    }

    return(xx)}

  x = lapply(x, function(x) cuti(x, y))
  x = bind_cols(x)
  return(x)

}
