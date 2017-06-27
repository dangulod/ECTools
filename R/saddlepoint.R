# Classes ----

setClass(Class = "saddlepoint",
         slots = c(
           call = "language",
           n_days = "numeric",
           data = "numeric",
           transform = "numeric",
           moments_o = "numeric",
           moments_t = "numeric",
           quantiles = "data.frame"),
         package = "ECTools")

#### MEJORA DE ESTAS FUNCIONES ####

get_prob4m = function(k,
                      s = s,
                      n_days = n_days) {

  d = n_days * k4m(k, s)
  d1 = n_days * k4m_1st(k, s)
  d2 = n_days * k4m_2nd(k, s)
  exp_term = d - d1 * s + 0.5 * s * s * d2

  cnd_tmp = pnorm(-sqrt(s * s * d2))

  if (cnd_tmp %in% 0 | is.nan(cnd_tmp)) {
    prob_tmp = 0
  } else {
    prob_tmp = exp(exp_term) * cnd_tmp
  }

  if (s < 0) {
    get_prob4m = prob_tmp
  } else {
    get_prob4m = 1 - prob_tmp
  }

  return(get_prob4m)
}

find_s4m = function(k,
                    s = s,                             # NO SE USA
                    n_days = n_days,
                    prob_root = prob_root,
                    s_tmp = s_tmp                      # NO SE USA
) {

  n_stdevs = 10
  stdev = sqrt(k[2])                                    # NO SE USA
  s_inf = -300
  s_sup = 345

  prob_inf = get_prob4m(k, s_inf, n_days) - prob_root
  prob_sup = get_prob4m(k, s_sup, n_days) - prob_root

  while((abs(prob_inf) + abs(prob_sup)) > 0.000001) {

    s_tmp_new = 0.5 * (s_inf + s_sup)
    prob_tmp = get_prob4m(k, s_tmp_new, n_days) - prob_root

    if (prob_tmp > 0) {
      s_sup = s_tmp_new
      prob_sup = prob_tmp
    } else if (prob_tmp < 0) {
      s_inf = s_tmp_new
      prob_inf = prob_tmp
    } else {
      s_tmp = s_tmp_new
    }
  }

  s_tmp = s_tmp_new
  find_s4m = prob_tmp

  return(s_tmp)

}

#' fit_saddlepoint
#'
#' @title Saddle-point aproximation
#'
#' @param x vector
#' @param N_SIMUL (Optional) Number of simulations. Default 10000
#' @param n_days (Optional) Number of periods. Default 1
#' @param nm number of momements to use. Must be 2, 4 or 6
#'
#' @export
#'
fit_saddlepoint = function(x = x,
                           N_SIMUL = 10000,
                           n_days = 1,
                           nm = 6) {

  call = match.call()

  if(!is.numeric(x)) stop("x must be numeric")

  if (!nm %in% c(2, 4, 6)) {

    stop("Number of momemts must be 2,4 or 6")

  }

  # añadir vector de momentos

  k = c(k1(x),
        k2(x),
        ifelse(nm > 2, k3(x), 0),
        ifelse(nm > 2, k4(x), 0),
        ifelse(nm > 4, k5(x), 0),
        ifelse(nm > 4, k6(x), 0))

  names(k) = paste("k", 1:6, sep = "_")

  object = new(Class = "saddlepoint",
               call = call,
               n_days = n_days,
               data = x,
               moments_o = k)

  aleat = rnorm(N_SIMUL)

  s_tmp = 0
  loss = rep(0,N_SIMUL)
  prob_tmp = pnorm(aleat)

  for (i in 1:N_SIMUL) {

    prob_at_root = find_s4m(k, s_tmp, n_days, prob_tmp[i], s_tmp)
    loss[i] = n_days * k4m_1st(k, prob_at_root)

  }

  object@transform = loss

  return(object)
}

# Methods ----

plot.saddlepoint = function(x = x) {

  object = x
  ggplot() +
    geom_density(aes(object@data,
                     color = "observed"),
                 size = 1) +
    geom_density(aes(object@transform,
                     color = "estimated"),
                 size = 1) +
    xlab(NULL) +
    theme_bw()

}

ks_test.saddlepoint = function(object, ...) {

  ks.test(object@data,
          object@transform,
          ...)

}

show.saddlepoint = function(object) {

  cat("An object of class \"saddlepoint\"\n")
  cat("\nCall:\n", deparse(object@call), "\n\n",sep="")
  cat("Available slots:\n")
  print(slotNames(object))

}



#' @export
setMethod("plot", signature(x = "saddlepoint", y = "missing"), plot.saddlepoint)
setMethod("show", "saddlepoint", show.saddlepoint)
setMethod("print", "saddlepoint", function(x, ...) str(x))
setMethod("ks_test", "saddlepoint", ks_test.saddlepoint)
setMethod("qqPlot", "saddlepoint", qqPlot.saddlepoint)
