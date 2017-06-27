#' @rdname skewn
#'
#' @name Skew-Normal Distribution
#'
#' @title  Density  function,  distribution  function,  quantiles  and  random  number  generation  for  the  skew normal (SN) and the extended skew-normal (ESN) distribution.
#'
#' @param x vector of quantiles. Missing values (NA’s) andInf’s are allowed.
#' @param p vector of probabilities. Missing values (NAs) are allowed
#' @param xi vector of location parameters.
#' @param omega vector of scale parameters; must be positive
#' @param alpha vector of slant parameter(s);+/- Infis allowed. Withpsn, it must be of length 1 ifengine="T.Owen". With qsn, it must be of length 1
#' @param tau a single value representing the ‘hidden mean’ parameter of theESNdistribution; tau=0 (default) corresponds to a SN distribution.
#' @param dp a vector of length 3 (in the SN case) or 4 (in the ESN case), whose components represent the individual parameters described above. If dp is specified, the individual parameters cannot be set.
#' @param n a positive integer representing the sample size.
#' @param tol a scalar value which regulates the accuracy of the result of qsn, measured on the probability scale.
#' @param log logical flag used in dsn (defaul tFALSE). When TRUE, the logarithm of the density values is returned.
#' @param engine a character string which selects the computing engine; this is either "T.Owen" or "biv.nt.prob", the latter from package mnormt. If tau != 0 or length(alpha)>1, "biv.nt.prob" must be used.
#' If this argument is missing, a default selection rule is applied.
#' @param solver a character string which selects the numerical method used for solving the quan-tile equation; possible options are "NR" (default) and "RFB", described in the ‘Details’ section
#' @param ... additional parameters passed to T.Owen
#'
#' @return density (dsn), probability (psn), quantile (qsn) or random sample (rsn) from the skew-normal dis-tribution with given xi,omega and alpha parameters or from the extended skew-normal if tau!=0
#'
#' @details
#'
#' psn and qsn make use of function T.Owen or biv.nt.prob
#'
#' In qsn, the choice solver="NR" selects the Newton-Raphson method for solving the quantile equation, while optionsolver="RFB" alternates a step of regula falsi with one of bisection.
#' The "NR" method is generally more efficient, but "RFB" is occasionally required in some problematic cases.
#'
#' @export
#'
#' @references Azzalini, A. (1985). A class of distributions which includes the normal ones.Scand. J. Statist.12,171-178.
#'
#' Azzalini, A. with the collaboration of Capitanio, A. (2014).The Skew-Normal and Related Families.Cambridge University Press, IMS Monographs series.
#'
#' @examples
#'
#' pdf <- dsn(seq(-3, 3, by=0.1), alpha=3)
#' cdf <- psn(seq(-3, 3, by=0.1), alpha=3)
#' q <- qsn(seq(0.1, 0.9, by=0.1), alpha=-2)
#' r <- rsn(100, 5, 2, 5)
#' qsn(1/10^(1:4), 0, 1, 5, 3, solver="RFB")
#'
dsn = function(x, xi=0, omega=1, alpha=0, tau=0, dp=NULL, log=FALSE)
{
  if(!is.null(dp)) {
    if(!missing(alpha))
      stop("You cannot set both 'dp' and component parameters")
    xi = dp[1]
    omega = dp[2]
    alpha = dp[3]
    tau = if(length(dp)>3) dp[4] else 0
  }
  z = (x-xi)/omega
  logN = (-log(sqrt(2*pi)) -logb(omega) -z^2/2)
  if(abs(alpha) < Inf)
    logS = pnorm(tau * sqrt(1+alpha^2) + alpha*z, log.p=TRUE)
  else
    logS  = log(as.numeric(sign(alpha)*z + tau > 0))
  logPDF = as.numeric(logN + logS - pnorm(tau, log.p=TRUE))
  logPDF = replace(logPDF, abs(x) == Inf, -Inf)
  logPDF = replace(logPDF, omega <= 0, NaN)
  if(log) logPDF else exp(logPDF)
}

#' @rdname skewn
#'
#' @export
#'
psn = function(x, xi=0, omega=1, alpha=0, tau=0, dp=NULL, engine, ...)
{
  if(!is.null(dp)) {
    if(!missing(alpha))
      stop("You cannot set both 'dp' and component parameters")
    xi = dp[1]
    omega = dp[2]
    alpha = dp[3]
    tau = if(length(dp)>3) dp[4] else 0L
  }
  z = as.numeric((x-xi)/omega)
  nz = length(z)
  na = length(alpha)
  if(missing(engine)) engine =
    if(na == 1 & nz > 3 & all(alpha*z > -5) & (tau == 0L))
      "T.Owen" else "biv.nt.prob"
  if(engine == "T.Owen") {
    if(tau != 0 | na > 1)
      stop("engine='T.Owen' not compatible with other arguments")
    p = pnorm(z) - 2 * T.Owen(z, alpha, ...)
  }
  else{ #  engine="biv.nt.prob"
    p = numeric(nz)
    alpha = cbind(z, alpha)[,2]
    delta = delta.etc(alpha)
    p.tau = pnorm(tau)
    for(k in seq_len(nz)) {
      if(abs(alpha[k]) == Inf){
        p[k] = if(alpha[k] > 0)
          (pnorm(pmax(z[k],-tau)) - pnorm(-tau))/p.tau
        else
          1- (pnorm(tau) - pnorm(pmin(z[k], tau)))/p.tau
      }
      else { # SNbook: formula (2.48), p.40
        R = matrix(c(1, -delta[k], -delta[k], 1), 2, 2)
        p[k]= mnormt::biv.nt.prob(0, rep(-Inf,2), c(z[k], tau), c(0, 0), R)/p.tau
      }
    }}
  p = pmin(1, pmax(0, as.numeric(p)))
  replace(p, omega <= 0, NaN)
}

#' @rdname skewn
#'
#' @export
#'
qsn = function(p, xi = 0, omega = 1, alpha = 0, tau=0, dp=NULL, tol = 1e-08,
               solver="NR", ...)
{ if(!is.null(dp)) {
  if(!missing(alpha))
    stop("You cannot set both 'dp' and component parameters")
  xi = dp[1]
  omega = dp[2]
  alpha = dp[3]
  tau = if(length(dp) > 3) dp[4] else 0
}
  max.q = sqrt(qchisq(p,1)) + tau
  min.q = -sqrt(qchisq(1-p,1)) + tau
  if(tau == 0) {
    if(alpha == Inf)  return(as.numeric(xi + omega * max.q))
    if(alpha == -Inf) return(as.numeric(xi + omega * min.q))
  }
  na = is.na(p) | (p < 0) | (p > 1)
  zero = (p == 0)
  one = (p == 1)
  p = replace(p, (na | zero | one), 0.5)
  dp0 = c(0, 1, alpha, tau)
  if(solver == "NR") {
    dp0 = c(0, 1, alpha, tau)
    cum = sn.cumulants(dp=dp0, n=4)
    g1 = cum[3]/cum[2]^(3/2)
    g2 = cum[4]/cum[2]^2
    x = qnorm(p)
    x = (x + (x^2 - 1) * g1/6 + x * (x^2 - 3) * g2/24 -
           x * (2 * x^2 - 5) * g1^2/36)
    x = cum[1] + sqrt(cum[2]) * x
    px = psn(x, dp=dp0, ...)
    max.err = 1
    while (max.err > tol) { # cat("qsn:", x, "\n")
      # cat('x, px:', format(c(x,px)),"\n")
      x1 = x - (px - p)/dsn(x, dp=dp0)
      # x1 = pmin(x1,max.q)
      # x1 = pmax(x1,min.q)
      x = x1
      px = psn(x, dp=dp0, ...)
      max.err = max(abs(px-p))
      if(is.na(max.err)) stop('failed convergence, try with solver="RFB"')
    }
    x = replace(x, na, NA)
    x = replace(x, zero, -Inf)
    x = replace(x, one, Inf)
    q = as.numeric(xi + omega * x)
  } else { if(solver == "RFB") {
    abs.alpha = abs(alpha)
    if(alpha < 0) p = (1-p)
    x = xa = xb = xc = fa = fb = fc = rep(NA, length(p))
    nc = rep(TRUE, length(p)) # not converged (yet)
    nc[(na| zero| one)] = FALSE
    fc[!nc] = 0
    xa[nc] = qnorm(p[nc])
    xb[nc] = sqrt(qchisq(p[nc], 1)) + abs(tau)
    fa[nc] = psn(xa[nc], 0, 1, abs.alpha, tau, ...) - p[nc]
    fb[nc] = psn(xb[nc], 0, 1, abs.alpha, tau, ...) - p[nc]
    regula.falsi = FALSE
    while (sum(nc) > 0) { # alternate regula falsi/bisection
      xc[nc] = if(regula.falsi)
        xb[nc] - fb[nc] * (xb[nc] - xa[nc])/(fb[nc] - fa[nc])    else
          (xb[nc] + xa[nc])/2
      fc[nc] = psn(xc[nc], 0, 1, abs.alpha, tau, ...) - p[nc]
      pos = (fc[nc] > 0)
      xa[nc][!pos] = xc[nc][!pos]
      fa[nc][!pos] = fc[nc][!pos]
      xb[nc][pos] = xc[nc][pos]
      fb[nc][pos] = fc[nc][pos]
      x[nc] = xc[nc]
      nc[(abs(fc) < tol)] = FALSE
      regula.falsi = !regula.falsi
    }
    # x = replace(x, na, NA)
    x = replace(x, zero, -Inf)
    x = replace(x, one, Inf)
    Sign = function(x) sign(x) + as.numeric(x==0)
    q = as.numeric(xi + omega * Sign(alpha)* x)
  } else stop("unknown solver")}
  names(q) = names(p)
  return(q)
}

#' @rdname skewn
#'
#' @export
#'
rsn = function(n=1, xi=0, omega=1, alpha=0, tau=0, dp=NULL)
{
  if(!is.null(dp)) {
    if(!missing(alpha))
      stop("You cannot set both 'dp' and component parameters")
    xi = uname(dp[1])
    omega = unname(dp[2])
    alpha = unname(dp[3])
    tau = if(length(dp)>3) dp[4] else 0
  }
  if(tau == 0) {
    u1 = rnorm(n)
    u2 = rnorm(n)
    id = (u2 > alpha*u1)
    u1[id] = (-u1[id])
    z = u1
  }
  else { # for ESN use transformation method
    delta = alpha/sqrt(1+alpha^2)
    truncN = qnorm(runif(n, min=pnorm(-tau), max=1))
    z = delta * truncN + sqrt(1-delta^2) * rnorm(n)
  }
  y = xi+omega*z
  return(y)
}


#' @rdname skewt
#'
#' @name Skew-t Distribution
#'
#' @title Density function,  distribution function,  quantiles and random number generation for the skew-t (ST) distribution
#'
#' @param x vector of quantiles. Missing values (NAs) are allowed.
#' @param p vector of probabililities
#' @param xi vector of location parameters.
#' @param omega vector of scale parameters; must be positive.
#' @param alpha vector of slant parameters. Withpstandqst, it must be of length 1.
#' @param nu a single positive value representing the degrees of freedom; it can be non-integer. Default value isnu=Infwhich corresponds to the skew-normal distribution.
#' @param dp a vector of length 4, whose elements represent location, scale (positive), slant and degrees of freedom, respectively.  Ifdpis specified, the individual parameters cannot be set.
#' @param n a positive integer representing the sample size.
#' @param log logical; if TRUE, densities are given as log-densities
#' @param tol a scalar value which regulates the accuracy of the result of qsn, measured on the probability scale.
#' @param method an integer value between 0 and 4 which selects the computing method; see ‘Details’ below for the meaning of these values. If method=0 (default value), an automatic choice is
#' made among the four actual computing methods, which depends on the other arguments.
#' @param ... aditional parameters passed tointegrateor pst
#'
#' @return Density (dst), probability (pst), quantiles (qst) and random sample (rst) from the skew-t distribution with given xi, omega, alpha and nu parameters.
#'
#' @details For evaluation of pst, and so indirectly of qst, four different methods are employed. Method 1 consists in using pmst with dimensiond=1. Method 2 applies integrate to the density function dst.
#' Method 3 again uses integrate too but with a different integrand, as given in Section 4.2 of Azzalini & Capitanio (2003), full version of the paper. Method 4 consists in the recursive procedure of
#' Jamalizadeh, Khosravi and Balakrishnan (2009), which is recalled in Complement 4.3 on Azza-lini & Capitanio (2014); the recursion over nu starts from the explicit expression fornu=1 given by psc.
#' Of these, Method 1 and 4 are only suitable for integer values of nu. Method 4 becomes pro-gressively less efficient as nu increases, because its value corresponds to the number of nested calls,
#' but the decay of efficiency is slower for larger values oflength(x). If the default argument value method=0 is retained, an automatic choice among the above four methods is made, which dependson
#' the values of nu, alpha, length(x). The numerical accuracy of methods 1, 2 and 3 can be regulated via the ... argument, while method 4 is conceptually exact, up to machine precision.
#'
#' If qst is called withnu>1e4, computation is transferred to qsn.
#'
#' @references Azzalini, A. and Capitanio, A. (2003).  Distributions generated by perturbation of symmetry with emphasis on a multivariate skew-tdistribution.J.Roy. Statist. Soc. B65, 367–389. Full version of
#' the paper athttp://arXiv.org/abs/0911.2342.
#'
#' Azzalini, A. with the collaboration of Capitanio, A. (2014).The Skew-normal and Related Families. Cambridge University Press, IMS Monographs series.
#'
#' Jamalizadeh, A., Khosravi, M., and Balakrishnan, N. (2009). Recurrence relations for distributions of a skew-$t$ and a linear combination of order statistics from a bivariate-$t$.Comp. Statist. Data An.53, 847–852
#'
#' @export
#'
#' @examples
#'
#' pdf <- dst(seq(-4, 4, by=0.1), alpha=3, nu=5)
#' rnd <- rst(100, 5, 2, -5, 8)
#' q <- qst(c(0.25, 0.50, 0.75), alpha=3, nu=5)
#' pst(q, alpha=3, nu=5)  # must give back c(0.25, 0.50, 0.75)
#' p1 <- pst(x=seq(-3,3, by=1), dp=c(0,1,pi, 3.5))
#' p2 <- pst(x=seq(-3,3, by=1), dp=c(0,1,pi, 3.5), method=2, rel.tol=1e-9)
#'
dst =  function (x, xi=0, omega=1, alpha=0, nu=Inf, dp=NULL, log=FALSE)
{
  if(!is.null(dp)) {
    if(!missing(alpha))
      stop("You cannot set both component parameters and dp")
    xi = dp[1]
    omega = dp[2]
    alpha = dp[3]
    nu = dp[4]
  }
  if (nu == Inf) return(dsn(x, xi, omega, alpha, log=log))
  if (nu == 1) return(dsc(x, xi, omega, alpha, log=log))
  z   = (x - xi)/omega
  pdf = dt(z, df=nu, log=log)
  cdf = pt(alpha*z*sqrt((nu+1)/(z^2+nu)), df=nu+1, log.p=log)
  if(log)
    logb(2) + pdf + cdf -logb(omega)
  else
    2 * pdf * cdf / omega
}


#' @rdname skewt
#'
#' @export
#'
rst = function (n=1, xi = 0, omega = 1, alpha = 0, nu=Inf, dp=NULL)
{
  if(!is.null(dp)) {
    if(!missing(alpha))
      stop("You cannot set both component parameters and dp")
    xi = unname(dp[1])
    omega = unname(dp[2])
    alpha = unname(dp[3])
    nu = unname(dp[4])
  }
  z = rsn(n, 0, omega, alpha)
  if(nu < Inf) {
    v = rchisq(n,nu)/nu
    y = z/sqrt(v) + xi
  }
  else y = z+xi
  return(y)
}

#' @rdname skewt
#'
#' @export
#'
pst = function (x, xi=0, omega=1, alpha=0, nu=Inf, dp=NULL, method=0, ...)
{
  if(!is.null(dp)) {
    if(!missing(alpha))
      stop("You cannot set both component parameters and dp")
    xi = dp[1]
    omega = dp[2]
    alpha = dp[3]
    nu = dp[4]
  }
  if(length(alpha) > 1) stop("'alpha' must be a single value")
  if(length(nu) > 1) stop("'nu' must be a single value")
  if (nu <= 0) stop("nu must be non-negative")
  if (nu == Inf) return(psn(x, xi, omega, alpha))
  if (nu == 1) return(psc(x, xi, omega, alpha))
  int.nu = (round(nu) == nu)
  if((method == 1 | method ==4) & !int.nu)
    stop("selected method does not work for non-integer nu")
  ok = !(is.na(x) | (x==Inf) | (x==-Inf))
  z = ((x-xi)/omega)[ok]
  if(abs(alpha) == Inf) {
    z0 = replace(z, alpha*z < 0, 0)
    p = pf(z0^2, 1, nu)
    return(if(alpha>0) p else (1-p))
  }
  fp = function(v, alpha, nu, t.value)
    psn(sqrt(v) * t.value, 0, 1, alpha) * dchisq(v * nu, nu) * nu
  if(method == 4 || (method ==0  && int.nu &&
                     (nu < (8.2 + 3.55* log(log(length(z)+1))))))
    p = pst_int(z, 0, 1, alpha, nu)  # "method 4"
  else  {
    p = numeric(length(z))
    for (i in seq_len(length(z))) {
      if(abs(z[i]) == Inf)
        p[i] = (1+sign(z[i]))/2
      else {
        if(method==1 | method == 0)
          p[i] = pmst(z[i], 0, matrix(1,1,1), alpha, nu, ...) # method 1
        else {
          # upper = if(absalpha> 1) 5/absalpha + 25/(absalpha*nu) else 5+25/nu
          upper = 10 + 50/nu
          if(method==2 || (method==0 & (z[i] < upper) ))
            p[i] = integrate(dst, -Inf, z[i], dp=c(0,1,alpha, nu), ...)$value
          # method 2
          else
            p[i] = integrate(fp, 0, Inf, alpha, nu, z[i], ...)$value
          # method 3
        }}
    }}
  pr = rep(NA, length(x))
  pr[x==Inf] = 1
  pr[x==-Inf] = 0
  pr[ok] = p
  return(pmax(0,pmin(1,pr)))
}


pst_int = function (x, xi=0, omega=1, alpha=0, nu=Inf)
{# Jamalizadeh, A. and Khosravi, M. and Balakrishnan, N. (2009)
  if(nu != round(nu) | nu < 1) stop("nu not integer or not positive")
  z = (x-xi)/omega
  if(nu == 1)
    atan(z)/pi + acos(alpha/sqrt((1+alpha^2)*(1+z^2)))/pi
  else { if(nu==2)
    0.5 - atan(alpha)/pi + (0.5 + atan(z*alpha/sqrt(2+z^2))/pi)*z/sqrt(2+z^2)
    else
      (pst_int(sqrt((nu-2)/nu)*z, 0, 1, alpha, nu-2) +
         pst_int(sqrt(nu-1)*alpha*z/sqrt(nu+z^2), 0, 1, 0, nu-1) * z *
         exp(lgamma((nu-1)/2) +(nu/2-1)*log(nu)-0.5*log(pi)-lgamma(nu/2)
             -0.5*(nu-1)*log(nu+z^2)))
  }
}

#' @rdname skewt
#'
#' @export
#'
qst = function (p, xi = 0, omega = 1, alpha = 0, nu=Inf, tol = 1e-8,
                 dp = NULL, method=0, ...)
{
  if(!is.null(dp)) {
    if(!missing(alpha))
      stop("You cannot set both component parameters and dp")
    xi = dp[1]
    omega = dp[2]
    alpha = dp[3]
    nu = dp[4]
  }
  if(length(alpha) > 1) stop("'alpha' must be a single value")
  if(length(nu) > 1) stop("'nu' must be a single value")
  if (nu <= 0) stop("nu must be non-negative")
  if (nu > 1e4) return(qsn(p, xi, omega, alpha))
  if (nu == 1) return(qsc(p, xi, omega, alpha))
  if (alpha == Inf)
    return(xi + omega * sqrt(qf(p, 1, nu)))
  if (alpha == -Inf)
    return(xi - omega * sqrt(qf(1 - p, 1, nu)))
  na = is.na(p) | (p < 0) | (p > 1)
  abs.alpha = abs(alpha)
  if(alpha < 0) p = (1-p)
  zero = (p == 0)
  one = (p == 1)
  x = xa = xb = xc = fa = fb = fc = rep(NA, length(p))
  nc = rep(TRUE, length(p)) # not converged (yet)
  nc[(na| zero| one)] = FALSE
  fc[!nc] = 0
  xa[nc] = qt(p[nc], nu)
  xb[nc] = sqrt(qf(p[nc], 1, nu))
  fa[nc] = pst(xa[nc], 0, 1, abs.alpha, nu, method=method, ...) - p[nc]
  fb[nc] = pst(xb[nc], 0, 1, abs.alpha, nu, method=method, ...) - p[nc]
  regula.falsi = FALSE
  while (sum(nc) > 0) { # alternate regula falsi/bisection
    xc[nc] = if(regula.falsi)
      xb[nc] - fb[nc] * (xb[nc] - xa[nc])/(fb[nc] - fa[nc])    else
        (xb[nc] + xa[nc])/2
    fc[nc] = pst(xc[nc], 0, 1, abs.alpha, nu, method=method, ...) - p[nc]
    pos = (fc[nc] > 0)
    xa[nc][!pos] = xc[nc][!pos]
    fa[nc][!pos] = fc[nc][!pos]
    xb[nc][pos] = xc[nc][pos]
    fb[nc][pos] = fc[nc][pos]
    x[nc] = xc[nc]
    nc[(abs(fc) < tol)] = FALSE
    regula.falsi = !regula.falsi
  }
  # x = replace(x, na, NA)
  x = replace(x, zero, -Inf)
  x = replace(x, one, Inf)
  Sign = function(x) sign(x) + as.numeric(x==0)
  q = as.numeric(xi + omega * Sign(alpha)* x)
  names(q) = names(p)
  return(q)
}



delta.etc = function(alpha, Omega=NULL)
{
  inf = which(abs(alpha) == Inf)
  if(is.null(Omega)){ # case d=1
    delta = alpha/sqrt(1+alpha^2)
    delta[inf] = sign(alpha[inf])
    return(delta)
  }
  else { # d>1
    if(any(dim(Omega) != rep(length(alpha),2))) stop("dimension mismatch")
    Ocor = cov2cor(Omega)
    if(length(inf) == 0) { # d>1, standard case
      Ocor.alpha = as.vector(Ocor %*% alpha)
      alpha.sq = sum(alpha * Ocor.alpha)
      delta = Ocor.alpha/sqrt(1+alpha.sq)
      alpha. = sqrt(alpha.sq)
      delta. = sqrt(alpha.sq/(1+alpha.sq))
    }
    else { # d>1, case with some abs(alpha)=Inf
      if(length(inf) > 1)
        warning("Several abs(alpha)==Inf, I handle them as 'equal-rate Inf'")
      k = rep(0,length(alpha))
      k[inf] = sign(alpha[inf])
      Ocor.k = as.vector(Ocor %*% k)
      delta = Ocor.k/sqrt(sum(k * Ocor.k))
      delta. = 1
      alpha. = Inf
    }
    return(
      list(delta=delta, alpha.star=alpha., delta.star=delta., Omega.cor=Ocor))
  }
}


sn.cumulants <- function(xi = 0, omega = 1, alpha = 0, tau=0,
                         dp=NULL, n=4)
{
  cumulants.half.norm <- function(n=4){
    n <- max(n,2)
    n <- as.integer(2*ceiling(n/2))
    half.n  <-  as.integer(n/2)
    m <- 0:(half.n-1)
    a <- sqrt(2/pi)/(gamma(m+1)*2^m*(2*m+1))
    signs <- rep(c(1, -1), half.n)[seq_len(half.n)]
    a <- as.vector(rbind(signs*a, rep(0,half.n)))
    coeff <- rep(a[1],n)
    for (k in 2:n) {
      ind <- seq_len(k-1)
      coeff[k] <- a[k] - sum(ind*coeff[ind]*a[rev(ind)]/k)
    }
    kappa <- coeff*gamma(seq_len(n)+1)
    kappa[2] <- 1 + kappa[2]
    return(kappa)
  }
  if(!is.null(dp)) {
    if(!missing(alpha))
      stop("You cannot set both component parameters and dp")
    dp <- c(dp,0)[1:4]
    dp <- matrix(dp, 1, ncol=length(dp))
  }
  else  dp <- cbind(xi,omega,alpha,tau)
  delta <- ifelse(abs(dp[,3])<Inf, dp[,3]/sqrt(1+dp[,3]^2), sign(dp[,3]))
  tau <- dp[,4]
  if(all(tau==0)) {
    kv <- cumulants.half.norm(n)
    if(length(kv)>n) kv <- kv[-(n+1)]
    kv[2] <- kv[2] - 1
    kappa <- outer(delta,1:n,"^") * matrix(rep(kv,nrow(dp)),ncol=n,byrow=TRUE)
  }
  else{ # ESN
    if(n>4){
      warning("n>4 not allowed with ESN distribution")
      n <- min(n, 4)
    }
    kappa <- matrix(0, nrow=length(delta), ncol=0)
    for (k in 1:n) kappa <- cbind(kappa, zeta(k,tau)*delta^k)
  }
  kappa[,2] <- kappa[,2] + 1
  kappa <- kappa * outer(dp[,2],(1:n),"^")
  kappa[,1] <- kappa[,1] + dp[,1]
  kappa[,,drop=TRUE]
}


T.Owen <- function(h, a, jmax=50, cut.point=8)
{
  T.int <-function(h, a, jmax, cut.point)
  {
    fui <- function(h,i) (h^(2*i))/((2^i)*gamma(i+1))
    seriesL <- seriesH <- NULL
    i  <- 0:jmax
    low<- (h <= cut.point)
    hL <- h[low]
    hH <- h[!low]
    L  <- length(hL)
    if (L > 0) {
      b    <- outer(hL, i, fui)
      cumb <- apply(b, 1, cumsum)
      b1   <- exp(-0.5*hL^2) * t(cumb)
      matr <- matrix(1, jmax+1, L) - t(b1)
      jk   <- rep(c(1,-1), jmax)[1:(jmax+1)]/(2*i+1)
      matr <- t(matr*jk) %*%  a^(2*i+1)
      seriesL  <- (atan(a) - as.vector(matr))/(2*pi)
    }
    if (length(hH) > 0)  seriesH <-
      atan(a)*exp(-0.5*(hH^2)*a/atan(a)) * (1+0.00868*(hH*a)^4)/(2*pi)
    series <- c(seriesL, seriesH)
    id <- c((1:length(h))[low],(1:length(h))[!low])
    series[id] <- series  # re-sets in original order
    series
  }
  if(!is.vector(a) | length(a)>1) stop("'a' must be a vector of length 1")
  if(!is.vector(h)) stop("'h' must be a vector")
  aa <- abs(a)
  ah <- abs(h)
  if(is.na(aa)) stop("parameter 'a' is NA")
  if(aa==Inf) return(sign(a)*0.5*pnorm(-ah)) # sign(a): 16.07.2007
  if(aa==0)   return(rep(0,length(h)))
  na  <- is.na(h)
  inf <- (ah == Inf)
  ah  <- replace(ah,(na|inf),0)
  if(aa <= 1)
    owen <- T.int(ah,aa,jmax,cut.point)
  else
    owen<- (0.5*pnorm(ah) + pnorm(aa*ah)*(0.5-pnorm(ah))
            - T.int(aa*ah,(1/aa),jmax,cut.point))
  owen <- replace(owen,na,NA)
  owen <- replace(owen,inf,0)
  return(owen*sign(a))
}


pmst <- function(x, xi=rep(0,length(alpha)), Omega, alpha, nu=Inf, dp=NULL, ...)
{
  if(!(missing(alpha) & missing(Omega)) && !is.null(dp))
    stop("You cannot set both component parameters and dp")
  if(!is.null(dp)){
    if(!is.null(dp$xi)) xi <- dp$xi     else
      if(!is.null(dp$beta)) xi <- as.vector(dp$beta)
      Omega <- dp$Omega
      alpha <- dp$alpha
      nu <- dp$nu
  }
  if(!is.vector(x)) stop("x must be a vector")
  if(any(abs(alpha) == Inf)) stop("Inf's in alpha are not allowed")
  if(nu == Inf) return(pmsn(x, xi, Omega, alpha))
  d <- length(alpha)
  Omega<- matrix(Omega,d,d)
  omega<- sqrt(diag(Omega))
  Ocor <- cov2cor(Omega)
  O.alpha <- as.vector(Ocor %*% alpha)
  delta <- O.alpha/sqrt(1 + sum(alpha*O.alpha))
  Obig <- matrix(rbind(c(1, -delta), cbind(-delta, Ocor)), d+1, d+1)
  if(nu == as.integer(nu)) {
    z0 <- c(0,(x-xi)/omega)
    if(nu < .Machine$integer.max)
      p <- 2 * mnormt::pmt(z0, mean=rep(0,d+1), S=Obig, df=nu, ...)
    else
      p <- 2 * mnormt::pmnorm(z0, mean=rep(0,d+1), varcov=Obig, ...)
  }
  else {# for fractional nu, use formula in Azzalini & Capitanio (2003),
    # full-length paper, last paragraph of Section 4.2[Distr.function])
    z <- (x-xi)/omega
    fp <- function(v, Ocor, alpha, nu, t.value) {
      pv <-  numeric(length(v))
      for(k in seq_len(length(v))) pv[k] <- (dchisq(v[k] * nu, nu) * nu *
                                               pmsn(sqrt(v[k]) * t.value, rep(0,d), Ocor, alpha) )
      pv}
    p <- integrate(fp, 0, Inf, Ocor, alpha, nu, z, ...)$value
  }
  p
}


pmsn <- function(x, xi=rep(0,length(alpha)), Omega, alpha, tau=0,
                 dp=NULL, ...)
{
  if(!(missing(alpha) & missing(Omega)) && !is.null(dp))
    stop("You cannot set both component parameters and dp")
  if(!is.null(dp)){
    xi <- dp$xi
    Omega <- dp$Omega
    alpha <- dp$alpha
    tau <- if(is.null(dp$tau)) 0 else dp$tau
  }
  if(any(abs(alpha) == Inf)) stop("Inf's in alpha are not allowed")
  d <- length(alpha)
  Omega <- matrix(Omega, d, d)
  omega <- sqrt(diag(Omega))
  delta_etc <- delta.etc(alpha, Omega)
  delta <- delta_etc$delta
  Ocor <- delta_etc$Omega.cor
  Obig <- matrix(rbind(c(1,-delta), cbind(-delta,Ocor)), d+1, d+1)
  x <- if (is.vector(x)) matrix(x, 1, d) else data.matrix(x)
  if (is.vector(xi)) xi <- outer(rep(1, nrow(x)), as.vector(matrix(xi,1,d)))
  z0 <- cbind(tau, t(t(x - xi))/omega)
  mnormt::pmnorm(z0, mean=rep(0,d+1), varcov=Obig, ...)/pnorm(tau)
}
