RATIOSM_OUT = function(A = A, B = B, C = C, D = D, PD = PD, LGD = LGD, weight = weight, EAD = EAD, CORR = CORR, UL = UL) {

  rc_out = RcOut(A = A, B = B, C = C, D = D, PD = PD, LGD = LGD, weight = weight, EAD = EAD, CORR = CORR, UL = UL)
  ratio = sum(EAD * rc_out) / sum(EAD)

  return(ratio)
}

fit = function(A = A, B = B, C = C, D = D, PD = PD, LGD = LGD, weight = weight, EAD = EAD, CORR = CORR, UL = UL, CONTRI = CONTRI) {

  r_out = RATIOSM_OUT(A = A, B = B, C = C, D = D, PD = PD, LGD = LGD, weight = weight, EAD = EAD, CORR = CORR, UL = UL)

  r_in = sum(CONTRI * weight) / sum(EAD)

  return(abs(r_in - r_out))
}

#' CeR optimization
#'
#' @param A A
#' @param B B
#' @param C C
#' @param D D
#' @param PD PD
#' @param LGD LGD
#' @param weight N_client
#' @param EAD EAD
#' @param CORR CORR
#' @param UL UL
#' @param CONTRI Conttribution
#' @param variation range of parameter variation
#' @param maxiter Numeric, the maximum number of iterations to run
#' @param parallel a logical argument specifying if parallel computing
#' should be used (TRUE) or not (FALSE, default) for evaluating the
#' fitness function. This argument could also be used to specify the
#' number of cores to employ; by default, this is taken from detectCores.
#' Finally, the functionality of parallelization depends on system OS: on
#' Windows only 'snow' type functionality is available, while on
#' Unix/Linux/Mac OSX both 'snow' and 'multicore' (default) functionalities
#' are available
#' @param tol the desired accuracy.
#' @param ... additional named or unnamed arguments to be passed to ga
#'
#' @export
#'
cer_optim = function(A = A,
                     B = B,
                     C = C,
                     D = D,
                     PD = PD,
                     LGD = LGD,
                     weight = weight,
                     EAD = EAD,
                     CORR = CORR,
                     UL = UL,
                     CONTRI = CONTRI,
                     variation = 0.05,
                     maxiter = 1e4,
                     parallel = F,
                     tol = 1e-9,
                     ...) {

  require(GA)

  r_in = sum(CONTRI * weight) / sum(EAD)
  a_int = c(A * (1 - variation) , A * (1 + variation))
  min = c(B, C, D) * (1 - variation)
  max = c(B, C, D) * (1 + variation)

  fitness = function(x) {

    B_o = x[1]
    C_o = x[2]
    D_o = x[3]

    A_opt = optimize(fit,
                     B = B,
                     C = C,
                     D = D,
                     PD = PD,
                     LGD = LGD,
                     weight = weight,
                     EAD = EAD,
                     CORR = CORR,
                     UL = UL,
                     CONTRI = CONTRI,
                     interval = a_int,
                     tol = tol)

    A_o = A_opt$minimum

    residual = resid(A = A_o,
                     B = B_o,
                     C = C_o,
                     D = D_o,
                     PD = PD,
                     LGD = LGD,
                     weight = weight,
                     EAD = EAD,
                     CORR = CORR,
                     UL = UL,
                     RcIn = (CONTRI / EAD))

    return(-residual)
  }

  GA = ga(type = "real-valued",
          fitness = fitness,
          min = min,
          max = max,
          maxiter = maxiter,
          optim = T,
          parallel = parallel,
          ...)

  b = GA@solution[1, 1]
  c = GA@solution[1, 2]
  d = GA@solution[1, 3]

  a = optimize(fit,
               B = GA@solution[1, 1],
               C = GA@solution[1, 2],
               D = GA@solution[1, 3],
               PD = PD,
               LGD = LGD,
               weight = weight,
               EAD = EAD,
               CORR = CORR,
               UL = UL,
               CONTRI = CONTRI,
               interval = a_int,
               tol = tol)

  param = NULL

  param$sol = c(a$minimum, b, c, d)
  names(param$sol) = LETTERS[1:4]
  param$GA = GA

  return(param)
}
