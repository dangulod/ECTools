#' Contribucion de riesgo del CeR
#'
#' @param param vector with parameters c(A,B,C,D)
#' @param x dataframe con los datos
#'
#'
RC_OUT = function(param = param, x = x) {

  # x dataframe con los datos
  # param vector with parameters c(A,B,C,D)

  A = param[1]
  B = param[2]
  C = param[3]
  D = param[4]

  AVPD = sum(x$PD * x$NCLIENT) / (sum(x$NCLIENT))
  SUMEAD=sum(x$EAD * x$NCLIENT)
  AVEAD= SUMEAD / sum(x$NCLIENT)
  AVCORR = (sum(x$CORR* x$NCLIENT)/sum(x$NCLIENT)) ^ 2

  RC_OUT = x$UL * (A + B * (x$EAD / SUMEAD) * log(x$EAD / AVEAD) + C * log(x$PD/AVPD) + D * x$CORR * x$CORR * log(x$CORR * x$CORR / AVCORR)) - x$PD * x$LGD

  # RC_OUT = sum(x$EAD * RC_OUT) / sum(x$EAD)

  return(RC_OUT)
}


#' Residuos
#'
#'
RES = function(param = param, x = x) {

  # x dataframe con los datos
  # param = parametros a calibrar

  RC_OUT_1 =  RC_OUT(param = param, x = x)

  RES = sqrt(median(x$NCLIENT * (RC_OUT_1 - x$RC_IN) ^ 2) / sum(x$NCLIENT))

  return(RES)
}

# SM dataframe
# A
# B
# C
# D
# var vector (with the upper and the lower boundaries)

#' CeROpt
#'
#' @return CeROpt
#'
CeROpt = function(SM = SM, A = A, B = B, C = C, D = D, var = c(0.95, 1.05)) {

  if(!is.data.frame(SM)) {
    stop("SM is not a dataframe")
  }

  NlcOptim = getFromNamespace("NlcOptim", "NlcOptim")

  param = c(A, B, C, D)
  min = as.matrix(param * var[1])
  max = as.matrix(param * var[2])

  fitness = function(param = param) {
    fitness = RES(param, SM)
    return(fitness)
  }

  RATIOSM_IN = sum(SM$CON_SM * SM$NCLIENT) / sum(SM$EAD)

  confun = function(param = param){
    f = NULL
    f = (sum(SM$EAD * RC_OUT(x = SM, param = param)) / sum(SM$EAD)) - RATIOSM_IN
    return(list(ceq=f,c=NULL))
  }

  Opt = NlcOptim(X = param,
                 objfun = fitness,
                 confun = confun,
                 maxIter = 100, lb = min, ub = max)
  CeROpt = as.numeric(Opt$p)
  names(CeROpt) = c("A", "B", "C", "D")
  return(CeROpt)
}
