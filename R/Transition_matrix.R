# class -----

#' @rdname TranMat
#'
#' @name Credit migration
#'
#' @title Functions to apply to the credit migration or transition matrix
#'
#' @usage
#'
#' mat = read_mat(file = "transition_matrix.xlsx")
#'
#' @param object TranMat object
#' @param horizon Number of year for the computation
#'
#' @details Credit migration, or a transition matrix, indicates changes in the quality of settled credit at a
#' particular company. Transition matrices are the main input in various applications of risk
#' management.
#'
#' Credit rating is a process where any credit rating observation can form one of several state
#' ratings. In this research, it is assumed that the credit rating process follows the Markov chain
#' process. This means that the probability placed on one state can only be determined by
#' knowing the state from its previous observation. The assumption of Markov chain in the
#' credit rating process implies that the credit transition is more time invariant or time
#' homogenous, where the transition probability remains the same towards time and constant
#' during the predetermined horizon.
#'
#' @return readmat reads the matrix from a excel file and puts it in the right format
#'
#' the function 'tran' compute the transtitio matrix to horizon years
#'
#' the function 'abs_pd' compute the absolute probability matrix of default to horizon years
#'
#' the function 'con_pd' compute the conditional probability matrix of default to horizon years
#'
#' the function 'cum_pd' compute the cumulative probability matrix of default to horizon years
#'
#' @references http://www.bis.org/ifc/publ/ifcb31u.pdf
#'
#' @examples
#' mat = read_mat(file = file)
#' tran(mat, horizon =  10)
#' abs_pd(mat, horizon =  10)
#' con_pd(mat, horizon =  10)
#' cum_pd(mat, horizon =  10)
#'
#' @export
#'
setClass(Class = "TranMat",
         slots = c(
           x = "matrix"),
         package = "ECTools")

# read_mat ----

#' @rdname TranMat
#'
#' @export
read_mat = function(file = file) {

  require(readxl)
  require(dplyr)
  require(tibble)

  TM = read_xlsx(file)

  rat = c("AAA", "AA", "A", "BBB", "BB", "B", "CCC", "Default")

  TM = TM %>%
    filter(X__1 %in% rat) %>%
    as.data.frame() %>%
    column_to_rownames(var = "X__1") %>%
    select_(.dots = rat) %>%
    as.matrix()

  mat = new(Class = "TranMat",
      x = TM)

  return(mat)
}

# method ----

tran.mat = function(object = object, horizon = horizon) {

  rat = c("AAA", "AA", "A", "BBB", "BB", "B", "CCC", "Default")
  if(!(identical(row.names(object@x), rat) & identical(colnames(object@x), rat))) {stop("matrix is not a transition matrix")}

  if(!is.matrix(object@x)) {stop("matrix must be a matrix")}
  if(!horizon > 1) {stop("horizon must be greater than 1")}

  tran_mat = object@x

  for (i in 2:horizon) {

    tran_mat  = tran_mat %*% object@x

  }
  return(tran_mat)

}

# cumulative pd
cum_pd.mat = function(object = object, horizon = horizon) {

  cum_pd_mat = as.matrix(object@x[, "Default"])

  for (i in 2:horizon) {

    tran_mat = tran(object = object, horizon = i)
    cum_pd_mat = cbind(cum_pd_mat, tran_mat[, "Default"])

  }
  return(cum_pd_mat)

}


abs_pd.mat = function(object = object, horizon = horizon) {

  cum_pd_mat = cum_pd(object = object, horizon = horizon)

  abs_pd_mat = cum_pd_mat

  for(i in 2:horizon) {

    abs_pd_mat[,i] = cum_pd_mat[,i] - cum_pd_mat[,i - 1]

  }

  return(abs_pd_mat)
}


con_pd.mat = function(object = object, horizon = horizon) {

  cum_pd_mat = cum_pd(object = object, horizon = horizon)

  con_pd_mat = cum_pd_mat

  for(i in 2:horizon) {

    con_pd_mat[,i] = (cum_pd_mat[,i] - cum_pd_mat[,i - 1]) / (1 - cum_pd_mat[,i - 1])

  }

  return(con_pd_mat)
}


#' @rdname TranMat
#'
#' @export
tran = function(object, horizon = horizon) UseMethod("tran", object)

#' @rdname TranMat
#'
#' @export
abs_pd = function(object, horizon = horizon) UseMethod("abs_pd", object)

#' @rdname TranMat
#'
#' @export
con_pd = function(object, horizon = horizon) UseMethod("con_pd", object)

#' @rdname TranMat
#'
#' @export
cum_pd = function(object, horizon = horizon) UseMethod("cum_pd", object)

setMethod("tran", "TranMat", tran.mat)
setMethod("abs_pd", "TranMat", abs_pd.mat)
setMethod("con_pd", "TranMat", con_pd.mat)
setMethod("cum_pd", "TranMat", cum_pd.mat)
