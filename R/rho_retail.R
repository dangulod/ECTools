#' the coefficient of correlation for retail exposures defined
#'
#' @param segment segment
#' @param PD Probability of default
#'
#' @usage rho_retail(segment, PD)
#'
#' @return the coefficient of correlation for retail exposures defined
#' @export
#'
#' @description
#'
#'
#' According to the Capital Requirements Reglament (575/2013), Article 154. The coefficient
#' of correlation for corporates, institutions and central governments and central banks shall
#' be calculated according to the following formula:
#'
#' R = 0.12 * ((1 - exp(-50 * PD)) / (1 - exp(-50))) + 0.24 * (1 - ((1 - exp(-50 * PD)) / (1 - exp(-50))))
#'
#' For all exposures to large financial sector entities, the co-efficient of correlation of
#' is multiplied by 1.25. For all exposures to unregulated financial entities, the coefficients of
#' correlation are multiplied by 1.25.
#'
#' According to the Capital Requirements Reglament (575/2013), Article 154. The coefficient
#' of correlation is defined as:
#'
#' R = 0.03 * ((1 - exp(-35 * PD)) / (1 - exp(-35))) + 0.16 * (1 - ((1 - exp(-35 * PD)) / (1 - exp(-35))))
#'
#' However,
#'
#' For retail exposures secured by immovable property collateral a coefficient of correlation R of
#' 0,15 shall replace the figure produced by the correlation formula above.
#'
#' For qualifying revolving retail exposures, a coefficient of
#' correlation R of 0.04 shall replace the figure produced by the correlation formula above.
#'
#' @references
#'
#' CRR 575/2013/EU, Article 153 "Risk weighted exposure amounts for exposures to corporates, institutions and central governments and central banks"
#' CRR 575/2013/EU, Article 154 "Risk weighted exposure amounts for retail exposures"
#'
#' An Explanatory Note on the Basel II IRB Risk Weight Functions (page 14)
#' http://www.bis.org/bcbs/irbriskweight.pdf
#'
#' @examples rho_retail("mortage", 0.14)
#' rho_retail("credit cards", 0.14)
#' rho_retail("loans", 0.14)
#'
rho_retail = function(segment = segment, PD = PD, multiplier = F)
{
  if (grepl(pattern = "hipoteca|mortage", x = tolower(segment)))
  {
    # Hipotecario Particulares
    return(0.15)
  } else if (grepl(pattern = "tarjeta|card|revolving", x = tolower(segment)))
  {
    # Tarjetas de Crédito
    return(0.04)
  } else if (grepl(pattern = "empresa|corpora|institution|government|bank|banco|soberano|público|public", x = tolower(segment)))
  {
    # Empresas
    # Corporativa GBM
    # Sector Público
    # Soberano ME
    # Bancos
    rho = 0.12 * ((1 - exp(-50 * PD)) / (1 - exp(-50))) + 0.24 * (1 - ((1 - exp(-50 * PD)) / (1 - exp(-50))))
    if(multiplier)
    {
      return(rho * 1.25)
    } else
    {
      return(rho)
    }
  } else
  {
    # Otros Créditos
    # Consumo Particulares
    # Pymes
    # Otros activos
    return(0.03 * ((1 - exp(-35 * PD)) / (1 - exp(-35))) + 0.16 * (1 - ((1 - exp(-35 * PD)) / (1 - exp(-35)))))
  }
}


#' Boostraping of the second central moment
#'
#' @param x a numeric vector
#' @param p fraction of numbers to select
#' @param n_iter number of iterations
#' @param q quantile to return from the distribution generated
#' @param rnd logical, if TRUE the number are randomly selected
#'
#' @return In each iteration a subset from the vector x is selected
#' for variance computation, in this way a distribution is obtained.
#' Then, the value returned will be the quantile selected from the variance
#' distribution
#'
#' @export
#'
k2_boost = function(x = x, p = p, n_iter = 100, q = 0.9995, rnd = T) {

  if (p <= 0 | p > 1) stop("p must b between the range (0,1]")

  x = na.omit(x)
  p = ceiling(p * length(x))
  k2_boost = rep(0, n_iter)

  if (rnd == T) {

    for (i in 1:n_iter) k2_boost[i] = k2(sample(x, p))
    # for (i in 1:n_iter) k2_boost[i] = k2(x[runif(p, 1, length(x))])

  } else {

    for (i in 1:n_iter) {

      start = floor(runif(1, 1, length(x) - p))
      end = start + p

      k2_boost[i] = k2(x[start:end])

    }
  }

  quant = as.numeric(quantile(k2_boost, q))

  return(quant)
}
