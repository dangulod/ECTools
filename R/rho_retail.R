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
#' Accoring to the Capital Requirements Reglament (575/2013), Article 154. The coefficient
#' of correlation is defined as:
#'
#' (1 - exp(-35 * PD)) / (1 - exp(-35))) + 0.16 * (1 - ((1 - exp(-35 * PD)) / (1 - exp(-35))))
#'
#' However,
#'
#' For retail exposures secured by immovable property collateral a coefficient of correlation R of
#' 0,15 shall replace the figure produced by the correlation formula above.
#'
#' For qualifying revolving retail exposures, a coefficient of
#' correlation R of 0,04 shall replace the figure produced by the correlation formula above.
#'
#' @references
#'
#' CRR 575/2013/EU, Article 154 "Risk weighted exposure amounts for retail exposures"
#'
#' An Explanatory Note on the Basel II IRB Risk Weight Functions (page 14)
#' http://www.bis.org/bcbs/irbriskweight.pdf
#'
#' @examples rho_retail("mortage", 0.14)
#' rho_retail("credit cards", 0.14)
#' rho_retail("loans", 0.14)
#'
rho_retail = function(segment = segment, PD = PD)
{
  if (grepl(pattern = "hipoteca|mortage", x = tolower(segment)))
  {
    return(0.15)
  } else if (grepl(pattern = "tarjeta|card|revolving", x = tolower(segment)))
  {
    return(0.04)
  } else
  {
    return(0.03 * ((1 - exp(-35 * PD)) / (1 - exp(-35))) + 0.16 * (1 - ((1 - exp(-35 * PD)) / (1 - exp(-35)))))
  }
}

