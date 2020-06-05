#' @title CSSEM Polynomial Method
#'
#' @description
#' A function to calculate CSEM for Scale Scores using Polynomial Method based on Lord's CSEM
#'
#' @param numOfItem a numeric number indicating number of items
#' @param convTable a data frame or matrix containing conversion table of raw score to scale score
#' @param K a numeric number indicating highest degree of polynomial regression, 10 in default
#'
#' @return a data frame containing CSSEM using Polynomial Method
#'
#' @author {Huan Liu, University of Iowa, \email{huan-liu-1@@uiowa.edu}}
#'
#' @export

CSSEMPolynomial <- function(numOfItem, convTable, K = 10){

  # csem Lord
  csemLordDat <- CSEMLord(numOfItem)

  # merge with converstion table
  cssemDat <- merge(csemLordDat, convTable, by = "rawScore")

  # change variable name
  names(cssemDat)[names(cssemDat) == 'csemLord'] <- 'csem'

  # apply polynomial method
  cssemPolynomial <- PolynomialMethod(cssemDat, K)

  # return results
  return(list("RSquared" = cssemPolynomial$RSquared, "CSSEMPolynomial" = cssemPolynomial$CSSEMPoly))

}












