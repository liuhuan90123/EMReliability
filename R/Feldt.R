#' @title Feldt
#'
#' @description
#' A function to calculate Feldt reliability in CTT
#'
#' @param dat a data frame or matrix containing raw data
#' @return a reliability number
#'
#' @author {Huan Liu, University of Iowa, \email{huan-liu-1@@uiowa.edu}}
#'
#' @export

Feldt <- function(dat){

  # sum of total variance
  varianceTotal <- sum(cov(dat))

  # sum of each variable variance
  varianceItemSum <- sum(diag(cov(dat)))

  # sum of covariance
  covarianceSum <- sum(rowSums(cov(dat))^2)

  feldt <- varianceTotal^2 / (varianceTotal^2 - covarianceSum) * (varianceTotal - varianceItemSum) / varianceTotal

  return(feldt)

}

