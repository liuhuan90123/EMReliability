#' @title Cronbach's Alpha
#'
#' @description
#' A function to calculate Cronbach's Alpha reliability in CTT
#'
#' @param dat a data frame or matrix containing raw data
#' @return a reliability number
#'
#' @author {Huan Liu, University of Iowa, \email{huan-liu-1@@uiowa.edu}}
#'
#' @export


CronbachAlpha <- function(dat){

  # num of culomns
  numColumn <- ncol(dat)

  # sum of each variable variance
  varianceItemSum <- sum(diag(cov(dat)))

  # sum of total variance
  varianceTotal <- sum(cov(dat))

  cronbachAlpha <- numColumn / (numColumn - 1) * (varianceTotal - varianceItemSum)/varianceTotal

  return(cronbachAlpha)

}
