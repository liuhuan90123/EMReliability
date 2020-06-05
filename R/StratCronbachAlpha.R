#' @title Stratified Cronbach's Alpha
#'
#' @description
#' A function to calculate Stratified Cronbach's Alpha reliability in CTT
#'
#' @param dat a data frame or matrix containing raw data
#' @param strat a vector containing number of items for each strat
#' @return a reliability number
#'
#' @author {Huan Liu, University of Iowa, \email{huan-liu-1@@uiowa.edu}}
#'
#' @export


StratCronbachAlpha <- function(dat, strat){
  # num of strats
  numStrat <- length(strat)

  # vector of reliability and variance
  reli <- c()
  vari <- c()

  # Cronbach's alpha and variance for strat #1
  for (k in 1:1) {
    reli[k] <- CronbachAlpha(dat[,1:strat[k]])
    vari[k] <- sum(cov(dat[,1:strat[k]]))
  }

  # Cronbach's alpha and variance for strat #2 to #numStrat
  for (k in 2:numStrat) {
    reli[k] <- CronbachAlpha(dat[,(sum(strat[1:k-1])+1):sum(strat[1:k])])
    vari[k] <- sum(cov(dat[,(sum(strat[1:k-1])+1):sum(strat[1:k])]))
  }

  # sum of total variance: variance of composite score
  varianceTotal <- sum(cov(dat))

  # varianceCompo
  stratCronbachAlpha <- 1-sum(vari*(1-reli))/varianceTotal

  return(stratCronbachAlpha)
}


