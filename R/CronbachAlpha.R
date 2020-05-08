#' @title Cronbach's Alpha and Stratified Cronbach's Alpha
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


CronbachAlpha <- function(dat, strat){

  if(missing(strat)){
    # num of culomns
    numColumn <- ncol(dat)

    # sum of each variable variance
    varianceItemSum <- sum(diag(cov(dat)))

    # sum of total variance
    varianceTotal <- sum(cov(dat))

    cronbachAlpha <- numColumn / (numColumn - 1) * (varianceTotal - varianceItemSum)/varianceTotal

    return(cronbachAlpha)

  }else{
    # num of strats
    numStrat <- length(strat)

    # Cronbach's alpha and variance for strat #1
    for (k in 1:1) {
      reli[k] <- CronbachAlpha(dat[,1:strat[k]])
      vari[k] <- sum(cov(dat[,1:strat[k]]))
    }

    # Cronbach's alpha and variance for strat #2 to #numStrat
    for (k in 2:ns) {
      reli[k] <- CronbachAlpha(dat[,(sum(strat[1:k-1])+1):sum(strat[1:k])])
      vari[k] <- sum(cov(dat[,(sum(strat[1:k-1])+1):sum(strat[1:k])]))
    }

    # sum of total variance: variance of composite score
    varianceTotal <- sum(cov(dat))

    # varianceCompo
    cronbachAlpha <- 1-sum(vari*(1-reli))/varianceTotal

    return(cronbachAlpha)

  }
}
