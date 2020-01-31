# Cronbach's Alpha

CronbachAlpha <- function(dat){

  numColumn <- ncol(dat)
  varianceItemSum <- sum(diag(cov(dat)))
  varianceTotal <- sum(cov(dat))
  cronbachAlpha <- numColumn / (numColumn - 1) * (varianceTotal - varianceItemSum)/varianceTotal
  cronbachAlpha

}
