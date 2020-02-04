### Cronbach's Alpha

CronbachAlpha <- function(dat){
  # num of culomns
  numColumn <- ncol(dat)
  # sum of each variable variance
  varianceItemSum <- sum(diag(cov(dat)))
  # sum of total variance
  varianceTotal <- sum(cov(dat))
  cronbachAlpha <- numColumn / (numColumn - 1) * (varianceTotal - varianceItemSum)/varianceTotal
  cronbachAlpha

}
