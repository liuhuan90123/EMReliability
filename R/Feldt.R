### Feldt


Feldt <- function(dat){
  # sum of total variance
  varianceTotal <- sum(cov(dat))
  # sum of each variable variance
  varianceItemSum <- sum(diag(cov(dat)))
  # sum of covariance
  covarianceSum <- sum(rowSums(cov(dat))^2)
  feldt <- varianceTotal^2 / (varianceTotal^2 - covarianceSum) * (varianceTotal - varianceItemSum) / varianceTotal
  feldt

}

