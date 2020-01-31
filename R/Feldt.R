### Feldt


Feldt <- function(dat){

  varianceTotal <- sum(cov(dat))
  varianceItemSum <- sum(diag(cov(dat)))
  covarianceSum <- sum(rowSums(cov(dat))^2)
  feldt <- varianceTotal^2 / (varianceTotal^2 - covarianceSum) * (varianceTotal - varianceItemSum) / varianceTotal
  feldt

}

