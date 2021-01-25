

KR20 <- function(dat){


  numOfItem <- ncol(dat)
  numOfPerson <- nrow(dat)

  varSum <- var(rowSums(dat))
  pqSum <- sum(apply(dat, 2, function(x) sum(x)*(numOfPerson-sum(x))/numOfPerson^2))

  KR20 <- numOfItem/(numOfItem - 1) * (1 - pqSum/varSum)

  return(KR20)

}




KR21 <- function(dat){


  numOfItem <- ncol(dat)
  numOfPerson <- nrow(dat)

  varSum <- var(rowSums(dat))
  pqSum <- numOfItem * sum(dat)/(numOfItem*numOfPerson) * (1 - sum(dat)/(numOfItem*numOfPerson))

  KR21 <- numOfItem/(numOfItem - 1) * (1 - pqSum/varSum)

  return(KR21)

}


KR20(rawData_B)
KR21(rawData_B)
CronbachAlpha(rawData_B)


KR20(rawData_A)
KR21(rawData_A)
CronbachAlpha(rawData_A)
