### CSSEM Polynomial Method

PolynomialMethod <- function(cssemDat, K){


  rSquaredDat <- as.data.frame(matrix(nrow = K, ncol = 1))

  # for loop to iterate different k
  for (k in 1:K){

    # k <- 9

    # fit model with k and get coefficients
    modelK <- lm(roundedSS ~ poly(rawScore, k, raw=TRUE), cssemDat)
    regCoef <- summary(modelK)$coefficients[, 1]
    regCoef

    rSquaredDat[k, 1]<- summary(modelK)$r.squared


    if(is.na(regCoef[k+1])){
      warning(paste("The maximum k accepted is", k-1, sep = " "))
      break
    }

    # calculate transformation coefficients fx: from 1 to K
    cssemDat$fx <- 0
    i <- k

    while(i > 1){

      cssemDat$fx <- cssemDat$fx +  regCoef[i+1] * (i * cssemDat$rawScore^(i-1))
      i <- i-1

    }

    cssemDat$fx <- cssemDat$fx + regCoef[i+1]

    # calculate cssem using polynomial method
    cssemDat$cssemPoly <- cssemDat$fx * cssemDat$csemLord

    names(cssemDat)[names(cssemDat) == 'cssemPoly'] <- paste("cssemPolyk", k, sep = "")

  }

  print("R Squared summary")
  print(as.data.frame(rSquaredDat[1:k,]))
  print("CSSEM Polynomial Method")
  print(cssemDat)


}
