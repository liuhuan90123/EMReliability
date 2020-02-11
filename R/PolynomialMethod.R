### CSSEM Polynomial Method

PolynomialMethod <- function(cssemDat, K){

  for (k in 1:K){

    # fit model with k and get coefficients
    modelK <- lm(roundedSS ~ poly(rawScore, k, raw=TRUE), cssemDat)
    regCoef <- summary(modelK)$coefficients[, 1]
    regCoef

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

  cssemDat

}
