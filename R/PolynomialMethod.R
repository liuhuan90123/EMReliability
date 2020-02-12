#' @title Polynomial Method for CSSEM
#'
#' @description
#' A function to implement polynomial method in calculating CSSEM
#'
#' @param cssemDat a ata frame or matrix containing conversion table of raw score to scale score, and csem of raw score
#' @param K degree of polynomial regression
#'
#' @return a data frame containing CSSEM using Polynomial Method using different k values
#'
#' @author {Huan Liu, University of Iowa, \email{huan-liu-1@@uiowa.edu}}
#'
#' @export

PolynomialMethod <- function(cssemDat, K){

  # create data frame to store r square
  rSquaredDat <- as.data.frame(matrix(nrow = K, ncol = 1))

  # for loop to iterate different k
  for (k in 1:K){

    # k <- 9

    # fit model with k
    modelK <- lm(roundedSS ~ poly(rawScore, k, raw=TRUE), cssemDat)

    # extract regression coefficients
    regCoef <- summary(modelK)$coefficients[, 1]

    # extract r square coefficient
    rSquaredDat[k, 1]<- summary(modelK)$r.squared

    # check whether regression coefficient of highest order is missing
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

    # rename variable with indicator k
    names(cssemDat)[names(cssemDat) == 'cssemPoly'] <- paste("cssemPolyk", k, sep = "")

  }

  return(list("R Squared summary" = as.data.frame(rSquaredDat[1:k,]),
         "CSSEM Polynomial Method" = cssemDat))

}
