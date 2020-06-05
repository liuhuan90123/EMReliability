#' @title Polynomial Method for CSSEM
#'
#' @description
#' A function to implement polynomial method in calculating CSSEM
#'
#' @param cssemDat a data frame or matrix containing conversion table of raw score to scale score, and csem of raw score
#' @param K a numeric number indicating highest degree of polynomial regression, 10 in default
#'
#' @return a data frame containing CSSEM using Polynomial Method using different k values
#'
#' @author {Huan Liu, University of Iowa, \email{huan-liu-1@@uiowa.edu}}
#'
#' @export

PolynomialMethod <- function(cssemDat, K = 10, gra = TRUE){

  # create data frame to store r square
  rSquaredDat <- as.data.frame(matrix(nrow = K, ncol = 1))

  # for loop to iterate different k
  for (k in 1:K){

    # fit model with k
    modelK <- lm(roundedSS ~ poly(rawScore, k, raw=TRUE), cssemDat)

    # extract regression coefficients
    regCoef <- summary(modelK)$coefficients[, 1]

    # extract r square coefficient
    rSquaredDat[k, 1]<- summary(modelK)$r.squared

    # check whether regression coefficient of highest order is missing
    if(is.na(regCoef[k+1])){

      message(paste("The maximum k accepted is", k-1, ".", sep = " "))

      break

    }

    # graph
    if(gra == TRUE){

      prd <- data.frame(rawScore = seq(from = range(cssemDat$rawScore)[1], to = range(cssemDat$rawScore)[2], length.out = 100))
      prd$predictedSS <- predict(modelK, newdata = prd, se.modelK = TRUE)

      p <- ggplot(prd, aes(x = rawScore, y = predictedSS)) +
        theme_bw() +
        geom_line() +
        geom_point(data = cssemDat, aes(x = rawScore, y = roundedSS)) +
        scale_y_continuous(name = "Scale Score") +
        scale_x_continuous(name = "Raw Score") +
        ggtitle(paste("Polynomial Method Fitted Line with k = ", k, sep = "")) +
        theme(plot.title = element_text(hjust = 0.5))

      print(p)
    }

    # calculate transformation coefficients fx: from 1 to K
    cssemDat$fx <- 0
    i <- k

    while(i > 1){

      cssemDat$fx <- cssemDat$fx +  regCoef[i+1] * (i * cssemDat$rawScore^(i-1))
      i <- i-1

    }

    cssemDat$fx <- cssemDat$fx + regCoef[i+1]

    # calculate SS: from 1 to K
    cssemDat$SS <- 0
    i <- k

    while(i > 0){

      cssemDat$SS <- cssemDat$SS +  regCoef[i+1] * cssemDat$rawScore^(i)
      i <- i-1

    }

    cssemDat$SS <- cssemDat$SS + regCoef[i+1]
    cssemDat$SS <- round(cssemDat$SS)

    # calculate cssem using polynomial method
    cssemDat$cssemPoly <- cssemDat$fx * (cssemDat$csem)

    # check negative values
    negCount <- sum(cssemDat$fx < 0)

    if(negCount > 0){
      message(paste("Negative transformation coefficient exists when k = ", k,
                    ". The corresponding CSSEM will not be available.", sep = " "))

      cssemDat$cssemPoly[cssemDat$fx<0] <- NA
    }

    # rename variable with indicator k
    names(cssemDat)[names(cssemDat) == 'cssemPoly'] <- paste("cssemPolyk", k, sep = "")
    names(cssemDat)[names(cssemDat) == 'fx'] <- paste("fxk", k, sep = "")
    names(cssemDat)[names(cssemDat) == 'SS'] <- paste("SSk", k, sep = "")

  }

  return(list("RSquared" = as.matrix(rSquaredDat[1:k,]), "CSSEMPoly" = cssemDat))

}






