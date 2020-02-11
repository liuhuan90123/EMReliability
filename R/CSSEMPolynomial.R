### CSSEM Polynomial Method



# read conversion table from cvs file
convTable <- read.csv("TestData/ConversionTableFormX.csv")
convTable$roundedSS <- round(convTable$unroundedSS)
convTableSub <- convTable[,c("rawScore", "roundedSS")]



CSSEMPolynomial <- function(numOfItem, convTable){


  # numOfItem <- 40

  csemLordDat <- CSEMLord(numOfItem)

  cssemDat <- merge(csemLordDat, convTableSub, by = "rawScore")

  ## fit model with k = 3 and get coefficients
  m1 <- lm(roundedSS ~ 1 + rawScore + I(rawScore^2) + I(rawScore^3), cssemDat)
  summary(m1)

  # first derivative
  f <- expression(8.573e-04 * x^3 -5.481e-02* x^2 + 1.578e+00 * x + 9.922e+01)
  D(f, "x")

  # apply formula
  cssemDat$fx <- 0.0008573 * (3 * cssemDat$rawScore^2) - 0.05481 * (2 * cssemDat$rawScore) + 1.578
  cssemDat$cssemPoly <- cssemDat$fx * cssemDat$csemLord

  # return polynomial cssem based on Lord csem
  return(cssemDat)

}

CSSEMPolynomial(40, convTableSub)




# read conversion table from cvs file
convTable <- read.csv("TestData/ConversionTableFormX.csv")
convTable$roundedSS <- round(convTable$unroundedSS)
convTableSub <- convTable[,c("rawScore", "roundedSS")]



numOfItem <- 40
csemLordDat <- CSEMLord(numOfItem)
cssemDat <- merge(csemLordDat, convTableSub, by = "rawScore")


# cssemPolyMethod <- matrix(nrow = numOfItem + 1, ncol = 2)

PolynomialMethod <- function(cssemDat, K){

  for (k in 1:10){

    # k <- 3 # test

    # fit model with k and get coefficients
    modelK <- lm(roundedSS ~ poly(rawScore, k, raw=TRUE), cssemDat)
    regCoef <- summary(modelK)$coefficients[, 1]
    regCoef

    # calculate transformation coefficients fx: from k to 1 k = 3
    # cssemDat$fx <- regCoef[k+1] * (k * cssemDat$rawScore^(k-1)) + regCoef[k] * ((k-1) * cssemDat$rawScore^(k-2)) + regCoef[(k-1)]

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


cssemDatLong <- reshape(cssemDat,
                        direction = "long",
                        varying = list(names(cssemDat)[5:14]),
                        v.names = "cssem",
                        idvar = c("rawScore", "roundedSS"),
                        timevar = "Kvalue",
                        times = cssemPolyk1:cssemPolyk10)





