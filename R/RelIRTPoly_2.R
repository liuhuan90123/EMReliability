#' @title IRT Reliability Polynomial Method
#'
#' @description
#' A function to calculate reliability for Scale Scores using Polynomial Method based on IRT CSEM
#'
#' @param itemPara a text file with parameters of sequence b and a, a is on the 1.702 metric
#' @param convTable a data frame or matrix containing conversion table of raw score to scale score, only raw score and scale score
#' @param K a numeric number indicating highest degree of polynomial regression
#' @param estType estimation method, MLE or EAP
#' @param rawData raw data matrix
#' @return a data frame containing CSSEM using Polynomial Method
#'
#' @author {Huan Liu, University of Iowa, \email{huan-liu-1@@uiowa.edu}}
#'
#' @export


RelIRTPoly_2 <- function(itemPara, convTable, K, estType){

  if (estType == "MLE"){

    # cssem using Polynomial method
    cssemMLEPolyDat <- CSSEMIRTPoly(itemPara, convTable, K, "MLE")$CSSEMPolyMLE

    # weight
    cssemMLEPolyDat$wt <- sapply(cssemMLEPolyDat$theta, FUN = function(x) dnorm(x)) / sum(sapply(cssemMLEPolyDat$theta, FUN = function(x) dnorm(x)))

    # SS variance
    SSVar <- sum(cssemMLEPolyDat$wt * (cssemMLEPolyDat$roundedSS - weighted.mean(cssemMLEPolyDat$roundedSS, cssemMLEPolyDat$wt))^2)

    # for loop to calculate reliability for scale score
    RelMLEPolyDat <- as.data.frame(matrix(nrow = K, ncol = 1))

    for(i in 1:K){

      errorVar <- eval(parse(text=(paste("sum(cssemMLEPolyDat$cssemPolyk", i, "^2 * cssemMLEPolyDat$wt)", sep = ""))))

      RelMLEPolyDat[i,1] <- 1 - errorVar/SSVar

    }

    RelMLEPolyDat$kValue <- 1:K

    # select reliability values not equal to 1, and larger than 0
    RelMLEPolyDat <- RelMLEPolyDat[RelMLEPolyDat$V1 > 0 & RelMLEPolyDat$V1 != 1,]

    # return results
    return(list("kValue" = RelMLEPolyDat$kValue, "RelMLEPoly" = RelMLEPolyDat$V1))


  }else if (estType == "EAP"){

    # cssem using Polynomial method
    cssemEAPPolyDat <- CSSEMIRTPoly(itemPara, convTable, K, "EAP")$CSSEMPolyEAP

    # weight
    cssemEAPPolyDat$wt <- sapply(cssemEAPPolyDat$theta, FUN = function(x) dnorm(x)) / sum(sapply(cssemEAPPolyDat$theta, FUN = function(x) dnorm(x)))

    # SS variance
    SSVar <- sum(cssemEAPPolyDat$wt * (cssemEAPPolyDat$roundedSS - weighted.mean(cssemEAPPolyDat$roundedSS, cssemEAPPolyDat$wt))^2)

    # for loop to calculate reliability for scale score
    RelEAPPolyDat <- as.data.frame(matrix(nrow = K, ncol = 1))

    for(i in 1:K){

      errorVar <- eval(parse(text=(paste("sum(cssemEAPPolyDat$cssemPolyk", i, "^2 * cssemEAPPolyDat$wt)", sep = ""))))

      RelEAPPolyDat[i,1] <- 1 - errorVar/SSVar

    }

    RelEAPPolyDat$kValue <- 1:K

    # select reliability values not equal to 1, and larger than 0
    RelEAPPolyDat <- RelEAPPolyDat[RelEAPPolyDat$V1 > 0 & RelEAPPolyDat$V1 != 1,]

    # return results
    return(list("kValue" = RelEAPPolyDat$kValue, "RelEAPPoly" = RelEAPPolyDat$V1))

  }else{

    warning("RelIRTPolyfunction only supports MLE and EAP estimation method!")

  }

}



# item parameters
itemPara <- itemPara_A

# theta and weights
theta <- NormalQuadraPoints(41)$nodes
weights <- NormalQuadraPoints(41)$weights

# CSEM
itemParaCSEM <- as.data.frame(CSEMIRT(theta, itemPara, "EAP"))


# fit conversion table using polynomial method to get f(theta) and f'(theta)









### CSSEM IRT Polynomial MLE New -----------------------------------------------

### method two --------------------------------------------------------------------


# theta
theta <- convTable_B$theta

# CSEM MLE
itemParaCSEM <- as.data.frame(CSEMIRT(theta, itemPara_B, "EAP"))

# merge data
itemParaCSEM <- merge(itemParaCSEM, convTable_B_Poly, by = "theta")

# change name to fit Polynomial Method function
names(itemParaCSEM) <- c("rawScore", "csem", "roundedSS")

# call PM function
cssemPolyMLE <- PolynomialMethod(itemParaCSEM, 20)



cssemDat <- itemParaCSEM
K <- 20

# create data frame to store r square
rSquaredDat <- as.data.frame(matrix(nrow = K, ncol = 1))
regCoefDat <- list()

# for loop to iterate different k
for (k in 1:K){

  # fit model with k
  modelK <- lm(roundedSS ~ poly(rawScore, k, raw=TRUE), cssemDat)

  # extract regression coefficients
  regCoef <- summary(modelK)$coefficients[, 1]
  regCoefDat[[k]] <- summary(modelK)$coefficients

  # extract r square coefficient
  rSquaredDat[k, 1]<- summary(modelK)$r.squared

  # check whether regression coefficient of highest order is missing
  if(is.na(regCoef[k+1])){

    message(paste("The maximum k accepted is", k-1, sep = " "))

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
  cssemDat$cssemPoly <- cssemDat$fx * cssemDat$csem

  # rename variable with indicator k
  names(cssemDat)[names(cssemDat) == 'cssemPoly'] <- paste("cssemPolyk", k, sep = "")

}

# return(list("RSquared" = as.matrix(rSquaredDat[1:k,]), "CSSEMPoly" = cssemDat))

regCoefDat[[4]]  ## get coefficients for k = 4



# # SS formula
#
# SS = 118.263733790 + 3.909028972 * x -0.360085771 * x^2 -0.036053623 * x^3 + 0.009014156 * x^4
#
# f=expression(118.263733790 + 3.909028972 * x -0.360085771 * x^2 -0.036053623 * x^3 + 0.009014156 * x^4)
# D(f,'x')
#
#
# # CSSEM formula
#
# CSSEM = 3.909028972 - 0.360085771 * (2 * x) - 0.036053623 * (3 * x^2) + 0.009014156 * (4 * x^3)

# SS formula



# names(itemParaCSEM) <- c("theta", 'csemMLE',"scaleScore")


# theta
theta_New <- NormalQuadraPoints(41)$nodes
weights_New <- NormalQuadraPoints(41)$weights

# CSEM MLE
itemParaCSEM_New <- as.data.frame(CSEMIRT(theta_New, itemPara_B, "EAP"))

# itemParaCSEM_New <- within(itemParaCSEM_New,{
#   scaleScoreNew =  118.35478428 + 3.92088450 * theta -0.47021240 * theta^2 -0.03740998 * theta^3 + 0.01328610 * theta^4
#   cssemNew = 3.92088450 -0.47021240 * (2 * csemMLE) -0.03740998  * (3 * csemMLE^2) + 0.01328610 * (4 * csemMLE^3)
#   roundedSS = round(scaleScoreNew)
# })

# form B
itemParaCSEM_New <- within(itemParaCSEM_New,{
  scaleScoreNew =  118.35478428 + 3.92088450 * theta -0.47021240 * theta^2 -0.03740998 * theta^3 + 0.01328610 * theta^4
  cssemNew = 3.92088450 -0.47021240 * (2 * csemEAP) -0.03740998  * (3 * csemEAP^2) + 0.01328610 * (4 * csemEAP^3)
  roundedSS = scaleScoreNew
})

# form A
# itemParaCSEM_New <- within(itemParaCSEM_New,{
#   scaleScoreNew =  118.263733790 + 3.909028972 * theta -0.360085771 * theta^2 -0.036053623 * theta^3 + 0.009014156 * theta^4
#   cssemNew = 3.909028972 -0.360085771 * (2 * csemMLE) -0.036053623  * (3 * csemMLE^2) + 0.009014156 * (4 * csemMLE^3)
#   roundedSS = round(scaleScoreNew)
# })

### read new posterior distribution rates

postDist <- read.table("TestData/PosteriorDistribution.txt", sep = " ")

plot(postDist$V1, postDist$V2)


itemParaCSEM_New$weights_new <- postDist$V2


# SS variance
SSVar <- sum(itemParaCSEM_New$weights_new * (itemParaCSEM_New$roundedSS - weighted.mean(itemParaCSEM_New$roundedSS, itemParaCSEM_New$weights_new))^2)

# error variance
errorVar <- sum(itemParaCSEM_New$cssemNew * itemParaCSEM_New$weights_new)

# reliability
RelMLEPolyNew <- 1 - errorVar/SSVar
RelMLEPolyNew





sum(itemParaCSEM_New$weights_new)




# change variable name
names(cssemPolyMLE$CSSEMPoly)[names(cssemPolyMLE$CSSEMPoly) == 'rawScore'] <- 'theta'
names(cssemPolyMLE$CSSEMPoly)[names(cssemPolyMLE$CSSEMPoly) == 'csem'] <- 'csemMLE'

# return results
return(list("RSquared" = cssemPolyMLE$RSquared, "CSSEMPolyMLE" = cssemPolyMLE$CSSEMPoly))





