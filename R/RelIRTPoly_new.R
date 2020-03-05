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


RelIRTPoly_new <- function(itemPara, convTable, K, estType){

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


