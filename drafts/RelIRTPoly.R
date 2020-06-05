#' @title IRT Reliability Polynomial Method
#'
#' @description
#' A function to calculate reliability for Scale Scores using Polynomial Method based on IRT CSEM
#'
#' @param itemPara a text file with parameters of sequence b and a, a is on the 1.702 metric
#' @param convTable a data frame or matrix containing conversion table of raw score to scale score
#' @param K a numeric number indicating highest degree of polynomial regression
#' @param estType estimation method, MLE or EAP
#' @param rawData raw data matrix
#' @return a data frame containing CSSEM using Polynomial Method
#'
#' @author {Huan Liu, University of Iowa, \email{huan-liu-1@@uiowa.edu}}
#'
#' @export



RelIRTPoly <- function(itemPara, convTable, K, estType, rawData){

  # raw score frequence
  rawFreq <- as.data.frame(table(rowSums(rawData)))
  names(rawFreq) <- c("rawScore", "freq")


  if (estType == "MLE"){

    # cssem using Polynomial method
    cssemMLEPolyDat <- CSSEMIRTPoly(itemPara, convTable, K, "MLE")$CSSEMPolyMLE

    # raw score
    cssemMLEPolyDat$rawScore <- c(0:nrow(itemPara))

    # merge data
    cssemMLEPolyDat <- merge(cssemMLEPolyDat, rawFreq, by = "rawScore")

    # weight
    cssemMLEPolyDat$wt <- cssemMLEPolyDat$freq / sum(cssemMLEPolyDat$freq)

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

    # raw score
    cssemEAPPolyDat$rawScore <- c(0:nrow(itemPara))

    # merge data
    cssemEAPPolyDat <- merge(cssemEAPPolyDat, rawFreq, by = "rawScore")

    # weight
    cssemEAPPolyDat$wt <- cssemEAPPolyDat$freq / sum(cssemEAPPolyDat$freq)

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


