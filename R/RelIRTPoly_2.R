#' @title IRT Reliability for Scale Score using Polynomial Method
#'
#' @description
#' A function to calculate reliability for Scale Scores using Polynomial Method based on IRT CSEM
#'
#' @param itemPara a text file with parameters of sequence b, a, and c on the 1.702 metric
#' @param convTable a data frame or matrix containing conversion table of raw score,rawScore", and scale score, "roundedSS"
#' @param K a numeric number indicating selected degree of polynomial regression for SS, default highest degree for cssem is 20
#' @param estType estimation method, MLE or EAP
#' @return a list containing k value and corresponding reliability coefficient
#'
#' @author {Huan Liu, University of Iowa, \email{huan-liu-1@@uiowa.edu}}
#'
#' @export


RelIRTPoly_2 <- function(itemPara, convTable, ks, estType){

  if (estType == "MLE"){

    # default K
    K <- 20

    # theta and weights
    theta <- NormalQuadraPoints(41)$nodes

    # CSEM
    itemParaCSEM <- as.data.frame(CSEMIRT(theta, itemPara, "MLE"))

    rSquaredDat <- as.data.frame(matrix(nrow = K, ncol = 1))

    # for loop to iterate different k
    for (k in 1:K){

      # fit model with k
      modelK <- lm(roundedSS ~ poly(theta, k, raw=TRUE), convTable)

      # extract regression coefficients
      regCoef <- summary(modelK)$coefficients[, 1]

      # extract r square coefficient
      rSquaredDat[k, 1]<- summary(modelK)$r.squared

      # check whether regression coefficient of highest order is missing
      if(is.na(regCoef[k+1])){

        message(paste("The maximum k accepted is", k-1, sep = " "))

        break

      }

      # calculate SS: from 1 to K
      itemParaCSEM$SS <- 0
      i <- k

      while(i > 0){

        itemParaCSEM$SS <- itemParaCSEM$SS +  regCoef[i+1] * itemParaCSEM$theta^(i)
        i <- i-1

      }

      itemParaCSEM$SS <- itemParaCSEM$SS + regCoef[i+1]
      itemParaCSEM$SS <- round(itemParaCSEM$SS)

      # rename variable with indicator k
      names(itemParaCSEM)[names(itemParaCSEM) == 'SS'] <- paste("SSk", k, sep = "")

    }

    # call PM method, select ks th SS
    cssemDat <- itemParaCSEM[,c("theta", "csemMLE",paste("SSk", ks, sep = ""))]
    names(cssemDat) <- c("rawScore", "csemMLE", "roundedSS")
    cssemMLEPolyDat <- PolynomialMethod(cssemDat, K)$CSSEMPoly

    # weight
    cssemMLEPolyDat$wt <- NormalQuadraPoints(41)$weights

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

    # default K
    K <- 20

    # theta and weights
    theta <- NormalQuadraPoints(41)$nodes

    # CSEM
    itemParaCSEM <- as.data.frame(CSEMIRT(theta, itemPara, "EAP"))

    rSquaredDat <- as.data.frame(matrix(nrow = K, ncol = 1))

    # for loop to iterate different k
    for (k in 1:K){

      # fit model with k
      modelK <- lm(roundedSS ~ poly(theta, k, raw=TRUE), convTable)

      # extract regression coefficients
      regCoef <- summary(modelK)$coefficients[, 1]

      # extract r square coefficient
      rSquaredDat[k, 1]<- summary(modelK)$r.squared

      # check whether regression coefficient of highest order is missing
      if(is.na(regCoef[k+1])){

        message(paste("The maximum k accepted is", k-1, sep = " "))

        break

      }

      # calculate SS: from 1 to K
      itemParaCSEM$SS <- 0
      i <- k

      while(i > 0){

        itemParaCSEM$SS <- itemParaCSEM$SS +  regCoef[i+1] * itemParaCSEM$theta^(i)
        i <- i-1

      }

      itemParaCSEM$SS <- itemParaCSEM$SS + regCoef[i+1]
      itemParaCSEM$SS <- round(itemParaCSEM$SS)

      # rename variable with indicator k
      names(itemParaCSEM)[names(itemParaCSEM) == 'SS'] <- paste("SSk", k, sep = "")

    }

    # call PM method, select ks th SS
    cssemDat <- itemParaCSEM[,c("theta", "csemEAP",paste("SSk", ks, sep = ""))]
    names(cssemDat) <- c("rawScore", "csemEAP", "roundedSS")
    cssemEAPPolyDat <- PolynomialMethod(cssemDat, K)$CSSEMPoly

    # weight
    cssemEAPPolyDat$wt <- NormalQuadraPoints(41)$weights

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


