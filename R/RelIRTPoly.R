#' @title IRT Reliability for Scale Score using Polynomial Method
#'
#' @description
#' A function to calculate reliability for Scale Scores using Polynomial Method based on IRT CSEM
#'
#' @param itemPara a data frame or matrix with parameters of sequence b, a and c on the 1.702 metric
#' @param convTable a data frame or matrix containing conversion table of raw score,"rawScore", and scale score, "roundedSS"
#' @param K highest degree for cssem, with default value 10
#' @param estType estimation method, MLE or EAP
#' @return a list containing k value and corresponding reliability coefficient
#'
#' @author {Huan Liu, University of Iowa, \email{huan-liu-1@@uiowa.edu}}
#'
#' @export


RelIRTPoly <- function(itemPara, convTable, K = 10, estType){

  itemPara <- itemPara_A # test
  convTable <- convTable_A
  K <- 10
  # estType <- "EAP"


  itemPara <- itemPara_B # test
  convTable <- convTable_B
  K <- 10



  if (estType == "MLE"){
    # theta and weights
    theta <- NormalQuadraPoints(41)$nodes

    # CSEM
    itemParaCSEM <- as.data.frame(CSEMIRT(theta, itemPara, "MLE"))
    # weight
    itemParaCSEM$wt <- NormalQuadraPoints(41)$weights

    rSquaredDat <- as.data.frame(matrix(nrow = K, ncol = 1))
    relDat <- as.data.frame(matrix(nrow = K, ncol = 1))

    # for loop to iterate different k
    for (k in 1:K){

      k <- 4 # test

      k <- 7 # test

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
      # itemParaCSEM$SS <- round(itemParaCSEM$SS)

      # calculate transformation coefficients fx: from 1 to K
      itemParaCSEM$fx <- 0
      i <- k
      while(i > 1){
        itemParaCSEM$fx <- itemParaCSEM$fx +  regCoef[i+1] * (i * itemParaCSEM$theta^(i-1))
        i <- i-1
      }

      itemParaCSEM$fx <- itemParaCSEM$fx + regCoef[i+1]

      # calculate cssem using polynomial method
      itemParaCSEM$cssemPoly <- itemParaCSEM$fx * itemParaCSEM$csemMLE

      # check negative values
      negCount <- sum(itemParaCSEM$fx < 0)
      if(negCount > 0){
        message(paste("Negative transformation coefficient exists when k = ", k,
                      ". The corresponding CSSEM will be corrected to 0 to calculate reliability.", sep = " "))

        itemParaCSEM$cssemPoly[itemParaCSEM$fx<0] <- 0
      }

      SSVar <- sum(itemParaCSEM$wt * (itemParaCSEM$SS - weighted.mean(itemParaCSEM$SS, itemParaCSEM$wt))^2)
      errVar <- sum(itemParaCSEM$cssemPoly^2 * itemParaCSEM$wt)
      SSVar
      errVar

      # relDat[k,1] <- 1 - errVar/SSVar
      relDat[k,1] <- 1 - errVar/(SSVar+errVar)

    }

    relDat$kValue <- 1:K

    # select reliability values not equal to 1, and larger than 0
    relDat <- relDat[relDat$V1 > 0 & relDat$V1 != 1,]
    relDat

    # return results
    return(list("kValue" = relDat$kValue, "RelMLEPoly" = relDat$V1))


  }else if (estType == "EAP"){

    k <- 4 # test

    k <- 7 # test


    # theta and weights
    theta <- NormalQuadraPoints(41)$nodes

    # CSEM
    itemParaCSEM <- as.data.frame(CSEMIRT(theta, itemPara, "EAP"))
    # weight
    itemParaCSEM$wt <- NormalQuadraPoints(41)$weights

    rSquaredDat <- as.data.frame(matrix(nrow = K, ncol = 1))
    relDat <- as.data.frame(matrix(nrow = K, ncol = 1))

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
      # itemParaCSEM$SS <- round(itemParaCSEM$SS)

      # calculate transformation coefficients fx: from 1 to K
      itemParaCSEM$fx <- 0
      i <- k
      while(i > 1){
        itemParaCSEM$fx <- itemParaCSEM$fx +  regCoef[i+1] * (i * itemParaCSEM$theta^(i-1))
        i <- i-1
      }

      itemParaCSEM$fx <- itemParaCSEM$fx + regCoef[i+1]

      # calculate cssem using polynomial method
      itemParaCSEM$cssemPoly <- itemParaCSEM$fx * itemParaCSEM$csemEAP

      # check negative values
      negCount <- sum(itemParaCSEM$fx < 0)
      if(negCount > 0){
        message(paste("Negative transformation coefficient exists when k = ", k,
                      ". The corresponding CSSEM will be corrected to 0 to calculate reliability.", sep = " "))

        itemParaCSEM$cssemPoly[itemParaCSEM$fx<0] <- 0
      }

      SSVar <- sum(itemParaCSEM$wt * (itemParaCSEM$SS - weighted.mean(itemParaCSEM$SS, itemParaCSEM$wt))^2)
      errVar <- sum(itemParaCSEM$cssemPoly^2 * itemParaCSEM$wt)

      SSVar
      errVar


      # relDat[k,1] <- 1 - errVar/SSVar
      relDat[k,1] <- 1 - errVar/(SSVar+errVar)

    }

    relDat$kValue <- 1:K

    # select reliability values not equal to 1, and larger than 0
    relDat <- relDat[relDat$V1 > 0 & relDat$V1 != 1,]
    relDat

    # return results
    return(list("kValue" = relDat$kValue, "RelEAPPoly" = relDat$V1))



  }else{

    warning("RelIRTPolyfunction only supports MLE and EAP estimation method!")

  }
}









