#' @title TestRelIRT
#'
#' @description
#' A function to calculate test reliability of IRT
#'
#' @param itemPara a data frame or matrix with parameters of sequence b, a and c on the 1.702 metric
#' @return test reliability of IRT
#'
#' @author {Huan Liu, University of Iowa, \email{huan-liu-1@@uiowa.edu}}
#' @export


TestRelIRT <- function(itemPara, convTable){



  itemPara <- itemPara_A
  convTable <- convTable_A



  if (ncol(itemPara) == 3){
    #  item parameters should be on the 1.702 metric
    names(itemPara) <- c("b", "a", "c")
  }

  if (ncol(itemPara) == 2){
    #  item parameters should be on the 1.702 metric
    names(itemPara) <- c("b", "a")
    itemPara$c <- 0
  }

  if (ncol(itemPara) == 1){
    #  item parameters should be on the 1.702 metric
    names(itemPara) <- c("b")
    itemPara$a <- 1
    itemPara$c <- 0
  }

  # num of items
  numOfItem <- nrow(itemPara)

  # num of quadratures
  numOfQuad <- 11

  # weights and nodes
  quadPoints <- NormalQuadraPoints(numOfQuad)

  # replicate item parameter and theta
  itemParaRep <-itemPara[rep(seq_len(numOfItem), each = numOfQuad),]
  itemParaRep$theta <- rep(quadPoints$nodes, each = 1, length.out = numOfQuad*numOfItem)

  # calculate information by theta
  itemParaRep <- within(itemParaRep, {
    P = c + (1 - c) / (1 + exp(-1.702 * a * (theta - b)))
    Q = 1 - P
    PQ = P * Q
    info = 1.702**2 * a**2 * P * Q
  })

  ## true score variance approach -------------------------------------------------------

  # true score variance

  # sum probability by theta
  itemParaAggr <- aggregate(itemParaRep, by=list(Category=itemParaRep$theta), FUN=sum)

  # add weights for each theta
  itemParaAggr$weights <- quadPoints$weights

  # calculate true score variance
  varianceTrue <- sum((itemParaAggr$P)^2 * itemParaAggr$weights) - (sum(itemParaAggr$P * itemParaAggr$weights))^2

  # error variance

  # claculate error variance
   <- sum(itemParaAggr$PQ * itemParaAggr$weights)


  varianceTrue/ (varianceTrue + varianceError2)



  ## observed score variance approach ------------------------------------------------------------

  # order by theta
  itemParaRep <- itemParaRep[order(itemParaRep$theta),]

  # define matrix of marginal distribution of theta
  fxTheta <- matrix(NA, nrow = numOfQuad, ncol = numOfItem + 1) # 41 num of quadratures, 41 num of raw sxores

  # for loop to calculate fxTheta
  for (i in 1:numOfQuad){

    probs <- matrix(c(itemParaRep[(1 + numOfItem * (i - 1)):(numOfItem * i),]$P),
                    nrow = numOfItem, ncol = 1, byrow = FALSE)

    fxTheta[i, ] <- LordWingersky(probs)$probability

  }

  # transform to data frame
  fxTheta <- as.data.frame(fxTheta)

  # add quadrature weights
  fxTheta$weights <- quadPoints$weights
  fxTheta$theta <- quadPoints$nodes

  # weighted mean
  fxTheta$weightedMean <- apply(fxTheta[1:(1+numOfItem)], 1, function(x)  weighted.mean(c(0:numOfItem), x))
  fxTheta$weightedMeanSS <- apply(fxTheta[1:(1+numOfItem)], 1, function(x)  weighted.mean(convTable$roundedSS, x))


  # weighted mean for raw score and scale score
  weightedMean <- sum(fxTheta$weightedMean  * fxTheta$weights)
  weightedMeanSS <- sum(fxTheta$weightedMeanSS  * fxTheta$weights)

  # variance of True X
  varianceTrueX <- sum(fxTheta$weights * (fxTheta$weightedMean - weightedMean)^2)
  varianceTrueXSS <- sum(fxTheta$weights * (fxTheta$weightedMeanSS - weightedMeanSS)^2)

  # variance of error
  for(i in 1:numOfQuad){
    fxTheta[i,"varRaw"] <- sum( fxTheta[i,1:(1+numOfItem)] * (c(0:numOfItem) - fxTheta[i,"weightedMean"])^2)
    fxTheta[i,"varSS"] <- sum( fxTheta[i,1:(1+numOfItem)] * (convTable$roundedSS - fxTheta[i,"weightedMeanSS"])^2)
  }


  # error variance
  varianceError <- sum(fxTheta$varRaw * fxTheta$weights)
  varianceErrorSS <- sum(fxTheta$varSS * fxTheta$weights)

  # variance of Obs X
  varianceObsX <- varianceError + varianceTrueX
  varianceObsXSS <- varianceErrorSS + varianceTrueXSS


  # variance of Obs X Approach 2

  # calculate weighted distribution
  fxThetaWeighted <- apply(fxTheta[,1:(1 + numOfItem)], 2, function(x) x * fxTheta[,"weights"])

  # sum weighted distribution
  fxDist <- as.data.frame(matrix(colSums(fxThetaWeighted[,1:(1 + numOfItem)]), nrow = (1 + numOfItem), ncol = 1))
  fxDist$X <- c(0:numOfItem)
  fxDist$roundedSS <- convTable$roundedSS
  names(fxDist) <- c("wts", "X","roundedSS")


  # weighted mean of Obs X and Y
  weightedMean <- sum(fxDist$X * fxDist$wts)/sum(fxDist$wts)
  weightedMeanSS <- sum(fxDist$roundedSS * fxDist$wts)/sum(fxDist$wts)

  # variance of Obs X
  varianceObsX <- sum(fxDist$wts * (fxDist$X - weightedMean)^2)
  varianceObsXSS <- sum(fxDist$wts * (fxDist$roundedSS - weightedMeanSS)^2)


  # IRT test reliability
  TestRelIRT <- 1 - varianceError/(varianceError + varianceTrueX)
  TestRelIRTSS <- 1 - varianceErrorSS/(varianceErrorSS + varianceTrueXSS)















  VarandRel <- as.data.frame(matrix(c(varianceError, varianceTrueX, varianceObsX, TestRelIRT,
                                      varianceErrorSS, varianceTrueXSS, varianceObsXSS, TestRelIRTSS
  )))
  rownames(VarandRel) <- c("Overall error variance for raw scores",
                           "True score variance for raw scores",
                           "Observed score variance for raw scores",
                           "Reliability for raw scores",
                           "Overall error variance for scale scores",
                           "True score variance for scale scores",
                           "Observed score variance for scale scores",
                           "Reliability for scale scores")
  colnames(VarandRel) <- "coefficient"
  VarandRel


  conditionalSEMs <- fxTheta[,c("theta", "weights","weightedMean", "varRaw", "weightedMeanSS", "varSS")]
  names(conditionalSEMs) <- c("Theta", "weights","Ex_Raw", "Raw_Variance", "Ex_Scale", "Scale_Variance")
  conditionalSEMs$Raw_CSEM <- sqrt(conditionalSEMs$Raw_Variance)
  conditionalSEMs$Scale_CSEM <- sqrt(conditionalSEMs$Scale_Variance)
  rownames(conditionalSEMs) <- 1:nrow(conditionalSEMs)



  return(list("Fitted Frequencies" = fxDist,
              "Variance and Reliability" = VarandRel,
              "Conditional SEMs" = conditionalSEMs[,c("Theta", "weights",
                                                      "Ex_Raw", "Ex_Scale", "Raw_CSEM", "Scale_CSEM")]))


}



