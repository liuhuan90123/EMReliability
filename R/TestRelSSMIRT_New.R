#' @title TestRelSSMIRT
#'
#' @description
#' A function to calculate test reliability of SS MIRT with 2PL model
#'
#' @param itemPara a data frame or matrix with parameters of sequence b, a1, a2,...,ai on the 1.702 metric
#' @param strat a vector containing number of items for each strat
#' @param cormat a correlation matrix for factors
#' @return test reliability of SS MIRT with 2PL model
#'
#' @author {Huan Liu, University of Iowa, \email{huan-liu-1@@uiowa.edu}}
#' @export

# library(LaplacesDemon)
# library(profvis)
# source("R/NormalQuadraPoints.R")
source("R/LordWingersky.R")
library(mvtnorm)


# TestRelMIRT(itemPara_SS_A, strat, cormat_A, convTable_A, numOfQuad)

# Form A
itemPara_SS_A <- read.table("TestData/SpanishLit_prm_A_SS.txt")[,c(7:10)]
cormat_A <- matrix(c(1, 0.9067069, 0.6994119,
                     0.9067069, 1, 0.4891160,
                     0.6994119,0.4891160,1), nrow = 3)
strat <- c(13, 12, 6)
# item parameter transformation
names(itemPara_SS_A) <- c("b", "a1","a2","a3")
itemPara_SS_A$a <- c(itemPara_SS_A$a1[1:13], itemPara_SS_A$a2[14:25], itemPara_SS_A$a3[26:31])
itemPara_SS_A[,"b"] <- -itemPara_SS_A[,"b"]/itemPara_SS_A[,"a"]
itemPara_SS_A[,"a"] <- itemPara_SS_A[,"a"]/1.702
itemPara_SS_A$a1[1:13] <- itemPara_SS_A$a[1:13]
itemPara_SS_A$a2[14:25] <- itemPara_SS_A$a[14:25]
itemPara_SS_A$a3[26:31] <- itemPara_SS_A$a[26:31]


# form B
itemPara_SS_B <- read.table("TestData/SpanishLit_prm_B_SS.txt")[,c(7:10)]
cormat_B <- matrix(c(1, 0.9722234, 0.5602197,
                     0.9722234, 1, 0.4795721,
                     0.5602197,0.4795721,1), nrow = 3)
strat <- c(13, 12, 6)
# item parameter transformation
names(itemPara_SS_B) <- c("b", "a1","a2","a3")
itemPara_SS_B$a <- c(itemPara_SS_B$a1[1:13], itemPara_SS_B$a2[14:25], itemPara_SS_B$a3[26:31])
itemPara_SS_B[,"b"] <- -itemPara_SS_B[,"b"]/itemPara_SS_B[,"a"]
itemPara_SS_B[,"a"] <- itemPara_SS_B[,"a"]/1.702
itemPara_SS_B$a1[1:13] <- itemPara_SS_B$a[1:13]
itemPara_SS_B$a2[14:25] <- itemPara_SS_B$a[14:25]
itemPara_SS_B$a3[26:31] <- itemPara_SS_B$a[26:31]


# read conversion tables
convTable_A <- read.csv("TestData/conversion_table_Form A.csv")
convTable_A <- convTable_A[1:32, c("RawScore", "roundedSS")]
names(convTable_A) <- c("y","roundedSS")

convTable_B <- read.csv("TestData/conversion_table_Form B.csv")
convTable_B <- convTable_B[1:32, c("RawScore", "roundedSS")]
names(convTable_B) <- c("y","roundedSS")

# test
itemPara <- itemPara_SS_B
cormat <- cormat_B
convTable <- convTable_B






TestRelMIRT <- function(itemPara, strat, cormat,convTable, numOfQuad = 11){ #modelType,

  # num of quadratures
  numOfQuad <- 11

  # number of factors
  numOfFactors <- length(strat)

  # set nodes and weights
  nodes <- seq(-4, 4, length.out = numOfQuad)

  # if(modelType = "bf"){
  #   numOfFactors <- numOfFactors + 1
  # }

  nodesM <- as.matrix(switch(numOfFactors,
                             {expand.grid(nodes)},
                             {expand.grid(nodes, nodes)},
                             {expand.grid(nodes, nodes, nodes)},
                             {expand.grid(nodes, nodes, nodes, nodes)},
                             {expand.grid(nodes, nodes, nodes, nodes, nodes)}))


  weightsUnwtd <- dmvnorm(nodesM, c(rep(0,numOfFactors)), cormat, log=FALSE) # mvtnorm
  nodesM <- as.data.frame(nodesM)
  nodesM$weightsWtd <- weightsUnwtd / sum(weightsUnwtd)


  # item parameters and quadrature points to P/Q

  # num of items
  numOfItem <- nrow(itemPara)

  # fxtheta distribution function
  FxTheta <- function(itemPara){

    # itemPara <- itemPara[1:strat[1],c("b", "a")]#test
    NormalQuadraPoints <- function(n){

      # set nodes ranging from -4 to 4
      nodes <- seq(-4, 4, length.out = n)

      # unnormalized weights
      weightsUnwtd <- sapply(nodes, FUN = function(x) dnorm(x))

      # normalized weightes
      weightsWtd <- weightsUnwtd / sum(weightsUnwtd)

      # return nodes and normalized weights
      return(list("nodes" = nodes, "weights" = weightsWtd))

    }

    # names item parameter
    names(itemPara) <- c("b", "a")

    # num of items
    numOfItem <- nrow(itemPara)

    # weights and nodes
    quadPoints <- NormalQuadraPoints(numOfQuad)

    # replicate item parameter and theta
    itemParaRep <-itemPara[rep(seq_len(numOfItem), each = numOfQuad),]
    itemParaRep$theta <- rep(quadPoints$nodes, each = 1, length.out = numOfQuad*numOfItem)

    # calculate information by theta
    itemParaRep <- within(itemParaRep, {
      P = 0 + (1 - 0) / (1 + exp(-1.702 * a * (theta - b)))
      Q = 1 - P
      PQ = P * Q
      info = 1.702**2 * a**2 * P * Q
    })

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
    fxTheta
  }

  # fxtheta distribution
  fxTheta1 <- FxTheta(itemPara[1:strat[1],c("b", "a")])
  fxTheta2 <- FxTheta(itemPara[(strat[1]+1):(strat[1]+strat[2]),c("b", "a")])
  fxTheta3 <- FxTheta(itemPara[(strat[1]+strat[2]+1):(strat[1]+strat[2]+strat[3]),c("b", "a")])

  names(fxTheta1) <- c(0:strat[1])
  names(fxTheta2) <- c(0:strat[2])
  names(fxTheta3) <- c(0:strat[3])

  fxTheta1$Var1 <- seq(-4, 4, length.out = numOfQuad)
  nodesM <- merge(x = nodesM, y = fxTheta1, by = "Var1", all.x = TRUE)

  fxTheta2$Var2 <- seq(-4, 4, length.out = numOfQuad)
  nodesM <- merge(x = nodesM, y = fxTheta2, by = "Var2", all.x = TRUE)

  fxTheta3$Var3 <- seq(-4, 4, length.out = numOfQuad)
  nodesM <- merge(x = nodesM, y = fxTheta3, by = "Var3", all.x = TRUE)


  for(i in 1:nrow(nodesM)){

    # fx distribution
    fx1 <- t(nodesM[i, 5:(5+strat[1])])
    fx2 <- t(nodesM[i, (5+strat[1]+1):(5+strat[1]+strat[2]+1)])
    fx3 <- t(nodesM[i, (5+strat[1]+strat[2]+1+1):(5+strat[1]+strat[2]+1+strat[3]+1)])

    rownames(fx1) <- 0:strat[1]
    rownames(fx2) <- 0:strat[2]
    rownames(fx3) <- 0:strat[3]

    xSum <- expand.grid(rownames(fx1), rownames(fx2), rownames(fx3))
    names(xSum) <- c("x1", "x2", "x3")

    fxSum <- expand.grid(fx1, fx2, fx3)
    names(fxSum) <- c("fx1", "fx2", "fx3")

    fxThetaSum <- cbind(fxSum, xSum)

    fxThetaSum$x1 <- as.numeric(as.character(fxThetaSum$x1))
    fxThetaSum$x2 <- as.numeric(as.character(fxThetaSum$x2))
    fxThetaSum$x3 <- as.numeric(as.character(fxThetaSum$x3))

    # fy distribution
    fxThetaSum$y <- fxThetaSum$x1 + fxThetaSum$x2 + fxThetaSum$x3
    fxThetaSum$wty <- fxThetaSum$fx1 * fxThetaSum$fx2 * fxThetaSum$fx3
    fy <- fxThetaSum[,c("y", "wty")]
    fyDist <- aggregate(fy$wty, by=list(Category=fy$y), FUN=sum)
    names(fyDist) <- c("y", "wts")

    # weighted mean of Obs Y (true y) and variance of Obs Y
    weightedMean <- sum(fyDist$y * fyDist$wts)/sum(fyDist$wts)
    varianceY <- sum(fyDist$wts * (fyDist$y - weightedMean)^2)

    # save results
    nodesM[i,"weightedMean"] <- weightedMean
    nodesM[i,"varianceY"] <- varianceY

    # SS

    fySSDist <- merge(fyDist, convTable, by = "y")

    # weighted mean of Obs Y (true y) and variance of Obs Y
    weightedMeanSS <- sum(fySSDist$roundedSS * fySSDist$wts)/sum(fySSDist$wts)
    varianceYSS <- sum(fySSDist$wts * (fySSDist$roundedSS - weightedMeanSS)^2)

    # store results
    nodesM[i,"weightedMeanSS"] <- weightedMeanSS
    nodesM[i,"varianceYSS"] <- varianceYSS

    nodesM[i,c(43:74)] <- t(fySSDist$wts)

  }

  # variance of error
  varianceError <- sum(nodesM$weightsWtd * nodesM$varianceY)
  varianceErrorSS <- sum(nodesM$weightsWtd * nodesM$varianceYSS)


  # weighted mean of True Y
  weightedMean <- sum(nodesM$weightedMean * nodesM$weightsWtd)/sum(nodesM$weightsWtd)
  weightedMeanSS <- sum(nodesM$weightedMeanSS * nodesM$weightsWtd)/sum(nodesM$weightsWtd)


  # variance of True Y
  varianceTrueY <- sum(nodesM$weightsWtd * (nodesM$weightedMean - weightedMean)^2)
  varianceTrueYSS <- sum(nodesM$weightsWtd * (nodesM$weightedMeanSS - weightedMeanSS)^2)

  # variance of Obs Y
  varianceObsY <- varianceError + varianceTrueY
  varianceObsYSS <- varianceErrorSS + varianceTrueYSS


  # variance of Obs Y Approach 2

  # sum of observed score variance
  fyThetaWeighted <- apply(nodesM[,43:(43 + numOfItem)], 2, function(x) x * nodesM[,"weightsWtd"])

  # sum weighted distribution
  fyObsDist <- as.data.frame(matrix(colSums(fyThetaWeighted[,1:(1 + numOfItem)]), nrow = (1 + numOfItem), ncol = 1))
  fyObsDist$y <- c(0:numOfItem) # test
  fyObsDist$roundedSS <- convTable$roundedSS
  names(fyObsDist) <- c("wts", "y","roundedSS")

  # weighted mean of Obs Y
  weightedMeanY <- sum(fyObsDist$y * fyObsDist$wts)/sum(fyObsDist$wts)
  weightedMeanSSY <- sum(fyObsDist$roundedSS * fyObsDist$wts)/sum(fyObsDist$wts)

  # variance of Obs Y
  varianceObsY2 <- sum(fyObsDist$wts * (fyObsDist$y - weightedMeanY)^2)
  varianceObsYSS2 <- sum(fyObsDist$wts * (fyObsDist$roundedSS - weightedMeanSSY)^2)

  fyObsDist <- fyObsDist[,c(2,3,1)]
  names(fyObsDist) <- c("Raw", "Scale", "Fitted_Freq")


  # MIRT test reliability
  TestRelSSMIRT <- 1 - varianceError/(varianceError + varianceTrueY)
  TestRelSSMIRTSS <- 1 - varianceErrorSS/(varianceErrorSS + varianceTrueYSS)


  VarandRel <- as.data.frame(matrix(c(varianceError, varianceTrueY, varianceObsY, TestRelSSMIRT,
              varianceErrorSS, varianceTrueYSS, varianceObsYSS, TestRelSSMIRTSS
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


  conditionalSEMs <- nodesM[,c("Var1", "Var2", "Var3", "weightsWtd","weightedMean", "varianceY", "weightedMeanSS", "varianceYSS")]
  names(conditionalSEMs) <- c("Theta1", "Theta2", "Theta3", "weights","Ex_Raw", "Raw_Variance", "Ex_Scale", "Scale_Variance")
  conditionalSEMs$Raw_CSEM <- sqrt(conditionalSEMs$Raw_Variance)
  conditionalSEMs$Scale_CSEM <- sqrt(conditionalSEMs$Scale_Variance)
  conditionalSEMs <- conditionalSEMs[with(conditionalSEMs, order(Theta1, Theta2, Theta3)), ]
  rownames(conditionalSEMs) <- 1:nrow(conditionalSEMs)

  return(list("Fitted Frequencies" = fyObsDist,
              "Variance and Reliability" = VarandRel,
              "Conditional SEMs" = conditionalSEMs[,c("Theta1", "Theta2", "Theta3", "weights",
                                                      "Ex_Raw", "Ex_Scale", "Raw_CSEM", "Scale_CSEM")]))
}

TestRelMIRT(itemPara_SS_A, strat, cormat_A, convTable_A, numOfQuad)
TestRelMIRT(itemPara_SS_B, strat, cormat_B, convTable_B, numOfQuad)


# profvis({
#   TestRelMIRT(itemPara_SS_A, strat, cormat_A, convTable_A, numOfQuad)
# })
#
# profvis({
#   TestRelMIRT(itemPara_SS_A, strat, cormat_A, convTable_A, numOfQuad=41)
# })

library(xlsx)
write.xlsx(conditionalSEMs, "conditionalSEMs_SS_B.xlsx")
