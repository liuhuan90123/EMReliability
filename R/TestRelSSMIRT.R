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

TestRelSSMIRT <- function(itemPara, strat, cormat){

  # num of items
  numOfItem <- nrow(itemPara)

  # num of quadratures
  numOfQuad <- 15

  # number of factors
  numOfFactors <- length(strat)

  # set nodes and weights
  nodes <- seq(-5, 5, length.out = numOfQuad)
  nodesM <- as.matrix(expand.grid(nodes,nodes,nodes))
  weightsUnwtd <- dmvn(nodesM, c(0,0,0), cormat, log=FALSE)
  nodesM <- as.data.frame(nodesM)
  nodesM$weightsWtd <- weightsUnwtd / sum(weightsUnwtd)

  # fxtheta distribution function
  FxTheta <- function(itemPara){

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

    # reverse column sequence
    fxTheta <- fxTheta[, c(ncol(fxTheta):1)]

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

  # for loop
  tau <- c()
  errvar <- c()
  fyDistMat <- matrix(NA,numOfQuad^3,32)

  n <- 0
  for (k in 1:numOfQuad){
    for (j in 1:numOfQuad){
      for (i in 1:numOfQuad){
        # index
        n <- n+1

        # fx distribution
        fx1 <- t(fxTheta1[i,])
        fx2 <- t(fxTheta2[j,])
        fx3 <- t(fxTheta3[k,])

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

        # store results
        tau[n] <- weightedMean
        errvar[n] <- sum(fyDist$wts * (fyDist$y - weightedMean)^2)
        fyDistMat[n,] <- t(fyDist$wts)
      }
    }
  }

  nodesM$tau <- tau
  nodesM$errvar <- errvar
  nodesM[,7:38] <- fyDistMat

  # sum of error variance
  varianceError <- sum(nodesM$weightsWtd*nodesM$errvar)

  # sum of observed score variance
  fyThetaWeighted <- apply(nodesM[,7:(7 + numOfItem)], 2, function(x) x * nodesM[,"weightsWtd"])

  # sum weighted distribution
  fyObsDist <- as.data.frame(matrix(colSums(fyThetaWeighted[,1:(1 + numOfItem)]), nrow = (1 + numOfItem), ncol = 1))
  fyObsDist$y <- c(numOfItem:0)
  names(fyObsDist) <- c("wts", "y")

  # weighted mean of Obs Y
  weightedMean <- sum(fyObsDist$y * fyObsDist$wts)/sum(fyObsDist$wts)

  # variance of Obs Y
  varianceObsY <- sum(fyObsDist$wts * (fyObsDist$y - weightedMean)^2)

  # MIRT test reliability
  TestRelSSMIRT <- 1 - varianceError/varianceObsY

  return(TestRelSSMIRT)

}














