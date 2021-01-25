itemPara_A <- read.table("TestData/ItemParaFormX.txt")
names(itemPara_A) <- c("b", "a")
itemPara_A[,"a"] <- itemPara_A[,"a"]/1.702






UIRTtoCTT <- function(itemPara){
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
  numOfQuad <- 41

  # weights and nodes
  quadPoints <- NormalQuadraPoints(numOfQuad)

  # replicate item parameter and theta
  itemParaRep <-itemPara[rep(seq_len(numOfItem), each = numOfQuad),]
  itemParaRep$theta <- rep(quadPoints$nodes, each = 1, length.out = numOfQuad*numOfItem)
  # itemParaRep$wt<- rep(quadPoints$weights, each = 1, length.out = numOfQuad*numOfItem)
  itemParaRep$wt <- quadPoints$weights


  # calculate information by theta
  itemParaRep <- within(itemParaRep, {
    P = c + (1 - c) / (1 + exp(-1.702 * a * (theta - b)))
    Q = 1 - P
    PQ = P * Q
    info = 1.702**2 * a**2 * P * Q
  })

  ### Test level ---------------------------------------------------------------------

  ## true score variance approach

  # true score variance

  # sum probability by theta: tautheta
  itemParaAggr <- aggregate(itemParaRep, by=list(Category=itemParaRep$theta), FUN=sum) # P as tautheta  : itemParaAggr$P

  # add weights for each theta
  itemParaAggr$weights <- quadPoints$weights

  # calculate true score variance
  varianceTrue <- sum((itemParaAggr$P)^2 * itemParaAggr$weights) - (sum(itemParaAggr$P * itemParaAggr$weights))^2 # muT

  # tautheta and mutau
  itemParaRep$tautheta <- rep(itemParaAggr$P, each = 1, length.out = numOfQuad*numOfItem)
  itemParaRep$mutau <- sum(itemParaAggr$P * itemParaAggr$weights)

  # error variance

  # claculate error variance
  varianceError <- sum(itemParaAggr$PQ * itemParaAggr$weights)


  ## item level ------------------

  # sum probability by theta

  itemParaRep$PQwtd <- itemParaRep$wt * itemParaRep$PQ
  sumItemVar <- sum(aggregate(itemParaRep$PQwtd, by=list(Category=itemParaRep$b), FUN=sum)$x)                # sum of item variance

  itemParaRep$Pwtd <- itemParaRep$wt * itemParaRep$P                                                         # covariance and lamda
  itemParaAggrIP <- aggregate(itemParaRep, by=list(Category=itemParaRep$b), FUN=sum)
  itemParaRep$mui <- rep(itemParaAggrIP$Pwtd, each = numOfQuad, length.out = numOfQuad*numOfItem)  # ui

  itemParaRep <- within(itemParaRep,{
    covi = ((P -mui)*(tautheta - mutau) + PQ) * wt

  })

  itemParaAggrCov <- aggregate(itemParaRep, by=list(Category=itemParaRep$b), FUN=sum)                       # covariance and lamda


  ## observed score variance approach ---------------------------------------------------------------------  # varianceObsX

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

  # error variance 2
  varianceError2 <- sum(apply(fxTheta, 1, function(x) sum(x * (c(numOfItem:0) - weighted.mean(c(numOfItem:0), x))^2))*quadPoints$weights)

  # observed score variance

  # add quadrature weights
  fxTheta$weights <- quadPoints$weights

  # calculate weighted distribution
  fxThetaWeighted <- apply(fxTheta[,1:(1 + numOfItem)], 2, function(x) x * fxTheta[,"weights"])

  # sum weighted distribution
  fxDist <- as.data.frame(matrix(colSums(fxThetaWeighted[,1:(1 + numOfItem)]), nrow = (1 + numOfItem), ncol = 1))
  fxDist$X <- c(numOfItem:0)
  names(fxDist) <- c("wts", "X")

  # weighted mean of Obs X
  weightedMean <- sum(fxDist$X * fxDist$wts)/sum(fxDist$wts)

  # variance of Obs X
  varianceObsX <- sum(fxDist$wts * (fxDist$X - weightedMean)^2)

  # test reliability
  #  = 1 - varianceError / (varianceError + varianceTrue)

  lamda <- itemParaAggrCov$covi/varianceObsX


  cronbachAlphaUIRT <- numOfItem / (numOfItem-1) * (varianceObsX - sumItemVar)/varianceObsX

  FeldtRajuUIRT <- 1/(1-sum(lamda^2)) * (1-sumItemVar/varianceObsX)

  TestRelIRT <- 1 - varianceError2 / varianceObsX


  # return coefficient
  return(list("cronbachAlphaUIRT" = cronbachAlphaUIRT, "FeldtRajuUIRT" = FeldtRajuUIRT, "TestRelIRT" = TestRelIRT,
              "varianceObsX" = varianceObsX, "sumItemVar" = sumItemVar))


}


UIRTtoCTT(itemPara_A)


StratUIRTtoCTT <- function(itemPara, strat){

  # variance of X
  varianceObsX <- UIRTtoCTT(itemPara)$varianceObsX

  # num of strats
  numStrat <- length(strat)

  # vector of reliability and variance
  cronbachAlphaG <- c()
  varG <- c()

  # Cronbach's alpha and variance for strat #1
  for (k in 1:1) {
    cronbachAlphaG[k] <- UIRTtoCTT(itemPara[1:strat[k],])$cronbachAlphaUIRT
    varG[k] <- UIRTtoCTT(itemPara[1:strat[k],])$varianceObsX
  }

  # Cronbach's alpha and variance for strat #2 to #numStrat
  for (k in 2:numStrat) {
    cronbachAlphaG[k] <- UIRTtoCTT(itemPara[(sum(strat[1:k-1])+1):sum(strat[1:k]),])$cronbachAlphaUIRT
    varG[k] <- UIRTtoCTT(itemPara[(sum(strat[1:k-1])+1):sum(strat[1:k]),])$varianceObsX
  }

  stratUIRTtoCTT <- 1-(sum(varG*(1-cronbachAlphaG)))/varianceObsX

  return(stratUIRTtoCTT)

}


stratCronbachAlphaUIRT <- StratUIRTtoCTT(itemPara_A, strat = c(10, 30))
stratCronbachAlphaUIRT


# cronbachAlphaUIRT
# 0.9012882
# FeldtRajuUIRT
# 0.9022233
# stratCronbachAlphaUIRT
# 0.9023846
# TestRelIRT
# 0.878756






