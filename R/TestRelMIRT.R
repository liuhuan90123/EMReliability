# Test Reliability

###### simple structure   D method --------------------------

FX <- function(itemPara, newWts){
  # transform item parameters to the logistic metric
  names(itemPara) <- c("b", "a")

  # num of items
  numOfItem <- nrow(itemPara)

  # num of quadratures
  numOfQuad <- 15

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

  ## true score variance

  # sum probability by theta
  itemParaAggr <- aggregate(itemParaRep, by=list(Category=itemParaRep$theta), FUN=sum)

  # add weights for each theta
  # itemParaAggr$weights <- quadPoints$weights
  itemParaAggr$weights <- newWts

  # calculate true score variance
  varianceTrue <- sum((itemParaAggr$P)^2 * itemParaAggr$weights) - (sum(itemParaAggr$P * itemParaAggr$weights))^2


  ## error variance

  # claculate error variance
  varianceError <- sum(itemParaAggr$PQ * itemParaAggr$weights)


  ## observed score variance

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
  # fxTheta <- fxTheta[, c(ncol(fxTheta):1)]

  # transform to data frame
  fxTheta <- as.data.frame(fxTheta)

  fxTheta

  # # add quadrature weights
  # fxTheta$weights <- quadPoints$weights
  #
  # # calculate weighted distribution
  # # fxThetaWeighted <- apply(fxTheta[,1:41], 2, function(x) x * fxTheta[,"weights"])
  #
  # fxThetaWeighted <- apply(fxTheta[,1:(1 + numOfItem)], 2, function(x) x * fxTheta[,"weights"])
  #
  # # sum weighted distribution
  # #fxDist <- as.data.frame(matrix(colSums(fxThetaWeighted[,1:41]), nrow = 41, ncol = 1))
  # fxDist <- as.data.frame(matrix(colSums(fxThetaWeighted[,1:(1 + numOfItem)]), nrow = (1 + numOfItem), ncol = 1))
  # fxDist$X <- c(numOfItem:0)
  # names(fxDist) <- c("wts", "X")
  # fxDist
}




#####simple structure MIRT ##########

# item parameters

# read item parameters from txt file
itemPara_A_SS <- read.table("TestData/SpanishLit_prm_A_SS.txt")[,c(7:10)]

# Form A
# read correlations
# cor <- c(0.9067069, # 1&2
#          0.6994119, # 1&3
#          0.4891160) # 2&3

cor <- matrix(c(1, 0.9067069, 0.6994119,
                0.9067069, 1, 0.4891160,
                0.6994119,0.4891160,1), nrow = 3)

# form B
# itemPara_A_SS <- read.table("TestData/SpanishLit_prm_B_SS.txt")[,c(7:10)]
#
# cor <- matrix(c(1, 0.97, 0.56,
#                 0.97, 1, 0.48,
#                 0.56,0.48,1), nrow = 3)

names(itemPara_A_SS) <- c("b", "a1","a2","a3")
itemPara_A_SS$a <- c(itemPara_A_SS$a1[1:13], itemPara_A_SS$a2[14:25], itemPara_A_SS$a3[26:31])
itemPara_A_SS[,"b"] <- -itemPara_A_SS[,"b"]/itemPara_A_SS[,"a"]
itemPara_A_SS[,"a"] <- itemPara_A_SS[,"a"]/1.702

itemPara_A_SS$a1[1:13] <- itemPara_A_SS$a[1:13]
itemPara_A_SS$a2[14:25] <- itemPara_A_SS$a[14:25]
itemPara_A_SS$a3[26:31] <- itemPara_A_SS$a[26:31]

itemPara <- itemPara_A_SS


# num of items
numOfItem <- nrow(itemPara)

# num of quadratures
numOfQuad <- 15

# weights and nodes
# quadPoints <- NormalQuadraPoints(numOfQuad)

# install.packages("mvtnorm")
# install.packages("pbivnorm")
# install.packages("LaplacesDemon")
library(LaplacesDemon)
library(mvtnorm)
library(pbivnorm)
# set nodes ranging from -5 to 5
n <- 15
nodes <- seq(-5, 5, length.out = n)


# Three dimensions with correlation matrix--------------------------------------------------------------

nodes3 <- as.matrix(expand.grid(nodes,nodes,nodes))

weightsUnwtd <- dmvn(nodes3, c(0,0,0), cor, log=FALSE) # 41^3
nodes3 <- as.data.frame(nodes3)
nodes3$weightsWtd <- weightsUnwtd / sum(weightsUnwtd)


### weights ---
newWts1 <- as.numeric(unlist(aggregate(nodes3$weightsWtd, by=list(Category=nodes3$Var1), FUN=sum)["x"]))
newWts2 <- as.numeric(unlist(aggregate(nodes3$weightsWtd, by=list(Category=nodes3$Var2), FUN=sum)["x"]))
newWts3 <- as.numeric(unlist(aggregate(nodes3$weightsWtd, by=list(Category=nodes3$Var3), FUN=sum)["x"]))



### fxtheta distribution ###########

itemPara1 <- itemPara_A_SS[1:13,c("b", "a")]
itemPara2 <- itemPara_A_SS[14:25,c("b", "a")]
itemPara3 <- itemPara_A_SS[26:31,c("b", "a")]


fxTheta1 <- FX(itemPara1, newWts = newWts1)
fxTheta2 <- FX(itemPara2,newWts2)
fxTheta3 <- FX(itemPara3,newWts3)


names(fxTheta1) <- c(0:13)
names(fxTheta2) <- c(0:12)
names(fxTheta3) <- c(0:6)




tau <- c()
errvar <- c()
fyDistMat <- matrix(,,32)


for (k in 1:15){
  for (j in 1:15){
    for (i in 1:15){
      # extract
      # i <- 1
      # j <- 1
      # k <- 1

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


      fxThetaSum$y <- fxThetaSum$x1 + fxThetaSum$x2 + fxThetaSum$x3
      fxThetaSum$wty <- fxThetaSum$fx1 * fxThetaSum$fx2 * fxThetaSum$fx3

      fy <- fxThetaSum[,c("y", "wty")]

      fyDist <- aggregate(fy$wty, by=list(Category=fy$y), FUN=sum)

      names(fyDist) <- c("y", "wts")

      # weighted mean of Obs Y # true y
      weightedMean <- sum(fyDist$y * fyDist$wts)/sum(fyDist$wts)
      weightedMean

      # variance of Obs X
      varianceY <- sum(fyDist$wts * (fyDist$y - weightedMean)^2)
      varianceY

      tau <- c(tau,weightedMean)
      errvar <- c(errvar, varianceY)
      fyDistMat <- rbind(fyDistMat, t(fyDist$wts))

    }
  }
}

nodes3$tau <- tau
nodes3$errvar <- errvar
nodes3[,7:38] <- fyDistMat[-1,]

# sum of error variance
Var <- sum(nodes3$weightsWtd*nodes3$errvar)


# sum of observed score variance

numOfItem <- 31

fyThetaWeighted <- apply(nodes3[,7:(7 + numOfItem)], 2, function(x) x * nodes3[,"weightsWtd"])

# sum weighted distribution
# fxDist <- as.data.frame(matrix(colSums(fxThetaWeighted[,1:32]), nrow = 32, ncol = 1))
fyObsDist <- as.data.frame(matrix(colSums(fyThetaWeighted[,1:(1 + numOfItem)]), nrow = (1 + numOfItem), ncol = 1))
fyObsDist$X <- c(numOfItem:0)
names(fyObsDist) <- c("wts", "y")
fyObsDist

# weighted mean of Obs Y
weightedMean <- sum(fyObsDist$y * fyObsDist$wts)/sum(fyObsDist$wts)
weightedMean


# variance of Obs Y
varianceObsY <- sum(fyObsDist$wts * (fyObsDist$y - weightedMean)^2)
varianceObsY


1 - Var/varianceObsY







##### Test reliability Bifactor Full -------------------------------------------


FX_BF <- function(itemPara, newWts){

  # itemPara <- itemPara1
  # newWts <- nodes2[,c(3)]
  # transform item parameters to the logistic metric
  names(itemPara) <- c("b", "ag", "ai")

  # num of items
  numOfItem <- nrow(itemPara)

  # num of quadratures
  numOfQuad <- 10^2

  # weights and nodes
  # quadPoints <- NormalQuadraPoints(numOfQuad)

  quadPoints <- nodes2[,c(1,2)]

  # replicate item parameter and theta
  itemParaRep <-itemPara[rep(seq_len(numOfItem), each = numOfQuad),]
  itemParaRep$thetag <- rep(quadPoints[,c(1)], each = 1, length.out = numOfQuad*numOfItem)
  itemParaRep$thetai <- rep(quadPoints[,c(2)], each = 1, length.out = numOfQuad*numOfItem)

  # calculate information by theta
  itemParaRep <- within(itemParaRep, {
    P = 0 + (1 - 0) / (1 + exp(ag*thetag + ai*thetai -b))
    Q = 1 - P
    PQ = P * Q
    # info = 1.702**2 * a**2 * P * Q
  })

  ## true score variance

  # sum probability by theta
  # itemParaAggr <- aggregate(itemParaRep, by=list(Category=itemParaRep$theta), FUN=sum)

  itemParaAggr <- aggregate(.~thetag+thetai, itemParaRep, sum)

  # add weights for each theta
  # itemParaAggr$weights <- quadPoints$weights
  itemParaAggr$weights <- newWts

  # calculate true score variance
  varianceTrue <- sum((itemParaAggr$P)^2 * itemParaAggr$weights) - (sum(itemParaAggr$P * itemParaAggr$weights))^2


  ## error variance

  # claculate error variance
  varianceError <- sum(itemParaAggr$PQ * itemParaAggr$weights)


  ## observed score variance

  # order by theta
  itemParaRep <- itemParaRep[order(itemParaRep$thetag,itemParaRep$thetai),]

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

# item parameters

# read item parameters from txt file
itemPara_BF <- read.table("TestData/SpanishLit_prm_A_BF.txt")[,c(7:11)]
names(itemPara_BF) <- c("b", "ag","a1","a2", "a3") # a1 is primary
itemPara_BF$ai <- c(itemPara_BF$a1[1:13], itemPara_BF$a2[14:25], itemPara_BF$a3[26:31])

# num of items
numOfItem <- nrow(itemPara)

# num of quadratures
numOfQuad <- 10

library(LaplacesDemon)
library(mvtnorm)
library(pbivnorm)

# set nodes ranging from -5 to 5
nodes <- seq(-5, 5, length.out = numOfQuad)


nodes2 <- as.matrix(expand.grid(nodes,nodes))
weightsUnwtd <- dmvn(nodes2, c(0,0), diag(2), log=FALSE) # 41^3
nodes2 <- as.data.frame(nodes2)
nodes2$weightsWtd <- weightsUnwtd / sum(weightsUnwtd)

itemPara1 <- itemPara_BF[1:13,c("b", "ag", "ai")]
itemPara2 <- itemPara_BF[14:25,c("b", "ag", "ai")]
itemPara3 <- itemPara_BF[26:31,c("b", "ag", "ai")]

########################## FX function ------------------------

newWts <- nodes2[,c(3)]
fxTheta1 <- FX_BF(itemPara1, newWts)
fxTheta2 <- FX_BF(itemPara2,newWts)
fxTheta3 <- FX_BF(itemPara3,newWts)


names(fxTheta1) <- c(0:13)
names(fxTheta2) <- c(0:12)
names(fxTheta3) <- c(0:6)

library(tidyverse)
library(profvis)
library(data.table)

tau <- c()
errvar <- c()
fyDistMat <- matrix(NA,100^3,32)


n <- 0
for (k in 1:10^2){
  for (j in 1:10^2){
    for (i in 1:10^2){

      i <- 1
      j <- 1
      k <- 1
      n <- n+1
      fx1 <- t(fxTheta1[i,])
      fx2 <- t(fxTheta2[j,])
      fx3 <- t(fxTheta3[k,])

      xSum <- expand.grid(rownames(fx1), rownames(fx2), rownames(fx3))
      # names(xSum) <- c("x1", "x2", "x3")

      fxSum <- expand.grid(fx1, fx2, fx3)
      # names(fxSum) <- c("fx1", "fx2", "fx3")

      # fxThetaSum <- cbind(fxSum, xSum)

      xSum[,1] <- as.numeric(as.character(xSum[,1]))
      xSum[,2]<- as.numeric(as.character(xSum[,2]))
      xSum[,3] <- as.numeric(as.character(xSum[,3]))

      fxSum$y <- rowSums(xSum)

      fxSum$wty <- fxSum[,1] * fxSum[,2] * fxSum[,3]

      fy <- fxSum[,c("y", "wty")]

      fyDist <- setDT(fy)[,.(A = sum(wty)), by = 'y']

      names(fyDist) <- c("y", "wts")
      weightedMean <- sum(fyDist$y * fyDist$wts)/sum(fyDist$wts)

      tau[n] <- weightedMean
      errvar[n] <- sum(fyDist$wts * (fyDist$y - weightedMean)^2)
      fyDistMat[n,] <- t(fyDist$wts)

      # print(n)

    }
  }
}


nodes
nodes2 <- as.matrix(expand.grid(nodes,nodes))

# replicate item parameter and theta


Repeat1 <- function(d, n) {
  return(do.call("rbind", replicate(n, d, simplify = FALSE)))
}

nodes2rep1 <- Repeat1(nodes2, 10^4)

node2rep22 <- nodes2[rep(seq_len(nrow(nodes2)), each = 10^2),]
nodes2rep2 <- Repeat1(node2rep22, 10^2)

nodes2rep3 <- nodes2[rep(seq_len(nrow(nodes2)), each = 10^4),]

weights6 <- expand.grid(newWts,newWts,newWts)

nodes6 <- cbind(nodes2rep1, nodes2rep2, nodes2rep3,weights6)

names(nodes6) <- c("thetag1", "theta1", "thetag2", "theta2", "thetag3", "theta3", "wt1", "wt2", "wt3")


nodes6$wtSum <- nodes6$wt1 *nodes6$wt2 *nodes6$wt3

sum(nodes6$wtSum)







nodes6$tau <- tau
nodes6$errvar <- errvar
nodes6[,13:44] <- fyDistMat

# sum of error variance
nodes6$errvar_wt <- nodes6$wtSum * nodes6$errvar
var <- sum(nodes6$errvar_wt)
var

# Var <- sum(nodes6$wtSum*nodes6$errvar)


# sum of observed score variance

numOfItem <- 31

fyThetaWeighted <- apply(nodes6[,13:(13 + numOfItem)], 2, function(x) x * nodes6[,"wtSum"])

# sum weighted distribution
# fxDist <- as.data.frame(matrix(colSums(fxThetaWeighted[,1:32]), nrow = 32, ncol = 1))
fyObsDist <- as.data.frame(matrix(colSums(fyThetaWeighted[,1:(1 + numOfItem)]), nrow = (1 + numOfItem), ncol = 1))
fyObsDist$X <- c(numOfItem:0)

# fyObsDist$X <- c(0:numOfItem)
names(fyObsDist) <- c("wts", "y")
fyObsDist

# weighted mean of Obs Y
weightedMean <- sum(fyObsDist$y * fyObsDist$wts)/sum(fyObsDist$wts)
weightedMean


# variance of Obs Y
varianceObsY <- sum(fyObsDist$wts * (fyObsDist$y - weightedMean)^2)
varianceObsY


1 - var/varianceObsY


# 1 - varianceObsY/var















######## BI-Factor Full Approach 2 ------------------------------------


# read item parameters from txt file
itemPara_BF <- read.table("TestData/SpanishLit_prm_A_BF.txt")[,c(7:11)]
names(itemPara_BF) <- c("b", "ag","a1","a2", "a3") # a1 is primary
itemPara_BF$ai <- c(itemPara_BF$a1[1:13], itemPara_BF$a2[14:25], itemPara_BF$a3[26:31])

# num of items
numOfItem <- nrow(itemPara)

# num of quadratures
numOfQuad <- 10

library(LaplacesDemon)
library(mvtnorm)
library(pbivnorm)

# set nodes ranging from -5 to 5
nodes <- seq(-5, 5, length.out = numOfQuad)


nodes4 <- as.matrix(expand.grid(nodes,nodes,nodes, nodes))
weightsUnwtd <- dmvn(nodes2, c(0,0,0,0), diag(4), log=FALSE) # 41^3
nodes4 <- as.data.frame(nodes4)
nodes4$weightsWtd <- weightsUnwtd / sum(weightsUnwtd)

itemPara1 <- itemPara_BF[1:13,c("b", "ag", "ai")]
itemPara2 <- itemPara_BF[14:25,c("b", "ag", "ai")]
itemPara3 <- itemPara_BF[26:31,c("b", "ag", "ai")]







FX_BF <- function(itemPara, newWts, nodes[,c(thetag, theta1)]){

  # itemPara <- itemPara1
  # newWts <- nodes2[,c(3)]
  # transform item parameters to the logistic metric
  names(itemPara) <- c("b", "ag", "ai")

  # num of items
  numOfItem <- nrow(itemPara)

  # num of quadratures
  numOfQuad <- 10^2

  # weights and nodes
  # quadPoints <- NormalQuadraPoints(numOfQuad)

  quadPoints <- nodes

  # replicate item parameter and theta
  itemParaRep <-itemPara[rep(seq_len(numOfItem), each = numOfQuad),]
  itemParaRep$thetag <- rep(quadPoints[,c(1)], each = 1, length.out = numOfQuad*numOfItem)
  itemParaRep$thetai <- rep(quadPoints[,c(2)], each = 1, length.out = numOfQuad*numOfItem)

  # calculate information by theta
  itemParaRep <- within(itemParaRep, {
    P = 0 + (1 - 0) / (1 + exp(ag*thetag + ai*thetai -b))
    Q = 1 - P
    PQ = P * Q
    # info = 1.702**2 * a**2 * P * Q
  })

  ## true score variance

  # sum probability by theta
  # itemParaAggr <- aggregate(itemParaRep, by=list(Category=itemParaRep$theta), FUN=sum)

  itemParaAggr <- aggregate(.~thetag+thetai, itemParaRep, sum)

  # add weights for each theta
  # itemParaAggr$weights <- quadPoints$weights
  itemParaAggr$weights <- newWts

  # calculate true score variance
  varianceTrue <- sum((itemParaAggr$P)^2 * itemParaAggr$weights) - (sum(itemParaAggr$P * itemParaAggr$weights))^2


  ## error variance

  # claculate error variance
  varianceError <- sum(itemParaAggr$PQ * itemParaAggr$weights)


  ## observed score variance

  # order by theta
  itemParaRep <- itemParaRep[order(itemParaRep$thetag,itemParaRep$thetai),]

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





















