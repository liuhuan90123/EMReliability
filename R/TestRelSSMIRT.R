



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





# FX(itemPara_A)


#####simple structure MIRT ##########

# item parameters

# read item parameters from txt file
itemPara_A_SS <- read.table("TestData/SpanishLit_prm_A_SS.txt")[,c(7:10)]
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



# read correlations

# cor <- read.table("TestData/SpanishLit_prm_A_SS_cor.txt", header = F)
# cor <- c(0.9067069, # 1&2
#          0.6994119, # 1&3
#          0.4891160) # 2&3


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


# # Two dimensions with no correlation
# nodesSS <- expand.grid(nodes,nodes)
#
# # unnormalized weights
# # weightsUnwtd <- apply(nodesSS, 1, FUN = function(x) dmvnorm(x))
# weightsUnwtd <- pbivnorm(nodesSS, cor[1])
# # M <- cor[1]
# # dmvnorm(nodesSS)
#
# # normalized weightes
# weightsWtd <- weightsUnwtd / sum(weightsUnwtd)
#
# # return nodes and normalized weights
# # return(list("nodes" = nodesSS, "weights" = weightsWtd))
# nodesSS$weights <- weightsWtd
#
# # pbivnorm(z[i], z[i], rho=rel)
# dmvnorm(c(3.50,	-2.50), sigma = M)


# Three dimensions with correlation matrix--------------------------------------------------------------


nodes3 <- as.matrix(expand.grid(nodes,nodes,nodes))


         #0.9067069, # 1&2
         #0.6994119, # 1&3
         #0.4891160) # 2&3

weightsUnwtd <- dmvn(nodes3, c(0,0,0), cor, log=FALSE) # 41^3
nodes3 <- as.data.frame(nodes3)
nodes3$weightsWtd <- weightsUnwtd / sum(weightsUnwtd)


### weights ---



newWts1 <- as.numeric(unlist(aggregate(nodes3$weightsWtd, by=list(Category=nodes3$Var1), FUN=sum)["x"]))

newWts2 <- as.numeric(unlist(aggregate(nodes3$weightsWtd, by=list(Category=nodes3$Var2), FUN=sum)["x"]))

newWts3 <- as.numeric(unlist(aggregate(nodes3$weightsWtd, by=list(Category=nodes3$Var3), FUN=sum)["x"]))







# sum(nodes3$weightsWtd)
# sum()

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



# Form B : 0.721
# Form A : 0.737

CSEMIRT(-1.880281, itemPara_A_SS[1:13, c("b", "a")], "MLE")




