# library(LaplacesDemon)
library(profvis)
# source("R/NormalQuadraPoints.R")
# source("R/LordWingersky.R")
library(mvtnorm)

if(!require(mvtnorm)) install.packages("mvtnorm")
library(mvtnorm)

######## BI-Factor Full Approach 2 ------------------------------------


# read item parameters conversion tables from txt file
# Form A
itemPara_BF <- read.table("TestData/SpanishLit_prm_A_BF.txt")[,c(7:11)]
convTable_A <- read.csv("TestData/conversion_table_Form A.csv")
convTable_A <- convTable_A[1:32, c("RawScore", "roundedSS")]
names(convTable_A) <- c("y","roundedSS")
convTable <- convTable_A

# Form B
# itemPara_BF <- read.table("TestData/SpanishLit_prm_B_BF.txt")[,c(7:11)]
# convTable_B <- read.csv("TestData/conversion_table_Form B.csv")
# convTable_B <- convTable_B[1:32, c("RawScore", "roundedSS")]
# names(convTable_B) <- c("y","roundedSS")
# convTable <- convTable_B


strat <- c(13, 12, 6)

names(itemPara_BF) <- c("b", "ag","a1","a2", "a3") # ag is primary
itemPara_BF$ai <- c(itemPara_BF$a1[1:13], itemPara_BF$a2[14:25], itemPara_BF$a3[26:31])

# num of items
numOfItem <- nrow(itemPara_BF)

# num of quadratures
numOfQuad <- 11


# set nodes ranging from -4 to 4
nodes <- seq(-4, 4, length.out = numOfQuad)
nodesM <- as.matrix(expand.grid(nodes,nodes,nodes, nodes))
weightsUnwtd <- dmvnorm(nodesM, c(0,0,0,0), diag(4), log=FALSE) # 41^3
nodesM <- as.data.frame(nodesM)
nodesM$weightsWtd <- weightsUnwtd / sum(weightsUnwtd)

itemPara1 <- itemPara_BF[1:13, c("b", "ag", "ai")]
itemPara2 <- itemPara_BF[14:25, c("b", "ag", "ai")]
itemPara3 <- itemPara_BF[26:31, c("b", "ag", "ai")]

# itemPara <- itemPara1

FX_BF <- function(itemPara){

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

  # transform item parameters to the logistic metric
  names(itemPara) <- c("b", "ag", "ai")

  # num of items
  numOfItem <- nrow(itemPara)

  # num of quadratures
  # numOfQuad <- numOfQuad^2

  # weights and nodes
  quadPoints <- expand.grid(NormalQuadraPoints(numOfQuad)$nodes, NormalQuadraPoints(numOfQuad)$nodes)

  # quadPoints <- nodes

  # replicate item parameter and theta
  itemParaRep <-itemPara[rep(seq_len(numOfItem), each = numOfQuad^2),]
  itemParaRep$thetag <- rep(quadPoints[,c(1)], each = 1, length.out = numOfQuad^2*numOfItem)
  itemParaRep$thetai <- rep(quadPoints[,c(2)], each = 1, length.out = numOfQuad^2*numOfItem)

  # calculate information by theta
  itemParaRep <- within(itemParaRep, {
    P = 0 + (1 - 0) / (1 + exp(-(ag*thetag + ai*thetai + b)))
    Q = 1 - P
    PQ = P * Q
    # info = 1.702**2 * a**2 * P * Q
  })

  # order by theta
  itemParaRep <- itemParaRep[order(itemParaRep$thetag,itemParaRep$thetai),]

  # define matrix of marginal distribution of theta
  fxTheta <- matrix(NA, nrow = numOfQuad^2, ncol = numOfItem + 1) # 41 num of quadratures, 41 num of raw sxores

  # for loop to calculate fxTheta
  for (i in 1:numOfQuad^2){

    probs <- matrix(c(itemParaRep[(1 + numOfItem * (i - 1)):(numOfItem * i),]$P),
                    nrow = numOfItem, ncol = 1, byrow = FALSE)

    fxTheta[i, ] <- LordWingersky(probs)$probability

  }

  # transform to data frame
  fxTheta <- as.data.frame(fxTheta)

  fxTheta

}


# fxtheta distribution
fxTheta1 <- FX_BF(itemPara1)
fxTheta2 <- FX_BF(itemPara2)
fxTheta3 <- FX_BF(itemPara3)

quadPoints <- expand.grid(seq(-4, 4, length.out = numOfQuad), seq(-4, 4, length.out = numOfQuad)) # Var1 Var2
fxTheta1[,c(15,16)] <- quadPoints
fxTheta2[,c(14,15)]<- quadPoints
fxTheta3[,c(8,9)]<- quadPoints

names(fxTheta1) <- c(0:strat[1],"theta2", "thetag")
names(fxTheta2) <- c(0:strat[2],"theta3", "thetag")
names(fxTheta3) <- c(0:strat[3],"theta4", "thetag")

names(nodesM) <- c("thetag", "theta2", "theta3", "theta4", "weightsWtd")


nodesM <- merge(x = nodesM, y = fxTheta1, by = c("theta2", "thetag"), all.x = TRUE)
nodesM <- merge(x = nodesM, y = fxTheta2, by = c("theta3", "thetag"), all.x = TRUE)
nodesM <- merge(x = nodesM, y = fxTheta3, by = c("theta4", "thetag"), all.x = TRUE)



for(i in 1:nrow(nodesM)){
  # i <- 1
  # fx distribution
  fx1 <- t(nodesM[i, 6:(6+strat[1])])
  fx2 <- t(nodesM[i, (6+strat[1]+1):(6+strat[1]+strat[2]+1)])
  fx3 <- t(nodesM[i, (6+strat[1]+strat[2]+1+1):(6+strat[1]+strat[2]+1+strat[3]+1)])

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

  nodesM[i,c(44:75)] <- t(fySSDist$wts)

}


nodesM <- nodesM[with(nodesM, order(thetag, theta2, theta3, theta4)), ]


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
fyThetaWeighted <- apply(nodesM[,44:(44 + numOfItem)], 2, function(x) x * nodesM[,"weightsWtd"])

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





# MIRT test reliability
TestRelBFMIRT <- 1 - varianceError/(varianceError + varianceTrueY)
TestRelBFMIRTSS <- 1 - varianceErrorSS/(varianceErrorSS + varianceTrueYSS)

varianceError
varianceTrueY
varianceObsY
TestRelBFMIRT

varianceErrorSS
varianceTrueYSS
varianceObsYSS
TestRelBFMIRTSS




VarandRel <- as.data.frame(matrix(c(varianceError, varianceTrueY, varianceObsY, TestRelBFMIRT,
                                    varianceErrorSS, varianceTrueYSS, varianceObsYSS, TestRelBFMIRTSS
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

conditionalSEMs <- nodesM[,c("thetag", "theta2", "theta3", "theta4", "weightsWtd","weightedMean", "varianceY", "weightedMeanSS", "varianceYSS")]
names(conditionalSEMs) <- c("Thetag", "Theta2", "Theta3", "Theta4", "weights","Ex_Raw", "Raw_Variance", "Ex_Scale", "Scale_Variance")
conditionalSEMs$Raw_CSEM <- sqrt(conditionalSEMs$Raw_Variance)
conditionalSEMs$Scale_CSEM <- sqrt(conditionalSEMs$Scale_Variance)
conditionalSEMs <- conditionalSEMs[with(conditionalSEMs, order(Thetag, Theta2, Theta3, Theta4)), ]
rownames(conditionalSEMs) <- 1:nrow(conditionalSEMs)
conditionalSEMs

# return(list("Variance and Reliability" = VarandRel,
#             "Conditional SEMs" = conditionalSEMs[,c("Theta1", "Theta2", "Theta3", "weights",
#                                                     "Ex_Raw", "Ex_Scale", "Raw_CSEM", "Scale_CSEM")]))


# Form A
# > VarandRel
#                                           coefficient
# Overall error variance for raw scores      5.2842179
# True score variance for raw scores        15.3932347
# Observed score variance for raw scores    20.6774526
# Reliability for raw scores                 0.7444454
# Overall error variance for scale scores    4.4141167
# True score variance for scale scores      11.7296405
# Observed score variance for scale scores  16.1437571
# Reliability for scale scores               0.7265744

# Form B
# > VarandRel
#                                          coefficient
# Overall error variance for raw scores      5.2806124
# True score variance for raw scores        14.0073056
# Observed score variance for raw scores    19.2879179
# Reliability for raw scores                 0.7262218
# Overall error variance for scale scores    4.5397396
# True score variance for scale scores      11.0183305
# Observed score variance for scale scores  15.5580701
# Reliability for scale scores               0.7082068


# library(xlsx)
# write.xlsx(conditionalSEMs, "conditionalSEMs_BF_A.xlsx")


