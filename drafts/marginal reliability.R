### Marginal reliability ---

# read item parameters from txt file
itemPara_A <- read.table("TestData/ItemParaFormX.txt")
names(itemPara_A) <- c("b", "a")
itemPara_A[,"a"] <- itemPara_A[,"a"]/1.702
itemPara_B <- read.table("TestData/ItemParaFormY.txt")
names(itemPara_B) <- c("b", "a")
itemPara_B[,"a"] <- itemPara_B[,"a"]/1.702

## Theoretical approach --


MarginalRelIRT_2 <- function(itemPara, estType){

  # weights and nodes
  quadPoints <- NormalQuadraPoints(41)

  if (estType == "MLE"){

    # calculate info
    itemParaInfo <- as.data.frame(Info(quadPoints$nodes, itemPara, "MLE"))

    # add weights for each theta
    itemParaInfo$weights <- quadPoints$weights

    # weighted information
    itemParaInfo$infoWeighted <- itemParaInfo$weights * (1/itemParaInfo$infoMLE)

    # marginal reliability MLE

    ### Information approach
    ## Green & NiceW1
    marginalRelMLEGreen <- 1 - sum(itemParaInfo$infoWeighted)
    marginalRelMLEGreen
    # 0.851

    ## Cheng & Nice2
    marginalRelMLEcheng <- (sum(itemParaInfo$infoMLE*itemParaInfo$weights)-1)/sum(itemParaInfo$infoMLE*itemParaInfo$weights) # upper bound
    marginalRelMLEcheng
    # 0.872

    # return coefficient
    return(marginalRelMLE)


  }else if (estType == "MAP"){

    # calculate info
    itemParaInfo <- as.data.frame(Info(quadPoints$nodes, itemPara, "MAP"))

    # add weights for each theta
    itemParaInfo$weights <- quadPoints$weights

    # inverse of information
    itemParaInfo$infoMAPInv <- 1 / itemParaInfo$infoMAP

    # marginal reliability MAP

    ## Green & NiceW1
    marginalRelMAPGreen <- 1 - sum(itemParaInfo$infoMAPInv * itemParaInfo$weights)
    marginalRelMAPGreen
    # 0.874

    ##
    marginalRelMLE <- 1 - sum(itemParaInfo$infoWeighted)
    marginalRelMLE <- 1/(sum(itemParaInfo$infoWeighted) + 1)

    # return coefficient
    return(marginalRelMAP)

  }else{

    warning("MarginalRelMLE function only supports MLE and MAP estimation method!")

  }

}


itemPara <- itemPara_A

# weights and nodes
quadPoints <- NormalQuadraPoints(41)


# calculate info
itemParaInfo <- as.data.frame(Info(quadPoints$nodes, itemPara, "MLE"))
# itemParaInfo$infoMAP <- itemParaInfo$infoMLE + 1

# add weights for each theta
itemParaInfo$weights <- quadPoints$weights

# inverse information
itemParaInfo$infoMLEInv <- 1/itemParaInfo$infoMLE
# itemParaInfo$infoMAPInv <- 1/itemParaInfo$infoMAP


### Information approach

# marginal reliability MLE

## Green & NiceW1
marginalRelMLEGreen <- 1 - sum(itemParaInfo$infoMLEInv * itemParaInfo$weights)
marginalRelMLEGreen
# 0.851

## Nice2 # upper bound
marginalRelMLENice2 <- (sum(itemParaInfo$infoMLE*itemParaInfo$weights)-1)/sum(itemParaInfo$infoMLE*itemParaInfo$weights)
marginalRelMLENice2
# 0.872



# marginal reliability MAP

## Green & NiceW1
marginalRelMAPGreen <- 1 - sum(1/(itemParaInfo$infoMLE +1) * itemParaInfo$weights)
marginalRelMAPGreen
# 0.874

## Anderson
sum(itemParaInfo$infoMLE/(itemParaInfo$infoMLE+1) * itemParaInfo$weights)



## NiceW2 & Cheng # upper bound
marginalRelMAPNice <- sum(itemParaInfo$infoMLE*itemParaInfo$weights)/(sum(itemParaInfo$infoMLE*itemParaInfo$weights) + 1)
marginalRelMAPNice
# 0.887









# marginal reliability MLE
MarginalRelMLE_A <- MarginalRelIRT_2(itemPara_A, "MLE")
MarginalRelMLE_A

MarginalRelMLE_B <- MarginalRelIRT_2(itemPara_B, "MLE")
MarginalRelMLE_B

# marginal reliability MAP
MarginalRelMAP_A <- MarginalRelIRT_2(itemPara_A, "MAP")
MarginalRelMAP_A

MarginalRelMAP_B <- MarginalRelIRT_2(itemPara_B, "MAP")
MarginalRelMAP_B







Info <- function(theta, itemPara, estType){

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

  # replicate item parameter and theta
  itemParaRep <-itemPara[rep(seq_len(nrow(itemPara)), each = length(theta)),]
  itemParaRep$theta <- rep(theta, each = 1, length.out = length(theta)*nrow(itemPara))

  # calculate information by theta 3PL
  itemParaRep <- within(itemParaRep, {
    P = c + (1 - c) / (1 + exp(-1.702 * a * (theta - b)))
    Q = 1 - P
    PQ = P * Q
    info = 1.702**2 * a**2 * (Q/P) * (P-c)^2 / (1-c)^2
  })

  # sum information by theta
  itemParaInfo <- aggregate(itemParaRep$info, by=list(Category=itemParaRep$theta), FUN=sum)
  names(itemParaInfo) <- c("theta", "infoMLE")

  # calculate info for MAP
  itemParaInfo$infoMAP <- itemParaInfo$infoMLE + 1 # cobtributed by population

  # return info for each theta
  if (estType == "MLE"){

    return(list("theta" = itemParaInfo$theta, "infoMLE" = itemParaInfo$infoMLE))

  }else if (estType == "MAP"){

    return(list("theta" = itemParaInfo$theta, "infoMAP" = itemParaInfo$infoMAP))

  }else{

    warning("Info function only supports MLE and MAP estimation method!")

  }

}


## Empirical approach --

# Form A
thetaSEMLE_A <- read.table("TestData/UShistory_X_MLE-sco.txt")[,c(4,5)]
summary(thetaSEMLE_A)
thetaSEMLE_A[thetaSEMLE_A == 99.99] <- NA
thetaSEMLE_A <- na.omit(thetaSEMLE_A)
marginalRelMLE2 <- (var(thetaSEMLE_A$V4)-mean(thetaSEMLE_A$V5^2))/var(thetaSEMLE_A$V4)
marginalRelMLE2

thetaSEEAP_A <- read.table("TestData/UShistory_X_EAP-sco.txt")[,c(3,4)]
# marginalRelEAP2 <- (var(thetaSEEAP_A$V3)-mean(thetaSEEAP_A$V4^2))/var(thetaSEEAP_A$V3)
# marginalRelEAP2

marginalRelEAP2 <- var(thetaSEEAP_A$V3)/(var(thetaSEEAP_A$V3)+mean(thetaSEEAP_A$V4^2))
marginalRelEAP2

thetaSEMAP_A <- read.table("TestData/UShistory_X_MAP-sco.txt")[,c(4,5)]
# marginalRelMAP2 <- (var(thetaSEMAP_A$V4)-mean(thetaSEMAP_A$V5^2))/var(thetaSEMAP_A$V4)
# marginalRelMAP2
marginalRelMAP2 <- var(thetaSEMAP_A$V4)/(var(thetaSEMAP_A$V4)+mean(thetaSEMAP_A$V5^2))
marginalRelMAP2


# Form B
thetaSEMLE_B <- read.table("TestData/UShistory_Y_MLE-sco.txt")[,c(4,5)]
summary(thetaSEMLE_B)
thetaSEMLE_B[thetaSEMLE_B == 99.99] <- NA
thetaSEMLE_B <- na.omit(thetaSEMLE_B)
marginalRelMLE2 <- (var(thetaSEMLE_B$V4)-mean(thetaSEMLE_B$V5^2))/var(thetaSEMLE_B$V4)
marginalRelMLE2


thetaSEEAP_B <- read.table("TestData/UShistory_Y_EAP-sco.txt")[,c(3,4)]
# marginalRelEAP2 <- (var(thetaSEEAP_A$V3)-mean(thetaSEEAP_A$V4^2))/var(thetaSEEAP_A$V3)
# marginalRelEAP2

marginalRelEAP2 <- var(thetaSEEAP_B$V3)/(var(thetaSEEAP_B$V3)+mean(thetaSEEAP_B$V4^2))
marginalRelEAP2

thetaSEMAP_B <- read.table("TestData/UShistory_Y_MAP-sco.txt")[,c(4,5)]
# marginalRelMAP2 <- (var(thetaSEMAP_A$V4)-mean(thetaSEMAP_A$V5^2))/var(thetaSEMAP_A$V4)
# marginalRelMAP2
marginalRelMAP2 <- var(thetaSEMAP_B$V4)/(var(thetaSEMAP_B$V4)+mean(thetaSEMAP_B$V5^2))
marginalRelMAP2





## verification part 1, information and SE

# MLE
quadPoints <- thetaSEMLE_A$V4
itemParaInfo <- as.data.frame(Info(quadPoints, itemPara_A, "MLE"))
itemParaInfo$infoMLEInv <- sqrt(1 / itemParaInfo$infoMLE)
itemParaInfo


# MAP
quadPoints <- thetaSEMAP_A$V4
itemParaInfo <- as.data.frame(Info(quadPoints, itemPara_A, "MAP"))
itemParaInfo$infoMAPInv <- sqrt(1 / itemParaInfo$infoMAP)
itemParaInfo


## simulation

# thetaq
set.seed(202)
N <- 10000
thetaq <- rnorm(N,0,1)
var(thetaq)
# thetaq <- rep(c(-1.323946,  0.547630,  0.064133,  1.051180, -0.648138, -0.909544,  0.330588, -0.527889,  0.347623), each = 100)
# thetaq

itemPara_A$c <- 0
itemPara_B$c <- 0

prob3PL <- function(D, theta, a, b, c){
  c + (1 - c) / (1 + exp(-D * a * (theta - b)))
}

gen3PL <- function(itemPara, I, thetaq){

  J <- nrow(itemPara)
  resp <- matrix(nrow = I, ncol = J)
  prob <- matrix(nrow = I, ncol = J)
  gamma <- matrix(runif(n = I*J, min = 0, max = 1), nrow = I, ncol = J)
  # theta <- rnorm(I, 0, 1)

  for (i in 1:I){
    for(j in 1:J){
      prob[i, j] <- prob3PL(D = 1.702, thetaq[i], itemPara[j,"a"], itemPara[j, "b"], itemPara[j, "c"])
    }
  }
  resp <- ifelse(prob > gamma, 1, 0)
  return(resp)
}

respEAP <- gen3PL(itemPara_B, N, thetaq)
respEAP

# write.table(respEAP, "resp_simu_B.txt", col.names = F, row.names = F)

EAPTheta <- function(itemPara, rawData, numOfQua =41){

  # thetaq
  thetaq <- NormalQuadraPoints(numOfQua)$nodes

  itemParaRep <- itemPara[rep(seq_len(nrow(itemPara)), each = numOfQua),]
  itemParaRep$thetaq <- rep(thetaq, each = 1, length.out = numOfQua*nrow(itemPara))

  itemParaRep <- within(itemParaRep, {

    P = c + (1 - c) / (1 + exp(-1.702 * a * (thetaq - b)))
    Q = 1 - P
    PQ = P * Q
    info = 1.702**2 * a**2 * (Q/P) * (P-c)^2 / (1-c)^2
  })

  itemParaRep <- itemParaRep[order(itemParaRep$thetaq), ]

  summaryEAP <- matrix(nrow = nrow(rawData), ncol = 2)

  for (i in 1:nrow(rawData)){ # change to apply

    # i =1

    itemParaRep$resp <- rep(rawData[i,], each = 1, length.out = numOfQua*nrow(itemPara))
    itemParaRep$resp <- as.numeric(itemParaRep$resp)
    itemParaRep$respRev <- 1 - itemParaRep$resp

    itemParaRep$T <- itemParaRep$P * itemParaRep$resp + itemParaRep$Q * itemParaRep$respRev

    # sum T by thetaq
    itemParaT <- aggregate(itemParaRep$T, by=list(Category=itemParaRep$thetaq), FUN=prod)

    names(itemParaT) <- c("thetaq", "T")

    itemParaT$wt <- NormalQuadraPoints(numOfQua)$weights
    # itemParaT$wt <- newWts1

    nume <- sum(itemParaT$thetaq * itemParaT$T * itemParaT$wt)
    deno <- sum(itemParaT$T * itemParaT$wt)

    thetaEAP <- nume/deno
    summaryEAP[i, 1] <- thetaEAP

    SDEAP <- sqrt(sum((itemParaT$thetaq - thetaEAP)^2 * itemParaT$T * itemParaT$wt) / deno)
    summaryEAP[i, 2] <- SDEAP

  }
  return(summaryEAP)
}

EAPSD <- EAPTheta(itemPara_B, respEAP)
EAPSD
1 - mean(EAPSD[,2])

var(EAPSD[,1])/ (var(EAPSD[,1]) + mean(EAPSD[,2]))


## simulated data


# Form A
thetaSEMLE_A <- read.table("TestData/UShistory_X_MLE_simu-sco.txt")[,c(4,5)]
summary(thetaSEMLE_A)
thetaSEMLE_A[thetaSEMLE_A == 99.99] <- NA
thetaSEMLE_A <- na.omit(thetaSEMLE_A)
marginalRelMLE2 <- (var(thetaSEMLE_A$V4)-mean(thetaSEMLE_A$V5^2))/var(thetaSEMLE_A$V4)
marginalRelMLE2

thetaSEEAP_A <- read.table("TestData/UShistory_X_EAP_simu-sco.txt")[,c(3,4)]
# marginalRelEAP2 <- (var(thetaSEEAP_A$V3)-mean(thetaSEEAP_A$V4^2))/var(thetaSEEAP_A$V3)
# marginalRelEAP2

marginalRelEAP2 <- var(thetaSEEAP_A$V3)/(var(thetaSEEAP_A$V3)+mean(thetaSEEAP_A$V4^2))
marginalRelEAP2

thetaSEMAP_A <- read.table("TestData/UShistory_X_MAP_simu-sco.txt")[,c(4,5)]
# marginalRelMAP2 <- (var(thetaSEMAP_A$V4)-mean(thetaSEMAP_A$V5^2))/var(thetaSEMAP_A$V4)
# marginalRelMAP2
marginalRelMAP2 <- var(thetaSEMAP_A$V4)/(var(thetaSEMAP_A$V4)+mean(thetaSEMAP_A$V5^2))
marginalRelMAP2



# Form B
thetaSEMLE_B <- read.table("TestData/UShistory_Y_MLE_simu-sco.txt")[,c(4,5)]
summary(thetaSEMLE_B)
thetaSEMLE_B[thetaSEMLE_B == 99.99] <- NA
thetaSEMLE_B <- na.omit(thetaSEMLE_B)
marginalRelMLE2 <- (var(thetaSEMLE_B$V4)-mean(thetaSEMLE_B$V5^2))/var(thetaSEMLE_B$V4)
marginalRelMLE2

thetaSEEAP_B <- read.table("TestData/UShistory_Y_EAP_simu-sco.txt")[,c(3,4)]
# marginalRelEAP2 <- (var(thetaSEEAP_A$V3)-mean(thetaSEEAP_A$V4^2))/var(thetaSEEAP_A$V3)
# marginalRelEAP2

marginalRelEAP2 <- var(thetaSEEAP_B$V3)/(var(thetaSEEAP_B$V3)+mean(thetaSEEAP_B$V4^2))
marginalRelEAP2

thetaSEMAP_B <- read.table("TestData/UShistory_Y_MAP_simu-sco.txt")[,c(4,5)]
# marginalRelMAP2 <- (var(thetaSEMAP_A$V4)-mean(thetaSEMAP_A$V5^2))/var(thetaSEMAP_A$V4)
# marginalRelMAP2
marginalRelMAP2 <- var(thetaSEMAP_B$V4)/(var(thetaSEMAP_B$V4)+mean(thetaSEMAP_B$V5^2))
marginalRelMAP2
