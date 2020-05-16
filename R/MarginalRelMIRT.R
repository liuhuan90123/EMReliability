# Marginal Reliability

#### simple structure MLE&EAP  P method ----------

## MLE ---------------------------------------

# number of factors
numOfFactors <- 3

# Form A
scoSS_MLE <- read.table("TestData/SpanishLit_sco_A_SS_MLE.txt")[,c(4:9, 10, 12, 15)]
cor <- c(0.91, # 1&2
         0.70, # 1&3
         0.49) # 2&3

# Form B
# scoSS_MLE <- read.table("TestData/SpanishLit_sco_B_SS_MLE.txt")[,c(4:9, 10, 12, 15)]
# cor <- c(0.97, # 1&2
#          0.56, # 1&3
#          0.48) # 2&3

# change variable name
names(scoSS_MLE) <- c("theta1", "theta2", "theta3", "se1", "se2", "se3", "var1", "var2", "var3")

# delete observations with missing values, 99.99 for flexMIRT output
scoSS_MLE[scoSS_MLE == 99.99] <- NA
scoSS_MLE <- na.omit(scoSS_MLE)

# composite theta
# scoSS_MLE$thetaC <- scoSS_MLE$theta1 + scoSS_MLE$theta2 + scoSS_MLE$theta3

# composite error variance
scoSS_MLE<- transform( scoSS_MLE,
                   varC = var1 + var2 + var3 + 2 *cor[1] * se1 * se2 + 2 *cor[2] * se1 * se3 + 2 *cor[3] * se2 * se3
)

# average of error variance
ErrorVarAvg <- mean(scoSS_MLE$varC)

# marginal reliability approach
MarginalRelMLE_SS <- ErrorVarAvg /(ErrorVarAvg + numOfFactors) # var(e)/(var(e) + var(theta))

# coefficients
MarginalRelMLE_SS


##### approach 2
# v1 <- mean(scoSS_MLE$var1)
# v2 <- mean(scoSS_MLE$var2)
# v3 <- mean(scoSS_MLE$var3)
#
# se1 <- sqrt(v1)
# se2 <- sqrt(v2)
# se3 <- sqrt(v3)
#
# vc <- v1+v2+v3+2 *cor[1] * se1 * se2 + 2 *cor[2] * se1 * se3 + 2 *cor[3] * se2 * se3
# vc
# vc/(vc+3)


####### EAP ---------------------------------------

# Form A
# scoSS_EAP <- read.table("TestData/SpanishLit_sco_A_SS_EAP.txt")[,c(3:14)]
# cor <- c(0.91, # 1&2
#          0.70, # 1&3
#          0.49) # 2&3

# Form B
scoSS_EAP <- read.table("TestData/SpanishLit_sco_B_SS_EAP.txt")[,c(3:14)]
cor <- c(0.97, # 1&2
         0.56, # 1&3
         0.48) # 2&3

# change variable name
names(scoSS_EAP) <- c("theta1", "theta2", "theta3", "se1", "se2", "se3", "var11", "var21", "var22", "var31","var32","var33")


# composite error variance
scoSS_EAP<- transform( scoSS_EAP,
                   varC = var11 + var22  + var33 +  2*(var21 + var31 + var32)
)

# average of error variance
ErrorVarAvg <- mean(scoSS_EAP$varC)

# marginal reliability approach
MarginalRelEAP_SS <- 1 - ErrorVarAvg /(2*(sum(cor)) + numOfFactors) # var(e)/(var(e) + var(theta))

# coefficients
MarginalRelEAP_SS



# CSEMIRT(0.330588, itemPara_A_SS[1:13, c("b", "a")], "EAP")
#
#
#
# se <- c(0.421373,	0.421900,	0.684679) #11,22,33
#
# var <- c(0.177556,	0.177999,	0.081154)#11,22,33
# cov <- c(0.153253, 0.034056,	0.468785)


### test of EAP SD: does not match flexMIRT output too

# Form B
itemPara_SS <- read.table("TestData/SpanishLit_prm_B_SS.txt")[,c(7:10)]

names(itemPara_SS) <- c("c","a1","a2","a3")
itemPara_SS$a <- c(itemPara_SS$a1[1:13], itemPara_SS$a2[14:25], itemPara_SS$a3[26:31])
itemPara_SS[,"b"] <- -itemPara_SS[,"c"]/itemPara_SS[,"a"]
itemPara_SS[,"a"] <- itemPara_SS[,"a"]/1.702

itemPara_SS$a1[1:13] <- itemPara_SS$a[1:13]
itemPara_SS$a2[14:25] <- itemPara_SS$a[14:25]
itemPara_SS$a3[26:31] <- itemPara_SS$a[26:31]

itemPara1 <- itemPara_SS[1:13,c("b", "a")]
itemPara2 <- itemPara_SS[14:25,c("b", "a")]
itemPara3 <- itemPara_SS[26:31,c("b", "a")]


# read raw data
# rawData <- read.table("TestData/FormA_31_3000.txt")
rawData <- read.table("TestData/FormB_31_3000.txt")




rawData1 <- rawData[,c(1:13)]
rawData1 <- rawData[,c(14:25)]
rawData1 <- rawData[,c(26:31)]

c <- 0

EAPTheta(itemPara1, rawData1, 15)



########### bi factor general D method ---------------------------------------


### bifactor General


# item parameters

# read item parameters from txt file
itemPara_BF_G <- read.table("TestData/SpanishLit_prm_A_BF.txt")[,c(7:11)]

# itemPara_BF_G <- read.table("SpanishLit-prm-B-BF.txt")[,c(7:11)]
# cor <- matrix(c(1, 0.9067069, 0.6994119,
#                 0.9067069, 1, 0.4891160,
#                 0.6994119,0.4891160,1), nrow = 3)

# form B
itemPara_BF_G <- read.table("TestData/SpanishLit_prm_B_BF.txt")[,c(7:11)]
#
# cor <- matrix(c(1, 0.97, 0.56,
#                 0.97, 1, 0.48,
#                 0.56,0.48,1), nrow = 3)



names(itemPara_BF_G) <- c("b", "ag","a1","a2", "a3") # a1 is primary
itemPara_BF_G$a <- c(itemPara_BF_G$a1[1:13], itemPara_BF_G$a2[14:25], itemPara_BF_G$a3[26:31])


itemPara_BF_G[,"b"] <- -itemPara_BF[,"b"]/(itemPara_BF[,"ag"] + itemPara_BF[,"a"]) #######?????????????????



itemPara_BF_G[,"a"] <- itemPara_BF_G[,"a"]/1.702

itemPara_BF_G$a1[1:13] <- itemPara_BF_G$a[1:13]
itemPara_BF_G$a2[14:25] <- itemPara_BF_G$a[14:25]
itemPara_BF_G$a3[26:31] <- itemPara_BF_G$a[26:31]


itemPara_BF_G$ag <- itemPara_BF_G$ag/1.702


itemPara_BF_G <- itemPara_BF_G[,c("b", "ag")]
names(itemPara_BF_G) <- c("b", "a")

library(EMReliability)
TestRelIRT(itemPara_BF_G)
MarginalRelIRT(itemPara_BF_G, estType = "MLE")
MarginalRelIRT(itemPara_BF_G, estType = "EAP")


########### bi factor full D method ---------------------------------------

# item parameters

# read item parameters from txt file
itemPara_BF <- read.table("TestData/SpanishLit_prm_B_BF.txt")[,c(7:11)]


names(itemPara_BF) <- c("b", "ag","a1","a2", "a3") # a1 is primary
itemPara_BF$a <- c(itemPara_BF$a1[1:13], itemPara_BF$a2[14:25], itemPara_BF$a3[26:31])


itemPara_BF[,"b"] <-  -itemPara_BF_G[,"b"]/(itemPara_BF_G[,"ag"] + itemPara_BF_G[,"a"])#######?????????????????



itemPara_BF[,"a"] <- itemPara_BF[,"a"]/1.702

itemPara_BF$a1[1:13] <- itemPara_BF$a[1:13]
itemPara_BF$a2[14:25] <- itemPara_BF$a[14:25]
itemPara_BF$a3[26:31] <- itemPara_BF$a[26:31]


itemPara_BF$ag <- itemPara_BF$ag/1.702


















### observed score variance

# rawvar1 <- var(rowSums(rawSS[,1:13]))
# rawvar2 <- var(rowSums(rawSS[,14:25]))
# rawvar3 <- var(rowSums(rawSS[,26:31]))
#
# rawvarC <- rawvar1 + rawvar2 + rawvar3 +
#   2 *cor[1] * sqrt(rawvar1) * sqrt(rawvar2) +
#   2 *cor[2] * sqrt(rawvar1)  * sqrt(rawvar3) +
#   2 *cor[3] * sqrt(rawvar2) * sqrt(rawvar3)

# # marginal reliability approach
# relmar <- varAvg/(varAvg+3) # var(e)/(var(e) + var(theta))

# # observed score approach 1 #var(e)/var(X_C)
# relObsVar1 <- 1-varAvg/rawvarC
#
# # observed score approach 1 #var(e)/var(X_C)
# relObsVar2 <- 1-varAvg/var(rowSums(rawSS))

# #coefficients
# relmar
# relObsVar1
# relObsVar2



EAPTheta <- function(itemPara, rawData, numOfQua =21){
  # itemPara <- itemPara_A
  # numOfQua <- 21
  # resp <- rawData

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

  for (i in 1:nrow(rawData)){

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





c <- 0

EAPTheta(itemPara1, rawData1, 15)


thetaEAP <- -1.77677885
SDEAP <- 0.6031666




sum((NormalQuadraPoints(21)$nodes - thetaEAP)^2 * NormalQuadraPoints(21)$weights) * thetaEAP












