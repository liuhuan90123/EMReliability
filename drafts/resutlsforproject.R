
### Results for project ------------------

# source("R/NormalQuadraPoints.R") # n set as 41
# source("R/Info.R")
# source("R/LordWingersky.R")
# source("R/CronbachAlpha.R") # raw data
# source("R/Feldt.R")
# source("R/MarginalRelMLE.R") # itemPara
# source("R/MarginalRelEAP.R") # itemPara
# source("R/TestRelIRT.R") # itemPara
# source("R/KolenRelIRT.R") # itemPara, convTable
# source("R/CSEMLord.R") # numer of item
# source("R/PolynomialMethod.R") # cssemDat, K=num of degree
# source("R/CSSEMPolynomial.R") # numOfItem, convTable, K
# source("R/CSEMLord.R") # number of item
# source("R/CSSEMMLEPoly.R") # itemPara, convTable, K
# source("R/CSSEMEAPPoly.R") # itemPara, convTable, K
# source("R/CSEMLord.R") # numer of item
# source("R/CSSEMBinomial.R") # numer of item
# source("R/CSSEMPolynomial.R") # numOfItem, convTable, K
# source("R/CSEMIRT.R") # itemPara
# source("R/CSSEMIRT.R") # itemPara, convTable

# load library
library(EMReliability)

# read raw data
rawData_A <- read.table("TestData/RawDataFormX.txt")
rawData_B <- read.table("TestData/RawDataFormY.txt")

# read item parameters from txt file
itemPara_A <- read.table("TestData/ItemParaFormX.txt")
names(itemPara_A) <- c("b", "a")
itemPara_A[,"a"] <- itemPara_A[,"a"]/1.702
itemPara_B <- read.table("TestData/ItemParaFormY.txt")
names(itemPara_B) <- c("b", "a")
itemPara_B[,"a"] <- itemPara_B[,"a"]/1.702

# read conversion tables
convTable_A <- read.csv("TestData/ConversionTableFormX.csv")
convTable_A$roundedSS <- round(convTable_A$unroundedSS)

convTable_B <- read.csv("TestData/ConversionTableFormY.csv")
convTable_B$roundedSS <- round(convTable_B$unroundedSS)

convTable_A_Poly <- convTable_A[,c("theta", "roundedSS")]
convTable_B_Poly <- convTable_B[,c("theta", "roundedSS")]


# test help functions ------------------------------------
NormalQuadraPoints(41)
LordWingersky(c(0.9,0.9,0.9))
Info(NormalQuadraPoints(41)$nodes, itemPara_A, "EAP")

# CronbachAlpha & GT
CronbachAlpha_A <- CronbachAlpha(rawData_A)
CronbachAlpha_A

CronbachAlpha_B <- CronbachAlpha(rawData_B)
CronbachAlpha_B

# Feldt
Feldt_A <- Feldt(rawData_A)
Feldt_A

Feldt_B <- Feldt(rawData_B)
Feldt_B

# test reliability IRT
TestRelIRT_A <- TestRelIRT(itemPara_A)
TestRelIRT_A

TestRelIRT_B <- TestRelIRT(itemPara_B)
TestRelIRT_B

# marginal reliability MLE
MarginalRelMLE_A <- MarginalRelIRT(itemPara_A, "MLE")
MarginalRelMLE_A

MarginalRelMLE_B <- MarginalRelIRT(itemPara_B, "MLE")
MarginalRelMLE_B

# marginal reliability EAP
MarginalRelEAP_A <- MarginalRelIRT(itemPara_A, "EAP")
MarginalRelEAP_A

MarginalRelEAP_B <- MarginalRelIRT(itemPara_B, "EAP")
MarginalRelEAP_B

# Kolen's method
KolenRelIRT_A <- KolenRelIRT(itemPara_A, convTable_A)
KolenRelIRT_A

KolenRelIRT_B <- KolenRelIRT(itemPara_B, convTable_B)
KolenRelIRT_B

source("R/RelIRTPoly_new.R") # itemPara
# Reliability for rounded SS using polynomial method
RelMLEPoly_A <- RelIRTPoly_new(itemPara_A, convTable_A_Poly, 20, "MLE")
RelMLEPoly_A

RelMLEPoly_B <- RelIRTPoly_new(itemPara_B, convTable_B_Poly, 20, "MLE")
RelMLEPoly_B

RelEAPPoly_A <- RelIRTPoly_new(itemPara_A, convTable_A_Poly, 20, "EAP")
RelEAPPoly_A
RelEAPPoly_B <- RelIRTPoly_new(itemPara_B, convTable_B_Poly, 20, "EAP")
RelEAPPoly_B


RelMLEPoly_A <- RelIRTPoly(itemPara_A, convTable_A_Poly, 20, "MLE", rawData_A)
RelMLEPoly_A

RelMLEPoly_B <- RelIRTPoly(itemPara_B, convTable_B_Poly, 20, "MLE", rawData_B)
RelMLEPoly_B

RelEAPPoly_A <- RelIRTPoly(itemPara_A, convTable_A_Poly, 20, "EAP", rawData_A)
RelEAPPoly_A
RelEAPPoly_B <- RelIRTPoly(itemPara_B, convTable_B_Poly, 20, "EAP", rawData_B)
RelEAPPoly_B




EAP <- merge(as.data.frame(RelEAPPoly_A), as.data.frame(RelEAPPoly_B), by = "kValue")
EAP <- merge(EAP, as.data.frame(RelMLEPoly_B),  by = "kValue")

write.csv(EAP, "EAP.csv")

write.csv(RelMLEPoly_A, "RelMLEPoly_A.csv")
write.csv(RelMLEPoly_B, "RelMLEPoly_B.csv")

### CSEM --------------------------------------------------------------

# CSEM Lord
csemLord <- CSEMLord(40)
csemLord

# CSEM MLE
csemMLE_A <- CSEMIRT(NormalQuadraPoints(41)$nodes, itemPara_A, "MLE")
csemMLE_A
csemMLE_B <- CSEMIRT(NormalQuadraPoints(41)$nodes, itemPara_B, "MLE")
csemMLE_B

# CSEM EAP
csemEAP_A <- CSEMIRT(NormalQuadraPoints(41)$nodes, itemPara_A, "EAP")
csemEAP_A
csemEAP_B <- CSEMIRT(NormalQuadraPoints(41)$nodes, itemPara_B, "EAP")
csemEAP_B

### CSSEM -------------------------------------------------------------
# CSSEM Binomial
csemBinomial_A <- CSSEMBinomial(40, convTable_A)
csemBinomial_A
csemBinomial_B <- CSSEMBinomial(40, convTable_B)
csemBinomial_B

# CSSEM Polynomial

# number of item
numOfItem <- 40

convTable_A_sub <- convTable_A[,c("rawScore", "roundedSS")]
convTable_B_sub <- convTable_B[,c("rawScore", "roundedSS")]

cssemPolynomial_A <- CSSEMPolynomial(numOfItem, convTable_A_sub, 20)
cssemPolynomial_A
cssemPolynomial_B <- CSSEMPolynomial(numOfItem, convTable_B_sub, 20)
cssemPolynomial_B

# CSSEM Kolen's Method
cssemKolen_A <- CSSEMKolen(itemPara_A, convTable_A)
cssemKolen_A
cssemKolen_B <- CSSEMKolen(itemPara_B, convTable_B)
cssemKolen_B

# CSSEM IRT MLE Polynomial

cssemMLEPoly_A <- CSSEMIRTPoly(itemPara_A, convTable_A_Poly, 20, "MLE")
cssemMLEPoly_A
cssemMLEPoly_B <- CSSEMIRTPoly(itemPara_B, convTable_B_Poly, 20, "MLE")
cssemMLEPoly_B

# CSSEM IRT EAP Polynomial

cssemEAPPoly_A <- CSSEMIRTPoly(itemPara_A, convTable_A_Poly, 20, "EAP")
cssemEAPPoly_A
cssemEAPPoly_B <- CSSEMIRTPoly(itemPara_B, convTable_B_Poly, 20, "EAP")
cssemEAPPoly_B




### CSSEM IRT Polynomial MLE New -----------------------------------------------




### method two --------------------------------------------------------------------


# theta
theta <- convTable_B$theta

# CSEM MLE
itemParaCSEM <- as.data.frame(CSEMIRT(theta, itemPara_B, "EAP"))

# merge data
itemParaCSEM <- merge(itemParaCSEM, convTable_B_Poly, by = "theta")

# change name to fit Polynomial Method function
names(itemParaCSEM) <- c("rawScore", "csem", "roundedSS")

# call PM function
cssemPolyMLE <- PolynomialMethod(itemParaCSEM, 20)



cssemDat <- itemParaCSEM
K <- 20

# create data frame to store r square
rSquaredDat <- as.data.frame(matrix(nrow = K, ncol = 1))
regCoefDat <- list()

# for loop to iterate different k
for (k in 1:K){

  # fit model with k
  modelK <- lm(roundedSS ~ poly(rawScore, k, raw=TRUE), cssemDat)

  # extract regression coefficients
  regCoef <- summary(modelK)$coefficients[, 1]
  regCoefDat[[k]] <- summary(modelK)$coefficients

  # extract r square coefficient
  rSquaredDat[k, 1]<- summary(modelK)$r.squared

  # check whether regression coefficient of highest order is missing
  if(is.na(regCoef[k+1])){

    message(paste("The maximum k accepted is", k-1, sep = " "))

    break

  }

  # calculate transformation coefficients fx: from 1 to K
  cssemDat$fx <- 0
  i <- k

  while(i > 1){

    cssemDat$fx <- cssemDat$fx +  regCoef[i+1] * (i * cssemDat$rawScore^(i-1))
    i <- i-1

  }

  cssemDat$fx <- cssemDat$fx + regCoef[i+1]

  # calculate cssem using polynomial method
  cssemDat$cssemPoly <- cssemDat$fx * cssemDat$csem

  # rename variable with indicator k
  names(cssemDat)[names(cssemDat) == 'cssemPoly'] <- paste("cssemPolyk", k, sep = "")

}

# return(list("RSquared" = as.matrix(rSquaredDat[1:k,]), "CSSEMPoly" = cssemDat))

regCoefDat[[4]]  ## get coefficients for k = 4



# # SS formula
#
# SS = 118.263733790 + 3.909028972 * x -0.360085771 * x^2 -0.036053623 * x^3 + 0.009014156 * x^4
#
# f=expression(118.263733790 + 3.909028972 * x -0.360085771 * x^2 -0.036053623 * x^3 + 0.009014156 * x^4)
# D(f,'x')
#
#
# # CSSEM formula
#
# CSSEM = 3.909028972 - 0.360085771 * (2 * x) - 0.036053623 * (3 * x^2) + 0.009014156 * (4 * x^3)

# SS formula



# names(itemParaCSEM) <- c("theta", 'csemMLE',"scaleScore")


# theta
theta_New <- NormalQuadraPoints(41)$nodes
weights_New <- NormalQuadraPoints(41)$weights

# CSEM MLE
itemParaCSEM_New <- as.data.frame(CSEMIRT(theta_New, itemPara_B, "EAP"))

# itemParaCSEM_New <- within(itemParaCSEM_New,{
#   scaleScoreNew =  118.35478428 + 3.92088450 * theta -0.47021240 * theta^2 -0.03740998 * theta^3 + 0.01328610 * theta^4
#   cssemNew = 3.92088450 -0.47021240 * (2 * csemMLE) -0.03740998  * (3 * csemMLE^2) + 0.01328610 * (4 * csemMLE^3)
#   roundedSS = round(scaleScoreNew)
# })

# form B
itemParaCSEM_New <- within(itemParaCSEM_New,{
  scaleScoreNew =  118.35478428 + 3.92088450 * theta -0.47021240 * theta^2 -0.03740998 * theta^3 + 0.01328610 * theta^4
  cssemNew = 3.92088450 -0.47021240 * (2 * csemEAP) -0.03740998  * (3 * csemEAP^2) + 0.01328610 * (4 * csemEAP^3)
  roundedSS = scaleScoreNew
})

# form A
# itemParaCSEM_New <- within(itemParaCSEM_New,{
#   scaleScoreNew =  118.263733790 + 3.909028972 * theta -0.360085771 * theta^2 -0.036053623 * theta^3 + 0.009014156 * theta^4
#   cssemNew = 3.909028972 -0.360085771 * (2 * csemMLE) -0.036053623  * (3 * csemMLE^2) + 0.009014156 * (4 * csemMLE^3)
#   roundedSS = round(scaleScoreNew)
# })

### read new posterior distribution rates

postDist <- read.table("TestData/PosteriorDistribution.txt", sep = " ")

plot(postDist$V1, postDist$V2)


itemParaCSEM_New$weights_new <- postDist$V2


# SS variance
SSVar <- sum(itemParaCSEM_New$weights_new * (itemParaCSEM_New$roundedSS - weighted.mean(itemParaCSEM_New$roundedSS, itemParaCSEM_New$weights_new))^2)

# error variance
errorVar <- sum(itemParaCSEM_New$cssemNew * itemParaCSEM_New$weights_new)

# reliability
RelMLEPolyNew <- 1 - errorVar/SSVar
RelMLEPolyNew





sum(itemParaCSEM_New$weights_new)




# change variable name
names(cssemPolyMLE$CSSEMPoly)[names(cssemPolyMLE$CSSEMPoly) == 'rawScore'] <- 'theta'
names(cssemPolyMLE$CSSEMPoly)[names(cssemPolyMLE$CSSEMPoly) == 'csem'] <- 'csemMLE'

# return results
return(list("RSquared" = cssemPolyMLE$RSquared, "CSSEMPolyMLE" = cssemPolyMLE$CSSEMPoly))















### Plot function
library(ggplot2)

# plot CSEM Lord

plot(csemLord$rawScore, csemLord$csemLord)

### Plot Kolen CSSEM  ---------------------------------------------------------------------------
cssemKolen_A <- as.data.frame(cssemKolen_A)
plot(cssemKolen_A$trueScaleScore, cssemKolen_A$cssemKolen)

png("CSSEM_KolenIRT_A.png",  width = 799, height = 596)

KA <- ggplot(cssemKolen_A, aes(x = trueScaleScore, y = cssemKolen)) +
  geom_point(size = 2) +
  scale_x_continuous(name = "True Scale Score", breaks  = seq(100, 130, 5)) +
  scale_y_continuous(name = "CSSEM Kolen's IRT Method", breaks  = seq(0, 3, 0.5),
                     limits = c(0,3)) +
  theme_bw()

print(KA)
dev.off()



### plot CSSEM Polynomial ---------------------------------------------------------------------------



#### Plot --------------------------------------
library(ggplot2)

### aggregate many to one ------------


cssemDatWide <- CSSEMPolynomial(40, convTableSub, K)$"CSSEM Polynomial Method"
cssemDat <- cssemDatWide # test
k <- 13 # test, accepted maximum + 1
cssemDat <- cssemDat[,c(3,5:(5+k-2))]
cssemDatAggre <- as.data.frame(apply(cssemDat[,c(-1)], 2, function(x) aggregate(x, by=list(Category=cssemDat$roundedSS), FUN=mean)))
cssemDatAggre <- cssemDatAggre[,c(1, seq(2, 24, 2))]



### plot all ks ---------------------------------------------------

k <- 13 # The maximum accepted K

cssemDatLong <- reshape(cssemDatAggre,
                        direction = "long",
                        varying = list(names(cssemDatAggre)[2:k]),
                        v.names = "cssempoly",
                        idvar = c("cssemPolyk1.Category"),
                        timevar = "Kvalue",
                        times = 1:(k-1))


library(ggplot2)
ggplot(cssemDatLong, aes(x = cssemPolyk1.Category, y = cssempoly, color = factor(Kvalue))) +
  geom_point() +
  scale_x_continuous(name = "Rounded Scale Score", breaks  = seq(100, 130, 5)) +
  scale_y_continuous(name = "CSSEM Polynomial Method") +
  geom_line() +
  theme_bw() +
  labs(colour="K value")



### range and confidence interval  with k from 1 to maximum accepted ----------

# moments
cssemDatAggre$max <- apply(cssemDatAggre[,c(-1)], 1, max)
cssemDatAggre$min <- apply(cssemDatAggre[,c(-1)], 1, min)
cssemDatAggre$mean <- apply(cssemDatAggre[,c(-1)], 1, mean)
cssemDatAggre$median <- apply(cssemDatAggre[,c(-1)], 1, median)
cssemDatAggre$sd <- apply(cssemDatAggre[,c(-1)], 1, sd)

# 95% confidence interval
cssemDatAggre_1 <- within(cssemDatAggre, {
  lower = mean - 1.96 * sd
  upper = mean + 1.96 * sd
})


### plot confidence interval --------------------------------------

# k = 1 to maximum accepted

ggplot(cssemDatAggre_1, aes(x = cssemPolyk1.Category, y = mean)) +
  geom_line(colour ="blue") +
  geom_point(colour="blue") +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill= "blue", alpha = 0.2) +
  scale_x_continuous(name = "Rounded Scale Score", breaks  = seq(100, 130, 5)) +
  scale_y_continuous(name = "CSSEM Polynomial Method") +
  theme_bw()

### range and confidence interval  with k from 3 to maximum accepted ---333333333333333333333333333333-------

# moments
cssemDatAggre$max <- apply(cssemDatAggre[,c(-1:-3)], 1, max)
cssemDatAggre$min <- apply(cssemDatAggre[,c(-1:-3)], 1, min)
cssemDatAggre$mean <- apply(cssemDatAggre[,c(-1:-3)], 1, mean)
cssemDatAggre$median <- apply(cssemDatAggre[,c(-1:-3)], 1, median)
cssemDatAggre$sd <- apply(cssemDatAggre[,c(-1:-3)], 1, sd)

# 95% confidence interval
cssemDatAggre <- within(cssemDatAggre, {
  lower = mean - 1.96 * sd
  upper = mean + 1.96 * sd
})


### plot confidence interval

# k = 1 to maximum accepted

ggplot(cssemDatAggre, aes(x = cssemPolyk1.Category, y = mean)) +
  geom_line(colour ="blue") +
  geom_point(colour="blue") +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill= "blue", alpha = 0.2) +
  scale_x_continuous(name = "Rounded Scale Score", breaks  = seq(100, 130, 5)) +
  scale_y_continuous(name = "CSSEM Polynomial Method") +
  theme_bw()






