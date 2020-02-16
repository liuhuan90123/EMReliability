
### Results for project ------------------


## IRT Test Reliability
# help functions
source("R/NormalQuadraPoints.R") # n set as 41

NormalQuadraPoints(41)

source("R/Info.R")

Info(NormalQuadraPoints(41)$nodes, itemPara_A, "MLE")

source("R/LordWingersky.R")

LordWingersky(c(0.9,0.9,0.9))


source("R/CronbachAlpha.R") # raw data
source("R/Feldt.R")

source("R/MarginalRelMLE.R") # itemPara
source("R/MarginalRelEAP.R") # itemPara
source("R/TestRelIRT.R") # itemPara


source("R/KolenRelIRT.R") # itemPara, convTable
source("R/CSEMLord.R") # numer of item
source("R/PolynomialMethod.R") # cssemDat, K=num of degree
source("R/CSSEMPolynomial.R") # numOfItem, convTable, K
source("R/CSEMLord.R") # number of item
source("R/CSSEMMLEPoly.R") # itemPara, convTable, K
source("R/CSSEMEAPPoly.R") # itemPara, convTable, K

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

# CronbachAlpha

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

# CSEM Lord
source("R/CSEMLord.R") # numer of item
csemLord <- CSEMLord(40)

# CSSEM Binomial
source("R/CSSEMBinomial.R") # numer of item
csemBinomial <- CSSEMBinomial(40, convTable_A)


# CSSEMPolynomial

source("R/CSSEMPolynomial.R") # numOfItem, convTable, K

# Note: 1. cssemDat should include rawScore, roundedSS, and Lordcsem
#       2. conversion table should include raw score and rounded scale score
# number of item
numOfItem <- 40

convTable_A_sub <- convTable_A[,c("rawScore", "roundedSS")]
convTable_B_sub <- convTable_B[,c("rawScore", "roundedSS")]

CSSEMPolynomial(numOfItem, convTable_A_sub, 15)
CSSEMPolynomial(numOfItem, convTable_B_sub, 15)



# CSSEMMLEPoly (change csemLord)

CSSEMMLEPoly(itemPara_A, convTable_A, 20)
CSSEMMLEPoly(itemPara_B, convTable_B, 20)

# CSSEMEAPPoly (change csemLord)

CSSEMEAPPoly(itemPara_A, convTable_A, 20)
CSSEMEAPPoly(itemPara_B, convTable_B, 20)



source("R/CSEMIRT.R") # itemPara


# CSEM IRT MLE&EAP
CSEMIRT(NormalQuadraPoints(41)$nodes, itemPara_A, "MLE")

CSEMIRT(NormalQuadraPoints(41)$nodes, itemPara_B, "EAP")


source("R/CSSEMIRT.R") # itemPara, convTable

CSSEMKolen(itemPara_A, convTable_A)

CSSEMKolen(itemPara_B, convTable_B)


##

convTable_A_Poly <- convTable_A[,c("theta", "roundedSS")]
convTable_B_Poly <- convTable_B[,c("theta", "roundedSS")]



CSSEMIRTPoly(itemPara, convTable_A_Poly, K, "EAP")

CSSEMIRTPoly(itemPara, convTable_A_Poly, K, "MLE")



### realiability for scale score


RelIRTPoly(itemPara_A, convTable_A_Poly, 20, "MLE", rawData_A)
RelIRTPoly(itemPara_B, convTable_B_Poly, 20, "MLE", rawData_B)

RelIRTPoly(itemPara_A, convTable_A_Poly, 20, "EAP", rawData_A)
RelIRTPoly(itemPara_B, convTable_B_Poly, 20, "EAP", rawData_B)

### Plot function


### Plot Kolen CSSEM  ---------------------------------------------------------------------------
cssemKolen <- as.data.frame(cssemKolen)


### true scale score ---------

cssemKolen$trueSS <- colSums(fxThetaTSS)[1:41]

names(cssemKolen) <- c("cssemKolen", "trueSS")

png("CSSEM_KolenIRT_A.png",  width = 799, height = 596)

library(ggplot2)

K <- ggplot(cssemKolen, aes(x = trueSS, y = cssemKolen)) +
  geom_point(size = 2) +
  scale_x_continuous(name = "True Scale Score", breaks  = seq(100,  130, 5)) +
  scale_y_continuous(name = "CSSEM_Kolen IRT Method") +
  theme_bw()

print(K)
dev.off()


write.csv(cssemKolen, "cssemKolen.csv")


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






