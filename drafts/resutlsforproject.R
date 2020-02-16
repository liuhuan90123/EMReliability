
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

# Reliability for rounded SS using polynomial method
RelMLEPoly_A <- RelIRTPoly(itemPara_A, convTable_A_Poly, 20, "MLE", rawData_A)
RelMLEPoly_A
RelMLEPoly_B <- RelIRTPoly(itemPara_B, convTable_B_Poly, 20, "MLE", rawData_B)
RelMLEPoly_B

RelEAPPoly_A <- RelIRTPoly(itemPara_A, convTable_A_Poly, 20, "EAP", rawData_A)
RelEAPPoly_A
RelEAPPoly_B <- RelIRTPoly(itemPara_B, convTable_B_Poly, 20, "EAP", rawData_B)
RelEAPPoly_B

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






