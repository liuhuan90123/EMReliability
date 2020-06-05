options(max.print=1000000)
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

sum(diag(cov(rawData_A)))



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


# criteria
# raw score
cor(rowSums(rawData_A),rowSums(rawData_B))

# SS(X)
rawScore_A <- as.data.frame(rowSums(rawData_A))
names(rawScore_A) <- "rawScore"
rawSS_A <- left_join(rawScore_A, convTable_A, by = "rawScore")

rawScore_B <- as.data.frame(rowSums(rawData_B))
names(rawScore_B) <- "rawScore"
rawSS_B <- left_join(rawScore_B, convTable_B, by = "rawScore")

cor(rawSS_A$roundedSS, rawSS_B$roundedSS)



# theta score
thetaMLE_A <- read.table("TestData/UShistory_X-sco.txt")[,4]
thetaMLE_B <- read.table("TestData/UShistory_Y-sco.txt")[,4]
cor(thetaMLE_A, thetaMLE_B)

thetaEAP_A <- read.table("TestData/UShistory_X_EAP-sco.txt")[,3]
thetaEAP_B <- read.table("TestData/UShistory_Y_EAP-sco.txt")[,3]
cor(thetaEAP_A, thetaEAP_B)

# SS(theta)

thetaToSS <- function(theta, convTable, itemPara){

  thetaMLE <- as.data.frame(theta)
  # default K
  K <- 20

  rSquaredDat <- as.data.frame(matrix(nrow = K, ncol = 1))

  # for loop to iterate different k
  for (k in 1:K){

    # fit model with k
    modelK <- lm(roundedSS ~ poly(theta, k, raw=TRUE), convTable)

    # extract regression coefficients
    regCoef <- summary(modelK)$coefficients[, 1]

    # extract r square coefficient
    rSquaredDat[k, 1]<- summary(modelK)$r.squared

    # check whether regression coefficient of highest order is missing
    if(is.na(regCoef[k+1])){

      message(paste("The maximum k accepted is", k-1, sep = " "))

      break

    }

    # calculate SS: from 1 to K
    thetaMLE$SS <- 0
    i <- k

    while(i > 0){

      thetaMLE$SS <- thetaMLE$SS +  regCoef[i+1] * thetaMLE$theta^(i)
      i <- i-1

    }

    thetaMLE$SS <- thetaMLE$SS + regCoef[i+1]
    thetaMLE$SS <- round(thetaMLE$SS)

    # rename variable with indicator k
    names(thetaMLE)[names(thetaMLE) == 'SS'] <- paste("SSk", k, sep = "")
  }
  thetaMLE
}

# MLE
thetaMLE_A <- read.table("TestData/UShistory_X-sco.txt")[,4]
thetaMLE_B <- read.table("TestData/UShistory_Y-sco.txt")[,4]

SS_MLE_A <- thetaToSS(thetaMLE_A, convTable_A, itemPara_A)
SS_MLE_B <- thetaToSS(thetaMLE_B, convTable_B, itemPara_B)

cor(SS_MLE_A$SSk4, SS_MLE_B$SSk7)


# EAP

thetaEAP_A <- read.table("TestData/UShistory_X_EAP-sco.txt")[,3]
thetaEAP_B <- read.table("TestData/UShistory_Y_EAP-sco.txt")[,3]

SS_EAP_A <- thetaToSS(thetaEAP_A, convTable_A, itemPara_A)
SS_EAP_B <- thetaToSS(thetaEAP_B, convTable_B, itemPara_B)

cor(SS_EAP_A$SSk4, SS_EAP_B$SSk7)




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

# read conversion tables
convTable_A <- read.csv("TestData/ConversionTableFormX.csv")
convTable_A$roundedSS <- round(convTable_A$unroundedSS)

convTable_B <- read.csv("TestData/ConversionTableFormY.csv")
convTable_B$roundedSS <- round(convTable_B$unroundedSS)


KolenRelIRT_A <- KolenRelIRT(itemPara_A, convTable_A)
KolenRelIRT_A

KolenRelIRT_B <- KolenRelIRT(itemPara_B, convTable_B)
KolenRelIRT_B

# Polynomial method reliability for EAP and MLE
# method 2
source("R/RelIRTPoly.R")
RelIRTPoly_MLE_A <- RelIRTPoly(itemPara_A, convTable_A, 10, "MLE")
RelIRTPoly_MLE_A
RelIRTPoly_MLE_B <- RelIRTPoly(itemPara_B, convTable_B, 10, "MLE")
RelIRTPoly_MLE_B
RelIRTPoly_EAP_A <- RelIRTPoly(itemPara_A, convTable_A, 10, "EAP")
RelIRTPoly_EAP_A
RelIRTPoly_EAP_B <- RelIRTPoly(itemPara_B, convTable_B, 10, "EAP")
RelIRTPoly_EAP_B



## EAP theta


# EAPTheta(itemPara_A, rawData_A)
# EAPTheta(itemPara_B, rawData_B)

# EAP <- merge(as.data.frame(RelEAPPoly_A), as.data.frame(RelEAPPoly_B), by = "kValue")
# EAP <- merge(EAP, as.data.frame(RelMLEPoly_B),  by = "kValue")
#
# write.csv(EAP, "EAP.csv")
#
# write.csv(RelMLEPoly_A, "RelMLEPoly_A.csv")
# write.csv(RelMLEPoly_B, "RelMLEPoly_B.csv")

### CSEM

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

### CSSEM

# CSSEM Binomial
cssemBinomial_A <- CSSEMBinomial(40, convTable_A)
cssemBinomial_A
cssemBinomial_B <- CSSEMBinomial(40, convTable_B)
cssemBinomial_B

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


### Plot CSEMS -------------------------------------------------------------------------------------------------------

library(ggplot2)


# 1 CSEM Lord ------------------------------------------
# plot CSEM Lord
plot(csemLord$rawScore, csemLord$csemLord)

# ggplot
csemLord <- as.data.frame(csemLord)

png("TestData/csemLord.png",  width = 799, height = 596)

L <- ggplot(csemLord, aes(x = rawScore, y = csemLord)) +
  geom_point(size = 2) +
  scale_x_continuous(name = "Raw Score", breaks  = seq(0, 40, 5)) +
  scale_y_continuous(name = "CSEM Lord Method", breaks  = seq(0, 4, 0.5),
                     limits = c(0,4)) +
  theme_bw()

print(L)
dev.off()


### data for excel
# install.packages("xlsx")
library(xlsx)
write.xlsx(csemLord, "csemLord.xlsx")




# plot CSEM MLE/EAP
plot(csemMLE_A$theta, csemMLE_A$csemMLE)
plot(csemMLE_B$theta, csemMLE_B$csemMLE)

plot(csemEAP_A$theta, csemEAP_A$csemEAP)
plot(csemEAP_B$theta, csemEAP_B$csemEAP)



# 2 CSEM MLE ------------------------------------------

csemMLE_A <- as.data.frame(csemMLE_A)

png("TestData/csemMLE_A.png",  width = 799, height = 596)

csemMLE_A_P <- ggplot(csemMLE_A, aes(x = theta, y = csemMLE)) +
  geom_point(size = 2) +
  scale_x_continuous(name = "Theta", breaks  = seq(-5, 5, 1)) +
  scale_y_continuous(name = "CSEM MLE", breaks  = seq(0, 3, 0.5),
                     limits = c(0,3)) +
  theme_bw()

print(csemMLE_A_P)
dev.off()


csemMLE_B <- as.data.frame(csemMLE_B)

png("TestData/csemMLE_B.png",  width = 799, height = 596)

csemMLE_B_P <- ggplot(csemMLE_B, aes(x = theta, y = csemMLE)) +
  geom_point(size = 2) +
  scale_x_continuous(name = "Theta", breaks  = seq(-5, 5, 1)) +
  scale_y_continuous(name = "CSEM MLE", breaks  = seq(0, 3, 0.5),
                     limits = c(0,3)) +
  theme_bw()

print(csemMLE_B_P)
dev.off()

# AB
csemMLE <- merge(csemMLE_A, csemMLE_B, "theta")
csemMLE <- reshape2::melt(csemMLE, id=c("theta"))
levels(csemMLE$variable)[levels(csemMLE$variable)=="csemMLE.x"] <- "Form A"
levels(csemMLE$variable)[levels(csemMLE$variable)=="csemMLE.y"] <- "Form B"
names(csemMLE) <- c("theta", "Form", "csemMLE")

png("TestData/csemMLE_AB.png",  width = 799, height = 596)

csemMLE_AB_P <- ggplot(csemMLE, aes(x = theta, y = csemMLE, group = Form, shape = Form)) +
  geom_point(size = 2) +
  scale_x_continuous(name = "Theta", breaks  = seq(-5, 5, 1)) +
  scale_y_continuous(name = "CSEM MLE", breaks  = seq(0, 3, 0.5),
                     limits = c(0,3)) +
  theme_bw() +
  theme(legend.title=element_blank())

print(csemMLE_AB_P)
dev.off()



# 3 CSEM EAP ------------------------------------------

csemEAP_A <- as.data.frame(csemEAP_A)

png("TestData/csemEAP_A.png",  width = 799, height = 596)

csemEAP_A_P <- ggplot(csemEAP_A, aes(x = theta, y = csemEAP)) +
  geom_point(size = 2) +
  scale_x_continuous(name = "Theta", breaks  = seq(-5, 5, 1)) +
  scale_y_continuous(name = "CSEM EAP", breaks  = seq(0, 1, 0.2),
                     limits = c(0,1)) +
  theme_bw()

print(csemEAP_A_P)
dev.off()


csemEAP_B <- as.data.frame(csemEAP_B)

png("TestData/csemEAP_B.png",  width = 799, height = 596)

csemEAP_B_P <- ggplot(csemEAP_B, aes(x = theta, y = csemEAP)) +
  geom_point(size = 2) +
  scale_x_continuous(name = "Theta", breaks  = seq(-5, 5, 1)) +
  scale_y_continuous(name = "CSEM EAP", breaks  = seq(0, 1, 0.2),
                     limits = c(0,1)) +
  theme_bw()

print(csemEAP_B_P)
dev.off()

# AB
csemEAP <- merge(csemEAP_A, csemEAP_B, "theta")
csemEAP <- reshape2::melt(csemEAP, id=c("theta"))
levels(csemEAP$variable)[levels(csemEAP$variable)=="csemEAP.x"] <- "Form A"
levels(csemEAP$variable)[levels(csemEAP$variable)=="csemEAP.y"] <- "Form B"
names(csemEAP) <- c("theta", "Form", "csemEAP")

png("TestData/csemEAP_AB.png",  width = 799, height = 596)

csemEAP_AB_P <- ggplot(csemEAP, aes(x = theta, y = csemEAP, group = Form, shape = Form)) +
  geom_point(size = 2) +
  scale_x_continuous(name = "Theta", breaks  = seq(-5, 5, 1)) +
  scale_y_continuous(name = "CSEM EAP", breaks  = seq(0, 1, 0.2),
                     limits = c(0,1)) +
  theme_bw()+
  theme(legend.title=element_blank())

print(csemEAP_AB_P)
dev.off()


### MLE and EAP based on posterior tehta distribution

# MLE
thetaSEMLE_A <- read.table("TestData/UShistory_X-sco.txt")[,c(4,5)]
thetaSEMLE_B <- read.table("TestData/UShistory_Y-sco.txt")[,c(4,5)]
# thetaSEMLE_A <- thetaSEMLE_A[order(thetaSEMLE_A$V3),]
# thetaSEMLE_A <- thetaSEMLE_A[seq(1,3000,75),]
# thetaSEMLE_B <- thetaSEMLE_B[order(thetaSEMLE_B$V3),]
# thetaSEMLE_B <- thetaSEMLE_B[seq(1,3000,75),]

# AB
ggplot() +
  geom_point(thetaSEMLE_A, mapping = aes(x = V4, y = V5, shape = "Form A"), size = 1.5) +
  geom_point(thetaSEMLE_B, mapping = aes(x = V4, y = V5, shape = "Form B"), size = 1.5) +
  scale_x_continuous(name = "Theta", breaks  = seq(-5, 5, 1)) +
  scale_y_continuous(name = "CSEM MLE", breaks  = seq(0, 2, 0.2),
                     limits = c(0,2)) +
  theme_bw() +
  theme(legend.title=element_blank())

#EAP
thetaSEEAP_A <- read.table("TestData/UShistory_X_EAP-sco.txt")[,c(3,4)]
thetaSEEAP_B <- read.table("TestData/UShistory_Y_EAP-sco.txt")[,c(3,4)]
# thetaSEEAP_A <- thetaSEEAP_A[order(thetaSEEAP_A$V3),]
# thetaSEEAP_A <- thetaSEEAP_A[seq(1,3000,75),]
# thetaSEEAP_B <- thetaSEEAP_B[order(thetaSEEAP_B$V3),]
# thetaSEEAP_B <- thetaSEEAP_B[seq(1,3000,75),]

# AB
png("TestData/csemEAP_AB_post.png",  width = 799, height = 596)
csemEAP_pos_AB_P <- ggplot() +
  geom_point(thetaSEEAP_A, mapping = aes(x = V3, y = V4,  shape = "Form A"), size = 1.5) +
  geom_point(thetaSEEAP_B, mapping = aes(x = V3, y = V4,  shape = "Form B"), size = 1.5) +
  scale_x_continuous(name = "Theta", breaks  = seq(-5, 5, 1)) +
  scale_y_continuous(name = "CSEM EAP", breaks  = seq(0, 1, 0.2),
                     limits = c(0,1)) +
  theme_bw() +
  theme(legend.title=element_blank())

print(csemEAP_pos_AB_P)
dev.off()




### 23 for each form ----------------------------------------------------------

png("TestData/CSSEM_Theta_A.png",  width = 799, height = 596)
SEM_Theta_A <- ggplot() +
  geom_point(csemMLE_A, mapping = aes(x = theta, y = csemMLE,  shape = "CSEM MLE Information"), size = 2) +
  geom_point(csemEAP_A, mapping = aes(x = theta, y = csemEAP,  shape = "CSEM EAP Information"), size = 2) +
  geom_point(thetaSEEAP_A, mapping = aes(x = V3, y = V4,  shape = "CSEM EAP Posterior"), size = 2) +
  geom_line(csemMLE_A, mapping = aes(x = theta, y = csemMLE,  shape = "CSEM MLE Information"), size = 0.6) +
  geom_line(csemEAP_A, mapping = aes(x = theta, y = csemEAP,  shape = "CSEM EAP Information"), size = 0.6) +
  geom_line(thetaSEEAP_A, mapping = aes(x = V3, y = V4,  shape = "CSEM EAP Posterior"), size = 0.6) +
  scale_x_continuous(name = "Theta", breaks  = seq(-5, 5, 1)) +
  scale_y_continuous(name = "CSEM", breaks  = seq(0, 2.5, 0.5),
                     limits = c(0,2.5)) +
  theme_bw() +
  theme(legend.title=element_blank())

print(SEM_Theta_A)
dev.off()

write.xlsx(csemMLE_A, "csemtheta_A.xlsx", sheetName="csemMLE_A")
write.xlsx(csemEAP_A, "csemtheta_A.xlsx", sheetName="csemEAP_A", append=TRUE)
write.xlsx(thetaSEEAP_A, "csemtheta_A.xlsx", sheetName="thetaSEEAP_A", append=TRUE)





png("TestData/CSSEM_Theta_B.png",  width = 799, height = 596)
SEM_Theta_B <- ggplot() +
  geom_point(csemMLE_B, mapping = aes(x = theta, y = csemMLE,  shape = "CSEM MLE Information"), size = 2) +
  geom_point(csemEAP_B, mapping = aes(x = theta, y = csemEAP,  shape = "CSEM EAP Information"), size = 2) +
  geom_point(thetaSEEAP_B, mapping = aes(x = V3, y = V4,  shape = "CSEM EAP Posterior"), size = 2) +
  geom_line(csemMLE_B, mapping = aes(x = theta, y = csemMLE,  shape = "CSEM MLE Information"), size = 0.6) +
  geom_line(csemEAP_B, mapping = aes(x = theta, y = csemEAP,  shape = "CSEM EAP Information"), size = 0.6) +
  geom_line(thetaSEEAP_B, mapping = aes(x = V3, y = V4,  shape = "CSEM EAP Posterior"), size = 0.6) +
  scale_x_continuous(name = "Theta", breaks  = seq(-5, 5, 1)) +
  scale_y_continuous(name = "CSEM", breaks  = seq(0, 2.5, 0.5),
                     limits = c(0,2.5)) +
  theme_bw() +
  theme(legend.title=element_blank())

print(SEM_Theta_B)
dev.off()

write.xlsx(csemMLE_B, "csemtheta_B.xlsx", sheetName="csemMLE_B")
write.xlsx(csemEAP_B, "csemtheta_B.xlsx", sheetName="csemEAP_B", append=TRUE)
write.xlsx(thetaSEEAP_B, "csemtheta_B.xlsx", sheetName="thetaSEEAP_B", append=TRUE)


# png("TestData/CSSEM_Theta_B.png",  width = 799, height = 596)
# SEM_Theta_B <- ggplot() +
#   geom_point(csemMLE_A, mapping = aes(x = theta, y = csemMLE,  shape = "CSEM MLE Information", colour = "CSEM MLE Information"), size = 2) +
#   geom_point(csemEAP_A, mapping = aes(x = theta, y = csemEAP,  shape = "CSEM EAP Information", colour = "CSEM EAP Information"), size = 2) +
#   geom_point(thetaSEEAP_A, mapping = aes(x = V3, y = V4,  shape = "CSEM EAP Posterior", colour = "CSEM EAP Posterior"), size = 2) +
#   geom_line(csemMLE_A, mapping = aes(x = theta, y = csemMLE,  shape = "CSEM MLE Information", colour = "CSEM MLE Information"), size = 0.6) +
#   geom_line(csemEAP_A, mapping = aes(x = theta, y = csemEAP,  shape = "CSEM EAP Information", colour = "CSEM EAP Information"), size = 0.6) +
#   geom_line(thetaSEEAP_A, mapping = aes(x = V3, y = V4,  shape = "CSEM EAP Posterior", colour = "CSEM EAP Posterior"), size = 0.6) +
#   scale_x_continuous(name = "Theta", breaks  = seq(-5, 5, 1)) +
#   scale_y_continuous(name = "CSEM", breaks  = seq(0, 2.5, 0.5),
#                      limits = c(0,2.5)) +
#   theme_bw() +
#   theme(legend.title=element_blank())
#
# print(SEM_Theta_B)
# dev.off()



# 4 CSSEM Binomial ------------------------------------------


# plot cssem Binomial
plot(cssemBinomial_A$roundedSS, cssemBinomial_A$cssemBinomial)
plot(cssemBinomial_B$roundedSS, cssemBinomial_B$cssemBinomial)

# many to one aggr
cssemBinomial_A <- as.data.frame(cssemBinomial_A)
cssemBinomial_B <- as.data.frame(cssemBinomial_B)

cssemBinomial_A_Aggre <- aggregate(cssemBinomial_A$cssemBinomial, by=list(Category=cssemBinomial_A$roundedSS), FUN=function(x){sqrt(sum(x^2)/length(x))})
cssemBinomial_B_Aggre <- aggregate(cssemBinomial_B$cssemBinomial, by=list(Category=cssemBinomial_B$roundedSS), FUN=function(x){sqrt(sum(x^2)/length(x))})

names(cssemBinomial_A_Aggre) <- names(cssemBinomial_B_Aggre) <- c("roundedSS", "cssemBinomial")

plot(cssemBinomial_A_Aggre$roundedSS, cssemBinomial_A_Aggre$cssemBinomial)
plot(cssemBinomial_B_Aggre$roundedSS, cssemBinomial_B_Aggre$cssemBinomial)


# ggplot


png("TestData/CSSEM_Binomial_A.png",  width = 799, height = 596)

CSSEM_Binomial_A_P <- ggplot(cssemBinomial_A_Aggre, aes(x = roundedSS, y = cssemBinomial)) +
  geom_point(size = 1.5) +
  scale_x_continuous(name = "Rounded Scale Score", breaks  = seq(100, 130, 5)) +
  scale_y_continuous(name = "CSSEM Binomial Method", breaks  = seq(0, 3, 0.5),
                     limits = c(0,3)) +
  theme_bw()

print(CSSEM_Binomial_A_P)
dev.off()


png("TestData/CSSEM_Binomial_B.png",  width = 799, height = 596)

CSSEM_Binomial_B_P <- ggplot(cssemBinomial_B_Aggre, aes(x = roundedSS, y = cssemBinomial)) +
  geom_point(size = 1.5) +
  scale_x_continuous(name = "Rounded Scale Score", breaks  = seq(100, 130, 5)) +
  scale_y_continuous(name = "CSSEM Binomial Method", breaks  = seq(0, 3, 0.5),
                     limits = c(0,3)) +
  theme_bw()

print(CSSEM_Binomial_B_P)
dev.off()

# AB
png("TestData/CSSEM_Binomial_AB.png",  width = 799, height = 596)
CSSEM_Binomial_AB <- ggplot() +
  geom_point(cssemBinomial_A_Aggre, mapping = aes(x = roundedSS, y = cssemBinomial,  shape = "Form A"), size = 2) +
  geom_point(cssemBinomial_B_Aggre, mapping = aes(x = roundedSS, y = cssemBinomial,  shape = "Form B"), size = 2) +
  scale_x_continuous(name = "Rounded Scale Score", breaks  = seq(100, 130, 5)) +
  scale_y_continuous(name = "CSSEM Binomial Method", breaks  = seq(0, 3, 0.5),
                     limits = c(0,3)) +
  theme_bw() +
  theme(legend.title=element_blank())

print(CSSEM_Binomial_AB)
dev.off()


# 5 CSSEM polynomial method ------------------------------------------
library(ggplot2)
library(xlsx)
# form A
cssemDat_A <- CSSEMPolynomial(40, convTable_A_sub, 20)$"CSSEMPolynomial"
write.xlsx(cssemDat_A, "slope_A_polynomial_k3.xlsx", sheetName="slope_A_polynomial_k3")

k <- 13 # test, accepted maximum + 1

# form B
cssemDat_B <- CSSEMPolynomial(40, convTable_B_sub, 20)$"CSSEMPolynomial"
write.xlsx(cssemDat_B, "slope_B_polynomial_k3.xlsx", sheetName="slope_B_polynomial_k3")
k <- 13 # test, accepted maximum + 1


# many to one aggr
cssemDat_A <- as.data.frame(cssemDat_A)
cssemDat_B <- as.data.frame(cssemDat_B)

cssemDat_A_Aggre <- aggregate(cssemDat_A$cssemPolyk3, by=list(Category=cssemDat_A$roundedSS), FUN=function(x){sqrt(sum(x^2)/length(x))})
cssemDat_B_Aggre <- aggregate(cssemDat_B$cssemPolyk3, by=list(Category=cssemDat_B$roundedSS), FUN=function(x){sqrt(sum(x^2)/length(x))})

names(cssemDat_A_Aggre) <- names(cssemDat_B_Aggre) <- c("roundedSS", "cssemPolyk3")

plot(cssemDat_A_Aggre$roundedSS, cssemDat_A_Aggre$cssemPolyk3)
plot(cssemDat_B_Aggre$roundedSS, cssemDat_B_Aggre$cssemPolyk3)



# AB
png("TestData/CSSEM_Polynomial_AB.png",  width = 799, height = 596)
CSSEM_Polynomial_AB <- ggplot() +
  geom_point(cssemDat_A_Aggre, mapping = aes(x = roundedSS, y = cssemPolyk3,  shape = "Form A"), size = 2) +
  geom_point(cssemDat_B_Aggre, mapping = aes(x = roundedSS, y = cssemPolyk3,  shape = "Form B"), size = 2) +
  scale_x_continuous(name = "Rounded Scale Score", breaks  = seq(100, 130, 5)) +
  scale_y_continuous(name = "CSSEM Polynomial Method", breaks  = seq(0, 3, 0.5),
                     limits = c(0,3)) +
  theme_bw() +
  theme(legend.title=element_blank())

print(CSSEM_Polynomial_AB)
dev.off()



# 6 CSSEM Kolen's method ------------------------------------------


plot(cssemKolen_A$trueScaleScore, cssemKolen_A$cssemKolen)
plot(cssemKolen_B$trueScaleScore, cssemKolen_B$cssemKolen)

# ggplot
cssemKolen_A <- as.data.frame(cssemKolen_A)

png("TestData/CSSEM_KolenIRT_A.png",  width = 799, height = 596)

KA <- ggplot(cssemKolen_A, aes(x = trueScaleScore, y = cssemKolen)) +
  geom_point(size = 2) +
  scale_x_continuous(name = "True Scale Score", breaks  = seq(100, 130, 5)) +
  scale_y_continuous(name = "CSSEM Kolen's IRT Method", breaks  = seq(0, 3, 0.5),
                     limits = c(0,3)) +
  theme_bw()

print(KA)
dev.off()

cssemKolen_B <- as.data.frame(cssemKolen_B)

png("TestData/CSSEM_KolenIRT_B.png",  width = 799, height = 596)

KB <- ggplot(cssemKolen_B, aes(x = trueScaleScore, y = cssemKolen)) +
  geom_point(size = 2) +
  scale_x_continuous(name = "True Scale Score", breaks  = seq(100, 130, 5)) +
  scale_y_continuous(name = "CSSEM Kolen's IRT Method", breaks  = seq(0, 3, 0.5),
                     limits = c(0,3)) +
  theme_bw()

print(KB)
dev.off()


# AB
png("TestData/CSSEM_KolenIRT_AB.png",  width = 799, height = 596)
CSSEM_KolenIRT_AB <- ggplot() +
  geom_point(cssemKolen_A, mapping = aes(x = trueScaleScore, y = cssemKolen,  shape = "Form A"), size = 2) +
  geom_point(cssemKolen_B, mapping = aes(x = trueScaleScore, y = cssemKolen,  shape = "Form B"), size = 2) +
  scale_x_continuous(name = "True Scale Score", breaks  = seq(100, 130, 5)) +
  scale_y_continuous(name = "CSSEM Kolen's IRT Method", breaks  = seq(0, 3, 0.5),
                     limits = c(0,3)) +
  theme_bw() +
  theme(legend.title=element_blank())

print(CSSEM_KolenIRT_AB)
dev.off()



##### aggregate cssem Kolen's CSSEM

cssemKolen_A <- as.data.frame(cssemKolen_A)
cssemKolen_A$roundedSS <- round(cssemKolen_A$trueScaleScore)

cssemKolen_B <- as.data.frame(cssemKolen_B)
cssemKolen_B$roundedSS <- round(cssemKolen_B$trueScaleScore)

# many to one aggr
cssemKolen_A_Aggre <- aggregate(cssemKolen_A$cssemKolen, by=list(Category=cssemKolen_A$roundedSS), FUN=function(x){sqrt(sum(x^2)/length(x))})
cssemKolen_B_Aggre <- aggregate(cssemKolen_B$cssemKolen, by=list(Category=cssemKolen_B$roundedSS), FUN=function(x){sqrt(sum(x^2)/length(x))})

names(cssemKolen_A_Aggre) <- names(cssemKolen_B_Aggre) <- c("roundedSS", "cssemKolen")

plot(cssemKolen_A_Aggre$roundedSS, cssemKolen_A_Aggre$cssemKolen)
plot(cssemKolen_B_Aggre$roundedSS, cssemKolen_B_Aggre$cssemKolen)


write.xlsx(cssemKolen_A_Aggre, "cssemKolen_A_Aggre.xlsx")
write.xlsx(cssemKolen_B_Aggre, "cssemKolen_B_Aggre.xlsx")











### 456 for each form ----------------------------------------------------------


# point + line + bw
png("TestData/CSSEM_SSX_A.png",  width = 799, height = 596)
SSEM_SSX_A <- ggplot() +
  geom_point(cssemKolen_A, mapping = aes(x = trueScaleScore, y = cssemKolen,  shape = "CSSEM Kolen's Method"), size = 2) +
  geom_point(cssemDat_A_Aggre, mapping = aes(x = roundedSS, y = cssemPolyk3,  shape = "CSSEM Polynomial Method"), size = 2) +
  geom_point(cssemBinomial_A_Aggre, mapping = aes(x = roundedSS, y = cssemBinomial,  shape = "CSSEM Binomial Method"), size = 2) +
  geom_line(cssemKolen_A, mapping = aes(x = trueScaleScore, y = cssemKolen,  shape = "CSSEM Kolen's Method"), size = 0.6) +
  geom_line(cssemDat_A_Aggre, mapping = aes(x = roundedSS, y = cssemPolyk3,  shape = "CSSEM Polynomial Method"), size = 0.6) +
  geom_line(cssemBinomial_A_Aggre, mapping = aes(x = roundedSS, y = cssemBinomial,  shape = "CSSEM Binomial Method"), size = 0.6) +
  scale_x_continuous(name = "Scale Score", breaks  = seq(100, 130, 5)) +
  scale_y_continuous(name = "CSSEM", breaks  = seq(0, 3, 0.5),
                     limits = c(0,3)) +
  theme_bw() +
  theme(legend.title=element_blank())

print(SSEM_SSX_A)
dev.off()


write.xlsx(cssemKolen_A, "csemSSX_A.xlsx", sheetName="cssemKolen_A")
write.xlsx(cssemDat_A_Aggre, "csemSSX_A.xlsx", sheetName="cssemDat_A_Aggre", append=TRUE)
write.xlsx(cssemBinomial_A_Aggre, "csemSSX_A.xlsx", sheetName="cssemBinomial_A_Aggre", append=TRUE)






png("TestData/CSSEM_SSX_B.png",  width = 799, height = 596)
SSEM_SSX_B <- ggplot() +
  geom_point(cssemKolen_B, mapping = aes(x = trueScaleScore, y = cssemKolen,  shape = "CSSEM Kolen's Method"), size = 2) +
  geom_point(cssemDat_B_Aggre, mapping = aes(x = roundedSS, y = cssemPolyk3,  shape = "CSSEM Polynomial Method"), size = 2) +
  geom_point(cssemBinomial_B_Aggre, mapping = aes(x = roundedSS, y = cssemBinomial,  shape = "CSSEM Binomial Method"), size = 2) +
  geom_line(cssemKolen_B, mapping = aes(x = trueScaleScore, y = cssemKolen,  shape = "CSSEM Kolen's Method"), size = 0.6) +
  geom_line(cssemDat_B_Aggre, mapping = aes(x = roundedSS, y = cssemPolyk3,  shape = "CSSEM Polynomial Method"), size = 0.6) +
  geom_line(cssemBinomial_B_Aggre, mapping = aes(x = roundedSS, y = cssemBinomial,  shape = "CSSEM Binomial Method"), size = 0.6) +
  scale_x_continuous(name = "Scale Score", breaks  = seq(100, 130, 5)) +
  scale_y_continuous(name = "CSSEM", breaks  = seq(0, 3, 0.5),
                     limits = c(0,3)) +
  theme_bw() +
  theme(legend.title=element_blank())

print(SSEM_SSX_B)
dev.off()


write.xlsx(cssemKolen_B, "csemSSX_B.xlsx", sheetName="cssemKolen_B")
write.xlsx(cssemDat_B_Aggre, "csemSSX_B.xlsx", sheetName="cssemDat_B_Aggre", append=TRUE)
write.xlsx(cssemBinomial_B_Aggre, "csemSSX_B.xlsx", sheetName="cssemBinomial_B_Aggre", append=TRUE)



# point  + bw
png("TestData/CSSEM_SSX_A_noline.png",  width = 799, height = 596)
SSEM_SSX_A_noline <- ggplot() +
  geom_point(cssemKolen_A, mapping = aes(x = trueScaleScore, y = cssemKolen,  shape = "CSSEM Kolen's Method"), size = 2) +
  geom_point(cssemDat_A_Aggre, mapping = aes(x = roundedSS, y = cssemPolyk3,  shape = "CSSEM Polynomial Method"), size = 2) +
  geom_point(cssemBinomial_A_Aggre, mapping = aes(x = roundedSS, y = cssemBinomial,  shape = "CSSEM Binomial Method"), size = 2) +
  # geom_line(cssemKolen_A, mapping = aes(x = trueScaleScore, y = cssemKolen,  shape = "CSSEM Kolen's Method"), size = 1) +
  # geom_line(cssemDat_A_Aggre, mapping = aes(x = roundedSS, y = cssemPolyk3,  shape = "CSSEM Polynomial Method"), size = 1) +
  # geom_line(cssemBinomial_A_Aggre, mapping = aes(x = roundedSS, y = cssemBinomial,  shape = "CSSEM Binomial Method"), size = 1) +
  scale_x_continuous(name = "Scale Score", breaks  = seq(100, 130, 5)) +
  scale_y_continuous(name = "CSSEM", breaks  = seq(0, 3, 0.5),
                     limits = c(0,3)) +
  theme_bw() +
  theme(legend.title=element_blank())

print(SSEM_SSX_A_noline)
dev.off()



# 7/8 CSSEM MLE/EAP polynomial method ------------------------------------------

# plot cssem IRT polynomial method

# CSSEM IRT MLE Polynomial
cssemDat_MLE_A <- CSSEMIRTPoly(itemPara_A, convTable_A_Poly, 20, "MLE")$"CSSEMPolyMLE"
write.xlsx(cssemDat_MLE_A, "slope_A_polynomial_k4.xlsx", sheetName="slope_A_polynomial_k4")

# k <- 8 # test, accepted maximum + 1 # when ploting, set k = 10 by judgement
cssemDat_MLE_B <- CSSEMIRTPoly(itemPara_B, convTable_B_Poly, 20, "MLE")$"CSSEMPolyMLE"
write.xlsx(cssemDat_MLE_B, "slope_B_polynomial_k7.xlsx", sheetName="slope_B_polynomial_k7")

# k <- 8 # test, accepted maximum + 1 # when ploting, set k = 10 by judgement


# many to one aggr
cssemDat_MLE_A <- as.data.frame(cssemDat_MLE_A)
cssemDat_MLE_B <- as.data.frame(cssemDat_MLE_B)

cssemDat_MLE_A_Aggre <- aggregate(cssemDat_MLE_A$cssemPolyk4, by=list(Category=cssemDat_MLE_A$roundedSS), FUN=function(x){sqrt(sum(x^2)/length(x))})
cssemDat_MLE_B_Aggre <- aggregate(cssemDat_MLE_B$cssemPolyk7, by=list(Category=cssemDat_MLE_B$roundedSS), FUN=function(x){sqrt(sum(x^2)/length(x))})

names(cssemDat_MLE_A_Aggre) <- names(cssemDat_MLE_B_Aggre) <- c("roundedSS", "cssemPoly")

plot(cssemDat_MLE_A_Aggre$roundedSS, cssemDat_MLE_A_Aggre$cssemPoly)
plot(cssemDat_MLE_B_Aggre$roundedSS, cssemDat_MLE_B_Aggre$cssemPoly)



# AB
png("TestData/CSSEM_Polynomial_MLE_AB.png",  width = 799, height = 596)
CSSEM_Polynomial_MLE_AB <- ggplot() +
  geom_point(cssemDat_MLE_A_Aggre, mapping = aes(x = roundedSS, y = cssemPoly,  shape = "Form A"), size = 2) +
  geom_point(cssemDat_MLE_B_Aggre, mapping = aes(x = roundedSS, y = cssemPoly,  shape = "Form B"), size = 2) +
  scale_x_continuous(name = "Rounded Scale Score", breaks  = seq(100, 130, 5)) +
  scale_y_continuous(name = "CSSEM MLE Polynomial Method", breaks  = seq(0, 3, 0.5),
                     limits = c(0,3)) +
  theme_bw() +
  theme(legend.title=element_blank())

print(CSSEM_Polynomial_MLE_AB)
dev.off()




# CSSEM IRT EAP Polynomial
cssemDat_EAP_A <- CSSEMIRTPoly(itemPara_A, convTable_A_Poly, 20, "EAP")$"CSSEMPolyEAP"
# k <- 10 # test, accepted maximum + 1 # when ploting, set k = 10 by judgement
cssemDat_EAP_B <- CSSEMIRTPoly(itemPara_B, convTable_B_Poly, 20, "EAP")$"CSSEMPolyEAP"
# k <- 5 # test, accepted maximum + 1 # when ploting, set k = 10 by judgement


# many to one aggr
cssemDat_EAP_A <- as.data.frame(cssemDat_EAP_A)
cssemDat_EAP_B <- as.data.frame(cssemDat_EAP_B)

cssemDat_EAP_A_Aggre <- aggregate(cssemDat_EAP_A$cssemPolyk4, by=list(Category=cssemDat_EAP_A$roundedSS), FUN=function(x){sqrt(sum(x^2)/length(x))})
cssemDat_EAP_B_Aggre <- aggregate(cssemDat_EAP_B$cssemPolyk7, by=list(Category=cssemDat_EAP_B$roundedSS), FUN=function(x){sqrt(sum(x^2)/length(x))})

names(cssemDat_EAP_A_Aggre) <- names(cssemDat_EAP_B_Aggre) <- c("roundedSS", "cssemPoly")

plot(cssemDat_EAP_A_Aggre$roundedSS, cssemDat_EAP_A_Aggre$cssemPoly)
plot(cssemDat_EAP_B_Aggre$roundedSS, cssemDat_EAP_B_Aggre$cssemPoly)



# AB
png("TestData/CSSEM_Polynomial_EAP_AB.png",  width = 799, height = 596)
CSSEM_Polynomial_EAP_AB <- ggplot() +
  geom_point(cssemDat_EAP_A_Aggre, mapping = aes(x = roundedSS, y = cssemPoly,  shape = "Form A"), size = 2) +
  geom_point(cssemDat_EAP_B_Aggre, mapping = aes(x = roundedSS, y = cssemPoly,  shape = "Form B"), size = 2) +
  scale_x_continuous(name = "Rounded Scale Score", breaks  = seq(100, 130, 5)) +
  scale_y_continuous(name = "CSSEM EAP Polynomial Method", breaks  = seq(0, 3, 0.5),
                     limits = c(0,3)) +
  theme_bw() +
  theme(legend.title=element_blank())

print(CSSEM_Polynomial_EAP_AB)
dev.off()



# 9 EAP Polynomial method for posterior sd
# form A: 4, form B: 7


# SS(theta)

thetaToSS <- function(theta, convTable, itemPara){

  thetaMLE <- as.data.frame(theta)
  # default K
  K <- 10

  rSquaredDat <- as.data.frame(matrix(nrow = K, ncol = 1))

  # for loop to iterate different k
  for (k in 1:K){
    # k <- 4

    # fit model with k
    modelK <- lm(roundedSS ~ poly(theta, k, raw=TRUE), convTable)

    # extract regression coefficients
    regCoef <- summary(modelK)$coefficients[, 1]

    # extract r square coefficient
    rSquaredDat[k, 1]<- summary(modelK)$r.squared

    # check whether regression coefficient of highest order is missing
    if(is.na(regCoef[k+1])){

      message(paste("The maximum k accepted is", k-1, sep = " "))

      break

    }

    # calculate SS: from 1 to K
    thetaMLE$SS <- 0
    i <- k

    while(i > 0){

      thetaMLE$SS <- thetaMLE$SS +  regCoef[i+1] * thetaMLE$theta^(i)
      i <- i-1

    }

    thetaMLE$SS <- thetaMLE$SS + regCoef[i+1]
    thetaMLE$SS <- round(thetaMLE$SS)


    # calculate cssem: from 1 to K
    thetaMLE$fx <- 0
    i <- k

    while(i > 1){
      thetaMLE$fx <- thetaMLE$fx +  regCoef[i+1] * (i * thetaMLE$theta^(i-1))
      i <- i-1

    }

    thetaMLE$fx <- thetaMLE$fx + regCoef[i+1]
    thetaMLE$cssemPoly <- thetaMLE$fx * (thetaMLE$csem)


    # rename variable with indicator k
    names(thetaMLE)[names(thetaMLE) == 'SS'] <- paste("SSk", k, sep = "")
    names(thetaMLE)[names(thetaMLE) == 'cssemPoly'] <- paste("cssemPolyk", k, sep = "")
  }
  thetaMLE
}



# MLE
# thetaMLE_A <- read.table("TestData/UShistory_X-sco.txt")[,4]
# thetaMLE_B <- read.table("TestData/UShistory_Y-sco.txt")[,4]
#
# SS_MLE_A <- thetaToSS(thetaMLE_A, convTable_A, itemPara_A)
# SS_MLE_B <- thetaToSS(thetaMLE_B, convTable_B, itemPara_B)
#
# cor(SS_MLE_A$SSk4, SS_MLE_B$SSk7)


# EAP ----SS

thetaEAP_A <- read.table("TestData/UShistory_X_EAP-sco.txt")[,c(3,4)]
thetaEAP_A <- as.data.frame(thetaEAP_A)
names(thetaEAP_A) <- c("theta", "csem")
head(thetaEAP_A)


SS_EAP_A <- thetaToSS(thetaEAP_A, convTable_A_Poly, itemPara_A)
head(SS_EAP_A)

# many to one aggr

thetaEAP_A_Aggre <- aggregate(SS_EAP_A$cssemPolyk4, by=list(Category=SS_EAP_A$SSk4), FUN=function(x){sqrt(sum(x^2)/length(x))})
names(thetaEAP_A_Aggre)  <- c("roundedSS", "cssemPoly")







# Form B

thetaEAP_B <- read.table("TestData/UShistory_Y_EAP-sco.txt")[,c(3,4)]
thetaEAP_B <- as.data.frame(thetaEAP_B)
names(thetaEAP_B) <- c("theta", "csem")
head(thetaEAP_B)


SS_EAP_B <- thetaToSS(thetaEAP_B, convTable_B_Poly, itemPara_B)
head(SS_EAP_B)

# many to one aggr

thetaEAP_B_Aggre <- aggregate(SS_EAP_B$cssemPolyk7, by=list(Category=SS_EAP_B$SSk7), FUN=function(x){sqrt(sum(x^2)/length(x))})
names(thetaEAP_B_Aggre)  <- c("roundedSS", "cssemPoly")



# k <- 4
#
# modelK <- lm(roundedSS ~ poly(theta, k, raw=TRUE), convTable)
#
# # extract regression coefficients
# regCoef <- summary(modelK)$coefficients[, 1]
#
# # thetaEAP_A <- within(thetaEAP_A,
# #                      {csem = csem^2})
#
# thetaEAP_A <- within(thetaEAP_A,
#                      {SS = 0.009014156* theta^4 +   -0.036053623* theta^3 +
#                        -0.360085771* theta^2 + 3.909028972* theta^1 + 118.263733790})
#
#
# thetaEAP_A <- within(thetaEAP_A,
#                      {cssem = sqrt((4 * 0.009014156* theta^3 +   3*-0.036053623* theta^2 +
#                        2*-0.360085771* theta^1 + 3.909028972) * csem^2)})
#
#
# ggplot(thetaEAP_A, aes(x = SS, y = cssem)) +
#   geom_point() +
#   scale_x_continuous(name = "Rounded Scale Score", breaks  = seq(100, 130, 5), limits = c(100,130)) +
#   scale_y_continuous(name = "CSSEM Polynomial Method", breaks  = seq(0, 4.5, 0.5), limits = c(0,4.5))







## 78 for each form --------------------------

# Form A

png("TestData/CSSEM_Polynomial_SSTheta_A.png",  width = 799, height = 596)
CSSEM_Polynomial_SSTheta_A <- ggplot() +
  geom_point(cssemDat_MLE_A_Aggre, mapping = aes(x = roundedSS, y = cssemPoly,  shape = "MLE"), size = 2) +
  geom_point(cssemDat_EAP_A_Aggre, mapping = aes(x = roundedSS, y = cssemPoly,  shape = "EAP Information"), size = 2) +
  geom_point(thetaEAP_A_Aggre, mapping = aes(x = roundedSS, y = cssemPoly,  shape = "EAP Posterior"), size = 2) +
  geom_line(cssemDat_MLE_A_Aggre, mapping = aes(x = roundedSS, y = cssemPoly,  shape = "MLE"), size = 0.6) +
  geom_line(cssemDat_EAP_A_Aggre, mapping = aes(x = roundedSS, y = cssemPoly,  shape = "EAP Information"), size = 0.6) +
  geom_line(thetaEAP_A_Aggre, mapping = aes(x = roundedSS, y = cssemPoly,  shape = "EAP Posterior"), size = 0.6) +
  scale_x_continuous(name = "Rounded Scale Score", breaks  = seq(100, 130, 5)) +
  scale_y_continuous(name = "CSSEM Polynomial Method", breaks  = seq(0, 3, 0.5),
                     limits = c(0,3)) +
  theme_bw() +
  theme(legend.title=element_blank())

print(CSSEM_Polynomial_SSTheta_A)
dev.off()


write.xlsx(cssemDat_MLE_A_Aggre, "csemSSTheta_A.xlsx", sheetName="cssemDat_MLE_A_Aggre")
write.xlsx(cssemDat_EAP_A_Aggre, "csemSSTheta_A.xlsx", sheetName="cssemDat_EAP_A_Aggre", append=TRUE)
write.xlsx(thetaEAP_A_Aggre, "csemSSTheta_A.xlsx", sheetName="thetaEAP_A_Aggre", append=TRUE)



# Form B
png("TestData/CSSEM_Polynomial_SSTheta_B.png",  width = 799, height = 596)
CSSEM_Polynomial_SSTheta_B <- ggplot() +
  geom_point(cssemDat_MLE_B_Aggre, mapping = aes(x = roundedSS, y = cssemPoly,  shape = "MLE"), size = 2) +
  geom_point(cssemDat_EAP_B_Aggre, mapping = aes(x = roundedSS, y = cssemPoly,  shape = "EAP"), size = 2) +
  geom_point(thetaEAP_B_Aggre, mapping = aes(x = roundedSS, y = cssemPoly,  shape = "EAP Posterior"), size = 2) +
  geom_line(cssemDat_MLE_B_Aggre, mapping = aes(x = roundedSS, y = cssemPoly,  shape = "MLE"), size = 0.6) +
  geom_line(cssemDat_EAP_B_Aggre, mapping = aes(x = roundedSS, y = cssemPoly,  shape = "EAP"), size = 0.6) +
  geom_line(thetaEAP_B_Aggre, mapping = aes(x = roundedSS, y = cssemPoly,  shape = "EAP Posterior"), size = 0.6) +
  scale_x_continuous(name = "Rounded Scale Score", breaks  = seq(100, 130, 5)) +
  scale_y_continuous(name = "CSSEM Polynomial Method", breaks  = seq(0, 3, 0.5),
                     limits = c(0,3)) +
  theme_bw() +
  theme(legend.title=element_blank())

print(CSSEM_Polynomial_SSTheta_B)
dev.off()

write.xlsx(cssemDat_MLE_B_Aggre, "csemSSTheta_B.xlsx", sheetName="cssemDat_MLE_B_Aggre")
write.xlsx(cssemDat_EAP_B_Aggre, "csemSSTheta_B.xlsx", sheetName="cssemDat_EAP_B_Aggre", append=TRUE)
write.xlsx(thetaEAP_B_Aggre, "csemSSTheta_B.xlsx", sheetName="thetaEAP_B_Aggre", append=TRUE)






# Form A k =4 roundedSS to theta

k <- 7

modelK <- lm(roundedSS ~ poly(theta, k, raw=TRUE), convTable_B_Poly)

# extract regression coefficients
regCoef <- summary(modelK)$coefficients[, 1]
regCoef


prd <- data.frame(theta = seq(from = -5, to = 5, length.out = 150))
prd$predictedSS <- predict(modelK, newdata = prd, se.modelK = TRUE)

write.xlsx(prd, "Fitted Line_B.xlsx", sheetName="Fitted Line_B")






### aggregate many to one ------------

cssemDat <- cssemDat[,c(3,5:(5+k-2))]
cssemDatAggre <- as.data.frame(apply(cssemDat[,c(-1)], 2, function(x) aggregate(x, by=list(Category=cssemDat$roundedSS), FUN=function(x){sqrt(sum(x^2)/length(x))})))
cssemDatAggre <- cssemDatAggre[,c(1, seq(2, 2*(k-1), 2))]


### plot all ks

cssemDatLong <- reshape(cssemDatAggre,
                        direction = "long",
                        varying = list(names(cssemDatAggre)[2:k]),
                        v.names = "cssempoly",
                        idvar = c("cssemPolyk1.Category"),
                        timevar = "Kvalue",
                        times = 1:(k-1))



png("TestData/CSSEM_poly_A.png",  width = 799, height = 596)

CSSEM_poly_A_P <- ggplot(cssemDatLong, aes(x = cssemPolyk1.Category, y = cssempoly, color = factor(Kvalue))) +
  geom_point() +
  scale_x_continuous(name = "Rounded Scale Score", breaks  = seq(100, 130, 5)) +
  scale_y_continuous(name = "CSSEM Polynomial Method") +
  geom_line() +
  theme_bw() +
  labs(colour="K value")

print(CSSEM_poly_A_P)
dev.off()




png("TestData/CSSEM_poly_B.png",  width = 799, height = 596)

CSSEM_poly_B_P <- ggplot(cssemDatLong, aes(x = cssemPolyk1.Category, y = cssempoly, color = factor(Kvalue))) +
  geom_point() +
  scale_x_continuous(name = "Rounded Scale Score", breaks  = seq(100, 130, 5)) +
  scale_y_continuous(name = "CSSEM Polynomial Method") +
  geom_line() +
  theme_bw() +
  labs(colour="K value")

print(CSSEM_poly_B_P)
dev.off()



### range and confidence interval  with k from 1 to maximum accepted ----------1111111111111111111111111---------

cssemDatAggre_1 <- cssemDatAggre


# moments
cssemDatAggre_1$max <- apply(cssemDatAggre_1[,c(-1)], 1, max)
cssemDatAggre_1$min <- apply(cssemDatAggre_1[,c(-1)], 1, min)
cssemDatAggre_1$mean <- apply(cssemDatAggre_1[,c(-1)], 1, mean)
cssemDatAggre_1$median <- apply(cssemDatAggre_1[,c(-1)], 1, median)
cssemDatAggre_1$sd <- apply(cssemDatAggre_1[,c(-1)], 1, sd)

# 95% confidence interval
cssemDatAggre_1 <- within(cssemDatAggre_1, {
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

cssemDatAggre_3 <- cssemDatAggre


# moments
cssemDatAggre_3$max <- apply(cssemDatAggre_3[,c(-1:-3)], 1, max)
cssemDatAggre_3$min <- apply(cssemDatAggre_3[,c(-1:-3)], 1, min)
cssemDatAggre_3$mean <- apply(cssemDatAggre_3[,c(-1:-3)], 1, mean)
cssemDatAggre_3$median <- apply(cssemDatAggre_3[,c(-1:-3)], 1, median)
cssemDatAggre_3$sd <- apply(cssemDatAggre_3[,c(-1:-3)], 1, sd)

# 95% confidence interval
cssemDatAggre_3 <- within(cssemDatAggre_3, {
  lower = mean - 1.96 * sd
  upper = mean + 1.96 * sd
})


### plot confidence interval

# k = 3 to maximum accepted

ggplot(cssemDatAggre_3, aes(x = cssemPolyk1.Category, y = mean)) +
  geom_line(colour ="blue") +
  geom_point(colour="blue") +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill= "blue", alpha = 0.2) +
  scale_x_continuous(name = "Rounded Scale Score", breaks  = seq(100, 130, 5)) +
  scale_y_continuous(name = "CSSEM Polynomial Method") +
  theme_bw()





### polynomial method model fit --------------------------



# csem Lord
csemLordDat <- CSEMLord(numOfItem)

# merge with converstion table
cssemDat <- merge(csemLordDat, convTable_A_sub, by = "rawScore")

# change variable name
names(cssemDat)[names(cssemDat) == 'csemLord'] <- 'csem'


library(ggplot2)

modelK <- lm(roundedSS ~ poly(rawScore, 3, raw=TRUE), cssemDat)

prd <- data.frame(rawScore = seq(from = range(cssemDat$rawScore)[1], to = range(cssemDat$rawScore)[2], length.out = 100))

prd$predictedSS <- predict(modelK, newdata = prd, se.modelK = TRUE)

ggplot(prd, aes(x = rawScore, y = predictedSS)) +
  theme_bw() +
  geom_line() +
  geom_point(data = cssemDat, aes(x = rawScore, y = roundedSS))






# fitted line

library(ggplot2)

K = 5
for (k in 1:K){

  modelK <- lm(roundedSS ~ poly(rawScore, k, raw=TRUE), cssemDat)
  prd <- data.frame(rawScore = seq(from = range(cssemDat$rawScore)[1], to = range(cssemDat$rawScore)[2], length.out = 100))
  prd$predictedSS <- predict(modelK, newdata = prd, se.modelK = TRUE)

  p <- ggplot(prd, aes(x = rawScore, y = predictedSS)) +
    theme_bw() +
    geom_line() +
    geom_point(data = cssemDat, aes(x = rawScore, y = roundedSS)) +
    scale_y_continuous(name = "Scale Score", breaks  = seq(100, 130, 5)) +
    scale_x_continuous(name = "Raw Score") +
    ggtitle(paste("k = ", k, sep = ""))

  print(p)

}
