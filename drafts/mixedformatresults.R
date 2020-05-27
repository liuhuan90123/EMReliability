# load library

# remove.packages("EMReliability")
# install_github("liuhuan91023/EMReliability")
# library(devtools)
# install_github("liuhuan90123/EMReliability")

library(EMReliability)

# read raw data
rawData_A <- read.table("TestData/FormA_31_3000.txt")
rawData_B <- read.table("TestData/FormB_31_3000.txt")

# read item parameters from txt file
itemPara_A_UIRT <- read.table("TestData/SpanishLit_prm_A_UIRT.txt")[,c(7,8)]
names(itemPara_A_UIRT) <- c("b", "a")
itemPara_A_UIRT[,"b"] <- -itemPara_A_UIRT[,"b"]/itemPara_A_UIRT[,"a"]
itemPara_A_UIRT[,"a"] <- itemPara_A_UIRT[,"a"]/1.702

itemPara_B_UIRT <- read.table("TestData/SpanishLit_prm_B_UIRT.txt")[,c(7,8)]
names(itemPara_B_UIRT) <- c("b", "a")
itemPara_B_UIRT[,"b"] <- -itemPara_B_UIRT[,"b"]/itemPara_B_UIRT[,"a"]
itemPara_B_UIRT[,"a"] <- itemPara_B_UIRT[,"a"]/1.702

# change extreme b value
itemPara_B_UIRT[3,1] <- -5


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

# Stratified Alpha
StratCronbachAlpha_A <- StratCronbachAlpha(rawData_A, strat = c(13, 12, 6))
StratCronbachAlpha_A
StratCronbachAlpha_B <- StratCronbachAlpha(rawData_B, strat = c(13, 12, 6))
StratCronbachAlpha_B

# Stratified Feldt CSEM

rawData_A <- read.table("TestData/FormA_31_3000.txt")
strat_A <-  c(13, 12, 6)
csemStratFeldt_A <- CSEMStratFeldt(rawData_A, strat_A)
plot(csemStratFeldt_A$rawScore, csemStratFeldt_A$csemStratFeldt)

rawData_B <- read.table("TestData/FormA_31_3000.txt")
strat_B <-  c(13, 12, 6)
csemStratFeldt_B <- CSEMStratFeldt(rawData_B, strat_B)
plot(csemStratFeldt_B$rawScore, csemStratFeldt_B$csemStratFeldt)



# UIRT --------------------

# test reliability IRT
TestRelIRT_A <- TestRelIRT(itemPara_A_UIRT)
TestRelIRT_A

TestRelIRT_B <- TestRelIRT(itemPara_B_UIRT)
TestRelIRT_B

# marginal reliability MLE
MarginalRelMLE_A <- MarginalRelIRT(itemPara_A_UIRT, "MLE")
MarginalRelMLE_A

MarginalRelMLE_B <- MarginalRelIRT(itemPara_B_UIRT, "MLE")
MarginalRelMLE_B

# marginal reliability EAP
MarginalRelEAP_A <- MarginalRelIRT(itemPara_A_UIRT, "EAP")
MarginalRelEAP_A

MarginalRelEAP_B <- MarginalRelIRT(itemPara_B_UIRT, "EAP")
MarginalRelEAP_B


# Kolen's method

# read conversion tables
convTable_A <- read.csv("TestData/conversion_table_Form A.csv")
convTable_A <- convTable_A[1:32, c("RawScore", "roundedSS")]

convTable_B <- read.csv("TestData/conversion_table_Form B.csv")
convTable_B <- convTable_B[1:32, c("RawScore", "roundedSS")]


KolenRelIRT_A_UIRT <- KolenRelIRT(itemPara_A_UIRT, convTable_A)
KolenRelIRT_A_UIRT

KolenRelIRT_B_UIRT <- KolenRelIRT(itemPara_B_UIRT, convTable_B)
KolenRelIRT_B_UIRT





# BI-Factor General

# read item parameters
# form A
# itemPara_BF_G <- read.table("TestData/SpanishLit_prm_A_BF.txt")[,c(7:11)]

# form B
itemPara_BF_G <- read.table("TestData/SpanishLit_prm_B_BF.txt")[,c(7:11)]

names(itemPara_BF_G) <- c("b", "ag","a1","a2", "a3") # a1 is primary
itemPara_BF_G$a <- c(itemPara_BF_G$a1[1:13], itemPara_BF_G$a2[14:25], itemPara_BF_G$a3[26:31])

itemPara_BF_G[,"b"] <- -itemPara_BF_G[,"b"]/(itemPara_BF_G[,"ag"] + itemPara_BF_G[,"a"]) #######?????????????????
itemPara_BF_G[,"a"] <- itemPara_BF_G[,"a"]/1.702

itemPara_BF_G$a1[1:13] <- itemPara_BF_G$a[1:13]
itemPara_BF_G$a2[14:25] <- itemPara_BF_G$a[14:25]
itemPara_BF_G$a3[26:31] <- itemPara_BF_G$a[26:31]
itemPara_BF_G$ag <- itemPara_BF_G$ag/1.702

itemPara_BF_G <- itemPara_BF_G[,c("b", "ag")]
names(itemPara_BF_G) <- c("b", "a")


TestRelIRT(itemPara_BF_G)
MarginalRelIRT(itemPara_BF_G, estType = "MLE")
MarginalRelIRT(itemPara_BF_G, estType = "EAP")

KolenRelIRT(itemPara_BF_G, convTable_B)

# Simple Structure 3 factors --------------------------------

# test reliability

# load packages
library(LaplacesDemon)
# library(mvtnorm)
# library(pbivnorm)

# item parameters and correlation matrix
# Form A
itemPara_SS_A <- read.table("TestData/SpanishLit_prm_A_SS.txt")[,c(7:10)]
cormat_A <- matrix(c(1, 0.9067069, 0.6994119,
                   0.9067069, 1, 0.4891160,
                   0.6994119,0.4891160,1), nrow = 3)
strat <- c(13, 12, 6)
# item parameter transformation
names(itemPara_SS_A) <- c("b", "a1","a2","a3")
itemPara_SS_A$a <- c(itemPara_SS_A$a1[1:13], itemPara_SS_A$a2[14:25], itemPara_SS_A$a3[26:31])
itemPara_SS_A[,"b"] <- -itemPara_SS_A[,"b"]/itemPara_SS_A[,"a"]
itemPara_SS_A[,"a"] <- itemPara_SS_A[,"a"]/1.702
itemPara_SS_A$a1[1:13] <- itemPara_SS_A$a[1:13]
itemPara_SS_A$a2[14:25] <- itemPara_SS_A$a[14:25]
itemPara_SS_A$a3[26:31] <- itemPara_SS_A$a[26:31]

# form B
itemPara_SS_B <- read.table("TestData/SpanishLit_prm_B_SS.txt")[,c(7:10)]
cormat_B <- matrix(c(1, 0.97, 0.56,
                0.97, 1, 0.48,
                0.56,0.48,1), nrow = 3)
strat <- c(13, 12, 6)
# item parameter transformation
names(itemPara_SS_B) <- c("b", "a1","a2","a3")
itemPara_SS_B$a <- c(itemPara_SS_B$a1[1:13], itemPara_SS_B$a2[14:25], itemPara_SS_B$a3[26:31])
itemPara_SS_B[,"b"] <- -itemPara_SS_B[,"b"]/itemPara_SS_B[,"a"]
itemPara_SS_B[,"a"] <- itemPara_SS_B[,"a"]/1.702
itemPara_SS_B$a1[1:13] <- itemPara_SS_B$a[1:13]
itemPara_SS_B$a2[14:25] <- itemPara_SS_B$a[14:25]
itemPara_SS_B$a3[26:31] <- itemPara_SS_B$a[26:31]

TestRelSSMIRT(itemPara_SS_A, strat, cormat_A)
TestRelSSMIRT(itemPara_SS_B, strat, cormat_B)



# marginal reliability MLE _ D method

# read correlations
corvec_A <- c(0.91, # 1&2
              0.70, # 1&3
              0.49) # 2&3


# read correlations
corvec_B <- c(0.97, # 1&2
              0.56, # 1&3
              0.48) # 2&3

MarginalRelSSMIRT_D_MLE(itemPara_SS_A, corvec_A, strat)
MarginalRelSSMIRT_D_MLE(itemPara_SS_B, corvec_B, strat)

# marginal reliability EAP





