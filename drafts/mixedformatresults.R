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
StratCronbachAlpha_A <- CronbachAlpha(rawData_A, strat = c(13, 12, 6))
StratCronbachAlpha_A
StratCronbachAlpha_B <- CronbachAlpha(rawData_B, strat = c(13, 12, 6))
StratCronbachAlpha_B

# Stratified Feldt CSEM

rawData_A <- read.table("TestData/FormA_31_3000.txt")
strat_A <-  c(13, 12, 6)
CSEMStratFeldt(rawData_A, strat_A)


rawData_B <- read.table("TestData/FormA_31_3000.txt")
strat_B <-  c(13, 12, 6)
CSEMStratFeldt(rawData_B, strat_B)




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


# Simple Structure 3 factors --------------------------------

# test reliability

# marginal reliability MLE

# marginal reliability EAP





