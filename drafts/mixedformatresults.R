# load library
library(EMReliability)

# read raw data
rawData_A <- read.table("TestData/FormA_31_3000.txt")
rawData_B <- read.table("TestData/FormB_31_3000.txt")

# read item parameters from txt file
itemPara_A <- read.table("TestData/SpanishLit_prm_A.txt")[,c(7,8)]
names(itemPara_A) <- c("b", "a")
itemPara_A[,"b"] <- -itemPara_A[,"b"]/itemPara_A[,"a"]
itemPara_A[,"a"] <- itemPara_A[,"a"]/1.702

itemPara_B <- read.table("TestData/SpanishLit_prm_B.txt")[,c(7,8)]
names(itemPara_B) <- c("b", "a")
itemPara_B[,"b"] <- -itemPara_B[,"b"]/itemPara_B[,"a"]
itemPara_B[,"a"] <- itemPara_B[,"a"]/1.702

# change extreme b value
itemPara_B[3,1] <- -5


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


# Stratified Alpha
CronbachAlpha(rawData_A, strat = c(13, 12, 6))
CronbachAlpha(rawData_B, strat = c(13, 12, 6))







