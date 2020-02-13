
### Results for project ------------------


## IRT Test Reliability

library(classify)
source("R/TestRelIRT.R") # itemPara
source("R/NormalQuadraPoints.R") # n set as 41
source("R/MarginalRelMLE.R") # itemPara
source("R/MarginalRelEAP.R") # itemPara
source("R/KolenRelIRT.R") # itemPara, convTable

# read item parameters from txt file
itemPara_A <- read.table("TestData/ItemParaFormX.txt")
itemPara_B <- read.table("TestData/ItemParaFormY.txt")

# read conversion tables
convTable_A <- read.csv("TestData/ConversionTableFormX.csv")
convTable_A$roundedSS <- round(convTable_A$unroundedSS)

convTable_B <- read.csv("TestData/ConversionTableFormY.csv")
convTable_B$roundedSS <- round(convTable_B$unroundedSS)


# test reliability IRT
TestRelIRT_A <- TestRelIRT(itemPara_A)
TestRelIRT_A

TestRelIRT_B <- TestRelIRT(itemPara_B)
TestRelIRT_B

# marginal reliability MLE
MarginalRelMLE_A <- MarginalRelMLE(itemPara_A)
MarginalRelMLE_A

MarginalRelMLE_B <- MarginalRelMLE(itemPara_B)
MarginalRelMLE_B

# marginal reliability EAP
MarginalRelEAP_A <- MarginalRelEAP(itemPara_A)
MarginalRelEAP_A

MarginalRelEAP_B <- MarginalRelEAP(itemPara_B)
MarginalRelEAP_B

# Kolen's method

KolenRelIRT_A <- KolenRelIRT(itemPara_A, convTable_A)
KolenRelIRT_A

KolenRelIRT_B <- KolenRelIRT(itemPara_B, convTable_B)
KolenRelIRT_B






