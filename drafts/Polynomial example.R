### Polynomial method

options(max.print=100000)

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

# conversion table sub
convTable_A_sub <- convTable_A[,c("rawScore", "roundedSS")]
convTable_B_sub <- convTable_B[,c("rawScore", "roundedSS")]

# conversion table poly
convTable_A_Poly <- convTable_A[,c("theta", "roundedSS")]
convTable_B_Poly <- convTable_B[,c("theta", "roundedSS")]



# CSSEM Polynomial Lord

cssemPolynomial_A <- CSSEMPolynomial(numOfItem = 40, convTable_A_sub, 20)
cssemPolynomial_A
cssemPolynomial_B <- CSSEMPolynomial(numOfItem = 40, convTable_B_sub, 20)
cssemPolynomial_B



# CSSEM Polynomial MLE

cssemDat_MLE_A <- CSSEMIRTPoly(itemPara_A, convTable_A_Poly, 20, "MLE")
cssemDat_MLE_A
cssemDat_MLE_B <- CSSEMIRTPoly(itemPara_B, convTable_B_Poly, 20, "MLE")
cssemDat_MLE_B

# CSSEM Polynomial MLE

cssemDat_MLE_A <- CSSEMIRTPoly(itemPara_A, convTable_A_Poly, 20, "EAP")
cssemDat_MLE_A
cssemDat_MLE_B <- CSSEMIRTPoly(itemPara_B, convTable_B_Poly, 20, "EAP")
cssemDat_MLE_B



