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

cssemPolynomial_A <- CSSEMPolynomial(numOfItem = 40, convTable_A_sub, 10)
cssemPolynomial_A
cssemPolynomial_B <- CSSEMPolynomial(numOfItem = 40, convTable_B_sub, 10)
cssemPolynomial_B

dat <- as.data.frame(cssemPolynomial_A$CSSEMPolynomial)
ggplot(dat, mapping = aes(x = rawScore, y = fxk3)) +
  geom_point()


dat2 <- as.data.frame(cssemPolynomial_A$CSSEMPolynomial)
ggplot(dat, mapping = aes(x = rawScore, y = fxk3)) +
  geom_point()

# CSSEM Polynomial MLE

cssemDat_MLE_A <- CSSEMIRTPoly(itemPara_A, convTable_A_Poly, 10, "MLE")
cssemDat_MLE_A
cssemDat_MLE_B <- CSSEMIRTPoly(itemPara_B, convTable_B_Poly, 10, "MLE")
cssemDat_MLE_B

rawFreq <- as.data.frame(table(rowSums(rawData_A)))
names(rawFreq) <- c("rawScore", "freq")
rawFreq$wt <- rawFreq$freq / sum(rawFreq$freq)



datN <- cssemDat_MLE_A$CSSEMPolyMLE
datN$wt <- c(0,0,0,rawFreq$wt)

SSVar <- sum(datN$wt * (datN$SSk4 - weighted.mean(datN$SSk4, datN$wt))^2)
errVar <- sum(datN$cssemPolyk4^2 * datN$wt)
1 - errVar/SSVar

# CSSEM Polynomial EAP

cssemDat_EAP_A <- CSSEMIRTPoly(itemPara_A, convTable_A_Poly, 10, "EAP")
cssemDat_EAP_A
cssemDat_EAP_B <- CSSEMIRTPoly(itemPara_B, convTable_B_Poly, 10, "EAP")
cssemDat_EAP_B



