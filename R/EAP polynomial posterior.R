# EAP polynomial posterior


#EAP
thetaSEEAP_A <- read.table("TestData/UShistory_X_EAP-sco.txt")[,c(3,4)]
names(thetaSEEAP_A) <- c("theta", "csemEAP")
itemParaCSEM <- thetaSEEAP_A

itemPara <- itemPara_A # test
convTable <- convTable_A
K <- 10
estType <- "EAP"



# 1 - mean(thetaSEEAP_B$V4^2)
thetaSEEAP_B <- read.table("TestData/UShistory_Y_EAP-sco.txt")[,c(3,4)]
names(thetaSEEAP_B) <- c("theta", "csemEAP")
itemParaCSEM <- thetaSEEAP_B

itemPara <- itemPara_B # test
convTable <- convTable_B
K <- 10
estType <- "EAP"





# # theta and weights
# theta <- NormalQuadraPoints(41)$nodes
#
# # CSEM
# itemParaCSEM <- as.data.frame(CSEMIRT(theta, itemPara, "EAP"))
# # weight
# itemParaCSEM$wt <- NormalQuadraPoints(41)$weights

rSquaredDat <- as.data.frame(matrix(nrow = K, ncol = 1))
relDat <- as.data.frame(matrix(nrow = K, ncol = 1))

# for loop to iterate different k
for (k in 1:K){

  k <- 4 # test

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
  itemParaCSEM$SS <- 0
  i <- k
  while(i > 0){
    itemParaCSEM$SS <- itemParaCSEM$SS +  regCoef[i+1] * itemParaCSEM$theta^(i)
    i <- i-1
  }

  itemParaCSEM$SS <- itemParaCSEM$SS + regCoef[i+1]
  itemParaCSEM$SS <- round(itemParaCSEM$SS)

  # calculate transformation coefficients fx: from 1 to K
  itemParaCSEM$fx <- 0
  i <- k
  while(i > 1){
    itemParaCSEM$fx <- itemParaCSEM$fx +  regCoef[i+1] * (i * itemParaCSEM$theta^(i-1))
    i <- i-1
  }

  itemParaCSEM$fx <- itemParaCSEM$fx + regCoef[i+1]

  # calculate cssem using polynomial method
  itemParaCSEM$cssemPoly <- itemParaCSEM$fx * itemParaCSEM$csemEAP

  # check negative values
  negCount <- sum(itemParaCSEM$fx < 0)
  if(negCount > 0){
    message(paste("Negative transformation coefficient exists when k = ", k,
                  ". The corresponding CSSEM will be corrected to 0 to calculate reliability.", sep = " "))

    itemParaCSEM$cssemPoly[itemParaCSEM$fx<0] <- 0
  }


  SSVar <- var(itemParaCSEM$SS)
  errVar <- mean(itemParaCSEM$cssemPoly^2)
  SSVar
  errVar

  1 - errVar/SSVar

  1 - errVar/(errVar + SSVar)


  SSVar <- sum(itemParaCSEM$wt * (itemParaCSEM$SS - weighted.mean(itemParaCSEM$SS, itemParaCSEM$wt))^2)
  errVar <- sum(itemParaCSEM$cssemPoly^2 * itemParaCSEM$wt)
  relDat[k,1] <- 1 - errVar/SSVar

}

relDat$kValue <- 1:K

# select reliability values not equal to 1, and larger than 0
relDat <- relDat[relDat$V1 > 0 & relDat$V1 != 1,]
relDat

# return results
return(list("kValue" = relDat$kValue, "RelEAPPoly" = relDat$V1))






# MLE polynomial posterior


# MLE
thetaSEMLE_A <- read.table("TestData/UShistory_X-sco.txt")[,c(4,5)]

thetaSEMLE_A <- thetaSEMLE_A[thetaSEMLE_A$V5 != 99.990000,]

# max(thetaSEMLE_A$V5)


names(thetaSEMLE_A) <- c("theta", "csemMLE")
itemParaCSEM <- thetaSEMLE_A

itemPara <- itemPara_A # test
convTable <- convTable_A
K <- 10
estType <- "MLE"


k <- 4 # test


thetaSEMLE_B <- read.table("TestData/UShistory_Y-sco.txt")[,c(4,5)]
thetaSEMLE_B <- thetaSEMLE_B[thetaSEMLE_B$V5 != 99.990000,]
names(thetaSEMLE_B) <- c("theta", "csemMLE")
itemParaCSEM <- thetaSEMLE_B

itemPara <- itemPara_B # test
convTable <- convTable_B
K <- 10
estType <- "MLE"



k <- 7 # test



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
itemParaCSEM$SS <- 0
i <- k
while(i > 0){
  itemParaCSEM$SS <- itemParaCSEM$SS +  regCoef[i+1] * itemParaCSEM$theta^(i)
  i <- i-1
}

itemParaCSEM$SS <- itemParaCSEM$SS + regCoef[i+1]
itemParaCSEM$SS <- round(itemParaCSEM$SS)

# calculate transformation coefficients fx: from 1 to K
itemParaCSEM$fx <- 0
i <- k
while(i > 1){
  itemParaCSEM$fx <- itemParaCSEM$fx +  regCoef[i+1] * (i * itemParaCSEM$theta^(i-1))
  i <- i-1
}

itemParaCSEM$fx <- itemParaCSEM$fx + regCoef[i+1]

# calculate cssem using polynomial method
itemParaCSEM$cssemPoly <- itemParaCSEM$fx * itemParaCSEM$csemMLE

# check negative values
negCount <- sum(itemParaCSEM$fx < 0)
if(negCount > 0){
  message(paste("Negative transformation coefficient exists when k = ", k,
                ". The corresponding CSSEM will be corrected to 0 to calculate reliability.", sep = " "))

  itemParaCSEM$cssemPoly[itemParaCSEM$fx<0] <- 0
}

SSVar <- var(itemParaCSEM$SS)
errVar <- mean(itemParaCSEM$cssemPoly^2)

SSVar
errVar





Ts <- 14.492
Es <- 1.675
EsInfo <- 1.682
S <- 12.197


r3 <- S/Ts
r3
r4 <- 1-Es/Ts
r4
r5 <- S/(S+Es)
r5

r8 <- 1-EsInfo/Ts
r8


varT <- 1

varThetahat <- 1.288

varE <- 0.150

r3 <- varT/varThetahat
r3
r4 <- 1-varE/varThetahat
r4
r5 <- varT/(varT+varE)
r5



library(ggplot2)





# SS form A

csemMIRTSS_A <- read.table("EM_SS_CSEM_m_A.OUT", header = T)
csemMIRTSS_A_raw <- csemMIRTSS_A[,c("Ex_Raw", "Raw_CSEM")]
names(csemMIRTSS_A_raw) <- c("theta", "csemEAP")
itemParaCSEM <- csemMIRTSS_A_raw

# itemPara <- itemPara_A # test
convTable <- convTable_A
names(convTable) <- c("rawScore", "roundedSS")
K <- 10


gra = T

# create data frame to store r square
rSquaredDat <- as.data.frame(matrix(nrow = K, ncol = 1))

# for loop to iterate different k
for (k in 1:K){

  # k <- 4
  # fit model with k
  modelK <- lm(roundedSS ~ poly(rawScore, k, raw=TRUE), convTable)

  # extract regression coefficients
  regCoef <- summary(modelK)$coefficients[, 1]

  # extract r square coefficient
  rSquaredDat[k, 1]<- summary(modelK)$r.squared

  # check whether regression coefficient of highest order is missing
  if(is.na(regCoef[k+1])){

    message(paste("The maximum k accepted is", k-1, ".", sep = " "))

    break

  }

  # graph
  if(gra == TRUE){

    prd <- data.frame(rawScore = seq(from = range(convTable$rawScore)[1], to = range(convTable$rawScore)[2], length.out = 100))
    prd$predictedSS <- predict(modelK, newdata = prd, se.modelK = TRUE)

    p <- ggplot(prd, aes(x = rawScore, y = predictedSS)) +
      theme_bw() +
      geom_line() +
      geom_point(data = convTable, aes(x = rawScore, y = roundedSS)) +
      scale_y_continuous(name = "Scale Score") +
      scale_x_continuous(name = "Raw Score") +
      ggtitle(paste("Polynomial Method Fitted Line with k = ", k, sep = "")) +
      theme(plot.title = element_text(hjust = 0.5))

    print(p)
  }


  # calculate SS: from 1 to K
  itemParaCSEM$SS <- 0
  i <- k
  while(i > 0){
    itemParaCSEM$SS <- itemParaCSEM$SS +  regCoef[i+1] * itemParaCSEM$theta^(i)
    i <- i-1
  }

  itemParaCSEM$SS <- itemParaCSEM$SS + regCoef[i+1]
  itemParaCSEM$roundedSS <- round(itemParaCSEM$SS)

  # calculate transformation coefficients fx: from 1 to K
  itemParaCSEM$fx <- 0
  i <- k
  while(i > 1){
    itemParaCSEM$fx <- itemParaCSEM$fx +  regCoef[i+1] * (i * itemParaCSEM$theta^(i-1))
    i <- i-1
  }

  itemParaCSEM$fx <- itemParaCSEM$fx + regCoef[i+1]

  # calculate cssem using polynomial method
  itemParaCSEM$cssemPoly <- itemParaCSEM$fx * itemParaCSEM$csemEAP

  # check negative values
  negCount <- sum(itemParaCSEM$fx < 0)
  if(negCount > 0){
    message(paste("Negative transformation coefficient exists when k = ", k,
                  ". The corresponding CSSEM will be corrected to 0 to calculate reliability.", sep = " "))

    itemParaCSEM$cssemPoly[itemParaCSEM$fx<0] <- 0
  }

  # rename variable with indicator k
  names(itemParaCSEM)[names(itemParaCSEM) == 'cssemPoly'] <- paste("cssemPolyk", k, sep = "")
  names(itemParaCSEM)[names(itemParaCSEM) == 'fx'] <- paste("fxk", k, sep = "")
  names(itemParaCSEM)[names(itemParaCSEM) == 'SS'] <- paste("SSk", k, sep = "")
  names(itemParaCSEM)[names(itemParaCSEM) == 'roundedSS'] <- paste("roundedSSk", k, sep = "")

}

write.xlsx(rSquaredDat, "rSquaredDat_SS_A.xlsx")
write.xlsx(itemParaCSEM, "itemParaCSEM_SS_A.xlsx")

itemParaCSEM_SS_A_Aggre <- aggregate(itemParaCSEM$cssemPolyk4, by=list(Category=itemParaCSEM$roundedSSk4), FUN=function(x){sqrt(sum(x^2)/length(x))})
names(itemParaCSEM_SS_A_Aggre) <- c("roundedSS", "cssempolynomialMIRT SS")
write.xlsx(itemParaCSEM_SS_A_Aggre, "itemParaCSEM_SS_A_Aggre.xlsx")







# SS form B

csemMIRTSS_B <- read.table("EM_SS_CSEM_m_B.OUT", header = T)
csemMIRTSS_B_raw <- csemMIRTSS_B[,c("Ex_Raw", "Raw_CSEM")]
names(csemMIRTSS_B_raw) <- c("theta", "csemEAP")
itemParaCSEM <- csemMIRTSS_B_raw

# itemPara <- itemPara_A # test
convTable <- convTable_B
names(convTable) <- c("rawScore", "roundedSS")
K <- 10


gra = T

# create data frame to store r square
rSquaredDat <- as.data.frame(matrix(nrow = K, ncol = 1))

# for loop to iterate different k
for (k in 1:K){

  # k <- 4
  # fit model with k
  modelK <- lm(roundedSS ~ poly(rawScore, k, raw=TRUE), convTable)

  # extract regression coefficients
  regCoef <- summary(modelK)$coefficients[, 1]

  # extract r square coefficient
  rSquaredDat[k, 1]<- summary(modelK)$r.squared

  # check whether regression coefficient of highest order is missing
  if(is.na(regCoef[k+1])){

    message(paste("The maximum k accepted is", k-1, ".", sep = " "))

    break

  }

  # graph
  if(gra == TRUE){

    prd <- data.frame(rawScore = seq(from = range(convTable$rawScore)[1], to = range(convTable$rawScore)[2], length.out = 100))
    prd$predictedSS <- predict(modelK, newdata = prd, se.modelK = TRUE)

    p <- ggplot(prd, aes(x = rawScore, y = predictedSS)) +
      theme_bw() +
      geom_line() +
      geom_point(data = convTable, aes(x = rawScore, y = roundedSS)) +
      scale_y_continuous(name = "Scale Score") +
      scale_x_continuous(name = "Raw Score") +
      ggtitle(paste("Polynomial Method Fitted Line with k = ", k, sep = "")) +
      theme(plot.title = element_text(hjust = 0.5))

    print(p)
  }


  # calculate SS: from 1 to K
  itemParaCSEM$SS <- 0
  i <- k
  while(i > 0){
    itemParaCSEM$SS <- itemParaCSEM$SS +  regCoef[i+1] * itemParaCSEM$theta^(i)
    i <- i-1
  }

  itemParaCSEM$SS <- itemParaCSEM$SS + regCoef[i+1]
  itemParaCSEM$roundedSS <- round(itemParaCSEM$SS)

  # calculate transformation coefficients fx: from 1 to K
  itemParaCSEM$fx <- 0
  i <- k
  while(i > 1){
    itemParaCSEM$fx <- itemParaCSEM$fx +  regCoef[i+1] * (i * itemParaCSEM$theta^(i-1))
    i <- i-1
  }

  itemParaCSEM$fx <- itemParaCSEM$fx + regCoef[i+1]

  # calculate cssem using polynomial method
  itemParaCSEM$cssemPoly <- itemParaCSEM$fx * itemParaCSEM$csemEAP

  # check negative values
  negCount <- sum(itemParaCSEM$fx < 0)
  if(negCount > 0){
    message(paste("Negative transformation coefficient exists when k = ", k,
                  ". The corresponding CSSEM will be corrected to 0 to calculate reliability.", sep = " "))

    itemParaCSEM$cssemPoly[itemParaCSEM$fx<0] <- 0
  }

  # rename variable with indicator k
  names(itemParaCSEM)[names(itemParaCSEM) == 'cssemPoly'] <- paste("cssemPolyk", k, sep = "")
  names(itemParaCSEM)[names(itemParaCSEM) == 'fx'] <- paste("fxk", k, sep = "")
  names(itemParaCSEM)[names(itemParaCSEM) == 'SS'] <- paste("SSk", k, sep = "")
  names(itemParaCSEM)[names(itemParaCSEM) == 'roundedSS'] <- paste("roundedSSk", k, sep = "")

}

write.xlsx(rSquaredDat, "rSquaredDat_SS_B.xlsx")
write.xlsx(itemParaCSEM, "itemParaCSEM_SS_B.xlsx")

itemParaCSEM_SS_B_Aggre <- aggregate(itemParaCSEM$cssemPolyk4, by=list(Category=itemParaCSEM$roundedSSk4), FUN=function(x){sqrt(sum(x^2)/length(x))})
names(itemParaCSEM_SS_B_Aggre) <- c("roundedSS", "cssempolynomialMIRT SS")
write.xlsx(itemParaCSEM_SS_B_Aggre, "itemParaCSEM_SS_B_Aggre.xlsx")









# BF form A

csemMIRTSS_A <- read.table("EMBF_CSEM_M_A.OUT", header = T)
csemMIRTSS_A_raw <- csemMIRTSS_A[,c("Ex_Raw", "Raw_CSEM")]
names(csemMIRTSS_A_raw) <- c("theta", "csemEAP")
itemParaCSEM <- csemMIRTSS_A_raw

# itemPara <- itemPara_A # test
convTable <- convTable_A
names(convTable) <- c("rawScore", "roundedSS")
K <- 10


gra = T

# create data frame to store r square
rSquaredDat <- as.data.frame(matrix(nrow = K, ncol = 1))

# for loop to iterate different k
for (k in 1:K){

  # k <- 4
  # fit model with k
  modelK <- lm(roundedSS ~ poly(rawScore, k, raw=TRUE), convTable)

  # extract regression coefficients
  regCoef <- summary(modelK)$coefficients[, 1]

  # extract r square coefficient
  rSquaredDat[k, 1]<- summary(modelK)$r.squared

  # check whether regression coefficient of highest order is missing
  if(is.na(regCoef[k+1])){

    message(paste("The maximum k accepted is", k-1, ".", sep = " "))

    break

  }

  # graph
  if(gra == TRUE){

    prd <- data.frame(rawScore = seq(from = range(convTable$rawScore)[1], to = range(convTable$rawScore)[2], length.out = 100))
    prd$predictedSS <- predict(modelK, newdata = prd, se.modelK = TRUE)

    p <- ggplot(prd, aes(x = rawScore, y = predictedSS)) +
      theme_bw() +
      geom_line() +
      geom_point(data = convTable, aes(x = rawScore, y = roundedSS)) +
      scale_y_continuous(name = "Scale Score") +
      scale_x_continuous(name = "Raw Score") +
      ggtitle(paste("Polynomial Method Fitted Line with k = ", k, sep = "")) +
      theme(plot.title = element_text(hjust = 0.5))

    print(p)
  }


  # calculate SS: from 1 to K
  itemParaCSEM$SS <- 0
  i <- k
  while(i > 0){
    itemParaCSEM$SS <- itemParaCSEM$SS +  regCoef[i+1] * itemParaCSEM$theta^(i)
    i <- i-1
  }

  itemParaCSEM$SS <- itemParaCSEM$SS + regCoef[i+1]
  itemParaCSEM$roundedSS <- round(itemParaCSEM$SS)

  # calculate transformation coefficients fx: from 1 to K
  itemParaCSEM$fx <- 0
  i <- k
  while(i > 1){
    itemParaCSEM$fx <- itemParaCSEM$fx +  regCoef[i+1] * (i * itemParaCSEM$theta^(i-1))
    i <- i-1
  }

  itemParaCSEM$fx <- itemParaCSEM$fx + regCoef[i+1]

  # calculate cssem using polynomial method
  itemParaCSEM$cssemPoly <- itemParaCSEM$fx * itemParaCSEM$csemEAP

  # check negative values
  negCount <- sum(itemParaCSEM$fx < 0)
  if(negCount > 0){
    message(paste("Negative transformation coefficient exists when k = ", k,
                  ". The corresponding CSSEM will be corrected to 0 to calculate reliability.", sep = " "))

    itemParaCSEM$cssemPoly[itemParaCSEM$fx<0] <- 0
  }

  # rename variable with indicator k
  names(itemParaCSEM)[names(itemParaCSEM) == 'cssemPoly'] <- paste("cssemPolyk", k, sep = "")
  names(itemParaCSEM)[names(itemParaCSEM) == 'fx'] <- paste("fxk", k, sep = "")
  names(itemParaCSEM)[names(itemParaCSEM) == 'SS'] <- paste("SSk", k, sep = "")
  names(itemParaCSEM)[names(itemParaCSEM) == 'roundedSS'] <- paste("roundedSSk", k, sep = "")

}

write.xlsx(rSquaredDat, "rSquaredDat_BF_A.xlsx")
write.xlsx(itemParaCSEM, "itemParaCSEM_BF_A.xlsx")

itemParaCSEM_SS_A_Aggre <- aggregate(itemParaCSEM$cssemPolyk4, by=list(Category=itemParaCSEM$roundedSSk4), FUN=function(x){sqrt(sum(x^2)/length(x))})
names(itemParaCSEM_SS_A_Aggre) <- c("roundedSS", "cssempolynomialMIRT BF")
write.xlsx(itemParaCSEM_SS_A_Aggre, "itemParaCSEM_BF_A_Aggre.xlsx")







# BF form B

csemMIRTSS_B <- read.table("EMBF_CSEM_M_B.OUT", header = T)
csemMIRTSS_B_raw <- csemMIRTSS_B[,c("Ex_Raw", "Raw_CSEM")]
names(csemMIRTSS_B_raw) <- c("theta", "csemEAP")
itemParaCSEM <- csemMIRTSS_B_raw

# itemPara <- itemPara_A # test
convTable <- convTable_B
names(convTable) <- c("rawScore", "roundedSS")
K <- 10


gra = T

# create data frame to store r square
rSquaredDat <- as.data.frame(matrix(nrow = K, ncol = 1))

# for loop to iterate different k
for (k in 1:K){

  # k <- 4
  # fit model with k
  modelK <- lm(roundedSS ~ poly(rawScore, k, raw=TRUE), convTable)

  # extract regression coefficients
  regCoef <- summary(modelK)$coefficients[, 1]

  # extract r square coefficient
  rSquaredDat[k, 1]<- summary(modelK)$r.squared

  # check whether regression coefficient of highest order is missing
  if(is.na(regCoef[k+1])){

    message(paste("The maximum k accepted is", k-1, ".", sep = " "))

    break

  }

  # graph
  if(gra == TRUE){

    prd <- data.frame(rawScore = seq(from = range(convTable$rawScore)[1], to = range(convTable$rawScore)[2], length.out = 100))
    prd$predictedSS <- predict(modelK, newdata = prd, se.modelK = TRUE)

    p <- ggplot(prd, aes(x = rawScore, y = predictedSS)) +
      theme_bw() +
      geom_line() +
      geom_point(data = convTable, aes(x = rawScore, y = roundedSS)) +
      scale_y_continuous(name = "Scale Score") +
      scale_x_continuous(name = "Raw Score") +
      ggtitle(paste("Polynomial Method Fitted Line with k = ", k, sep = "")) +
      theme(plot.title = element_text(hjust = 0.5))

    print(p)
  }


  # calculate SS: from 1 to K
  itemParaCSEM$SS <- 0
  i <- k
  while(i > 0){
    itemParaCSEM$SS <- itemParaCSEM$SS +  regCoef[i+1] * itemParaCSEM$theta^(i)
    i <- i-1
  }

  itemParaCSEM$SS <- itemParaCSEM$SS + regCoef[i+1]
  itemParaCSEM$roundedSS <- round(itemParaCSEM$SS)

  # calculate transformation coefficients fx: from 1 to K
  itemParaCSEM$fx <- 0
  i <- k
  while(i > 1){
    itemParaCSEM$fx <- itemParaCSEM$fx +  regCoef[i+1] * (i * itemParaCSEM$theta^(i-1))
    i <- i-1
  }

  itemParaCSEM$fx <- itemParaCSEM$fx + regCoef[i+1]

  # calculate cssem using polynomial method
  itemParaCSEM$cssemPoly <- itemParaCSEM$fx * itemParaCSEM$csemEAP

  # check negative values
  negCount <- sum(itemParaCSEM$fx < 0)
  if(negCount > 0){
    message(paste("Negative transformation coefficient exists when k = ", k,
                  ". The corresponding CSSEM will be corrected to 0 to calculate reliability.", sep = " "))

    itemParaCSEM$cssemPoly[itemParaCSEM$fx<0] <- 0
  }

  # rename variable with indicator k
  names(itemParaCSEM)[names(itemParaCSEM) == 'cssemPoly'] <- paste("cssemPolyk", k, sep = "")
  names(itemParaCSEM)[names(itemParaCSEM) == 'fx'] <- paste("fxk", k, sep = "")
  names(itemParaCSEM)[names(itemParaCSEM) == 'SS'] <- paste("SSk", k, sep = "")
  names(itemParaCSEM)[names(itemParaCSEM) == 'roundedSS'] <- paste("roundedSSk", k, sep = "")

}

write.xlsx(rSquaredDat, "rSquaredDat_BF_B.xlsx")
write.xlsx(itemParaCSEM, "itemParaCSEM_BF_B.xlsx")

itemParaCSEM_SS_B_Aggre <- aggregate(itemParaCSEM$cssemPolyk4, by=list(Category=itemParaCSEM$roundedSSk4), FUN=function(x){sqrt(sum(x^2)/length(x))})
names(itemParaCSEM_SS_B_Aggre) <- c("roundedSS", "cssempolynomialMIRT BF")
write.xlsx(itemParaCSEM_SS_B_Aggre, "itemParaCSEM_BF_B_Aggre.xlsx")


