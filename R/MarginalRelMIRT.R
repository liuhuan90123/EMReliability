# Marginal Reliability

#### simple structure MLE&EAP  P method ----------

## MLE ---------------------------------------

# number of factors
numOfFactors <- 3

# Form A
# scoSS_MLE <- read.table("TestData/SpanishLit_sco_A_SS_MLE.txt")[,c(4:9, 10, 12, 15)]
# cor <- c(0.91, # 1&2
#          0.70, # 1&3
#          0.49) # 2&3

# Form B
scoSS_MLE <- read.table("TestData/SpanishLit_sco_B_SS_MLE.txt")[,c(4:9, 10, 12, 15)]
cor <- c(0.97, # 1&2
         0.56, # 1&3
         0.48) # 2&3

# Form Test
# scoSS_MLE <- read.table("TestData/SpanishLit-sco-test.txt")[,c(4:9, 10, 12, 15)]
# cor <- c(0.91, # 1&2
#          0.71, # 1&3
#          0.51) # 2&3


# change variable name
names(scoSS_MLE) <- c("theta1", "theta2", "theta3", "se1", "se2", "se3", "var1", "var2", "var3")

# delete observations with missing values, 99.99 for flexMIRT output
scoSS_MLE[scoSS_MLE == 99.99] <- NA
scoSS_MLE <- na.omit(scoSS_MLE)

# composite error variance
scoSS_MLE<- transform( scoSS_MLE,
                   varC = var1 + var2 + var3 + 2 *cor[1] * se1 * se2 + 2 *cor[2] * se1 * se3 + 2 *cor[3] * se2 * se3
)

# average of error variance
ErrorVarAvg <- mean(scoSS_MLE$varC)

# marginal reliability approach
MarginalRelSSMIRT_MLE_P <- ErrorVarAvg /(ErrorVarAvg + numOfFactors) # var(e)/(var(e) + var(theta))

# coefficients
MarginalRelSSMIRT_MLE_P




####### EAP ---------------------------------------

# Form A
scoSS_EAP <- read.table("TestData/SpanishLit_sco_A_SS_EAP.txt")[,c(3:14)]
cor <- c(0.91, # 1&2
         0.70, # 1&3
         0.49) # 2&3

# Form B
# scoSS_EAP <- read.table("TestData/SpanishLit_sco_B_SS_EAP.txt")[,c(3:14)]
# cor <- c(0.97, # 1&2
#          0.56, # 1&3
#          0.48) # 2&3

# change variable name
names(scoSS_EAP) <- c("theta1", "theta2", "theta3", "se1", "se2", "se3", "var11", "var21", "var22", "var31","var32","var33")


# composite error variance
scoSS_EAP<- transform( scoSS_EAP,
                   varC = var11 + var22  + var33 +  2*(var21 + var31 + var32)
)

# average of error variance
ErrorVarAvg <- mean(scoSS_EAP$varC)

# marginal reliability approach
MarginalRelSSMIRT_EAP_P  <- 1 - ErrorVarAvg /(2*(sum(cor)) + numOfFactors) # var(e)/(var(e) + var(theta))

# coefficients
MarginalRelSSMIRT_EAP_P








########### bi factor full D method ---------------------------------------

# item parameters

# read item parameters from txt file
itemPara_BF <- read.table("TestData/SpanishLit_prm_B_BF.txt")[,c(7:11)]


names(itemPara_BF) <- c("b", "ag","a1","a2", "a3") # a1 is primary
itemPara_BF$a <- c(itemPara_BF$a1[1:13], itemPara_BF$a2[14:25], itemPara_BF$a3[26:31])


itemPara_BF[,"b"] <-  -itemPara_BF_G[,"b"]/(itemPara_BF_G[,"ag"] + itemPara_BF_G[,"a"])#######?????????????????



itemPara_BF[,"a"] <- itemPara_BF[,"a"]/1.702

itemPara_BF$a1[1:13] <- itemPara_BF$a[1:13]
itemPara_BF$a2[14:25] <- itemPara_BF$a[14:25]
itemPara_BF$a3[26:31] <- itemPara_BF$a[26:31]


itemPara_BF$ag <- itemPara_BF$ag/1.702































