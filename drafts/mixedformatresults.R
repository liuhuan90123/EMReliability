# load library

# remove.packages("EMReliability")
# install_github("liuhuan91023/EMReliability")
# library(devtools)
# install_github("liuhuan90123/EMReliability")

library(EMReliability)

# read raw data
rawData_A <- read.table("TestData/FormA_31_3000.txt")
rawData_B <- read.table("TestData/FormB_31_3000.txt")

datA1 <- rawData_A[,c(1:13)]
datA2 <- rawData_A[,c(14:25)]
datA3 <- rawData_A[,c(26:31)]

datB1 <- rawData_B[,c(1:13)]
datB2 <- rawData_B[,c(14:25)]
datB3 <- rawData_B[,c(26:31)]


write.table(datA1, "datA1.dat", row.names = F, col.names = F)
write.table(datA2, "datA2.dat", row.names = F, col.names = F)
write.table(datA3, "datA3.dat", row.names = F, col.names = F)

write.table(datB1, "datB1.dat", row.names = F, col.names = F)
write.table(datB2, "datB2.dat", row.names = F, col.names = F)
write.table(datB3, "datB3.dat", row.names = F, col.names = F)





# read item parameters from txt file
itemPara_A_UIRT <- read.table("TestData/SpanishLit_prm_A_UIRT.txt")[,c(7,8)]
names(itemPara_A_UIRT) <- c("b", "a")
itemPara_A_UIRT[,"b"] <- -itemPara_A_UIRT[,"b"]/itemPara_A_UIRT[,"a"]
itemPara_A_UIRT[,"a"] <- itemPara_A_UIRT[,"a"]/1.702

itemPara_B_UIRT <- read.table("TestData/SpanishLit_prm_B_UIRT.txt")[,c(7,8)]
names(itemPara_B_UIRT) <- c("b", "a")
itemPara_B_UIRT[,"b"] <- -itemPara_B_UIRT[,"b"]/itemPara_B_UIRT[,"a"]
itemPara_B_UIRT[,"a"] <- itemPara_B_UIRT[,"a"]/1.702


itemPara_A_UIRT$c <- 0
itemPara_A_UIRT <- itemPara_A_UIRT[,c("a", "b", "c")]

itemPara_B_UIRT$c <- 0
itemPara_B_UIRT <- itemPara_B_UIRT[,c("a", "b", "c")]

write.table(itemPara_A_UIRT, "itemPara_A_UIRT.txt", row.names = F, col.names = F)
write.table(itemPara_B_UIRT, "itemPara_B_UIRT.txt", row.names = F, col.names = F)

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
write.xlsx(csemStratFeldt_A, "csemStratFeldt_A.xlsx", sheetName="csemStratFeldt_A")

csemStratFeldt_A_Aggre <- aggregate(csemStratFeldt_A$csemStratFeldt, by=list(Category=csemStratFeldt_A$rawScore), FUN=function(x){sqrt(sum(x^2)/length(x))})
csemStratFeldt_A_Aggre

rawData_B <- read.table("TestData/FormA_31_3000.txt")
strat_B <-  c(13, 12, 6)
csemStratFeldt_B <- CSEMStratFeldt(rawData_B, strat_B)
plot(csemStratFeldt_B$rawScore, csemStratFeldt_B$csemStratFeldt)
write.xlsx(csemStratFeldt_B, "csemStratFeldt_B.xlsx", sheetName="csemStratFeldt_B")

csemStratFeldt_B_Aggre <- aggregate(csemStratFeldt_B$csemStratFeldt, by=list(Category=csemStratFeldt_B$rawScore), FUN=function(x){sqrt(sum(x^2)/length(x))})
csemStratFeldt_B_Aggre




names(csemStratFeldt_A_Aggre) <- names(csemStratFeldt_B_Aggre) <- c("raw Score", "csemStratFeldt")


write.xlsx(csemStratFeldt_A_Aggre, "csemStratFeldt_A_Aggre.xlsx")
write.xlsx(csemStratFeldt_B_Aggre, "csemStratFeldt_B_Aggre.xlsx")







### MIRT P method

## SS form A
csemMIRTSS_A <- read.table("EM_SS_CSEM_m_A.OUT", header = T)
csemMIRTSS_A_raw <- csemMIRTSS_A[,c("Ex_Raw", "Raw_CSEM")]

write.xlsx(csemMIRTSS_A_raw, "csemMIRTSS_A_raw.xlsx", sheetName="csemMIRTSS_A_raw")

csemMIRTSS_A_raw$Ex_Raw <- round(csemMIRTSS_A_raw$Ex_Raw)
csemMIRTSS_A_raw_Aggre <- aggregate(csemMIRTSS_A_raw$Raw_CSEM, by=list(Category=csemMIRTSS_A_raw$Ex_Raw), FUN=function(x){sqrt(sum(x^2)/length(x))})
csemMIRTSS_A_raw_Aggre

names(csemMIRTSS_A_raw_Aggre) <- c("raw Score", "csemMIRTSS_A_raw_Aggre")

write.xlsx(csemMIRTSS_A_raw_Aggre, "csemMIRTSS_A_raw_Aggre.xlsx", sheetName="csemMIRTSS_A_raw_Aggre")


## SS form B
csemMIRTSS_B <- read.table("EM_SS_CSEM_m_B.OUT", header = T)
csemMIRTSS_B_raw <- csemMIRTSS_B[,c("Ex_Raw", "Raw_CSEM")]

write.xlsx(csemMIRTSS_B_raw, "csemMIRTSS_B_raw.xlsx", sheetName="csemMIRTSS_B_raw")

csemMIRTSS_B_raw$Ex_Raw <- round(csemMIRTSS_B_raw$Ex_Raw)
csemMIRTSS_B_raw_Aggre <- aggregate(csemMIRTSS_B_raw$Raw_CSEM, by=list(Category=csemMIRTSS_B_raw$Ex_Raw), FUN=function(x){sqrt(sum(x^2)/length(x))})
csemMIRTSS_B_raw_Aggre

names(csemMIRTSS_B_raw_Aggre) <- c("raw Score", "csemMIRTSS_B_raw_Aggre")

write.xlsx(csemMIRTSS_B_raw_Aggre, "csemMIRTSS_B_raw_Aggre.xlsx", sheetName="csemMIRTSS_B_raw_Aggre")




## BF form A
csemMIRTBF_A <- read.table("EMBF_CSEM_M_A.OUT", header = T)
csemMIRTBF_A_raw <- csemMIRTBF_A[,c("Ex_Raw", "Raw_CSEM")]

write.xlsx(csemMIRTBF_A_raw, "csemMIRTBF_A_raw.xlsx", sheetName="csemMIRTBF_A_raw")

csemMIRTBF_A_raw$Ex_Raw <- round(csemMIRTBF_A_raw$Ex_Raw)
csemMIRTBF_A_raw_Aggre <- aggregate(csemMIRTBF_A_raw$Raw_CSEM, by=list(Category=csemMIRTBF_A_raw$Ex_Raw), FUN=function(x){sqrt(sum(x^2)/length(x))})
csemMIRTBF_A_raw_Aggre

names(csemMIRTBF_A_raw_Aggre) <- c("raw Score", "csemMIRTBF_A_raw_Aggre")

write.xlsx(csemMIRTBF_A_raw_Aggre, "csemMIRTBF_A_raw_Aggre.xlsx", sheetName="csemMIRTBF_A_raw_Aggre")


## BF form B
csemMIRTBF_B <- read.table("EMBF_CSEM_M_B.OUT", header = T)
csemMIRTBF_B_raw <- csemMIRTBF_B[,c("Ex_Raw", "Raw_CSEM")]

write.xlsx(csemMIRTBF_B_raw, "csemMIRTBF_B_raw.xlsx", sheetName="csemMIRTBF_B_raw")

csemMIRTBF_B_raw$Ex_Raw <- round(csemMIRTBF_B_raw$Ex_Raw)
csemMIRTBF_B_raw_Aggre <- aggregate(csemMIRTBF_B_raw$Raw_CSEM, by=list(Category=csemMIRTBF_B_raw$Ex_Raw), FUN=function(x){sqrt(sum(x^2)/length(x))})
csemMIRTBF_B_raw_Aggre

names(csemMIRTBF_B_raw_Aggre) <- c("raw Score", "csemMIRTBF_B_raw_Aggre")

write.xlsx(csemMIRTBF_B_raw_Aggre, "csemMIRTBF_B_raw_Aggre.xlsx", sheetName="csemMIRTBF_B_raw_Aggre")








### MIRT P method     Scale score

## SS form A
csemMIRTSS_A <- read.table("EM_SS_CSEM_m_A.OUT", header = T)
csemMIRTSS_A_ss <- csemMIRTSS_A[,c("Ex_Scale", "Scale_CSEM")]

write.xlsx(csemMIRTSS_A_ss, "csemMIRTSS_A_ss.xlsx", sheetName="csemMIRTSS_A_ss")

csemMIRTSS_A_ss$Ex_Scale <- round(csemMIRTSS_A_ss$Ex_Scale)
csemMIRTSS_A_ss_Aggre <- aggregate(csemMIRTSS_A_ss$Scale_CSEM, by=list(Category=csemMIRTSS_A_ss$Ex_Scale), FUN=function(x){sqrt(sum(x^2)/length(x))})
csemMIRTSS_A_ss_Aggre

names(csemMIRTSS_A_ss_Aggre) <- c("scale Score", "csemMIRTSS_A_ss_Aggre")

write.xlsx(csemMIRTSS_A_ss_Aggre, "csemMIRTSS_A_ss_Aggre.xlsx", sheetName="csemMIRTSS_A_ss_Aggre")


## SS form B
csemMIRTSS_B <- read.table("EM_SS_CSEM_m_B.OUT", header = T)
csemMIRTSS_B_ss <- csemMIRTSS_B[,c("Ex_Scale", "Scale_CSEM")]

write.xlsx(csemMIRTSS_B_ss, "csemMIRTSS_B_ss.xlsx", sheetName="csemMIRTSS_B_ss")

csemMIRTSS_B_ss$Ex_Scale <- round(csemMIRTSS_B_ss$Ex_Scale)
csemMIRTSS_B_ss_Aggre <- aggregate(csemMIRTSS_B_ss$Scale_CSEM, by=list(Category=csemMIRTSS_B_ss$Ex_Scale), FUN=function(x){sqrt(sum(x^2)/length(x))})
csemMIRTSS_B_ss_Aggre

names(csemMIRTSS_B_ss_Aggre) <- c("scale Score", "csemMIRTSS_B_ss_Aggre")

write.xlsx(csemMIRTSS_B_ss_Aggre, "csemMIRTSS_B_ss_Aggre.xlsx", sheetName="csemMIRTSS_B_ss_Aggre")




## BF form A
csemMIRTBF_A <- read.table("EMBF_CSEM_M_A.OUT", header = T)
csemMIRTBF_A_ss <- csemMIRTBF_A[,c("Ex_Scale", "Scale_CSEM")]

write.xlsx(csemMIRTBF_A_ss, "csemMIRTBF_A_ss.xlsx", sheetName="csemMIRTBF_A_ss")

csemMIRTBF_A_ss$Ex_Scale <- round(csemMIRTBF_A_ss$Ex_Scale)
csemMIRTBF_A_ss_Aggre <- aggregate(csemMIRTBF_A_ss$Scale_CSEM, by=list(Category=csemMIRTBF_A_ss$Ex_Scale), FUN=function(x){sqrt(sum(x^2)/length(x))})
csemMIRTBF_A_ss_Aggre

names(csemMIRTBF_A_ss_Aggre) <- c("scale Score", "csemMIRTBF_A_ss_Aggre")

write.xlsx(csemMIRTBF_A_ss_Aggre, "csemMIRTBF_A_ss_Aggre.xlsx", sheetName="csemMIRTBF_A_ss_Aggre")


## BF form B
csemMIRTBF_B <- read.table("EMBF_CSEM_M_B.OUT", header = T)
csemMIRTBF_B_ss <- csemMIRTBF_B[,c("Ex_Scale", "Scale_CSEM")]

write.xlsx(csemMIRTBF_B_ss, "csemMIRTBF_B_ss.xlsx", sheetName="csemMIRTBF_B_ss")

csemMIRTBF_B_ss$Ex_Scale <- round(csemMIRTBF_B_ss$Ex_Scale)
csemMIRTBF_B_ss_Aggre <- aggregate(csemMIRTBF_B_ss$Scale_CSEM, by=list(Category=csemMIRTBF_B_ss$Ex_Scale), FUN=function(x){sqrt(sum(x^2)/length(x))})
csemMIRTBF_B_ss_Aggre

names(csemMIRTBF_B_ss_Aggre) <- c("scale Score", "csemMIRTBF_B_ss_Aggre")

write.xlsx(csemMIRTBF_B_ss_Aggre, "csemMIRTBF_B_ss_Aggre.xlsx", sheetName="csemMIRTBF_B_ss_Aggre")

















# CSEM Lord
csemLord <- CSEMLord(31)
csemLord
plot(csemLord$rawScore, csemLord$csemLord)

library(xlsx)
write.xlsx(csemLord, "csemLord_31.xlsx", sheetName="csemLord_31")






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



# test reliability IRT
TestRelIRT_A <- TestRelIRT(itemPara_A_UIRT, convTable_A)
TestRelIRT_A

TestRelIRT_B <- TestRelIRT(itemPara_B_UIRT, convTable_B)
TestRelIRT_B



write.xlsx(TestRelIRT_A$`Conditional SEMs`, "testCSEM_single_A.xlsx")
write.xlsx(TestRelIRT_B$`Conditional SEMs`, "testCSEM_single_B.xlsx")


TestRelIRT_A_SS <- TestRelIRT_A$`Conditional SEMs`[,c("Ex_Scale", "Scale_CSEM")]
TestRelIRT_A_SS$roundedSS <- round(TestRelIRT_A_SS$Ex_Scale)
TestRelIRT_A_SS_Aggre <- aggregate(TestRelIRT_A_SS$Scale_CSEM, by=list(Category=TestRelIRT_A_SS$roundedSS), FUN=function(x){sqrt(sum(x^2)/length(x))})

TestRelIRT_B_SS <- TestRelIRT_B$`Conditional SEMs`[,c("Ex_Scale", "Scale_CSEM")]
TestRelIRT_B_SS$roundedSS <- round(TestRelIRT_B_SS$Ex_Scale)
TestRelIRT_B_SS_Aggre <- aggregate(TestRelIRT_B_SS$Scale_CSEM, by=list(Category=TestRelIRT_B_SS$roundedSS), FUN=function(x){sqrt(sum(x^2)/length(x))})


names(TestRelIRT_A_SS_Aggre) <- names(TestRelIRT_B_SS_Aggre) <- c("roundedSS", "cssemKolen")


write.xlsx(TestRelIRT_A_SS_Aggre, "TestRelIRT_A_SS_Aggre.xlsx")
write.xlsx(TestRelIRT_B_SS_Aggre, "TestRelIRT_B_SS_Aggre.xlsx")



TestRelIRT_A_SS <- TestRelIRT_A$`Conditional SEMs`[,c("Ex_Raw", "Raw_CSEM")]
TestRelIRT_A_SS$roundedSS <- round(TestRelIRT_A_SS$Ex_Raw)
TestRelIRT_A_SS_Aggre <- aggregate(TestRelIRT_A_SS$Raw_CSEM, by=list(Category=TestRelIRT_A_SS$roundedSS), FUN=function(x){sqrt(sum(x^2)/length(x))})

TestRelIRT_B_SS <- TestRelIRT_B$`Conditional SEMs`[,c("Ex_Raw", "Raw_CSEM")]
TestRelIRT_B_SS$roundedSS <- round(TestRelIRT_B_SS$Ex_Raw)
TestRelIRT_B_SS_Aggre <- aggregate(TestRelIRT_B_SS$Raw_CSEM, by=list(Category=TestRelIRT_B_SS$roundedSS), FUN=function(x){sqrt(sum(x^2)/length(x))})


names(TestRelIRT_A_SS_Aggre) <- names(TestRelIRT_B_SS_Aggre) <- c("Raw", "cssemKolen")


write.xlsx(TestRelIRT_A_SS_Aggre, "TestRelIRT_A_Raw_Aggre.xlsx")
write.xlsx(TestRelIRT_B_SS_Aggre, "TestRelIRT_B_Raw_Aggre.xlsx")







# KolenRelIRT_A_UIRT <- TestRelIRT(itemPara_A_UIRT, convTable_A)
# KolenRelIRT_A_UIRT
#
# KolenRelIRT_B_UIRT <- KolenRelIRT(itemPara_B_UIRT, convTable_B)
# KolenRelIRT_B_UIRT
#
#
# CSSEMKolen_A <- as.data.frame(CSSEMKolen(itemPara_A_UIRT, convTable_A))
# CSSEMKolen_A
# write.xlsx(CSSEMKolen_A, "CSSEM_Kolen_A_MUIRT.xlsx")
#
# CSSEMKolen_B <- as.data.frame(CSSEMKolen(itemPara_B_UIRT, convTable_B))
# CSSEMKolen_B
# write.xlsx(CSSEMKolen_B, "CSSEM_Kolen_B_MUIRT.xlsx")

# BI-Factor General

# read item parameters
# form A
itemPara_BF_G <- read.table("TestData/SpanishLit_prm_A_BF.txt")[,c(7:11)]

# form B
# itemPara_BF_G <- read.table("TestData/SpanishLit_prm_B_BF.txt")[,c(7:11)]

names(itemPara_BF_G) <- c("b", "ag","a1","a2", "a3") # a1 is primary
itemPara_BF_G$a <- c(itemPara_BF_G$a1[1:13], itemPara_BF_G$a2[14:25], itemPara_BF_G$a3[26:31])

itemPara_BF_G[,"b"] <- -itemPara_BF_G[,"b"]/(itemPara_BF_G[,"ag"]) #######+ itemPara_BF_G[,"a"]
itemPara_BF_G[,"a"] <- itemPara_BF_G[,"a"]/1.702

itemPara_BF_G$a1[1:13] <- itemPara_BF_G$a[1:13]
itemPara_BF_G$a2[14:25] <- itemPara_BF_G$a[14:25]
itemPara_BF_G$a3[26:31] <- itemPara_BF_G$a[26:31]
itemPara_BF_G$ag <- itemPara_BF_G$ag/1.702

itemPara_BF_G <- itemPara_BF_G[,c("b", "ag")]
names(itemPara_BF_G) <- c("b", "a")


TestRelIRT(itemPara_BF_G)

KolenRelIRT(itemPara_BF_G, convTable_A)
# KolenRelIRT(itemPara_BF_G, convTable_B)

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



# num of quadratures
numOfQuad <- 41

### bivariate normal distribution

nodes <- seq(-4, 4, length.out = numOfQuad)
nodesM <- as.matrix(expand.grid(nodes,nodes))
weightsUnwtd <- dmvnorm(nodesM, c(0,0), cormat, log=FALSE) # mvtnorm
nodesM <- as.data.frame(nodesM)
nodesM$weightsWtd <- weightsUnwtd #/ sum(weightsUnwtd)
nodesM$weightsWtd <- round(nodesM$weightsWtd, 10)

write.table(format(nodesM, digits=10), "thetatwo.txt", col.names = F, row.names = F, quote = F)


### trivariate normal distribution



# set nodes and weights
nodes <- seq(-4, 4, length.out = numOfQuad)
nodesM <- as.matrix(expand.grid(nodes,nodes,nodes))
weightsUnwtd <- dmvnorm(nodesM, c(0,0,0), cormat_A, log=FALSE) # mvtnorm
nodesM <- as.data.frame(nodesM)
nodesM$weightsWtd <- weightsUnwtd #/ sum(weightsUnwtd)


# nodesM <- nodesM[,c("Var3", "Var2", "Var1", "weightsWtd")]

nodesM$Var3 <- format(nodesM$Var3, scientific = FALSE)
nodesM$Var2 <- format(nodesM$Var2, scientific = FALSE)
nodesM$Var1 <- format(nodesM$Var1, scientific = FALSE)
nodesM$weightsWtd <- round(nodesM$weightsWtd, 10)

write.table(format(nodesM, digits=10), "thetathreeformA.txt", col.names = F, row.names = F, quote = F)




# set nodes and weights
nodes <- seq(-4, 4, length.out = numOfQuad)
nodesM <- as.matrix(expand.grid(nodes,nodes,nodes))
weightsUnwtd <- dmvnorm(nodesM, c(0,0,0), cormat_B, log=FALSE) # mvtnorm
nodesM <- as.data.frame(nodesM)
nodesM$weightsWtd <- weightsUnwtd #/ sum(weightsUnwtd)


# nodesM <- nodesM[,c("Var3", "Var2", "Var1", "weightsWtd")]

nodesM$Var3 <- format(nodesM$Var3, scientific = FALSE)
nodesM$Var2 <- format(nodesM$Var2, scientific = FALSE)
nodesM$Var1 <- format(nodesM$Var1, scientific = FALSE)
nodesM$weightsWtd <- round(nodesM$weightsWtd, 10)

write.table(format(nodesM, digits=10), "thetathreeformB.txt", col.names = F, row.names = F, quote = F)




### four factor normal

# set nodes and weights
nodes <- seq(-4, 4, length.out = numOfQuad)
nodesM <- as.matrix(expand.grid(nodes,nodes,nodes,nodes))
weightsUnwtd <- dmvnorm(nodesM, c(0,0,0,0), diag(4), log=FALSE) # mvtnorm
nodesM <- as.data.frame(nodesM)
nodesM$weightsWtd <- weightsUnwtd #/ sum(weightsUnwtd)


# nodesM <- nodesM[,c("Var3", "Var2", "Var1", "weightsWtd")]
nodesM$Var4 <- format(nodesM$Var4, scientific = FALSE)
nodesM$Var3 <- format(nodesM$Var3, scientific = FALSE)
nodesM$Var2 <- format(nodesM$Var2, scientific = FALSE)
nodesM$Var1 <- format(nodesM$Var1, scientific = FALSE)
nodesM$weightsWtd <- round(nodesM$weightsWtd, 10)

write.table(format(nodesM, digits=10), "thetafourformA.txt", col.names = F, row.names = F, quote = F)
