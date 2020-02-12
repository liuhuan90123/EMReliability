

### Test Realiability of IRT 2PL
MarginalRelEAP()
MarginalRelMLE(itemPara)

#library(classify)

# read item parameters from txt file
itemPara <- read.table("TestData/ItemParaFormX.txt")

names(itemPara) <- c("b", "a")
itemPara[,"a"] <- itemPara[,"a"]/1.702
itemPara$c <- 0

itemPara <- itemPara[c("a", "b", "c")]


write.csv(itemPara, "itemParaLee.csv", col.names = F, row.names = F)



TestRelIRT(itemPara)



### CSEM binomial

# read conversion table from cvs file
convTable <- read.csv("TestData/ConversionTableFormX.csv")
convTable$roundedSS <- round(convTable$unroundedSS)
convTableSub <- convTable[,c("rawScore", "roundedSS")]


CSSEMBinomial(40, convTableSub)



source("R/CSEMMLE.R")
theta <- c(seq(-5, 5, 0.1))
CSEMMLE(itemPara, theta)




source("R/CSEMEAP.R")
CSEMEAP(itemPara, theta)
