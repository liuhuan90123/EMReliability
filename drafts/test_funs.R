library(statmod)



# read item parameters from txt file

itemPara <- read.table("TestData/ItemParaFormX.txt")
MarginalRelEAP(itemPara)
MarginalRelMLE(itemPara)

### Test Realiability of IRT 2PL

library(statmod)
#library(classify)

# read item parameters from txt file
itemPara <- read.table("TestData/ItemParaFormX.txt")


TestRelIRT(itemPara)
