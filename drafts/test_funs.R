library(statmod)



# read item parameters from txt file

itemPara <- read.table("TestData/ItemParaFormX.txt")
### Test Realiability of IRT 2PL
MarginalRelEAP()
MarginalRelMLE(itemPara)

#library(classify)

# read item parameters from txt file
itemPara <- read.table("TestData/ItemParaFormX.txt")


TestRelIRT(itemPara)
