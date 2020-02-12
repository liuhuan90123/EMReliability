
### Results for project ------------------


## IRT Test Reliability

library(classify)
source("R/TestRelIRT.R") # itemPara

# read item parameters from txt file
itemPara <- read.table("TestData/ItemParaFormX.txt")
TestRelIRT_A <- TestRelIRT(itemPara)

itemPara <- read.table("TestData/ItemParaFormY.txt")
TestRelIRT_B <- TestRelIRT(itemPara)






