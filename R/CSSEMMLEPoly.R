### CSSEM IRT MLE Polynomial Method


# read item parameter
itemPara <- read.table("TestData/ItemParaFormX.txt")

# read conversion table from cvs file
convTable <- read.csv("TestData/ConversionTableFormX.csv")
convTable$roundedSS <- round(convTable$unroundedSS)
convTableSub <- convTable[,c("thetaScore", "roundedSS")]

source("R/PolynomialMethod.R")

library(statmod)

CSSEMMLEPoly <- function(itemPara, convTable){

  # transform item parameters to the 1.702 metric
  names(itemPara) <- c("b", "a")
  itemPara[,"a"] <- itemPara[,"a"]/1.702

  # weights and nodes
  quadPoints <- gauss.quad.prob(41, dist = "normal", mu = 0, sigma = 1)

  # replace nodes with theta score in coversion table
  quadPoints$nodes <- convTable$thetaScore

  # replicate item parameter and theta
  itemParaRep <-itemPara[rep(seq_len(nrow(itemPara)), each = 41),]
  itemParaRep$theta <- rep(quadPoints$nodes, each = 1, length.out = 41*nrow(itemPara))

  # calculate information by theta
  itemParaRep <- within(itemParaRep, {
    P = 0 + (1 - 0) / (1 + exp(-1.702 * a * (theta - b)))
    Q = 1 - P
    PQ = P * Q
    info = 1.702**2 * a**2 * P * Q
  })

  # sum information by theta
  itemParaInfo <- aggregate(itemParaRep$info, by=list(Category=itemParaRep$theta), FUN=sum)
  names(itemParaInfo) <- c("thetaScore", "infoSum")

  # merge data
  itemParaInfo <- merge(itemParaInfo, convTable, by = "thetaScore")

  # transform inforamtion to csem

  itemParaInfo$csemMLE <- sqrt(1/itemParaInfo$infoSum)

  # change name to fit Polynomial Method function
  names(itemParaInfo) <- c("rawScore", "infoSum", "roundedSS", "csemLord")

  PolynomialMethod(itemParaInfo, 20)

}


CSSEMMLEPoly(itemPara, convTableSub)



### select column without negative values?




