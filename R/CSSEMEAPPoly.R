### CSSEM IRT EAP Polynomial Method




# read item parameter
itemPara <- read.table("TestData/ItemParaFormX.txt")

# read conversion table from cvs file
convTable <- read.csv("TestData/ConversionTableFormX.csv")
convTable$roundedSS <- round(convTable$unroundedSS)
convTableSub <- convTable[,c("thetaScore", "roundedSS")]

source("R/PolynomialMethod.R")

library(statmod)

CSSEMEAPPoly <- function(itemPara, convTable, K){

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

  itemParaInfo$csemEAP <- sqrt(1/(itemParaInfo$infoSum+1))

  # change name to fit Polynomial Method function
  names(itemParaInfo) <- c("rawScore", "infoSum", "roundedSS", "csemLord")

  PolynomialMethod(itemParaInfo, K)

}

K <- 20
CSSEMEAPPoly(itemPara, convTableSub, K)



### select column without negative values?

### reliability

# rounded scale score reliability; k = 4

cssemEAPPolyDat <- CSSEMEAPPoly(itemPara, convTableSub, K)
cssemEAPPolyDat$rawScore <- c(0:40)


### raw score frequency

# read raw data
rawData <- read.table("TestData/RawDataFormX.txt", header = F, sep = " ")

rawFreq <- as.data.frame(table(rowSums(rawData)))
names(rawFreq) <- c("rawScore", "freq")



cssemEAPPolyDat <- merge(cssemEAPPolyDat, rawFreq, by = "rawScore")

# weight
cssemEAPPolyDat$wt <- cssemEAPPolyDat$freq / sum(cssemEAPPolyDat$freq)

library(SDMTools)

# select k = 4
EAP <- 1 - sum(cssemEAPPolyDat$cssemPolyk4^2 * cssemEAPPolyDat$wt)/wt.var(cssemEAPPolyDat$roundedSS, cssemEAPPolyDat$wt)
EAP







