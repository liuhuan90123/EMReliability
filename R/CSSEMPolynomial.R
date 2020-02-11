### CSSEM Polynomial Method

source("R/CSEMLord.R")
source("R/PolynomialMethod")

# Note: 1. cssemDat should include rawScore, roundedSS, and Lordcsem
#       2. conversion table should include raw score and rounded scale score


# read conversion table from cvs file
convTable <- read.csv("TestData/ConversionTableFormX.csv")
convTable$roundedSS <- round(convTable$unroundedSS)
convTableSub <- convTable[,c("rawScore", "roundedSS")]



CSSEMPolynomial <- function(numOfItem, convTable){

  csemLordDat <- CSEMLord(numOfItem)
  cssemDat <- merge(csemLordDat, convTableSub, by = "rawScore")
  PolynomialMethod(cssemDat)

}

CSSEMPolynomial(40, convTableSub)







### agrregate to unique scale score, cssem ------------------

### R square in the variable name,
### confidence interval for each rounded scale score
### try 20, 30, 50, 100, and 200 in the K list


cssemDatLong <- reshape(cssemDat,
                        direction = "long",
                        varying = list(names(cssemDat)[5:14]),
                        v.names = "cssem",
                        idvar = c("rawScore", "csemLord", "roundedSS", "fx"),
                        timevar = "Kvalue",
                        times = 1:10)


library(ggplot2)
ggplot(cssemDatLong, aes(x = rawScore, y = cssem, color = factor(Kvalue))) +
  geom_point() +
  scale_x_continuous(name = "Rounded Scale Score", breaks  = seq(100, 130, 5)) +
  scale_y_continuous(name = "CSSEM Polynomial Method") +
  theme_bw()



