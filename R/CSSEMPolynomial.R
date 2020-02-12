### CSSEM Polynomial Method

source("R/CSEMLord.R")
source("R/PolynomialMethod.R")

# Note: 1. cssemDat should include rawScore, roundedSS, and Lordcsem
#       2. conversion table should include raw score and rounded scale score


# read conversion table from cvs file
convTable <- read.csv("TestData/ConversionTableFormX.csv")
convTable$roundedSS <- round(convTable$unroundedSS)
convTableSub <- convTable[,c("rawScore", "roundedSS")]

numOfItem <- 40
K <- 15



CSSEMPolynomial <- function(numOfItem, convTable, K){

  csemLordDat <- CSEMLord(numOfItem)
  cssemDat <- merge(csemLordDat, convTableSub, by = "rawScore")
  PolynomialMethod(cssemDat, K)

}


CSSEMPolynomial(numOfItem, convTableSub, K)



#### Plot --------------------------------------
library(ggplot2)



### aggregate many to one ------------


cssemDatWide <- CSSEMPolynomial(40, convTableSub, K)



cssemDat <- cssemDatWide # test
k <- 13 # test, accepted maximum + 1
cssemDat <- cssemDat[,c(3,5:(5+k-2))]
cssemDatAggre <- as.data.frame(apply(cssemDat[,c(-1)], 2, function(x) aggregate(x, by=list(Category=cssemDat$roundedSS), FUN=mean)))

cssemDatAggre <- cssemDatAggre[,c(1, seq(2, 24, 2))]




### range and confidence interval  with k from 1 to maximum accepted ----------


cssemDatAggre$max <- apply(cssemDatAggre[,c(-1)], 1, max)
cssemDatAggre$min <- apply(cssemDatAggre[,c(-1)], 1, min)
cssemDatAggre$mean <- apply(cssemDatAggre[,c(-1)], 1, mean)
cssemDatAggre$median <- apply(cssemDatAggre[,c(-1)], 1, median)
cssemDatAggre$sd <- apply(cssemDatAggre[,c(-1)], 1, sd)

# 95% confidence interval
cssemDatAggre_1 <- within(cssemDatAggre, {
       lower = mean - 1.96 * sd
       upper = mean + 1.96 * sd
       })


### plot confidence interval

# k = 1 to maximum accepted

ggplot(cssemDatAggre_1, aes(x = cssemPolyk1.Category, y = mean)) +
  geom_line(colour ="blue") +
  geom_point(colour="blue") +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill= "blue", alpha = 0.2) +
  scale_x_continuous(name = "Rounded Scale Score", breaks  = seq(100, 130, 5)) +
  scale_y_continuous(name = "CSSEM Polynomial Method") +
  theme_bw()






### range and confidence interval  with k from 3 to maximum accepted ---33333333333-------


cssemDatAggre$max <- apply(cssemDatAggre[,c(-1:-3)], 1, max)
cssemDatAggre$min <- apply(cssemDatAggre[,c(-1:-3)], 1, min)
cssemDatAggre$mean <- apply(cssemDatAggre[,c(-1:-3)], 1, mean)
cssemDatAggre$median <- apply(cssemDatAggre[,c(-1:-3)], 1, median)
cssemDatAggre$sd <- apply(cssemDatAggre[,c(-1:-3)], 1, sd)

# 95% confidence interval
cssemDatAggre <- within(cssemDatAggre, {
  lower = mean - 1.96 * sd
  upper = mean + 1.96 * sd
})


### plot confidence interval

# k = 1 to maximum accepted

ggplot(cssemDatAggre, aes(x = cssemPolyk1.Category, y = mean)) +
  geom_line(colour ="blue") +
  geom_point(colour="blue") +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill= "blue", alpha = 0.2) +
  scale_x_continuous(name = "Rounded Scale Score", breaks  = seq(100, 130, 5)) +
  scale_y_continuous(name = "CSSEM Polynomial Method") +
  theme_bw()



# hist(as.numeric((cssemDatAggre[11,c(2:13)])))


### questions ---------------------------------------------

## maximum and minimum of k values
## R square matrix : check r squared change of k from 1 to 40
## unique scale scores
## range and confidence interval for each rounded scale score, print confidence interval


### plot ---------------------------------------------------

k <- 13 # The maximum accepted K

cssemDatLong <- reshape(cssemDatAggre,
                        direction = "long",
                        varying = list(names(cssemDatAggre)[2:k]),
                        v.names = "cssempoly",
                        idvar = c("cssemPolyk1.Category"),
                        timevar = "Kvalue",
                        times = 1:(k-1))


library(ggplot2)
ggplot(cssemDatLong, aes(x = cssemPolyk1.Category, y = cssempoly, color = factor(Kvalue))) +
  geom_point() +
  scale_x_continuous(name = "Rounded Scale Score", breaks  = seq(100, 130, 5)) +
  scale_y_continuous(name = "CSSEM Polynomial Method") +
  geom_line() +
  theme_bw() +
  labs(colour="K value")







