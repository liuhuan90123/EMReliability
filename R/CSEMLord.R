### CSEM Lord Method

library(SDMTools)
library(ggplot2)

# read raw data
rawData <- read.table("TestData/RawDataFormX.txt", header = F, sep = " ")

# read conversion table from cvs file
convTable <- read.csv("TestData/ConversionTableFormX.csv")

convTable$roundedSS <- round(convTable$unroundedSS)

convTableSub <- convTable[,c("rawScore", "roundedSS")]

csemDat <- matrix(nrow = 41, ncol = 2)

for (i in 0:40){

  pi = i/40

  binoDat <- as.data.frame(c(0:40))

  names(binoDat) <- "raw"

  binoDat$prob <- choose(40, binoDat$raw) * (pi)^binoDat$raw * (1 - pi)^(40-binoDat$raw)

  binoDat <- cbind(binoDat, convTableSub)

  binoDat$xprob <- binoDat$rawScore * binoDat$prob
  binoDat$x2prob <- binoDat$rawScore^2 * binoDat$prob
  binoDat$ssprob <- binoDat$roundedSS * binoDat$prob
  binoDat$ss2prob <- binoDat$roundedSS^2 * binoDat$prob

  seX <- sqrt(40/(40-1)) * sqrt(sum(binoDat$x2prob) - sum(binoDat$xprob)^2)
  seSS <- sqrt(40/(40-1)) * sqrt(sum(binoDat$ss2prob) - sum(binoDat$ssprob)^2)

  csemDat[i+1, 1] <- seX
  csemDat[i+1, 2] <- seSS

}


csemDat <- as.data.frame(csemDat)

names(csemDat) <- c("csemLord", "cssemBinomial")

csemDat$rawScore <- c(0:40)

# calculate raw score frequence

rawFreq <- as.data.frame(table(rowSums(rawData)))
names(rawFreq) <- c("rawScore", "freq")



convTable <- merge(convTable, rawFreq, by = "rawScore")
convTable$wt <- convTable$freq / sum(convTable$freq)
convTable <- merge(convTable, csemDat, by = "rawScore")



# raw scale reliability using Lord CSEM
LordRel <- 1 - sum(convTable$csemLord^2 * convTable$wt)/wt.var(convTable$rawScore, convTable$wt)

# Scale score reliability using binomial method
SSRelCTTBinomial <- 1 - sum(convTable$cssemBinomial^2 * convTable$wt)/wt.var(convTable$roundedSS, convTable$wt)


# plot for Lord CSEM
png("CSEMLORD.png",  width = 799, height = 596)

CSEMLORDpng <- ggplot(convTable, aes(x = rawScore, y = csemLord)) +
  geom_point() +
  scale_x_continuous(name = "Raw Score", breaks  = seq(0, 40, 5)) +
  scale_y_continuous(name = "CSEM Lord Method") +
  theme_bw()

print(CSEMLORDpng)
dev.off()


## calculate average CSSEM: average of variance, not standard error

convTable$varBinomial <- convTable$cssemBinomial^2
convTableBinomial <- aggregate(convTable$varBinomial, by=list(Category=convTable$roundedSS), FUN=mean)
names(convTableBinomial) <- c("roundedSS", "varBinomialAvg")
convTableBinomial$cssemBinomialAvg <- sqrt(convTableBinomial$varBinomialAvg)


# plot for Binomial CSSEM
png("CSSEMCTTBinomial.png",  width = 799, height = 596)

CSSEMCTTBinomialpng <- ggplot(convTableBinomial, aes(x = roundedSS, y = cssemBinomialAvg)) +
  geom_point() +
  scale_x_continuous(name = "Rounded Scale Score", breaks  = seq(100, 130, 5)) +
  scale_y_continuous(name = "CSSEM Binomial Method") +
  theme_bw()

print(CSSEMCTTBinomialpng)
dev.off()



## polynomial method

m1 <- lm(roundedSS ~ 1 + rawScore + I(rawScore^2) + I(rawScore^3), convTable)
summary(m1)

# first derivative
f <- expression(8.842e-04 * x^3 -5.671e-02* x^2 + 1.618e+00 * x + 9.899e+01)
D(f, "x")

# apply formula
convTable$fx <- 0.0008842 * (3 * convTable$rawScore^2) - 0.05671 * (2 * convTable$rawScore) + 1.618
convTable$cssemPoly <- convTable$fx * convTable$csemLord


### Scale score reliability with polynomial method
SSRelCTTPoly <- 1 - sum(convTable$cssemPoly^2 * convTable$wt)/wt.var(convTable$roundedSS, convTable$wt)



## calculate average CSSEM: average of variance, not standard error

convTable$varPoly <- convTable$cssemPoly^2
convTablePoly <- aggregate(convTable$varPoly, by=list(Category=convTable$roundedSS), FUN=mean)
names(convTablePoly) <- c("roundedSS", "varPolyAvg")
convTablePoly$cssemPolyAvg <- sqrt(convTablePoly$varPolyAvg)


# plot for Poly CSSEM
png("CSSEMCTTPoly.png",  width = 799, height = 596)

CSSEMCTTPolypng <- ggplot(convTablePoly, aes(x = roundedSS, y = cssemPolyAvg)) +
  geom_point() +
  scale_x_continuous(name = "Rounded Scale Score", breaks  = seq(100, 130, 5)) +
  scale_y_continuous(name = "CSSEM Poly Method") +
  theme_bw()

print(CSSEMCTTPolypng)
dev.off()

