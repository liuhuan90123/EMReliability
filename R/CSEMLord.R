### CSEM Lord Method

library(SDMTools)
library(ggplot2)

# read raw data
rawData <- read.table("TestData/RawDataFormX.txt", header = F, sep = " ")

# read conversion table from cvs file
convTable <- read.csv("TestData/ConversionTableFormX.csv")

convTable$roundedSS <- round(convTable$unroundedSS)


ss <- convTable[,c("rawScore", "roundedSS")]

csem_ss <- matrix(nrow = 41, ncol = 2)

for (i in 0:40){

  # pi =20
  pi = i/40

  bidat <- as.data.frame(c(0:40))

  names(bidat) <- "raw"

  bidat$pr <- choose(40, bidat$raw) * (pi)^bidat$raw * (1 - pi)^(40-bidat$raw)

  bidat <- cbind(bidat, ss)

  bidat$xpr <- bidat$rawScore * bidat$pr
  bidat$x2pr <- bidat$rawScore^2 * bidat$pr
  bidat$spr <- bidat$roundedSS * bidat$pr
  bidat$s2pr <- bidat$roundedSS^2 * bidat$pr

  sd_x <- sqrt(40/(40-1)) * sqrt(sum(bidat$x2pr) - sum(bidat$xpr)^2)
  sd_s <- sqrt(40/(40-1)) * sqrt(sum(bidat$s2pr) - sum(bidat$spr)^2)

  csem_ss[i+1, 1] <- sd_x
  csem_ss[i+1, 2] <- sd_s

}


csem_ss <- as.data.frame(csem_ss)

names(csem_ss) <- c("csemLord", "csemBinomial")

csem_ss$rawScore <- c(0:40)

# calculate raw score frequence

rawFreq <- as.data.frame(table(rowSums(rawData)))
names(rawFreq) <- c("rawScore", "freq")



convTable <- merge(convTable, rawFreq, by = "rawScore")
convTable$wt <- convTable$freq / sum(convTable$freq)
convTable <- merge(convTable, csem_ss, by = "rawScore")



# raw scale reliability using Lord CSEM

LordRel <- 1 - sum(convTable$csemLord^2 * convTable$wt)/wt.var(convTable$rawScore, convTable$wt)

# Scale score reliability using binomial method

SSRelCTTBinomial <- 1 - sum(convTable$csemBinomial^2 * convTable$wt)/wt.var(convTable$roundedSS, convTable$wt)


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

convTable$varBinomial <- convTable$csemBinomial^2
convTableBinomial <- aggregate(convTable$varBinomial, by=list(Category=convTable$roundedSS), FUN=mean)
names(convTableBinomial) <- c("roundedSS", "varBinomialAvg")
convTableBinomial$csemBinomialAvg <- sqrt(convTableBinomial$varBinomialAvg)


# plot for Binomial CSSEM
png("CSSEMCTTBinomial.png",  width = 799, height = 596)

CSSEMCTTBinomialpng <- ggplot(convTableBinomial, aes(x = roundedSS, y = csemBinomialAvg)) +
  geom_point() +
  scale_x_continuous(name = "Rounded Scale Score", breaks  = seq(100, 130, 5)) +
  scale_y_continuous(name = "CSSEM Binomial Method") +
  theme_bw()

print(CSSEMCTTBinomialpng)
dev.off()



## polynomial method

m1 <- lm(Rounded.SS ~ 1 + RawScore + I(RawScore^2) + I(RawScore^3), CSEM_raw_uniq_ord)
summary(m1)

# first derivative
f <- expression(8.842e-04 * x^3 -5.671e-02* x^2 + 1.618e+00 * x + 9.899e+01)
D(f, "x")

# apply formula
CSEM_raw_uniq_ord$fx <- 0.0008842 * (3 * CSEM_raw_uniq_ord$RawScore^2) - 0.05671 * (2 * CSEM_raw_uniq_ord$RawScore) + 1.618
CSEM_raw_uniq_ord$CSSEM_poly <- CSEM_raw_uniq_ord$fx * CSEM_raw_uniq_ord$CSEM_Lord


### Scale score reliability with polynomial method
CTT_SS_Reliability_poly_X <- 1 - sum(CSEM_raw_uniq_ord$CSSEM_poly^2 * CSEM_raw_uniq_ord$wt)/wt.var(CSEM_raw_uniq_ord$Rounded.SS, CSEM_raw_uniq_ord$wt)




CSEMLord <- function(rawData, ){



}



