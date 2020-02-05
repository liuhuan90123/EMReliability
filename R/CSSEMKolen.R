### CSSEM Kolen's method

library(statmod)

# transform item parameters to the logistic metric
names(itemPara) <- c("b", "a")
itemPara[,"a"] <- itemPara[,"a"]/1.701

# weights and nodes
quadPoints <- gauss.quad.prob(41, dist = "normal", mu = 0, sigma = 1)

# replicate item parameter and theta
itemParaRep <-itemPara[rep(seq_len(nrow(itemPara)), each = 41),]
itemParaRep$theta <- rep(quadPoints$nodes, each = 1, length.out = 41*nrow(itemPara))

# calculate information by theta
itemParaRep <- within(itemParaRep, {
  P = 0 + (1 - 0) / (1 + exp(-1.701 * a * (theta - b)))
  Q = 1 - P
  PQ = P * Q
  info = 1.701**2 * a**2 * P * Q
})


## observed score variance

# order by theta
itemParaRep <- itemParaRep[order(itemParaRep$theta),]

# num of quadratures
numQuad <- 41

# define matrix of marginal distribution of theta
fxTheta <- matrix(NA, nrow = numQuad, ncol = 41) # 41 num of quadratures, 41 num of items

# calculate marginal distribution of theta
for (i in 1:numQuad){

  probs <- matrix(c(itemParaRep[(1 + 40 * (i - 1)):(40 * i),]$P,
                    itemParaRep[(1 + 40 * (i - 1)):(40 * i),]$Q),
                  nrow = 40, ncol = 2, byrow = FALSE)
  cats <- c(rep(2, 40))
  fxTheta[i, ] <- wlord(probs,cats)

}

# transform to data frame
fxTheta <- as.data.frame(fxTheta)

# add quadrature weights
fxTheta$weights <- quadPoints$weights


### IRT: CSEM for scale scores: Kolen's method ---------------------------------------------------------

fx_CSSEM_IRT <- as.data.frame(t(fxTheta[,-42]))



# conversion table

# read conversion table from cvs file
# convTable <- read.csv("TestData/ConversionTableFormX.csv")
# convTable$roundedSS <- round(convTable$unroundedSS)





SScon <- read.xlsx("conversion_table_Form A.xlsx",1)
SScon$Rounded.SS <- round(SScon$Unrounded.SS,0)


fx_CSSEM_IRT$SS <- rev(SScon$Rounded.SS)


# true scale score

fx_CSSEM_IRT_n <- fx_CSSEM_IRT %>%
  mutate_each(funs(. * SS), starts_with("V"))

# # add quadrature weights
# fxTheta$weights <- quadPoints$weights
#
# # calculate weighted distribution
# fxThetaWeighted <- apply(fxTheta[,1:41], 2, function(x) x * fxTheta[,"weights"])




fx_CSSEM_IRT_n <- rbind(fx_CSSEM_IRT, colSums(fx_CSSEM_IRT_n))

# CSSEM condtional on theta

CSSEM_IRT <- matrix(NA,nrow = 15, ncol = 1)

for (i in 1:15){
  # i = 1
  CSSEM_IRT[i, 1] <- sqrt(sum((fx_CSSEM_IRT_n[c(1:41),16] - fx_CSSEM_IRT_n[42, i])^2 * fx_CSSEM_IRT_n[c(1:41),i]))
}

# avarage CSSEM across theta distribution

ErrorVarIRT <- sum(CSSEM_IRT^2 * quad_points$weights)

# variance of scale score

fx_prXi <- fxTheta %>%
  mutate_each(funs(. * weights), starts_with("V"))

meanSS <- sum(rev(SScon$Rounded.SS) * colSums(fx_prXi[,-42]))

SSVarIRT <- sum((rev(SScon$Rounded.SS) - meanSS)^2 * colSums(fx_prXi[,-42]))

# reliability

RelIRTSSKolen <- 1 - ErrorVarIRT / SSVarIRT
RelIRTSSKolen


CSSEM_IRT <- as.data.frame(CSSEM_IRT)
CSSEM_IRT$theta <- quad_points$nodes

##true scale score
CSSEM_IRT$TrueSS <- colSums(fx_CSSEM_IRT_n)[1:15]

# plot CSSEM using Kolen's method
png("CSSEM_KolenIRT_A.png",  width = 799, height = 596)


K <- ggplot(CSSEM_IRT, aes(x = TrueSS, y = V1)) +
  geom_point(size = 2) +
  scale_x_continuous(name = "True Scale Score", breaks  = seq(100,  130, 5)) +
  scale_y_continuous(name = "CSSEM_Kolen IRT Method") +
  theme_bw()

print(K)
dev.off()








