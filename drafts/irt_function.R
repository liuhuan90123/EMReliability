# ### IRT reliability ###
#
# ### packages ------------------------------------------------------------------
#
# if (!require(irtplay)) install.packages("irtplay")
# if (!require(statmod)) install.packages("statmod")
# if (!require(spatstat)) install.packages("spatstat")
# if (!require(classify)) install.packages("classify")
#
# library(classify)
# library(spatstat)
# library(statmod)
# library(xlsx)
# library(tidyverse)
# library(irtreliability)
# library(stats4)
# library(lattice)
# library(mirt)
# library(rpf)
# library(irtplay)
# library(ltm)
# library(dplyr)
#
#
# ### functions --------------------------------------
#
# # read.flexmirt , rpf
# # gauss.quad.prob , statmod
# # %>% , magrittr
# # wlord , classify
#
# ### data ----------------------------------------------------------------------
#
# rawdata <- read.table("FormX_40_3000.txt", header = F, sep = " ")
#
#
# ### item parameters ----------------------------------------------------------------------
#
# ItemPara <- matrix(read.flexmirt("UShistory_X-prm.txt")$G1$param, ncol = 2, byrow = TRUE)
# ItemPara[,1] <- ItemPara[,1]/1.701 # items are on the logistic metric
# ItemPara <- as.data.frame(ItemPara)
# names(ItemPara) <- c("a_2PL", "b_2PL")

#' IRT reliability Function
#'
#' This function allows you to calculate IRT reliability coefficients
#' @param ItemPara data frame containing items parameters on 1.7 logistic metric
#' @param rawdata data frame containing raw data, 0/1
#' @keywords irt, reliability
#' @export
#' @examples
#' IRTRel()



### function structure -------------------------------------------------------------------

IRTRel <- function(ItemPara, rawdata){

### weights and nodes ---------------------------------------------------------

quad_points <- gauss.quad.prob(15, dist = "normal", mu = 0, sigma = 1)


### replicate quadrature points for each theta ---------------------------------

ItemPara_rep <-ItemPara[rep(seq_len(nrow(ItemPara)), each = 15),]
ItemPara_rep$theta <- rep(quad_points$nodes, each = 1, length.out = 600)


### calculate P values ----------------------------------------------------------

ItemPara_rep <- ItemPara_rep %>%
  transform(P = 0 + (1 - 0) / (1 + exp(-1.701 * a_2PL * (theta - b_2PL)))) %>%
  transform(Q = 1 - P) %>%
  transform(I = 1.701**2 * a_2PL**2 * P * Q) %>%
  transform(PQ = P * Q)


### marginal reliability with information -----------------------------------------

ItemPara_rep$weights <- quad_points$weights
ItemPara_rep$I_weighted <- ItemPara_rep$weights * ItemPara_rep$I

ItemPara_infor <- ItemPara_rep %>%
  group_by(theta) %>%
  summarise(tot_I = sum(I_weighted),
            tot_PQ = sum(PQ),
            tot_P = sum(P))

MarginalReliability <- sum(ItemPara_infor$tot_I) / (sum(ItemPara_infor$tot_I) + 1)

### test reliability: true score variance ----------------------------------------

ItemPara_infor$weights <- quad_points$weights
VarTrue <- sum((ItemPara_infor$tot_P)^2 * ItemPara_infor$weights)-(sum(ItemPara_infor$tot_P * ItemPara_infor$weights))^2


### test reliability: error variance ------------------------------------------------

VarError <- sum(ItemPara_infor$tot_PQ * ItemPara_infor$weights)


### test reliability: observed score variance, model based  ----------------------------------------

ItemPara_rep <- ItemPara_rep[order(ItemPara_rep$theta),]

fx_theta <- matrix(NA, nrow = 15, ncol = 41)

for (i in 1:15){

  probs <- matrix(c(ItemPara_rep[(1 + 40 * (i - 1)):(40 * i),]$P, ItemPara_rep[(1 + 40 * (i - 1)):(40 * i),]$Q), nrow = 40, ncol = 2, byrow = FALSE)
  cats <- c(rep(2, 40))
  fx_theta[i, ] <- wlord(probs,cats)

}

fx_theta <- as.data.frame(fx_theta)
fx_theta$weights <- quad_points$weights


fx_theta_n <- fx_theta %>%
  mutate_each(funs(. * weights), starts_with("V"))

fx <- colSums(fx_theta_n[,1:41])
fx_dit <- as.data.frame(matrix(fx, nrow = 41, ncol = 1))
fx_dit$X <- c(0:40)

xm <- weighted.mean(fx_dit$X, fx_dit$V1)
VarObsXModel <- sum(fx_dit$V1 * (fx_dit$X - xm)^2)

TestReliabilityObsModel <- 1 - VarError / VarObsXModel
TestReliabilityObsModel


### test reliability: observed score variance, data based  ----------------------------------------

fx_rawdata <- rawdata
fx_rawdata$rs <- rowSums(fx_rawdata )
VarObsXData <- var(fx_rawdata $rs)


### calculation formulas -----------------------

TestReliabilityTrueModel <- VarTrue / VarObsXModel
TestReliabilityTrueData <- VarTrue / VarObsXData

TestReliabilityObsXModel <- 1 - VarError / VarObsXModel
TestReliabilityObsXData <- 1 - VarError / VarObsXData


### summary ----------------------------------------

MarginalReliability

TestReliabilityTrueModel
TestReliabilityObsXModel

TestReliabilityTrueData
TestReliabilityObsXData

}




