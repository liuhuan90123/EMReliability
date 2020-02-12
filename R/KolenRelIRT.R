#' @title IRT Reliability of Kolen's Method
#'
#' @description
#' A function to calculate IRT Reliability of Kolen's Method
#'
#' @param itemPara a text file with parameters of sequence b and a, a is on the 1.702 metric
#' @param convTable matrix or data frame contains raw score to scale score
#' @return a reliability number
#'
#' @author {Huan Liu, University of Iowa, \email{huan-liu-1@@uiowa.edu}}
#'
#' @export

# read item parameters from txt file
itemPara <- read.table("TestData/ItemParaFormX.txt")

# read conversion table
convTable <- read.csv("TestData/ConversionTableFormX.csv")
convTable$roundedSS <- round(convTable$unroundedSS)

library(classify)


CSSEMKolen <- function(itemPara, convTable){

  # transform item parameters to the 1.702 metric
  names(itemPara) <- c("b", "a")
  itemPara[,"a"] <- itemPara[,"a"]/1.702

  # number of quadrature
  numOfQuad <- 41

  # number of Items
  numOfItem <- nrow(itemPara)

  # weights and nodes
  quadPoints <- NormalQuadraPoints(41)

  # replicate item parameter and theta
  itemParaRep <-itemPara[rep(seq_len(numOfItem), each = numOfQuad),]
  itemParaRep$theta <- rep(quadPoints$nodes, each = 1, length.out = numOfQuad*numOfItem)

  # calculate information by theta
  itemParaRep <- within(itemParaRep, {
    P = 0 + (1 - 0) / (1 + exp(-1.702 * a * (theta - b)))
    Q = 1 - P
    PQ = P * Q
    info = 1.702**2 * a**2 * P * Q
  })

  # reorder matrix by theta
  itemParaRep <- itemParaRep[order(itemParaRep$theta),]

  # create matrix to store f(x|theta)
  fxTheta <- matrix(NA, nrow = numOfQuad, ncol = numOfItem + 1)

  # for loop to calculate fxTheta
  for (i in 1:numOfQuad){

    probs <- matrix(c(itemParaRep[(1 + numOfItem * (i - 1)):(numOfItem * i),]$P,
                      itemParaRep[(1 + numOfItem * (i - 1)):(numOfItem * i),]$Q),
                    nrow = numOfItem, ncol = 2, byrow = FALSE)

    cats <- c(rep(2, numOfItem))

    fxTheta[i, ] <- wlord(probs,cats)

  }

  # transform to data frame
  fxTheta <- as.data.frame(fxTheta)

  # transform data frame fxTheta
  fxThetaT <- as.data.frame(t(fxTheta))

  # reverse SS
  fxThetaT$SS <- rev(convTable$roundedSS)

  # true scale score
  fxThetaTSS <- as.data.frame(apply(fxThetaT[c(1:41)], 2, function(x) x * fxThetaT$SS))
  fxThetaTSS$SS <- rev(convTable$roundedSS)

  # merge data
  fxThetaTSS <- rbind(fxThetaT, colSums(fxThetaTSS))

  # CSSEM condtional on theta
  cssemKolen <- matrix(NA,nrow = 41, ncol = 1)

  for (i in 1:41){

    cssemKolen[i, 1] <- sqrt(sum((fxThetaTSS[c(1:41),42] - fxThetaTSS[42, i])^2 * fxThetaTSS[c(1:41),i]))

  }

  # error variance: avarage CSSEM across theta distribution
  errorVarKolen <- sum(cssemKolen^2 * quadPoints$weights)

  # variance of scale score
  fxPrXi <- as.data.frame(apply(fxTheta[c(1:41)], 2, function(x) x * quadPoints$weights))

  # mean of scale score
  meanSS <- sum(rev(convTable$roundedSS) * colSums(fxPrXi))

  # variance of scale score
  SSVarKolen <- sum((rev(convTable$roundedSS) - meanSS)^2 * colSums(fxPrXi))

  # reliability
  kolenRelIRT <- 1 - errorVarKolen / SSVarKolen
  kolenRelIRT

}
