#' @title CSSEM Kolen's Method
#'
#' @description
#' A function to calculate CSEM for Scale Scores in IRT using Kolen's method
#' True scale score
#'
#' @param numOfItem a numeric number indicating number of items
#' @param convTable a data frame or matrix containing conversion table of raw score to scale score
#'
#' @return a data frame containing CSSEM using Kolen's Method
#'
#' @author {Huan Liu, University of Iowa, \email{huan-liu-1@@uiowa.edu}}
#'
#' @export


CSSEMKolen <- function(itemPara, convTable){

  # item parameters on the 1.702 metric
  names(itemPara) <- c("b", "a")

  # number of quadrature
  numOfQuad <- 41

  # number of Items
  numOfItem <- nrow(itemPara)

  # weights and nodes
  quadPoints <- NormalQuadraPoints(numOfQuad)

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

    probs <- matrix(c(itemParaRep[(1 + numOfItem * (i - 1)):(numOfItem * i),]$P),
                    nrow = numOfItem, ncol = 1, byrow = FALSE)

    fxTheta[i, ] <- LordWingersky(probs)$probability

  }

  # reverse column sequence
  fxTheta <- fxTheta[, c(ncol(fxTheta):1)]

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
  cssemKolen <- matrix(NA, nrow = 41, ncol = 1)

  for (i in 1:41){

    cssemKolen[i, 1] <- sqrt(sum((fxThetaTSS[c(1:41),42] - fxThetaTSS[42, i])^2 * fxThetaTSS[c(1:41),i]))

  }

  # true scale score
  trueSS <- colSums(fxThetaTSS)[1:41]

  return(list("trueScaleScore" = as.numeric(trueSS), "cssemKolen" = as.numeric(cssemKolen)))

}





