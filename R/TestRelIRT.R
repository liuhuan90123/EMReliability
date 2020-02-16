#' @title TestRelIRT
#'
#' @description
#' A function to calculate test reliability of IRT
#'
#' @param itemPara a text file with parameters of sequence b and a, a is on the 1.702 metric
#' @return a reliability number
#'
#' @author {Huan Liu, University of Iowa, \email{huan-liu-1@@uiowa.edu}}
#' @export

TestRelIRT <- function(itemPara){

  # transform item parameters to the logistic metric
  names(itemPara) <- c("b", "a")

  # num of items
  numOfItem <- nrow(itemPara)

  # num of quadratures
  numOfQuad <- 41

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

  ## true score variance

  # sum probability by theta
  itemParaAggr <- aggregate(itemParaRep, by=list(Category=itemParaRep$theta), FUN=sum)

  # add weights for each theta
  itemParaAggr$weights <- quadPoints$weights

  # calculate true score variance
  varianceTrue <- sum((itemParaAggr$P)^2 * itemParaAggr$weights) - (sum(itemParaAggr$P * itemParaAggr$weights))^2


  ## error variance

  # claculate error variance
  varianceError <- sum(itemParaAggr$PQ * itemParaAggr$weights)


  ## observed score variance

  # order by theta
  itemParaRep <- itemParaRep[order(itemParaRep$theta),]

  # define matrix of marginal distribution of theta
  fxTheta <- matrix(NA, nrow = numOfQuad, ncol = numOfItem + 1) # 41 num of quadratures, 41 num of raw sxores

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

  # add quadrature weights
  fxTheta$weights <- quadPoints$weights

  # calculate weighted distribution
  fxThetaWeighted <- apply(fxTheta[,1:41], 2, function(x) x * fxTheta[,"weights"])

  # sum weighted distribution
  fxDist <- as.data.frame(matrix(colSums(fxThetaWeighted[,1:41]), nrow = 41, ncol = 1))
  fxDist$X <- c(numOfItem:0)
  names(fxDist) <- c("wts", "X")

  # weighted mean of Obs X
  weightedMean <- sum(fxDist$X * fxDist$wts)/sum(fxDist$wts)

  # variance of Obs X
  varianceObsX <- sum(fxDist$wts * (fxDist$X - weightedMean)^2)

  # test reliability
  TestRelIRT <- 1 - varianceError / varianceObsX # varianceTrue / varianceObsX

  # return coefficient
  return(TestRelIRT)
}



