#' @title TestRelIRT is a function
#'
#' @description
#' A function to calculate marginal reliability of 2PL IRT with MLE
#'
#' @param itemPara a text file
#' @return a reliability number
#'
#' @author {Huan Liu, University of Iowa, \email{huan-liu-1@@uiowa.edu}}
#' @export

TestRelIRT <- function(itemPara){

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

  # num of quadratures
  numQuad <- 41

  # define matrix of marginal distribution of theta
  fxTheta <- matrix(NA, nrow = numQuad, ncol = 41) # 41 num of quadratures, 41 num of raw sxores

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

  # calculate weighted distribution
  fxThetaWeighted <- apply(fxTheta[,1:41], 2, function(x) x * fxTheta[,"weights"])

  # sum weighted distribution
  fxDist <- as.data.frame(matrix(colSums(fxThetaWeighted[,1:41]), nrow = 41, ncol = 1))
  fxDist$X <- c(nrow(itemPara):0)
  names(fxDist) <- c("wts", "X")

  # weighted mean of Obs X
  weightedMean <- sum(fxDist$X * fxDist$wts)/sum(fxDist$wts)

  # variance of Obs X
  varianceObsX <- sum(fxDist$wts * (fxDist$X - weightedMean)^2)

  # test reliability formula

  TestRelIRT <- 1 - varianceError / varianceObsX # varianceTrue / varianceObsX
  TestRelIRT
}



