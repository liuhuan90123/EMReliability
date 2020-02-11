#' @title Marginal Reliability of IRT 2PL: MLE
#'
#' @description
#' A function to calculate marginal reliability of 2PL IRT with MLE
#'
#' @param itemPara a text file with parameters of sequence b and a, a is on the 1.7 metric
#' @return a reliability number
#'
#' @author {Huan Liu, University of Iowa, \email{huan-liu-1@@uiowa.edu}}
#'
#' @export

MarginalRelMLE <- function(itemPara){

  # transform item parameters to the 1.702 metric
  names(itemPara) <- c("b", "a")
  itemPara[,"a"] <- itemPara[,"a"]/1.702

  # weights and nodes
  quadPoints <- gauss.quad.prob(41, dist = "normal", mu = 0, sigma = 1)

  # replicate item parameter and theta
  itemParaRep <-itemPara[rep(seq_len(nrow(itemPara)), each = 41),]
  itemParaRep$theta <- rep(quadPoints$nodes, each = 1, length.out = 41*nrow(itemPara))

  # calculate information by theta
  itemParaRep <- within(itemParaRep, {
    P = 0 + (1 - 0) / (1 + exp(-1.702 * a * (theta - b)))
    Q = 1 - P
    PQ = P * Q
    info = 1.702**2 * a**2 * P * Q
  })

  # sum information by theta
  itemParaInfo <- aggregate(itemParaRep$info, by=list(Category=itemParaRep$theta), FUN=sum)
  names(itemParaInfo) <- c("theta", "infoSum")

  # add weights for each theta
  itemParaInfo$weights <- quadPoints$weights

  ## MLE method

  # weighted information
  itemParaInfo$infoWeighted <- itemParaInfo$weights * itemParaInfo$infoSum

  # marginal reliability MLE
  marginalRelMLE <- sum(itemParaInfo$infoWeighted)/(sum(itemParaInfo$infoWeighted) + 1)
  marginalRelMLE

}



