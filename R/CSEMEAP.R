#' @title CSEM of IRT 2PL: EAP
#'
#' @description
#' A function to calculate CSEM of 2PL IRT with EAP
#'
#' @param itemPara a text file with parameters of sequence b and a, a is on the 1.702 metric
#' @param theta a vector, matrix, or data frame containing theta values
#' @return a data frame containing IRT EAP CSEM
#'
#' @author {Huan Liu, University of Iowa, \email{huan-liu-1@@uiowa.edu}}
#'
#' @export

CSEMEAP<- function(itemPara, theta){

  # transform item parameters to the 1.702 metric
  names(itemPara) <- c("b", "a")
  itemPara[,"a"] <- itemPara[,"a"]/1.702

  # weights and nodes
  # quadPoints <- gauss.quad.prob(41, dist = "normal", mu = 0, sigma = 1)

  # replicate item parameter and theta
  itemParaRep <-itemPara[rep(seq_len(nrow(itemPara)), each = length(theta)),]
  itemParaRep$theta <- rep(theta, each = 1, length.out = length(theta)*nrow(itemPara))

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

  # calculate CSEM for each theta
  itemParaInfo$csemEAP <- sqrt(1/(itemParaInfo$infoSum+1))

  itemParaInfo

}

