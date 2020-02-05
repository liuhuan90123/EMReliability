#' @export


MarginalRelEAP <- function(itemPara){
  if (missing(itemPara)) {
    stop("You should have itemPara!!")
  }
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

  # sum information by theta
  itemParaInfo <- aggregate(itemParaRep$info, by=list(Category=itemParaRep$theta), FUN=sum)
  names(itemParaInfo) <- c("theta", "infoSum")

  # add weights for each theta
  itemParaInfo$weights <- quadPoints$weights

  ## EAP method

  # inverse of information
  itemParaInfo$infoSumInv <- 1 / (itemParaInfo$infoSum + 1) # cobtributed by population

  marginalRelEAP <- 1 - sum(itemParaInfo$infoSumInv * itemParaInfo$weights)
  marginalRelEAP

}







