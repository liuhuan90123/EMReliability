#' @title IRT Information for 2PL model
#'
#' @description
#' A function to calculate information for 2PL model
#'
#' @param theta ability theta vector
#' @param itemPara a text file with parameters of sequence b and a and c, a is on the 1.702 metric
#' @param estType estimation method, MLE or EAP
#'
#' @return a data frame containing info
#'
#' @author {Huan Liu, University of Iowa, \email{huan-liu-1@@uiowa.edu}}
#'
#' @export



Info <- function(theta, itemPara, estType){

  # itemPara <- itemPara_A # test
  # theta <- NormalQuadraPoints(41)$nodes

  if (ncol(itemPara) == 3){
    #  item parameters should be on the 1.702 metric
    names(itemPara) <- c("b", "a", "c")
  }

  if (ncol(itemPara) == 2){
    #  item parameters should be on the 1.702 metric
    names(itemPara) <- c("b", "a")
    itemPara$c <- 0
  }

  if (ncol(itemPara) == 1){
    #  item parameters should be on the 1.702 metric
    names(itemPara) <- c("b")
    itemPara$a <- 1
    itemPara$c <- 0
  }

  # replicate item parameter and theta
  itemParaRep <-itemPara[rep(seq_len(nrow(itemPara)), each = length(theta)),]
  itemParaRep$theta <- rep(theta, each = 1, length.out = length(theta)*nrow(itemPara))


  # # calculate information by theta 2PL
  # itemParaRep <- within(itemParaRep, {
  #   P = 0 + (1 - 0) / (1 + exp(-1.702 * a * (theta - b)))
  #   Q = 1 - P
  #   PQ = P * Q
  #   info = 1.702**2 * a**2 * P * Q
  # })

  # calculate information by theta 3PL
  itemParaRep <- within(itemParaRep, {
    P = c + (1 - c) / (1 + exp(-1.702 * a * (theta - b)))
    Q = 1 - P
    PQ = P * Q
    info = 1.702**2 * a**2 * (Q/P) * (P-c)^2 / (1-c)^2
  })

  # sum information by theta
  itemParaInfo <- aggregate(itemParaRep$info, by=list(Category=itemParaRep$theta), FUN=sum)
  names(itemParaInfo) <- c("theta", "infoMLE")

  # calculate info for EAP
  itemParaInfo$infoEAP <- itemParaInfo$infoMLE + 1 # cobtributed by population

  # return info for each theta
  if (estType == "MLE"){

    return(list("theta" = itemParaInfo$theta, "infoMLE" = itemParaInfo$infoMLE))

  }else if (estType == "EAP"){

    return(list("theta" = itemParaInfo$theta, "infoEAP" = itemParaInfo$infoEAP))

  }else{

    warning("Info function only supports MLE and EAP estimation method!")

  }

}

