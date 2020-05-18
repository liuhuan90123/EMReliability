#' @title Marginal Reliability of IRT with MLE or EAP
#'
#' @description
#' A function to calculate marginal reliability of unidimensional IRT with MLE or EAP
#'
#' @param itemPara a data frame or matrix with parameters of sequence b, a and c on the 1.702 metric
#' @param estType ability estimation method, MLE or EAP
#' @return Marginal Reliability of IRT
#'
#' @author {Huan Liu, University of Iowa, \email{huan-liu-1@@uiowa.edu}}
#'
#' @export

MarginalRelIRT <- function(itemPara, estType){

  # weights and nodes
  quadPoints <- NormalQuadraPoints(41)

  if (estType == "MLE"){

    # calculate info
    itemParaInfo <- as.data.frame(Info(quadPoints$nodes, itemPara, "MLE"))

    # add weights for each theta
    itemParaInfo$weights <- quadPoints$weights

    # weighted information
    itemParaInfo$infoWeighted <- itemParaInfo$weights * itemParaInfo$infoMLE

    # marginal reliability MLE
    marginalRelMLE <- sum(itemParaInfo$infoWeighted)/(sum(itemParaInfo$infoWeighted) + 1)

    # return coefficient
    return(marginalRelMLE)


  }else if (estType == "EAP"){

    # calculate info
    itemParaInfo <- as.data.frame(Info(quadPoints$nodes, itemPara, "EAP"))

    # add weights for each theta
    itemParaInfo$weights <- quadPoints$weights

    # inverse of information
    itemParaInfo$infoEAPInv <- 1 / itemParaInfo$infoEAP

    # marginal reliability EAP
    marginalRelEAP <- 1 - sum(itemParaInfo$infoEAPInv * itemParaInfo$weights)

    # return coefficient
    return(marginalRelEAP)

  }else{

    warning("MarginalRelMLE function only supports MLE and EAP estimation method!")

  }

}



