#' @title Marginal Reliability of IRT 2PL: MLE&EAP
#'
#' @description
#' A function to calculate marginal reliability of 2PL IRT with MLE&EAP
#'
#' @param itemPara a text file with parameters of sequence b and a, a is on the 1.702 metric
#' @param estType estimation method, MLE or EAP
#' @return a reliability number
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



