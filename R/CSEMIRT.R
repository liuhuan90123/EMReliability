#' @title CSEM of IRT 2PL: EAP&MLE
#'
#' @description
#' A function to calculate CSEM of 2PL IRT with EAP&MLE
#'
#' @param theta a vector, matrix, or data frame containing theta values
#' @param itemPara a text file with parameters of sequence b and a, a is on the 1.702 metric
#' @param estType estimation method, MLE or EAP
#'
#' @return a data frame containing IRT EAP&MLE CSEM
#'
#' @author {Huan Liu, University of Iowa, \email{huan-liu-1@@uiowa.edu}}
#'
#' @export


CSEMIRT <- function(theta, itemPara, estType){

  # return info for each theta
  if (estType == "MLE"){

    # calculate info
    itemParaInfo <- as.data.frame(Info(theta, itemPara, "MLE"))

    # calculate CSEM for each theta
    itemParaInfo$csemMLE <- sqrt(1/itemParaInfo$infoMLE)

    # return csem
    return(list("theta" = theta, "csemMLE" = itemParaInfo$csemMLE))


  }else if (estType == "EAP"){

    # calculate info
    itemParaInfo <- as.data.frame(Info(theta, itemPara, "EAP"))

    # calculate CSEM for each theta
    itemParaInfo$csemEAP <- sqrt(1/itemParaInfo$infoEAP)

    # return csem
    return(list("theta" = theta, "csemEAP" = itemParaInfo$csemEAP))

  }else{

    warning("csemIRT function only supports MLE and EAP estimation method!")

  }
}
