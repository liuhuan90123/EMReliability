#' @title CSSEM IRT Polynomial Method
#'
#' @description
#' A function to calculate CSEM for Scale Scores in IRT using Polynomial Method
#'
#' @param numOfItem a numeric number indicating number of items
#' @param convTable a data frame or matrix containing conversion table of raw score to scale score
#' @param K a numeric number indicating degree of polynomial regression
#' @param estType estimation method, MLE or EAP
#'
#' @return a list containing R square and CSSEM using Polynomial Method
#'
#' @author {Huan Liu, University of Iowa, \email{huan-liu-1@@uiowa.edu}}
#'
#' @export

CSSEMIRTPoly <- function(itemPara, convTable, K, estType){

  # theta
  theta <- convTable$theta


  # return info for each theta
  if (estType == "MLE"){

    # CSEM MLE
    itemParaCSEM <- as.data.frame(CSEMIRT(theta, itemPara, "MLE"))

    # merge data
    itemParaCSEM <- merge(itemParaCSEM, convTable, by = "theta")

    # change name to fit Polynomial Method function
    names(itemParaCSEM) <- c("rawScore", "csem", "roundedSS")

    # call PM function
    cssemPolyMLE <- PolynomialMethod(itemParaCSEM, K)

    # change variable name
    names(cssemPolyMLE$CSSEMPoly)[names(cssemPolyMLE$CSSEMPoly) == 'rawScore'] <- 'theta'
    names(cssemPolyMLE$CSSEMPoly)[names(cssemPolyMLE$CSSEMPoly) == 'csem'] <- 'csemMLE'

    # return results
    return(list("RSquared" = cssemPolyMLE$RSquared, "CSSEMPolyMLE" = cssemPolyMLE$CSSEMPoly))

  }else if (estType == "EAP"){

    # CSEM MLE
    itemParaCSEM <- as.data.frame(CSEMIRT(theta, itemPara, "EAP"))

    # merge data
    itemParaCSEM <- merge(itemParaCSEM, convTable, by = "theta")

    # change name to fit Polynomial Method function
    names(itemParaCSEM) <- c("rawScore", "csem", "roundedSS")

    # call PM function
    cssemPolyEAP <- PolynomialMethod(itemParaCSEM, K)

    # change variable name
    names(cssemPolyEAP$CSSEMPoly)[names(cssemPolyEAP$CSSEMPoly) == 'rawScore'] <- 'theta'
    names(cssemPolyEAP$CSSEMPoly)[names(cssemPolyEAP$CSSEMPoly) == 'csem'] <- 'csemEAP'

    # return results
    return(list("RSquared" = cssemPolyEAP$RSquared, "CSSEMPolyEAP" = cssemPolyEAP$CSSEMPoly))

  }else{

    warning("CSSEMIRTPoly function only supports MLE and EAP estimation method!")

  }

}






