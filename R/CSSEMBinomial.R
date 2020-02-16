#' @title CSSEM Binomial Method
#'
#' @description
#' A function to calculate Conditional Standard Error of Measurement for Scale Scores in CTT
#'
#' @param numOfItem a numeric number indicating number of items
#' @param convTable a data frame or matrix containing conversion table of raw score to scale score: variable name
#'                  should be rawScore & roundedSS
#'
#' @return a data frame containing CSSEM using Binomial Method
#'
#' @author {Huan Liu, University of Iowa, \email{huan-liu-1@@uiowa.edu}}
#'
#' @export


CSSEMBinomial <- function(numOfItem, convTable){

  # create matrix to store csem
  cssemDat <- matrix(nrow = numOfItem + 1, ncol = 1)

  # for loop to calculate prob for each true score
  for (i in 0:numOfItem){

    # take raw score as true score
    pi = i/numOfItem

    # create data frame for raw score
    binoDat <- as.data.frame(c(0:numOfItem))
    names(binoDat) <- "raw"

    # calculate prob using binomial
    binoDat$prob <- choose(numOfItem, binoDat$raw) * (pi)^binoDat$raw * (1 - pi)^(numOfItem-binoDat$raw)

    # merge prob with conversion table
    binoDat <- cbind(binoDat, convTable)

    # calculate prob*X and prob*X^2
    binoDat$ssprob <- binoDat$roundedSS * binoDat$prob
    binoDat$ss2prob <- binoDat$roundedSS^2 * binoDat$prob

    # calculate and store Lord's csem for this true score
    seSS <- sqrt(numOfItem/(numOfItem-1)) * sqrt(sum(binoDat$ss2prob) - sum(binoDat$ssprob)^2)
    cssemDat[i+1, 1] <- seSS

  }

  # return raw score and Binomial csem
  return(list("roundedSS" = convTable$roundedSS, "cssemBinomial" = as.numeric(cssemDat)))

}


