#' @title CSEM Lord Mehtod
#'
#' @description
#' A function to calculate Conditional Standard Error of Measurement in CTT
#'
#' @param numOfItem a numeric number indicating number of items
#' @return a data frame containing Lord CSEM
#'
#' @author {Huan Liu, University of Iowa, \email{huan-liu-1@@uiowa.edu}}
#'
#' @export


CSEMLord <- function(numOfItem){

  # create matrix to store csem
  csemDat <- matrix(nrow = numOfItem + 1, ncol = 1)


  # for loop to calculate prob for each true score
  for (i in 0:numOfItem){

    # take raw score as true score
    pi = i/numOfItem

    # create data frame for raw score
    binoDat <- as.data.frame(c(0:numOfItem))
    names(binoDat) <- "raw"

    # calculate prob using binomial
    binoDat$prob <- choose(numOfItem, binoDat$raw) * (pi)^binoDat$raw * (1 - pi)^(numOfItem-binoDat$raw)

    # calculate prob*X and prob*X^2
    binoDat$xprob <- binoDat$raw * binoDat$prob
    binoDat$x2prob <- binoDat$raw^2 * binoDat$prob

    # calculate and store Lord's csem for this true score
    seX <- sqrt(numOfItem/(numOfItem-1)) * sqrt(sum(binoDat$x2prob) - sum(binoDat$xprob)^2)
    csemDat[i+1, 1] <- seX

  }

  # format csem data frame
  csemDat <- as.data.frame(csemDat)
  names(csemDat) <- c("csemLord")
  csemDat$rawScore <- c(0:numOfItem)

  # return Lord csem
  return(csemDat)

}

