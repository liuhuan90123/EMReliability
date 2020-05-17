#' @title CSEM Strat Feldt Mehtod
#'
#' @description
#' A function to calculate Conditional Standard Error of Measurement for composite score in CTT
#'
#' @param dat a data frame or matrix containing raw data
#' @param strat a vector containing number of item for each strat
#' @return a data frame containing Feldt CSEM
#'
#' @author {Huan Liu, University of Iowa, \email{huan-liu-1@@uiowa.edu}}
#'
#' @export



CSEMStratFeldt <- function(dat, strat){

  # num of strats
  numStrat <- length(strat)

  # create matrix to store csem
  csemDat <- as.data.frame(rowSums(dat), nrow = nrow(dat), ncol = 1)

  # for loop to calculate prob for each true score
  for (i in 1:nrow(dat)){

    varK <- c()
    # error variance for strat #1
    for (k in 1:1) {
      varK[k] <- sum(dat[i, 1:strat[1]])*(strat[1]-sum(dat[i, 1:strat[1]]))/(strat[1]-1)
    }

    # error variance for strat #2 to #numStrat
    for (k in 2:numStrat) {
      varK[k] <- sum(dat[i,(sum(strat[1:k-1])+1):sum(strat[1:k])])*(strat[k]-sum(dat[i,(sum(strat[1:k-1])+1):sum(strat[1:k])]))/(strat[k]-1)
    }

    # CSEM Strat Feldt for each examinee
    csemStratFeldt <- sqrt(sum(varK))

    csemDat[i,"csemStratFeldt"] <- csemStratFeldt
  }

  # return raw score and Feldt strat csem
  names(csemDat) <- c("rawScore", "csemStratFeldt")
  return(csemDat)


}



