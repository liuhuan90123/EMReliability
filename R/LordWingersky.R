#' @title Lord Wingersky recursive formula
#'
#' @description
#' A function to calculate raw score distribution conditional on quadrature points
#'
#' @param probs probabilities that a theta value can correctly answer each item
#'
#' @return probability of raw score distribution for the specfic theta
#'
#' @author {Huan Liu, University of Iowa, \email{huan-liu-1@@uiowa.edu}}
#'
#' @export


LordWingersky <- function(probs){

  # check p matrix
  if(!is.matrix(probs))
    probs <- as.matrix(probs)

  # number of items
  numOfItems <- dim(probs)[1]

  # matrix for storing final probs
  obsX <- array(0, c(numOfItems+1, numOfItems))
  qprobs <- 1 - probs
  l <- 1
  obsX[1,l] <- qprobs[1,]
  obsX[2,l] <- probs[1,]
  obsXR <- obsX[1,l]
  obsXR2 <- obsX[2,l]

  for(i in 2:numOfItems){

    for(r in 1:(i + 1)){
      if(r == 1)
        obsXR <- obsXR * qprobs[i,]

      if(r > 1 && r < (i + 1))
        obsX[r, l+1] <- obsX[r, l] * qprobs[i, ] + obsX[r-1, l] * probs[i, ]

      if(r == (i + 1)){
        obsXR2 <- obsXR2 * probs[i,]
        obsX[r, l+1] <- obsXR2
      }
    }
    obsX[1, l+1] <- obsXR
    l <- l + 1
  }

  return(list("ObservedScore" = 0:numOfItems, "probability" = obsX[ ,numOfItems]))

}

