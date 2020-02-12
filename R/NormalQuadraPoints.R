#' @title normal quadrature points and normalized weights
#'
#' @description
#' A function to generate quadrature points and normalized weights
#'
#' @param n number of quadrature points evenly distributed from -5 to +5
#'
#' @return a list containing nodes and normalized weights
#'
#' @author {Huan Liu, University of Iowa, \email{huan-liu-1@@uiowa.edu}}
#'
#' @export


NormalQuadraPoints <- function(n){

  # set nodes ranging from -5 to 5
  nodes <- seq(-5, 5, length.out = n)

  # unnormalized weights
  weightsUnwtd <- sapply(nodes, FUN = function(x) dnorm(x))

  # normalized weightes
  weightsWtd <- weightsUnwtd / sum(weightsUnwtd)

  # return nodes and normalized weights
  return(list(nodes, weightsWtd))

}

NormalQuadraPoints(41)
