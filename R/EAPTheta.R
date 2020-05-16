#' @title EAPTheta
#'
#' @description
#' A function to calculate examinee ability using EAP with item parameters known
#'
#' @param itemPara a matrix or data frame with parameters of sequence b, a and c on the 1.702 metric
#' @param rawData a matrix or data frame of raw responses
#' @return a data frame for EAP theta and SD
#'
#' @author {Huan Liu, University of Iowa, \email{huan-liu-1@@uiowa.edu}}
#' @export



EAPTheta <- function(itemPara, rawData){

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

  # number of quadrature points
  numOfQua <- 41

  # theta q
  thetaq <- NormalQuadraPoints(numOfQua)$nodes

  # replicate item parameter and theta q
  itemParaRep <- itemPara[rep(seq_len(nrow(itemPara)), each = numOfQua),]
  itemParaRep$thetaq <- rep(thetaq, each = 1, length.out = numOfQua*nrow(itemPara))

  # calculate information by theta q
  itemParaRep <- within(itemParaRep, {

    P = c + (1 - c) / (1 + exp(-1.702 * a * (thetaq - b)))
    Q = 1 - P
    PQ = P * Q
    info = 1.702**2 * a**2 * (Q/P) * (P-c)^2 / (1-c)^2
  })

  # order by theta q
  itemParaRep <- itemParaRep[order(itemParaRep$thetaq), ]

  # for loop to calculate theta and SE using EAP
  summaryEAP <- matrix(nrow = nrow(rawData), ncol = 2)
  for (i in 1:nrow(rawData)){

    # calculate likelihood
    itemParaRep$resp <- rep(rawData[i,], each = 1, length.out = numOfQua*nrow(itemPara))
    itemParaRep$resp <- as.numeric(itemParaRep$resp)
    itemParaRep$respRev <- 1 - itemParaRep$resp
    itemParaRep$T <- itemParaRep$P * itemParaRep$resp + itemParaRep$Q * itemParaRep$respRev

    # sum T by thetaq
    itemParaT <- aggregate(itemParaRep$T, by=list(Category=itemParaRep$thetaq), FUN=prod)
    names(itemParaT) <- c("thetaq", "T")
    itemParaT$wt <- NormalQuadraPoints(numOfQua)$weights

    # numerator and denominator for fomula
    nume <- sum(itemParaT$thetaq * itemParaT$T * itemParaT$wt)
    deno <- sum(itemParaT$T * itemParaT$wt)

    # theta
    thetaEAP <- nume/deno
    summaryEAP[i, 1] <- thetaEAP

    # SD
    SDEAP <- sqrt(sum((itemParaT$thetaq - thetaEAP)^2 * itemParaT$T * itemParaT$wt) / deno)
    summaryEAP[i, 2] <- SDEAP

  }
  # return coefficient
  return(summaryEAP)
}




