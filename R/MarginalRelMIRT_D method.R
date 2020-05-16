
library(mvtnorm)

numOfFactors <- 3

# item parameters and correlation

# Form A
# itemPara_SS <- read.table("TestData/SpanishLit_prm_A_SS.txt")[,c(7:10)]
#
# names(itemPara_SS) <- c("c","a1","a2","a3")
# itemPara_SS$a <- c(itemPara_SS$a1[1:13], itemPara_SS$a2[14:25], itemPara_SS$a3[26:31])
# itemPara_SS[,"b"] <- -itemPara_SS[,"c"]/itemPara_SS[,"a"]
# itemPara_SS[,"a"] <- itemPara_SS[,"a"]/1.702
#
# itemPara_SS$a1[1:13] <- itemPara_SS$a[1:13]
# itemPara_SS$a2[14:25] <- itemPara_SS$a[14:25]
# itemPara_SS$a3[26:31] <- itemPara_SS$a[26:31]
#
# # read correlations
# cormat <- matrix(c(1, 0.91, 0.70,
#                 0.91, 1, 0.49,
#                 0.70,0.49,1), nrow = 3)
#
# corvec <- c(0.91, # 1&2
#          0.70, # 1&3
#          0.49) # 2&3


# Form B
itemPara_SS <- read.table("TestData/SpanishLit_prm_B_SS.txt")[,c(7:10)]

names(itemPara_SS) <- c("c","a1","a2","a3")
itemPara_SS$a <- c(itemPara_SS$a1[1:13], itemPara_SS$a2[14:25], itemPara_SS$a3[26:31])
itemPara_SS[,"b"] <- -itemPara_SS[,"c"]/itemPara_SS[,"a"]
itemPara_SS[,"a"] <- itemPara_SS[,"a"]/1.702

itemPara_SS$a1[1:13] <- itemPara_SS$a[1:13]
itemPara_SS$a2[14:25] <- itemPara_SS$a[14:25]
itemPara_SS$a3[26:31] <- itemPara_SS$a[26:31]

# read correlations
cormat <- matrix(c(1, 0.97, 0.56,
                0.97, 1, 0.48,
                0.56,0.48,1), nrow = 3)


corvec <- c(0.97, # 1&2
            0.56, # 1&3
            0.48) # 2&3



# quadrature points and weights for 3 factors

# num of quadratures
numOfQuad <- 15

# set nodes ranging from -5 to 5
nodes <- seq(-5, 5, length.out = numOfQuad)

# expand nodes with 3 dimentions
nodes3 <- as.matrix(expand.grid(nodes,nodes,nodes))


# weight
weightsUnwtd <- dmvn(nodes3, c(0,0,0), cormat, log=FALSE) # 41^3
nodes3 <- as.data.frame(nodes3)
nodes3$weightsWtd <- weightsUnwtd / sum(weightsUnwtd)


names(nodes3) <- c("theta1", "theta2", 'theta3',"weightsWtd")

### weights ---
# newWts1 <- as.numeric(unlist(aggregate(nodes3$weightsWtd, by=list(Category=nodes3$Var1), FUN=sum)["x"]))
# newWts2 <- as.numeric(unlist(aggregate(nodes3$weightsWtd, by=list(Category=nodes3$Var2), FUN=sum)["x"]))
# newWts3 <- as.numeric(unlist(aggregate(nodes3$weightsWtd, by=list(Category=nodes3$Var3), FUN=sum)["x"]))


# CSEM
nodes3$se1 <- apply(as.data.frame(nodes3[,"theta1"]), 1, CSEMIRT, itemPara_SS[1:13, c("b", "a")], "EAP")
nodes3$se2 <- apply(as.data.frame(nodes3[,"theta2"]), 1, CSEMIRT, itemPara_SS[14:25, c("b", "a")], "EAP")
nodes3$se3 <- apply(as.data.frame(nodes3[,"theta3"]), 1, CSEMIRT, itemPara_SS[26:31, c("b", "a")], "EAP")


nodes3 <- transform( nodes3,
                    varC = se1^2 + se2^2  + se3^2  + 2 *corvec[1] * se1 * se2 + 2 *corvec[2] * se1 * se3 + 2 *corvec[3] * se2 * se3
)

# average of error variance
ErrorVarAvg <- sum(nodes3$varC * nodes3$weightsWtd)
# ErrorVarAvg

# marginal reliability approach
MarginalRelEAP_SS <- ErrorVarAvg /(ErrorVarAvg + numOfFactors) # var(e)/(var(e) + var(theta))

# coefficients
MarginalRelEAP_SS


# v1 <- sum(nodes3$se1^2 * nodes3$weightsWtd)
# v2 <- sum(nodes3$se2^2 * nodes3$weightsWtd)
# v3 <- sum(nodes3$se3^2 * nodes3$weightsWtd)
#
# se1 <- sqrt(v1)
# se2 <- sqrt(v2)
# se3 <- sqrt(v3)
#
# vc <- v1+v2+v3+2 *corvec[1] * se1 * se2 + 2 *corvec[2] * se1 * se3 + 2 *corvec[3] * se2 * se3
# vc
# # Approach 2
# vc/(vc+3)
# #A 0.7079257
# # B 0.7126083






CSEMIRT <- function(theta, itemPara, estType){

  # return info for each theta
  if (estType == "MLE"){

    # calculate info
    itemParaInfo <- as.data.frame(Info(theta, itemPara, "MLE"))

    # calculate CSEM for each theta
    itemParaInfo$csemMLE <- sqrt(1/itemParaInfo$infoMLE)

    # return csem
    # return(list("theta" = theta, "csemMLE" = itemParaInfo$csemMLE))
    return(itemParaInfo$csemMLE)

  }else if (estType == "EAP"){

    # calculate info
    itemParaInfo <- as.data.frame(Info(theta, itemPara, "EAP"))

    # calculate CSEM for each theta
    itemParaInfo$csemEAP <- sqrt(1/itemParaInfo$infoEAP)

    # return csem
    # return(list("theta" = theta, "csemEAP" = itemParaInfo$csemEAP))
    return(itemParaInfo$csemEAP)

  }else{

    warning("csemIRT function only supports MLE and EAP estimation method!")

  }
}


