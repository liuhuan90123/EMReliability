g <- -1.8826439
c <- -2.0107269
a <- 1.6477878

a1.7 <- a/1.7
a1.7

b1.7 <- -c/a
b1.7

g1.7 <- 1/(1+exp(-g))
g1.7

library(LaplacesDemon)
cormat <- matrix(c(1, 0.92,
                  0.92, 1), nrow = 2)

# set nodes and weights
numOfQuad <- 41
nodes <- seq(-4, 4, length.out = numOfQuad)
nodesM <- as.matrix(expand.grid(nodes,nodes))
weightsUnwtd <- dmvn(nodesM, c(0,0), cormat, log=FALSE)
nodesM <- as.data.frame(nodesM)
nodesM$weightsWtd <- weightsUnwtd / sum(weightsUnwtd)




nodesM <- as.matrix(expand.grid(nodes,nodes,nodes))
weightsUnwtd <- dmvn(nodesM, c(0,0,0), cormat, log=FALSE)
nodesM <- as.data.frame(nodesM)
nodesM$weightsWtd <- weightsUnwtd / sum(weightsUnwtd)



options(scipen=999)
# Form A
cormat_A <- matrix(c(1, 0.9067069, 0.6994119,
                     0.9067069, 1, 0.4891160,
                     0.6994119,0.4891160,1), nrow = 3)
# set nodes and weights
nodes <- seq(-4, 4, length.out = 41)
nodesM <- as.matrix(expand.grid(nodes,nodes,nodes))
weightsUnwtd <- dmvn(nodesM, c(0,0,0), cormat_A, log=FALSE)
nodesM <- as.data.frame(nodesM)
nodesM$weightsWtd <- weightsUnwtd # / sum(weightsUnwtd)

nodesM <- nodesM[,c("Var3", "Var2", "Var1", "weightsWtd")]

nodesM$Var3 <- format(nodesM$Var3, scientific = FALSE)
nodesM$Var2 <- format(nodesM$Var2, scientific = FALSE)
nodesM$Var1 <- format(nodesM$Var1, scientific = FALSE)
nodesM$weightsWtd <- round(nodesM$weightsWtd, 10)

# options(scipen=10)

write.table(nodesM, "thetathreeformA.txt", col.names = F, row.names = F, quote = F)


# item parameters
itemPara_SS_A <- itemPara_SS_A[,c("a", "b")]
itemPara_SS_A$c <- 0



# Form B
cormat_B <- matrix(c(1, 0.9722234, 0.5602197,
                     0.9722234, 1, 0.4795721,
                     0.5602197,0.4795721,1), nrow = 3)
# set nodes and weights
nodes <- seq(-4, 4, length.out = 41)
nodesM <- as.matrix(expand.grid(nodes,nodes,nodes))
weightsUnwtd <- dmvn(nodesM, c(0,0,0), cormat_B, log=FALSE)
nodesM <- as.data.frame(nodesM)
nodesM$weightsWtd <- weightsUnwtd# / sum(weightsUnwtd)



nodesM <- nodesM[,c("Var3", "Var2", "Var1", "weightsWtd")]

nodesM$Var3 <- format(nodesM$Var3, scientific = FALSE)
nodesM$Var2 <- format(nodesM$Var2, scientific = FALSE)
nodesM$Var1 <- format(nodesM$Var1, scientific = FALSE)
# nodesM$weightsWtd <- format(nodesM$weightsWtd, scientific = FALSE)
nodesM$weightsWtd <- round(nodesM$weightsWtd, 10)

options(scipen=5)

write.table(nodesM, "thetathreeformB.txt", col.names = F, row.names = F, quote = F)


# item parameters
itemPara_SS_B <- itemPara_SS_B[,c("a", "b")]
itemPara_SS_B$c <- 0
itemPara_SS_B

library(MASS)
multhetaA <- mvrnorm(n = 10000, c(0,0,0), cormat_A, tol = 1e-6, empirical = FALSE)
write.table(multhetaA, "multhetaA.txt", col.names = F, row.names = F, quote = F)


multhetaB <- mvrnorm(n = 10000, c(0,0,0), cormat_B, tol = 1e-6, empirical = FALSE)
write.table(multhetaB, "multhetaB.txt", col.names = F, row.names = F, quote = F)
