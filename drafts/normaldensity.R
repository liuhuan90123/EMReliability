### normal density function

options(scipen = 999)

# function
normalDensity <- function(x, mean, sd){

  1 / sqrt(2 * pi * sd^2) * exp(-0.5 * ((x-mean)/sd)^2)

}


# Kole&Lee's examples
normalWtKolen <- read.table("TestData/NW_Kolen.txt")
names(normalWtKolen) <- c("nodes", "weightsKolen")
normalWtLee <- read.table("TestData/NW_Lee.txt")
names(normalWtLee) <- c("nodes", "weightsLee")

# merge data
normalWtCom <- merge(normalWtKolen, normalWtLee, by = "nodes")

# compare K&L
normalWtCom$diffKL <- normalWtCom$weightsKolen - normalWtCom$weightsLee
sum(normalWtCom$diffKL)

# calculate Huan
normalWtCom$weightsHuanUnwtd <- apply(normalWtCom["nodes"], 1, FUN=function(x) dnorm(x))

normalWtCom$weightsHuanWtd <- normalWtCom$weightsHuanUnwtd/sum(normalWtCom$weightsHuanUnwtd)

sum(normalWtCom$weightsHuanWtd)

# compare Huan to K&L
normalWtCom$diffKH <- normalWtCom$weightsKolen - normalWtCom$weightsHuanWtd
sum(normalWtCom$diffKH)
normalWtCom$diffLH <- normalWtCom$weightsLee - normalWtCom$weightsHuanWtd
sum(normalWtCom$diffLH)


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



