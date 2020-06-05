
itemPara <- itemPara_A
convTable <- convTable_A_Poly
K <- 10



# theta and weights
theta <- NormalQuadraPoints(41)$nodes
# theta <- rnorm(10000)

# CSEM
itemParaCSEM <- as.data.frame(CSEMIRT(theta, itemPara, "MLE"))
# weight
itemParaCSEM$wt <- NormalQuadraPoints(41)$weights

rSquaredDat <- as.data.frame(matrix(nrow = K, ncol = 1))
relDat <- as.data.frame(matrix(nrow = K, ncol = 1))

# for loop to iterate different k
for (k in 1:K){
  k <- 4

  # fit model with k
  modelK <- lm(roundedSS ~ poly(theta, k, raw=TRUE), convTable)

  # extract regression coefficients
  regCoef <- summary(modelK)$coefficients[, 1]

  # extract r square coefficient
  rSquaredDat[k, 1]<- summary(modelK)$r.squared

  # check whether regression coefficient of highest order is missing
  if(is.na(regCoef[k+1])){

    message(paste("The maximum k accepted is", k-1, sep = " "))

    break

  }
  # (Intercept) poly(theta, k, raw = TRUE)1 poly(theta, k, raw = TRUE)2 poly(theta, k, raw = TRUE)3
  # 118.263733790                 3.909028972                -0.360085771                -0.036053623
  # poly(theta, k, raw = TRUE)4
  # 0.009014156
  # x <- -3.75
  # 0.009014156 * x^4 + -0.036053623*x^3 +  -0.360085771* x^2 + 3.909028972*x + 118.263733790
  # 4 * 0.009014156 * x^3 + 3* -0.036053623*x^2 +  2*-0.360085771* x + 3.909028972
  # sqrt(0.630908^2 * 3.7522059)

  # calculate SS: from 1 to K
  itemParaCSEM$SS <- 0
  i <- k

  while(i > 0){

    itemParaCSEM$SS <- itemParaCSEM$SS +  regCoef[i+1] * itemParaCSEM$theta^(i)
    i <- i-1

  }

  itemParaCSEM$SS <- itemParaCSEM$SS + regCoef[i+1]
  itemParaCSEM$SS <- round(itemParaCSEM$SS)


  # calculate transformation coefficients fx: from 1 to K
  itemParaCSEM$fx <- 0
  i <- k

  while(i > 1){

    itemParaCSEM$fx <- itemParaCSEM$fx +  regCoef[i+1] * (i * itemParaCSEM$theta^(i-1))
    i <- i-1

  }

  itemParaCSEM$fx <- itemParaCSEM$fx + regCoef[i+1]


  # calculate cssem using polynomial method
  itemParaCSEM$cssemPoly <- sqrt(itemParaCSEM$fx * (itemParaCSEM$csemMLE)^2)

  # check negative values
  negCount <- sum(itemParaCSEM$fx < 0)

  if(negCount > 0){
    message(paste("Negative transformation coefficient exists when k = ", k,
                  ". The corresponding CSSEM will not be available.", sep = " "))
  }else{

    SSVar <- sum(itemParaCSEM$wt * (itemParaCSEM$SS - weighted.mean(itemParaCSEM$SS, itemParaCSEM$wt))^2)
    errVar <- sum(itemParaCSEM$cssemPoly^2 * itemParaCSEM$wt)
    relDat[k,1] <- 1 - errVar/SSVar
  }
  sd(itemParaCSEM$SS)

  # rename variable with indicator k
  names(itemParaCSEM)[names(itemParaCSEM) == 'fx'] <- paste("fxk", k, sep = "")
  names(itemParaCSEM)[names(itemParaCSEM) == 'SS'] <- paste("SSk", k, sep = "")
  names(itemParaCSEM)[names(itemParaCSEM) == 'cssemPoly'] <- paste("cssemPolyk", k, sep = "")

}

write.xlsx(itemParaCSEM, "Poly_A_k4_rel.xlsx")


