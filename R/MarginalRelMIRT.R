# Marginal Reliability

#### simple structure MLE&EAP  P method ----------

## MLE ---------------------------------------

# number of factors
numOfFactors <- 3

# Form A
# scoSS_MLE <- read.table("TestData/SpanishLit_sco_A_SS_MLE.txt")[,c(4:9, 10, 12, 15)]
# cor <- c(0.9067069, # 1&2
#          0.6994119, # 1&3
#          0.4891160) # 2&3

# Form B
scoSS_MLE <- read.table("TestData/SpanishLit_sco_B_SS_MLE.txt")[,c(4:9, 10, 12, 15)]
cor <- c(0.9722234, # 1&2
         0.5602197, # 1&3
         0.4795721) # 2&3

# Form Test
# scoSS_MLE <- read.table("TestData/SpanishLit-sco-test.txt")[,c(4:9, 10, 12, 15)]
# cor <- c(0.91, # 1&2
#          0.71, # 1&3
#          0.51) # 2&3


# change variable name
names(scoSS_MLE) <- c("theta1", "theta2", "theta3", "se1", "se2", "se3", "var1", "var2", "var3")

# delete observations with missing values, 99.99 for flexMIRT output
scoSS_MLE[scoSS_MLE == 99.99] <- NA
scoSS_MLE <- na.omit(scoSS_MLE)

# composite error variance
scoSS_MLE <- transform( scoSS_MLE,
                   varC = var1 + var2 + var3,# + 2 *cor[1] * se1 * se2 + 2 *cor[2] * se1 * se3 + 2 *cor[3] * se2 * se3,
                   thetaSum = theta1 + theta2 + theta3
)

var(scoSS_MLE$theta1)
mean(scoSS_MLE$var1)
var(scoSS_MLE$theta2)
mean(scoSS_MLE$var2)
var(scoSS_MLE$theta3)
mean(scoSS_MLE$var3)



### classification -----------------------

cutscore <- -2
theta <- scoSS_MLE$thetaSum
hist(theta)
sem <- sqrt(scoSS_MLE$varC)

os<-theta
nn<-length(os)  # nn, number of examinee
nc <- length(cutscore)  # number of cutscore

if(nn != length(sem)) stop("Ability and se of different length")

esacc<-matrix(NA,length(cutscore), nn, dimnames = list(paste("cut at",round(cutscore,3)), round(os,3)))
escon <-esacc


  j=1 # test
  cuts<-c(-Inf, cutscore[j], Inf)
  categ<-cut(os,cuts,labels=FALSE,right=FALSE) # cut function in r

  for(i in 1:nn) {
    esacc[j,i]<-(pnorm(cuts[categ[i]+1],os[i],sem[i])-pnorm(cuts[categ[i]],os[i],sem[i]))
    escon[j,i]<-((pnorm(cuts[2], os[i],sem[i]) - pnorm(cuts[1],os[i],sem[i]))^2	+ (pnorm(cuts[3], os[i],sem[i]) - pnorm(cuts[2],os[i],sem[i]))^2	 )
  }

  ans<- (list("Marginal" = cbind("Accuracy" = rowMeans(esacc), "Consistency" = rowMeans(escon)), "Conditional" = list("Accuracy" =t(esacc), "Consistency" = t(escon))))

  ans$Marginal


  # >   ans$Marginal
  # Accuracy Consistency
  # cut at -2 0.851334   0.7891051


  # >   ans$Marginal
  # Accuracy Consistency
  # cut at 0 0.8210774   0.7515717

  # >   ans$Marginal
  # Accuracy Consistency
  # cut at 2 0.8608774   0.8093206


  # >   ans$Marginal
  # Accuracy Consistency
  # cut at 4 0.9031811   0.8681039

  # >   ans$Marginal
  # Accuracy Consistency
  # cut at 6 0.9552208   0.9340635

# average of error variance
ErrorVarAvg <- mean(scoSS_MLE$varC)
ErrorVarAvg
# 3.733342

ThetaEstVar <- var(scoSS_MLE$thetaSum)
ThetaEstVar
# 8.928242

# TrueVar <- numOfFactors
TrueVar <- 2*(sum(cor)) + numOfFactors
TrueVar

r3 <- TrueVar/ThetaEstVar
r3

r4 <- 1 - ErrorVarAvg/ThetaEstVar
r4


r5 <- TrueVar/(TrueVar+ErrorVarAvg)
r5








VarTheta1 <- var(scoSS_MLE$theta1)
VarTheta2 <- var(scoSS_MLE$theta2)
VarTheta3 <- var(scoSS_MLE$theta3)

VarTheta1 + VarTheta2 + VarTheta3


set1 <- sqrt(VarTheta1)
set2 <- sqrt(VarTheta2)
set3 <- sqrt(VarTheta3)

ThetaEstVar2 <- VarTheta1 + VarTheta2 + VarTheta3 + 2 *cor[1] * set1 * set2 + 2 *cor[2] * set1 * set3 + 2 *cor[3] * set2 * set3
ThetaEstVar2
# 13.2846

# marginal reliability approach
MarginalRelSSMIRT_MLE_P <- ErrorVarAvg /(ErrorVarAvg + numOfFactors + 2*(sum(cor))) # var(e)/(var(e) + var(theta))


(2*(sum(cor)) + numOfFactors)  /(2*(sum(cor)) + numOfFactors + ErrorVarAvg)


MarginalRelSSMIRT_MLE_P <- (ThetaEstVar2 - ErrorVarAvg) /ThetaEstVar2


# coefficients
MarginalRelSSMIRT_MLE_P










####### EAP ---------------------------------------

# Form A
scoSS_EAP <- read.table("TestData/SpanishLit_sco_A_SS_EAP.txt")[,c(3:14)]
cor <- c(0.9067069, # 1&2
         0.6994119, # 1&3
         0.4891160) # 2&3

# Form B
# scoSS_EAP <- read.table("TestData/SpanishLit_sco_B_SS_EAP.txt")[,c(3:14)]
# cor <- c(0.9722234, # 1&2
#          0.5602197, # 1&3
#          0.4795721) # 2&3

# change variable name
names(scoSS_EAP) <- c("theta1", "theta2", "theta3", "se1", "se2", "se3", "var11", "var21", "var22", "var31","var32","var33")

var(scoSS_EAP$theta1)
mean(scoSS_EAP$var11)

var(scoSS_EAP$theta2)
mean(scoSS_EAP$var22)

var(scoSS_EAP$theta3)
mean(scoSS_EAP$var33)

# composite error variance
scoSS_EAP<- transform( scoSS_EAP,
                   varC = var11 + var22  + var33 +  2*(var21 + var31 + var32),
                   thetaSum = theta1 + theta2 + theta3
)

### classification -----------------------

cutscore <- 2
theta <- scoSS_EAP$thetaSum
sem <- sqrt(scoSS_EAP$varC)

os<-theta
nn<-length(os)  # nn, number of examinee
nc <- length(cutscore)  # number of cutscore

if(nn != length(sem)) stop("Ability and se of different length")

esacc<-matrix(NA,length(cutscore), nn, dimnames = list(paste("cut at",round(cutscore,3)), round(os,3)))
escon <-esacc


j=1 # test
cuts<-c(-Inf, cutscore[j], Inf)
categ<-cut(os,cuts,labels=FALSE,right=FALSE) # cut function in r

for(i in 1:nn) {
  esacc[j,i]<-(pnorm(cuts[categ[i]+1],os[i],sem[i])-pnorm(cuts[categ[i]],os[i],sem[i]))
  escon[j,i]<-((pnorm(cuts[2], os[i],sem[i]) - pnorm(cuts[1],os[i],sem[i]))^2	+ (pnorm(cuts[3], os[i],sem[i]) - pnorm(cuts[2],os[i],sem[i]))^2	 )
}

ans<- (list("Marginal" = cbind("Accuracy" = rowMeans(esacc), "Consistency" = rowMeans(escon)), "Conditional" = list("Accuracy" =t(esacc), "Consistency" = t(escon))))

ans$Marginal




# average of error variance
ErrorVarAvg <- mean(scoSS_EAP$varC)
ErrorVarAvg

# 1.925

# true score variance
TrueVar <- (2*(sum(cor)) + numOfFactors)
TrueVar
ObsVar <- var(scoSS_EAP$thetaSum)
ObsVar


r71 <- ObsVar/TrueVar
r71
# marginal reliability approach
MarginalRelSSMIRT_EAP_P  <- 1 - ErrorVarAvg /TrueVar # var(e)/(var(e) + var(theta))

# coefficients
MarginalRelSSMIRT_EAP_P

r73 <- ObsVar/(ObsVar + ErrorVarAvg)
r73


ObsVar + ErrorVarAvg - TrueVar








### Empirical approach

VarTheta1 <- var(scoSS_EAP$theta1)
VarTheta2 <- var(scoSS_EAP$theta2)
VarTheta3 <- var(scoSS_EAP$theta3)

VarTheta1
VarTheta2
VarTheta3

set1 <- sqrt(VarTheta1)
set2 <- sqrt(VarTheta2)
set3 <- sqrt(VarTheta3)

set1
set2
set3

ThetaEstVar <- VarTheta1 + VarTheta2 + VarTheta3 + 2 * ( cor[1] * set1 * set2 + cor[2] * set1 * set3 + cor[3] * set2 * set3 )
ThetaEstVar


(ThetaEstVar)/(ThetaEstVar+ ErrorVarAvg)


















########### bi factor full D method ---------------------------------------

# item parameters

# read item parameters from txt file
itemPara_BF <- read.table("TestData/SpanishLit_prm_B_BF.txt")[,c(7:11)]


names(itemPara_BF) <- c("b", "ag","a1","a2", "a3") # a1 is primary
itemPara_BF$a <- c(itemPara_BF$a1[1:13], itemPara_BF$a2[14:25], itemPara_BF$a3[26:31])


itemPara_BF[,"b"] <-  -itemPara_BF_G[,"b"]/(itemPara_BF_G[,"ag"] + itemPara_BF_G[,"a"])#######?????????????????



itemPara_BF[,"a"] <- itemPara_BF[,"a"]/1.702

itemPara_BF$a1[1:13] <- itemPara_BF$a[1:13]
itemPara_BF$a2[14:25] <- itemPara_BF$a[14:25]
itemPara_BF$a3[26:31] <- itemPara_BF$a[26:31]


itemPara_BF$ag <- itemPara_BF$ag/1.702































