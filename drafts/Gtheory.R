# install.packages("gtheory")


library(gtheory)


#Conduct a univariate G study.
#Compare to results on page 116 of Brennan (2001).
data(Brennan.3.2)
formula.Brennan.3.2 <- "Score ~ (1 | Person) + (1 | Task) +
(1 | Rater:Task) + (1 | Person:Task)"
gstudy(data = Brennan.3.2, formula = formula.Brennan.3.2)
#Conduct a multivariate G study.
#Compare to results on page 270 of Brennan (2001).
data(Rajaratnam.2)
formula.Rajaratnam.2 <- "Score ~ (1 | Person) + (1 | Item)"
gstudy(data = Rajaratnam.2, formula = formula.Rajaratnam.2, colname.strata = "Subtest",
       colname.objects = "Person")


Person <- as.factor(rep(1:6,each = 4))
Item <- as.factor(rep(1:4,times = 6))
Score <- c(9,9,7,4,9,8,4,6,8,8,6,2,
           9,8,6,3,10,9,8,7,6,4,5,1)
pi_dat <- data.frame(Person,Item,Score)




rawData_A <- read.table("TestData/RawDataFormX.txt")
rawData_B <- read.table("TestData/RawDataFormY.txt")


rawData_A <- read.table("TestData/FormA_31_3000.txt")
rawData_B <- read.table("TestData/FormB_31_3000.txt")



rawData <- rawData_A
rawData <- rawData_B



numOfItem <- ncol(rawData)
numOfPerson <- nrow(rawData)

SSp <- numOfItem * sum((rowSums(rawData)/numOfItem)^2) - numOfPerson * numOfItem * (sum(rawData)/(numOfItem*numOfPerson))^2
SSi <- numOfPerson * sum((colSums(rawData)/numOfPerson)^2) - numOfPerson * numOfItem * (sum(rawData)/(numOfItem*numOfPerson))^2
SSpi <- sum(rawData^2) - numOfItem * sum((rowSums(rawData)/numOfItem)^2) - numOfPerson * sum((colSums(rawData)/numOfPerson)^2) + numOfPerson * numOfItem * (sum(rawData)/(numOfItem*numOfPerson))^2

MSp <- SSp / (numOfPerson - 1)
MSi <- SSi / (numOfItem - 1)
MSpi <- SSpi / ((numOfPerson - 1)*(numOfItem - 1))


VARp <- (MSp - MSpi)/ numOfItem
VARi <- (MSi - MSpi)/ numOfPerson
VARpi <- MSpi


VARp
VARi
VARpi


numOfItemD <- 50

gSEM <- as.data.frame(matrix(1:numOfItemD, nrow = numOfItemD, ncol = 1))
names(gSEM) <- "numOfItemD"

gSEM$VARabso <- (VARi + VARpi)/gSEM$numOfItemD
gSEM$VARrela <- (VARpi)/gSEM$numOfItemD

gSEM$SEMabso <- sqrt(gSEM$VARabso)
gSEM$SEMrela <- sqrt(gSEM$VARrela)

gSEM$gcoef <- VARp/(VARp + gSEM$VARrela)
gSEM$dcoef <- VARp/(VARp + gSEM$VARabso)



library(ggplot2)



ggplot(gSEM, aes(numOfItemD)) +
  geom_text(aes(y = SEMabso, label="Delta"), parse=TRUE) +
  geom_text(aes(y = SEMrela, label="delta"), parse=TRUE) +
  ggtitle("SEM Plot") + xlab("D Study Sample Sizes") + ylab("SEM")


ggplot(gSEM, aes(numOfItemD)) +
  geom_text(aes(y = gcoef, label="rho"), parse=TRUE) +
  geom_text(aes(y = dcoef, label="phi"), parse=TRUE) +
  ggtitle("Cofficient Plot") + xlab("D Study Sample Sizes") + ylab("Cofficient")


library(xlsx)
write.xlsx(gSEM, "gSEM_B_comp.xlsx")






# new
Person <- as.factor(rep(1:numOfPerson,each = numOfItem))
Item <- as.factor(rep(1:numOfItem,times = numOfPerson))
Score <- as.vector(t(as.matrix(rawData)))
pi_dat <- data.frame(Person,Item,Score)


# summary(aov(Score~Person+Item, data
#             = pi_dat))

formula1 <- Score ~ (1|Person)+(1|Item)
g1 <- gstudy(data = pi_dat, formula1)
d1 <- dstudy(g1,colname.objects="Person",
             colname.scores="Score",
             data= pi_dat)


n_i <- c(1:50)
#relative error variance
rel_err_var <- g1$components[3,2]/n_i
#absolute error variance
abs_err_var <-
  g1$components[2,2]/n_i+g1$components[3,
                                       2]/n_i
#calculate generalizability coefficient
gen_coef <-
  g1$components[1,2]/(g1$components[1,2]
                      + rel_err_var)
#calculate dependability coefficient
dep_coef <-
  g1$components[1,2]/(g1$components[1,2]
                      + abs_err_var)
round(rel_err_var,2)
round(abs_err_var,2)
round(gen_coef,2)
round(dep_coef,2)



###  g theory stratified design ---------------


compEx <- read.xlsx("gtheorycompositedata.xlsx", 1)

compEx$M1 <- (compEx$i1 + compEx$i2)/2
compEx$M2 <- (compEx$i3 + compEx$i4 + compEx$i5 + compEx$i6)/4
compEx$M3 <- (compEx$i7 + compEx$i8)/2

compEx$M <- rowSums(compEx[,1:8])/8



rawData <- compEx[,1:8]


rawData <- compEx[,9:11]


numOfItem <- ncol(rawData)
numOfPerson <- nrow(rawData)

SSp <- numOfItem * sum((rowSums(rawData)/numOfItem)^2) - numOfPerson * numOfItem * (sum(rawData)/(numOfItem*numOfPerson))^2
SSi <- numOfPerson * sum((colSums(rawData)/numOfPerson)^2) - numOfPerson * numOfItem * (sum(rawData)/(numOfItem*numOfPerson))^2
SSpi <- sum(rawData^2) - numOfItem * sum((rowSums(rawData)/numOfItem)^2) - numOfPerson * sum((colSums(rawData)/numOfPerson)^2) + numOfPerson * numOfItem * (sum(rawData)/(numOfItem*numOfPerson))^2

MSp <- SSp / (numOfPerson - 1)
# 13.42857
MSi <- SSi / (numOfItem - 1)
MSpi <- SSpi / ((numOfPerson - 1)*(numOfItem - 1))



sum((compEx[1,1:8]-4.5)^2)/7

######   CTT formula comound binominal -------------------




8*4.5^2 - 2* 4.5^2 - 4 * 3.75^2 - 2*6^2



var(as.vector(as.matrix(compEx[1,1:8])))


1.714286^2

compEx[1,9:11]
4.5 + 3.75 +  6

SSi <- 1 * sum((c(4.5,3.75,6)/1)^2) - 1 * 3 * (14.25/(3*1))^2


var(compEx[1,1:8])


v1 <- c(4,  5)
v2 <- c(3 , 3 , 5 , 4)
v3 <- c(5 , 7)


var(v1)
var(v2)
var(v3)




