### create a R package

# guide
# https://hilaryparker.com/2014/04/29/writing-an-r-package-from-scratch/


library(devtools)
library(roxygen2)

#1
setwd("C:\\Users\\Huan\\OneDrive - University of Iowa\\")
create("EMRel")

setwd("./ReliabilityTest")
document()


setwd("..")
install("ReliabilityTest")

library(ReliabilityTest)

?IRTRel


install_github('ReliabilityTest','github_username')


install_github("liuhuan90123/EMReliability")

library(EMReliability)
lsf.str(EMReliability)

CronbachAlpha

