###################################################
############ SUMMARY STATISTICS ###################


library(readxl)
library(dplyr)
library(rethinking)
library(rstan)
library(tidyr)
library(dplyr)
library(grid)
library(gtools)
library(stringr)

# Please setwd to the spite folder
# setwd()

`%ni%` <- Negate(`%in%`)
dat<-read.csv("data/data.csv")

table(dat$`IN/OUTGROUP`, dat$SNTASKA)

m<-matrix(NA, 6, 5)
cnames <- c("Task A (15, 5) vs (10, 10)", "Task B (13, 18) vs (10, 10)",  "Task C (10, 5) vs (10, 10)", "Task D (15, 10) vs (10, 10)", "Percent")
colnames(m) <- cnames
rnames <- c("Altruism", "Aversion to AI",  "Egalitarianism", "Aversion to DI", "Spite", "Percent (10, 10)")
rownames(m) <- rnames
m[1,1] <- "(10, 10)"
m[1,2] <- "(13, 18)"
m[1,3] <- "(10, 10)"
m[1,4] <- "(15, 10)"
m[1,5] <-  paste0(round(mean(dat$TRUEALT, na.rm = T),2)*100, "%")
#AI
m[2,1] <- "(10, 10)"
m[2,2] <- " "
m[2,3] <- "(10, 10)"
m[2,4] <- " "
m[2,5] <-  paste0(round(mean(dat$TRUEALT, na.rm = T),2)*100, "%")
#EG
m[3,1] <- "(10, 10)"
m[3,2] <- "(10, 10)"
m[3,3] <- "(10, 10)"
m[3,4] <- "(10, 10)"
m[3,5] <-  paste0(round(mean(dat$EGAL, na.rm = T),2)*100, "%")
#DI
m[4,1] <- " "
m[4,2] <- "(10, 10)"
m[4,3] <- " "
m[4,4] <- "(10, 10)"
m[4,5] <-  paste0(round(mean(dat$disadv, na.rm = T),2)*100, "%")
#SPITE
m[5,1] <- "(15, 5)"
m[5,2] <- "(10, 10)"
m[5,3] <- "(10, 5)"
m[5,4] <- "(10, 10)"
m[5,5] <-  paste0(round(mean(dat$TRUESPITE, na.rm = T),2)*100, "%")


m[6, 4]<- paste0(round(mean(dat$`RESTASKA (10,10 vs 10, 15)` == 0, na.rm = T)*100), "%")
m[6, 3]<- paste0(round(mean(dat$`RESTASKB (10, 10 vs 10,5)` == 0, na.rm = T)*100), "%")
m[6, 2]<- paste0(round(mean(dat$`RESTASKC(10,10 vs 13,18)` == 0, na.rm = T)*100), "%")
m[6, 1]<- paste0(round(mean(dat$`RESTASKD(10,10 vs 15,5)` == 0, na.rm = T)*100), "%")

m[6, 5] = " "

xtable(m)