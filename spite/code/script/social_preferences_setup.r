########################################################################################
########################## Social Preferences ##########################################

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

vars <- c("TRUEALT", "TRUESPITE", "EGAL", "disadv", "adadv", "selfish")

#CONSTRUCT A

temp  <- (dat$AGE)
quant <- (quantcut(dat$AGE, 10))
dat$AGE2 <- as.integer(quant)
distances  <- (1/2)*(quantile(temp, seq(0.01, 1, by =.1), na.rm = TRUE) + quantile(temp, seq(0, .99, by =.1), na.rm = TRUE)) 
LA <- length(distances)

DmatA <- matrix(nrow=length(distances), ncol=length(distances))
for(j in 1:length(distances)){
  for(k in 1:length(distances)){
    DmatA[j,k] <- abs(distances[j]-distances[k])
  }
}


#CONSTRUCT D

temp  <- (dat$'25DROUGHT' + dat$`25FLOOD` + dat$`25INFEST` )
quant <- 
distances  <- min(temp, na.rm = TRUE):max(temp, na.rm = TRUE) 
LD <- length(distances)

DmatD <- matrix(nrow=length(distances), ncol=length(distances))
for(j in 1:length(distances)){
  for(k in 1:length(distances)){
    DmatD[j,k] <- abs(distances[j]-distances[k])
  }
}


####################################################################################
############################ BEFORE 25 #############################################

#RUN LOOP


for (k in vars){
  
  d <- data.frame(
    y =  dat[,k],
    village = dat$village,
    C = as.integer(interaction(dat$district, dat$`IN/OUTGROUP`)),
    DROUGHT = temp+1L,
    AGE = dat$AGE2,
    TA = dat$SNTASKA+1L,
    TB = dat$SNTASKB+1L,
    TC = dat$SNTASKC+1L,
    TD = dat$SNTASKD+1L,
    IN = -standardize((dat$LANDINHERETED)),
    G = dat$SEX,
    ET = dat$ETHNICITY
  )
  
  d$DROUGHT[d$DROUGHT == 98] <- NA
  d<-d[complete.cases(d),]
  names(d)[1] <- "y"
  dl  <- as.list(d)
  dl$N <- length(dl$y)
  
  dl  <- as.list(d)
  dl$N <- length(dl$y)
  dl$LD <- LD
  dl$LA <- LA
  dl$DmatDROUGHT <- DmatD/10
  dl$DmatAGE <-DmatA
  dl$ET <-as.integer(as.factor(dl$ET))
  dl$NET <- length(unique(dl$ET))
  
  saveRDS(dl, "dl.RDS")
  
  #Run Model  
#  base=stan("code/models/gp_aggregate_Age.stan", data = dl, cores = 4)
  base=stan("code/models/agg_gp_b25_full.stan", data = dl, cores = 4, chains = 4)

  
  saveRDS(base, paste0("data/_",k,"_b25_ag_full.stan"))
 # saveRDS(base2, paste("data/",k,"gp_aggeregate_Age_inter.RDS"))
  
}


###############################################################################
############################## TOTAL ##########################################



#RUN LOOP

#CONSTRUCT D

temp  <- (dat$TOTALDROUGHT + dat$TOTALFLOOD + dat$TOTALINFEST )
quant <- 
  distances  <- min(temp, na.rm = TRUE):max(temp, na.rm = TRUE) 
LD <- length(distances)

DmatD <- matrix(nrow=length(distances), ncol=length(distances))
for(j in 1:length(distances)){
  for(k in 1:length(distances)){
    DmatD[j,k] <- abs(distances[j]-distances[k])
  }
}


vars <- c("TRUEALT", "TRUESPITE", "EGAL", "disadv", "adadv", "selfish", "competitive")

for (k in vars){
  
  d <- data.frame(
    y =  dat[,k],
    village = dat$village,
    C = as.integer(interaction(dat$district, dat$`IN/OUTGROUP`)),
    DROUGHT = temp, 
    AGE = dat$AGE2,
    G = dat$SEX,
    IN = standardize(dat$LANDINHERETED),
    TA = dat$SNTASKA+1L,
    TB = dat$SNTASKB+1L,
    TC = dat$SNTASKC+1L,
    TD = dat$SNTASKD+1L,
    ET = dat$ETHNICITY
  )
  
  d$DROUGHT[d$DROUGHT == 98] <- NA
  d<-d[complete.cases(d),]
  names(d)[1] <- "y"
  dl  <- as.list(d)
  dl$N <- length(dl$y)
  dl  <- as.list(d)
  dl$N <- length(dl$y)
  dl$LD <- LD
  dl$LA <- LA
  dl$DmatDROUGHT <- DmatD/10
  dl$DmatAGE <-DmatA
  dl$ET <-as.integer(as.factor(dl$ET))
  dl$NET <- length(unique(dl$ET))
  saveRDS(dl, "dl.RDS")
  
  #Run Models
  base=stan("code/models/b25_ag_full.stan", data = dl, cores = 4)
  saveRDS(base, paste0("data/_",k,"_ALL_AGAGE.RDS"))
  post <- extract.samples(base)
  
  
  
}

# Predict

r<-summary(glm(dl$y ~ dl$AGE + dl$DROUGHT + dl$IN + dl$DROUGHT*dl$IN + (dl$DROUGHT | dl$C), family = "binomial"))

D = 10
I = 2
inv_logit(coef(r)[1,1] + coef(r)[3,1]*D + coef(r)[4,1]*I + coef(r)[5,1]*I*D)

##########################################################################################
############################### ONLY IRT #################################################



for (k in vars){
  
  d <- data.frame(
    y =  dat[,k],
    region = as.integer(dat$district),
    village = dat$village,
    AGE = dat$AGE2,
    C = as.integer(interaction(dat$district, dat$`IN/OUTGROUP`)),
    IN = -standardize((dat$LANDINHERETED)),
    G = dat$SEX,
    ET = dat$ETHNICITY,
    outgroup = as.integer(dat$`IN/OUTGROUP`+1),
    dgift = ifelse(dat$GIFTDROUGHT > 0 , 0, 1), # Initial coding double check
    dgov = ifelse(dat$GOVAIDDROUGHT > 0, 0, 1), # Initial coding double check
    ddebt = ifelse(dat$DEBTDROUGHT > 0, 1, 0),
    dmig = ifelse(dat$MIGRATEDROUGHT > 0, 1, 0),
    dwild = ifelse(dat$WILDFOODDROUGHT > 0, 1, 0),
    dliq =  ifelse(dat$LIQUIDATIONDROUGHT >0, 1,0),
    dliv =  ifelse(dat$LIVESTOCKDROUGHT >0, 1,0),
    fgift = ifelse(dat$GIFTFLOOD > 0 , 0, 1), # Initial coding double check
    fgov = ifelse(dat$GOVAIDFLOOD > 0, 0, 1), # Initial coding double check
    fdebt = ifelse(dat$DEBTFLOOD > 0, 1, 0),
    fmig = ifelse(dat$MIGRATEFLOOD > 0, 1, 0),
    fwild = ifelse(dat$WILDFOODFLOOD > 0, 1, 0),
    fliq =  ifelse(dat$LIQUIDATIONFLOOD >0, 1,0),
    igift = ifelse(dat$GIFTINFEST > 0 , 0, 1), # Initial coding double check
    igov = ifelse(dat$GOVAIDINFEST > 0, 0, 1), # Initial coding double check
    idebt = ifelse(dat$DEBTINFEST > 0, 1, 0),
    imig = ifelse(dat$MIGRATEINFEST > 0, 1, 0),
    iwild = ifelse(dat$WILDFOODINFEST > 0, 1, 0),
    iliq =  ifelse(dat$LIQUIDATIONINFEST >0, 1,0)
  )
  
  names(d)[1] <- "y"
  d<-d[complete.cases(d),]
  
  dl  <- as.list(d[,1:9])
  dl$Y <- as.matrix(d[,10:28])
  dl$M1 <- ncol(dl$Y)
  dl$N <- length(dl$y)
  dl$LA <- LA
  dl$DmatAGE <-DmatA
  dl$ET <-as.integer(as.factor(dl$ET))
  dl$NET <- length(unique(dl$ET))
  
  
  #Run Model
  #out1=stan("code/models/irt_ag.stan", data = dl, cores = 4)
  base=stan("code/models/irt_ag_age.stan", data = dl, cores = 4)
  saveRDS(base, paste0("data/_",k,"_IRT_AG.RDS"))
  
  
}


#######################################################################################
#################### SEVERITY OF EXPOSURE BINARY ######################################




#Most Sever Drought
a<-str_split(dat$SEVEREDROUGHT, " ")
mostD <-c()
for(i in 1:length(a)){
  mostD<-c(mostD, a[[i]][1])
}

mostD <- as.numeric(mostD)
mostD[mostD == 0] <- NA
mostD <- mostD-(2014-dat$AGE)
mostD <- ifelse(is.na(mostD) | mostD > 26, 26, mostD)
Don <- ifelse(mostD == 26, 0, 1)

temp  <- mostD
distances  <- 1:26  
LD <- length(distances)

DmatD <- matrix(nrow=length(distances), ncol=length(distances))
for(j in 1:length(distances)){
  for(k in 1:length(distances)){
    DmatD[j,k] <- abs(distances[j]-distances[k])
  }
}

#Most Sever FLOOD
a<-str_split(dat$SEVEREFLOOD, " ")
mostF <-c()
for(i in 1:length(a)){
  mostF<-c(mostF, a[[i]][1])
}

mostF <- as.numeric(mostF)
mostF[mostF == 0] <- NA
mostF <- mostF-(2014-dat$AGE)
mostF <- ifelse(is.na(mostF) | mostF > 26, 26, mostF)
Fon <- ifelse(mostD == 26, 0, 1)

temp  <- mostF
distances  <- 1:26 
LF <- length(distances)

DmatF <- matrix(nrow=length(distances), ncol=length(distances))
for(j in 1:length(distances)){
  for(k in 1:length(distances)){
    DmatF[j,k] <- abs(distances[j]-distances[k])
  }
}

#Most Sever IFNESTATION
a<-str_split(dat$SEVEREINFEST, " ")
mostI <-c()
for(i in 1:length(a)){
  mostI<-c(mostI, a[[i]][1])
}

mostI <- as.numeric(mostI)
mostI[mostI == 0] <- NA
mostI <- mostI-(2014-dat$AGE)
mostI <- ifelse(is.na(mostI) | mostI >26, 26, mostI)
#CONSTRUCT I
mostI[39]=26

Ion <- ifelse(mostI==26, 0, 1)

temp  <- mostI
distances  <- 1:26  
LI <- length(distances)

DmatI <- matrix(nrow=length(distances), ncol=length(distances))
for(j in 1:length(distances)){
  for(k in 1:length(distances)){
    DmatI[j,k] <- abs(distances[j]-distances[k])
  }
}

S = Ion + Fon + Don
S = ifelse(S > 0 , 2, 1)

#RUN LOOP
vars <- c("TRUEALT", "TRUESPITE", "EGAL", "disadv", "adadv", "selfish", "competitive" )



for (k in vars){
  d <- data.frame(
    y =  dat[,k],
    village = dat$village,
    region = dat$district,
    C = as.integer(interaction(dat$district, dat$`IN/OUTGROUP`)),
    DROUGHT = mostD,
    FLOOD = mostF,
    INFEST = mostI,
    AGE = dat$AGE2,
    S =S,
    outgroup = dat$`IN/OUTGROUP`,
    C = as.integer(interaction(dat$district, dat$`IN/OUTGROUP`)),
    IN = -standardize((dat$LANDINHERETED)),
    G = dat$SEX,
    ET = dat$ETHNICITY
    
  )
  d<-d[complete.cases(d),]
  names(d)[1] <- "y"
  dl  <- as.list(d)
  dl$N <- length(dl$y)
  
  dl  <- as.list(d)
  dl$N <- length(dl$y)
  dl$LA <- LA
  dl$DmatAGE <-DmatA
  dl$ET <-as.integer(as.factor(dl$ET))
  dl$NET <- length(unique(dl$ET))

  base=stan("code/models/sever_ag_age.stan", data = dl, cores = 4)
  saveRDS(base, paste0("data/_",k,"_years_before_25AG.RDS"))
  
}
