setwd("C:/Users/jeffrey_andrews/OneDrive/Documents/phd")
vars <- c("TRUEALT", "adadv",  "EGAL", "disadv",  "TRUESPITE")
names <- c("(Altruism)", "(Adv. aversion)", "(Egalitarianism)", "(Dis. aversion)", "(Spite)")

collist <- c("lightblue", "darkblue", "orange", "firebrick")
png("output/interaction.png", height = 350*5, width = 1200)

par(mfrow=c(5,3), mar = c(4, 4.5,1,1), oma = c(6,2,1,1))
parental_wealth <- c(-.1,.1)
parental_wealth <- c(-1,1)

cnt2 <- 1
for(var in vars){

  base2<-readRDS(paste0("data/_",var,"_gp_aggeregate_Age_inter.RDS"))
  post <- extract.samples(base2)
  lmu <- list()
  lPI <- list()
  cnt <- 1
  
  for(j in parental_wealth){
    for(i in 1:2){
      if(i == 1)  out = c(1,2,3)
      if(i == 2)  out = c(4,5,6)
      
      x.seq <- 1:15
      link = function(x) rbinom(4000, 1000, inv_logit(rowMeans(post$aV[,])*post$Sigma_V + (rowMeans(post$aX[,out])*post$Sigma_CC) +
                                                        post$village_mu +
                                                        rowMeans(post$bG[, out, ]) +
                                                        rowMeans(post$kD[, out,x]) +
                                                        (rowMeans(post$bIN[, out])*post$Sigma_IN + post$bIN_mu)*j +
                                                        rowMeans(post$kD2[, out,x])*j +
                                                        rowMeans(post$kA[,out,])))  
      mu<-sapply(x.seq, link)/1000
      lPI[[cnt]] <- apply(mu, 2, PI, .67)
      lmu[[cnt]] <- colMeans(mu)
      cnt <- cnt+1
    } 
  }  
  #Plot
  plot(lmu[[1]], type = "l", ylim = c(0, .4), col = "lightblue", xlab = ifelse(cnt2 == 5, "Number of shocks before 25", ""),
       ylab = paste("P", names[cnt2]), cex.lab = 1.5, cex.axis = 1.5)
  #rethinking::shade(lPI[[1]], x.seq, col = col.alpha(collist[1], .2))
  lines(x.seq, lmu[[2]], col = "darkblue")
  lines(x.seq, lmu[[3]], col = "orange")
  lines(x.seq, lmu[[4]], col = "firebrick")
  #rethinking::shade(lPI[[2]], x.seq, col = col.alpha("blue", .05))
  
  
  ######################################################
  ################## IRT ###############################
  
  base2<-readRDS(paste0("data/_",var,"_IRT_AG_interact.RDS"))
  post <- extract.samples(base2)
  lmu <- list()
  lPI <- list()
  cnt <- 1
  for(k in parental_wealth){
      for (j in 1:2){
        x.seq <- seq(from = -4, to = 4, length.out =100)
        if(j == 1)  out = c(1,2,3)
        if(j == 2)  out = c(4,5,6)
        link = function(x) rbinom(4000, 1000, inv_logit(rowMeans(post$av)*post$sigma_V +
                                                          post$mu +
                                                          rowMeans(post$aX[, ,j]) +
                                                          rowMeans(post$bz[,,j])*x  +
                                                         (rowMeans(post$bIN[, out])*post$Sigma_IN + post$bIN_mu) *k +
                                                          rowMeans(post$bz2[,,j])*x*k+
                                                          rowMeans(post$kA[,out,])+
                                                          rowMeans(post$bG[,out,])
                                                          ))    
        mu<-sapply(x.seq, link)/1000
        lPI[[cnt]] <- apply(mu, 2, PI, .67)
        lmu[[cnt]] <- colMeans(mu)
        cnt <- cnt+1
      }
    }
  
    #Plot
    plot(lmu[[1]]~ x.seq, type = "l", ylim = c(0, .4), col = "lightblue",
         xlab = ifelse(cnt2 == 5, "IRT exposure", ""), 
         ylab = paste0("P", names[cnt2]), cex.lab = 1.5, cex.axis = 1.5)
    #rethinking::shade(lPI[[1]], x.seq, col = col.alpha(collist[1], .2))
    lines(x.seq, lmu[[2]], col = "darkblue")
    lines(x.seq, lmu[[3]], col = "orange")
    lines(x.seq, lmu[[4]], col = "firebrick")
    #rethinking::shade(lPI[[2]], x.seq, col = col.alpha("blue", .2))
    
  
  #############################################################################
  ############ SEVER BEFORE 25 ################################################
  
  
  base2 <- readRDS(paste0("data/_",var,"_years_before_25AG_interact.RDS"))
  post<-extract.samples(base2)
  
  lmu <- list()
  lPI <- list()
  cnt <- 1
  for(j in parental_wealth){

    for(i in 1:2){
      if(i == 1)  out = c(1,2,3)
      if(i == 2)  out = c(4,5,6)
      
      x.seq <- 1:2
      link = function(x) rbinom(4000, 1000, inv_logit(rowMeans(post$aV[,])*post$sigma_V +
                                                        rowMeans(post$aX[,out])*post$sigma_X +
                                                        rowMeans(post$bz[, x,out]) +
                                                        rowMeans(post$bz2[, x,out])*j +
                                                        (rowMeans(post$bIN[, out])*post$Sigma_IN + post$bIN_mu) * j +
                                                        rowMeans(post$bG[, out,]) +
                                                        rowMeans(post$kA[, i, ])+ post$a) )  
      mu<-sapply(x.seq, link)/1000
      lPI[[cnt]] <- apply(mu, 2, PI, .67)
      lmu[[cnt]] <- (mu)
      cnt <- cnt+1
    } 
   }
if(cnt2 == 5) namest  =c("None Before 25", "None Before 25", "Before 25", "Before 25", "None Before 25", "None Before 25", "Before 25", "Before 25")
if(cnt2 != 5) namest  =c("", "", "", "", "", "", "", "")
  boxplot(lmu[[1]][, 1], lmu[[3]][, 1], lmu[[1]][, 2], lmu[[3]][, 2], lmu[[2]][, 1], lmu[[4]][, 1], lmu[[2]][, 2], lmu[[4]][, 2], ylim = c(0, .5),
          col = c( col.alpha(collist[1], .2), col.alpha(collist[3], .2), col.alpha(collist[1], .2), col.alpha(collist[3],  .2), col.alpha(collist[2], .2), col.alpha(collist[4], .2), col.alpha(collist[2], .2), col.alpha(collist[4],  .2)),
          ylab = paste0("P", names[cnt2]), names = namest, cex.lab = 1.5, cex.axis = 1.5, las = 2)
          
  cnt2 <- cnt2 +1
}
dev.off()