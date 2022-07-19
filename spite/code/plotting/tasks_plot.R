####################################################
################# PLOTING FOR PHD ##################

vars <- c("taskd", "taskc", "taskb", "taska")
names <- c("(15,5)", "(13, 18)", "(10, 5)", "(10, 15)")
#note that the reverse coding if for rhetorical ease in the paper!

png("output/tasks_agg_mar.png", height = 400*4, width =  1200)
cnt2 <- 1
par(mfrow=c(4,3), mar=c(5,5,1,1), oma = c(5,4,1,1))

for(k in vars){
  collist <- c("darkorange", "dodgerblue")
  
  ###########################
  ######## BEOFRE 25 ########
  
  ifelse(cnt2 == 4, "Number of shocks before 25", " ")
  
  base<-readRDS(paste0("data/_",k,"_gp_aggeregate_Age.RDS"))
  post <- extract.samples(base)
  
  lmu <- list()
  lPI <- list()
  cnt <- 1
  for(i in 1:2){
    if(i == 1)  out = c(1,2,3)
    if(i == 2)  out = c(4,5,6)
    
    x.seq <- 1:15
    link = function(x) rbinom(4000, 1000, inv_logit(rowMeans(post$aV[,])*post$Sigma_V + (rowMeans(post$aX[,out])*post$Sigma_CC) +
                                                      post$village_mu +
                                                      rowMeans(post$bG[, out, ]) +
                                                      rowMeans(post$kD[, out,x])))  
    mu<-sapply(x.seq, link)/1000
    lPI[[cnt]] <- apply(mu, 2, PI, .67)
    lmu[[cnt]] <- colMeans(mu)
    cnt <- cnt+1
  } 
  #Plot
  plot(lmu[[1]], type = "l", ylim = c(0, 1), col = collist[1], xlab = ifelse(cnt2 == 4, "Number of shocks before 25", " "),
       ylab = paste("P", names[cnt2]), cex.lab = 1.5, cex.axis = 1.5)
  rethinking::shade(lPI[[1]], x.seq, col = col.alpha(collist[1], .2))
  lines(x.seq, lmu[[2]], col = collist[2])
  rethinking::shade(lPI[[2]], x.seq, col = col.alpha(collist[2], .2))
  
  
  #######################
  ####### IRT ###########
  
  
  base<-readRDS(paste0("data/_",k,"_IRT_AG.RDS"))
  post <- extract.samples(base)
  lmu <- list()
  lPI <- list()
  cnt <- 1
  for(i in 1:3){
    for (j in 1:2){
      if(i == 1)  vil = c(2,3)
      if(i == 2)  vil = c(6,5)
      if(i == 3)  vil = c(1,4)
      x.seq <- seq(from = -4, to = 4, length.out =100)
      
      link = function(x) rbinom(4000, 1000, inv_logit(rowMeans(post$av)*post$sigma_V + post$mu + rowMeans(post$aX[, , j]) +
                                                        rowMeans(post$bz[,,j])*x ))    
      mu<-sapply(x.seq, link)/1000
      lPI[[cnt]] <- apply(mu, 2, PI, .67)
      lmu[[cnt]] <- colMeans(mu)
      cnt <- cnt+1
    }
  }
  #Plot
  plot(x.seq, lmu[[1]], type = "l", ylim = c(0, 1), col = collist[1], xlab = ifelse(cnt2 == 4, "Latent exposure score", " "),
       ylab = paste("P", names[cnt2]), cex.lab = 1.5,cex.axis = 1.5)
  rethinking::shade(lPI[[1]], x.seq, col = col.alpha(collist[1], .2))
  lines(x.seq, lmu[[2]], col = collist[2])
  rethinking::shade(lPI[[2]], x.seq, col = col.alpha(collist[2], .2))
  
  
  ################################
  ######## Severe before 25 ######
  
  base <- readRDS(paste0("data/_",k,"_years_before_25AG.RDS"))
  post<-extract.samples(base)
  lmu <- list()
  lPI <- list()
  cnt <- 1
  for(i in 1:2){
    if(i == 1)  out = c(1,2,3)
    if(i == 2)  out = c(4,5,6)
    
    x.seq <- 2:1
    link = function(x) rbinom(4000, 1000, inv_logit(rowMeans(post$aV[,vil])*post$sigma_V + rowMeans(post$aX[,out])*post$sigma_X +
                                                      rowMeans(post$bz[, x,out]) +
                                                      rowMeans(post$bG[, out,]) +
                                                      post$kA[, i, 3]+ post$a) )  
    mu<-sapply(x.seq, link)/1000
    lPI[[cnt]] <- apply(mu, 2, PI, .67)
    lmu[[cnt]] <- (mu)
    cnt <- cnt+1
  } 
  
  if(cnt2 ==4)namest= c("Before 25", "None before 25", "Before 25", "None before 25")
  if(cnt2 !=4)namest=c(" ", " ", " ", " ")
  boxplot(lmu[[1]][, 1], lmu[[1]][, 2], lmu[[2]][, 2], lmu[[2]][, 2], ylim = c(0, 1),
          col = c( col.alpha(collist[1], .5), col.alpha(collist[1], .5), col.alpha(collist[2], .5), col.alpha(collist[2], .5)),
          ylab = paste(" P", names[cnt2]),
          names=namest,
          cex.lab = 1.5, cex.axis = 1.5, las = 2)
  cnt2 <- cnt2 +1
}

mtext("Task A", side = 2, adj = .895, line = 0, outer = T)
mtext("Task B", side = 2, adj = .64, line = 0, outer = T)
mtext("Task C", side = 2, adj = .38, line = 0, outer = T)
mtext("Task D", side = 2, adj = .125, line = 0, outer = T)


dev.off()