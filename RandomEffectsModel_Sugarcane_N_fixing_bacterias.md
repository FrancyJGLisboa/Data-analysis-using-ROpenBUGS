######################################

MLA<-read.table(file="MLA.txt",header=T)


library(R2OpenBUGS)
library(mcmcplots)
library(coda)


Rmodel <- function () {
  
  for(i in 1:n){
    y[i] ~ dnorm(mu[i],tau.res)
    e.y[i] <- y[i] - mu[i]
    mu[i] <- alpha +
      theta.treat[Treat[i]] + 
      theta.time[Time[i]] + 
      theta.BFN[BFN[i]] + 
      theta.time.treat[Time[i],Treat[i]] +
      theta.time.BFN[Time[i],BFN[i]] + 
      theta.treat.BFN[Treat[i],BFN[i]]
  }
  
  
  for(i.tr in 1:n.treat) {             
    theta.treat[i.tr] ~ dnorm(0, tau.treat)
  }
  
  for(i.tm in 1:n.time) {             
    theta.time[i.tm] ~ dnorm(0, tau.time)
  }
  
  for(i.bf in 1:n.BFN) {         
    theta.BFN[i.bf] ~ dnorm(0, tau.BFN)
    
  }
  
  for(i.tm in 1:n.time) {  
    
    for(i.tr in 1:n.treat) { 
      
      theta.time.treat[i.tm,i.tr] ~ dnorm(0, tau.time.treat)
    }
  }
  
  for(i.tm in 1:n.time) {  
    
    for(i.bf in 1:n.BFN) { 
      
      theta.time.BFN[i.tm,i.bf] ~ dnorm(0, tau.time.BFN)
    }
  }
  
  
  
  
  
  for(i.tr in 1:n.treat) {  
    
    for(i.bf in 1:n.BFN) { 
      
      theta.treat.BFN[i.tr,i.bf] ~ dnorm(0, tau.treat.BFN)
    }
  }
  
  
  
  
  #Prior on variance components (finite population standard deviation)
  
  alpha ~ dnorm(0,0.00001)
  
  tau.treat <-pow(sigma.treat,-2)
  sigma.treat ~ dunif(0,100)
  
  tau.time <- pow(sigma.time,-2)
  sigma.time ~ dunif(0,100)
  
  tau.BFN <-pow(sigma.BFN,-2)
  sigma.BFN ~ dunif(0,100)
  
  tau.time.treat <- pow(sigma.time.treat,-2)
  sigma.time.treat ~ dunif(0,100)
  
  tau.time.BFN <-pow(sigma.time.BFN,-2)
  sigma.time.BFN ~ dunif(0,100)
  
  tau.treat.BFN <- 1/(sigma.treat.BFN*sigma.treat.BFN)
  sigma.treat.BFN ~ dunif(0,100)
  
  tau.res <- pow(sigma.res,-2)
  sigma.res ~ dunif(0,100)
  
  
  
  
  #finite population standard deviation (variance components)
  s.res <- sd(e.y[])
  s.theta.treat <- sd(theta.treat[]) 
  s.theta.time<- sd(theta.time[])
  s.theta.BFN <- sd(theta.BFN[])
  s.theta.time.treat <- sd(theta.time.treat[,])
  s.theta.time.BFN <- sd(theta.time.BFN[,])
  s.theta.treat.BFN <- sd(theta.treat.BFN[,])
  
  ###########################
  
  ########
  ## Derived quantities
  
  
  for (i.tr in 1:n.treat){
    g.treat[i.tr] <- theta.treat[i.tr]-mean(theta.treat[])
  }
  
  ## main time effects
  for (i.tm in 1:n.time){
    g.time[i.tm] <- theta.time[i.tm]-mean(theta.time[])
  }
  
  ## main BFN effects
  for (i.bf in 1:n.BFN){
    g.BFN[i.bf] <- theta.BFN[i.bf]-mean(theta.BFN[])
  }
  
  ## time - treat interaction
  for (i in 1:n.time){
    for(j in 1:n.treat){
      g.time.treat[i,j] <- theta.time.treat[i,j]-mean(theta.time.treat[,])
      
    }
  }
  
  ## time - BFN interaction
  for (i in 1:n.time){
    for(j in 1:n.BFN){
      g.time.BFN[i,j] <- theta.time.BFN[i,j]-mean(theta.time.BFN[,])
      
    }
  }
  
  ## treat - BFN interaction
  for (i in 1:n.treat){
    for(j in 1:n.BFN){
      g.treat.BFN[i,j] <- theta.treat.BFN[i,j]-mean(theta.treat.BFN[,])
      
    }
  }
}


help(bugs)
filename <- file.path("C:/Users/Francy Lisboa/Desktop/RDados/Rmodel.txt")
filename

write.model(Rmodel,file.path("C:/Users/Francy Lisboa/Desktop/RDados/Rmodel.txt"))



# R function
## Supply data for variables
bugs.in <- function(infile=MLA){
  y <- infile$Dens
  n <- length(y)
  Time <- as.numeric(ordered(infile$Time))
  BFN <- as.numeric(ordered(infile$BFN))
  Treat <- as.numeric(ordered(infile$Treat))
  
  n.time <- max(Time)
  n.BFN <- max(BFN)
  n.treat <- max(Treat)
  
  ## Input data and initial values
  ## Initial values (for all unknown coefficients to be estimated)
  ## Three parallel MCMC chains to facilitate convergence diagnostic
  bugs.dat <- list(n=n, n.time=n.time, n.BFN=n.BFN, n.treat=n.treat, y=y-mean(y), Time=Time, Treat=Treat, BFN=BFN )
  
  inits1 <- list(alpha=1, theta.treat = rep(0, n.treat), theta.time=rep(0, n.time),
                 theta.BFN = rep(0,n.BFN),
                 theta.time.treat=matrix(0, ncol=n.time, nrow=n.treat),                   
                 theta.time.BFN=matrix(0, ncol=n.time, nrow=n.BFN),
                 theta.treat.BFN=matrix(0, ncol=n.treat, nrow=n.BFN),
                 sigma.res =1, sigma.treat =1, sigma.time=1, sigma.BFN=1,
                 sigma.time.treat=1, sigma.time.BFN=1, sigma.treat.BFN =1)
  
  inits2 <- list(alpha=0, theta.treat = rep(1, n.treat), theta.time=rep(1, n.time),
                 theta.BFN = rep(1,n.BFN),
                 theta.time.treat=matrix(1, ncol=n.time, nrow=n.treat),                   
                 theta.time.BFN=matrix(1, ncol=n.time, nrow=n.BFN),
                 theta.treat.BFN=matrix(1, ncol=n.treat, nrow=n.BFN),
                 sigma.res =0.5, sigma.treat =0.5, sigma.time=0.5, sigma.BFN=0.5,
                 sigma.time.treat=0.5, sigma.time.BFN=0.5, sigma.treat.BFN =0.5)
  
  inits3 <- list(alpha=0.5, theta.treat = rep(0.5, n.treat), theta.time=rep(0.5, n.time),
                 theta.BFN = rep(0.5,n.BFN),
                 theta.time.treat=matrix(0.5, ncol=n.time, nrow=n.treat),                   
                 theta.time.BFN=matrix(0.5, ncol=n.time, nrow=n.BFN),
                 theta.treat.BFN=matrix(0.5, ncol=n.treat, nrow=n.BFN),
                 sigma.res =0, sigma.treat =0, sigma.time=0, sigma.BFN=0,
                 sigma.time.treat=0, sigma.time.BFN=0, sigma.treat.BFN =0)
  
  
  inits <- list (inits1, inits2, inits3)
  
  parameters <- c("g.treat","g.time","g.BFN","g.time.treat","g.time.BFN","g.treat.BFN",
                  "s.theta.treat","s.theta.time","s.theta.BFN","s.theta.time.treat",
                  "s.theta.time.BFN","s.theta.treat.BFN","s.res")
  
  return(list(para=parameters, data=bugs.dat, inits=inits))
}


## Create input and initial data
input.to.bugs <- bugs.in()


## Do the simulation
bugs.out.S <- bugs(input.to.bugs$data, input.to.bugs$inits,input.to.bugs$para,
                   model.file="Rmodel.txt",
                   n.chains=3, n.iter=1000, n.thin=1, n.burnin=100, debug=T)
print(bugs.out.S )



gelman.diag(bugs.out.S,multivariate=F)


help(caterplot)
caterplot(bugs.out.S, parms=c("s.theta.treat","s.theta.time","s.theta.BFN","s.theta.time.treat",
                              "s.theta.time.BFN","s.theta.treat.BFN","s.res"),
                               style="plain",col="black",pch=15,cex=1,
                              lwd=c(1,1),quantiles=list(c(0.025,0.975)))


caterplot(bugs.out.S, parms=c("g.time"),style="plain",
          col="black",pch=15,cex=1,
          lwd=c(1,1),quantiles=list(c(0.025,0.975)))

caterplot(bugs.out.S, parms=c("g.time.treat"),style="plain",
          col="black",pch=15,cex=1,
          lwd=c(1,1),quantiles=list(c(0.025,0.975)))
