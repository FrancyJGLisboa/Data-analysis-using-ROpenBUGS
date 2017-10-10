# Data-analysis-using-ROpenBUGS

Data analysis routines involving linear mixed models under bayesian perspective by using OpenBUGS from R.

The code is about  the use of a varying-intercept model (random-effect model) to estimate the relative contribution of covariates on of colembola abundance under different green manure crops (GMC) in Rio de Janeiro-Brazil. The experiment was in a Latin-Square split-plot layout with the five GMC distrubuted across rows and columns each under two tillage systems. Here we are closely following the ideas from Gelman and Hill (2007) where they state that all variability sources should be treated as random-effects. In other words, we are assuming that the theta effects of a given variability source are a sample from a bigger population of effects and their estimates are not independent from each other. Gelman and Hill (2007) use the terms varying interceps and varying slopes to overcome the yet unclear and debatable definitions of fixed and random effects. The model below is a varying-intercept model in the sense that 1) effects are left to vary within a given categorical source of variability; 2) the effects within a source of variability are assumed to come from a commom probability distribution; 3) the estimate of effects within a source of variability is influenced by the sample size of the source (partial pooling effects).

 
# Experiment
In the Latin-Square split-plot the green manure crops are organized in a way that each crop occurs only one time in a row and only one time in a column. Additionaly, each crop plot was subdivided into two subplots each representing one type of tillage (conventional, no-tillage). Therefore, the abundance of colembolas was taken from a total of 50 observations (5 rows x 5 columns x 2 tillage systems).
  
  
# ROpenBUGS
To run this analysis I opted by the R package ROpenBUGS which allows to run Bayesian models  in R using the open source version of the WinBUGS program, the OpenBUGS. The model specification in ROpenBUGS was:

# General parameters
n.row = 5
n.col = 5
n.plots = 50
n.green = 5
n.til = 2

# Building the model
cat("model{

    # Data model Likelihood
        
    for(i in 1:n.plots){  # loop over 50 observations
    
    # We assume the data comes from a normal distribution with mean mu and precision tau.res
    y[i] ~ dnorm(mu[i],tau.res)
    
    # the linear predictor then assumes that the target variable (i.e. its variability) is a product of the additive combination of the follwing variables: 
    
    mu[i] <- alpha + # intercept 
    theta.row[ROW[i]] + # row effects
    theta.col[COL[i]] + # column effects
    theta.green[TREAT[i]] + # Green Manure Effects
    theta.plot[PLOT[i]] + # Plot effects
    theta.til[CT[i]] + # Tillage effects
    theta.row.til[ROW[i],CT[i]] + # Row-Tillage interaction effects
    theta.col.til[COL[i],CT[i]] + # Column-Tillage interaction effects
    theta.green.til[TREAT[i],CT[i]] + # Green Manure-Tillage interaction effects
    theta.plot.til[PLOT[i],CT[i]]  # Plot-Tillage effects

    }
    
    # Parameter models
    # Using Gelman and Hill (2007)'s parlance each source of variability is a batch of effects.
    # Within each batch the effects are assumed to be from a Normal distribution with mean 0 and precision (inverse of the variance) tau.
    
    # Modelling row effect parameters
    for(row in 1:n.row) { theta.row[row] ~ dnorm(0, tau.row)}
    
    # Modelling column effect parameters
    for(col in 1:n.col) { theta.col[col] ~ dnorm(0, tau.col)}
    
    # Modelling green manure effect parameters
    for(green in 1:n.green) { theta.green[green] ~ dnorm(0, tau.green)}
    
    # Modelling plot effect parameters
    for(plot in 1:n.plots) { theta.plots[plots] ~ dnorm(0, tau.plots)}
   
    # Modelling tillage effect parameters
    for(til in 1:n.til) { theta.til[til] ~ dnorm(0, tau.till)}
    
    # Modelling Row-Tillage interaction effect parameters
    for(row in 1:n.row) { for( til in 1:n.til) { theta.rt[row, til] ~ dnorm(0, tau.rt)}}
        
    # Modelling Column-Tillage interaction effect parameters
    for(col in 1:n.col) { for( til in 1:n.til) { theta.ct[col, til] ~ dnorm(0, tau.ct)}}
    
    # Modelling Green Manure-Tillage interaction effect parameters
    for(green in 1:n.green) { for( til in 1:n.til) { theta.gt[green, til] ~ dnorm(0, tau.gt)}}
    
    # Modelling Plot-Tillage interaction effect parameters
    for(plots in 1:n.plots) { for( til in 1:n.til) { theta.pt[plots, til] ~ dnorm(0, tau.pt)}}
    
        
    # Prior on variance components
    
    alpha ~dnorm(0,0.001)
    tau.res <- tau.de
    
    tau.row <- 1/(sigma.row*sigma.row)
    sigma.row ~ dunif(0,100)
    
    tau.col <- 1/(sigma.col*sigma.col)
    sigma.col ~ dunif(0,100)
    
    tau.green <- 1/(sigma.green*sigma.green)
    sigma.green ~ dunif(0,100)
    
    tau.plots <- 1/(sigma.plots*sigma.plots)
    sigma.np ~ dunif(0,100)
    
    tau.til <- 1/(sigma.til*sigma.til)
    sigma.til ~ dunif(0,100)
    
    tau.ct <- 1/(sigma.ct*sigma.ct)
    sigma.ct ~ dunif(0,100)
    
    tau.rt <- 1/(sigma.rt*sigma.rt)
    sigma.rt ~ dunif(0,100)
    
    tau.gt <- 1/(sigma.gt*sigma.gt)
    sigma.gt ~ dunif(0,100)
    
    tau.pt <- 1/(sigma.pt*sigma.pt)
    sigma.pt ~ dunif(0,100)
    
    # finite population standard deviation
    for( i in 1:50) {
    e.y[i] <- y[i] - mu[i] # residues at the lowest hierarchical level.
    }
    s.y <- sd(e.y[])
    s.theta.row<- sd(theta.row[])
    s.theta.col <- sd(theta.col[])
    s.theta.green <- sd(theta.green[])
    s.theta.plots <- sd(theta.plots[])
    s.theta.til <- sd(theta.til[,])
    s.theta.ct <- sd(theta.ct[,])
    s.theta.rt <- sd(theta.rt[,])
    s.theta.gt <- sd(theta.gt[,])
    s.theta.pt <- sd(theta.pt[,])

    # green manure effects size (magnitude of effec from the estimated grand mean)
    for (i in 1:n.green){ g.green[i] <- theta.green[i]-mean(theta.green[])}
        
    # tillage effects
    for (i in 1:n.til){ g.til[i] <- theta.til[i]-mean(theta.til[])}
        
    # Green Manure - Tillage effect size
    for (i in 1:n.green){for(j in 1:n.til){ g.gt[i,j] <- theta.gt[i,j]-mean(theta.gt[,])
        
   }
 }
     # End of the model
    
    }",file="VARPART_green_manure_model.txt")
    
    
    
#####################################
# Bundle data
win.data <- list(y=as.numeric(y), n.row=n.row, n.col=n.col,n.plots=n.plots, n.treat=n.treat,n.ct=n.ct, ROW=ROW, COL=COL, TREAT=TREAT, CT=CT, PLOT=PLOT)


# Inits function
inits <- function(){ list( alpha = rnorm(1),
                           sigma.row= rlnorm(1),
                           sigma.col = rlnorm(1),
                           sigma.green = rlnorm(1),
                           sigma.plots = rlnorm(1),
                           sigma.til = rlnorm(1),
                           sigma.gt = rlnorm(1),
                           sigma.rt = rlnorm(1),
                           sigma.ct = rlnorm(1),
                           sigma.pt = rlnorm(1))}
# Parameters to estimate
params <- c("g.gree","g.til","g.gt","s.theta.green",
"s.theta.col","s.theta.row","s.theta.plots",
"s.theta.til","s.theta.gt","s.theta.rt",
"s.theta.ct","s.theta.pt","theta.green",
"theta.row","theta.col","theta.plots",
"theta.til","theta.gt","theta.rt",
"theta.ct","theta.pt")

params
# MCMC settings
ni <- 600000
nb <- 550000
nt <- 30
nc <- 3

## Run the MCMC simulation
data.model <- jags.model(file="VARPART_green_manure_model.txt", data=win.data ,inits = inits,
                         n.chains=nc, n.adapt=nb)

## Compile MCMC samples, monitor all parameters,
## hyperparameters, and replicated values
data.samples <- coda.samples(data.model,variable.names=params, n.iter=ni, thin=nt, n.adapt=nb)  

mcmcplot(data.samples)
gelman.diag(data.samples,multivariate=F)



## End of the code
