# Introduction

Here we are closely following the ideas from Gelman and Hill (2007) where they state that all variability sources should be
treated as random-effects. In other words, we are assuming that the theta effects of a given variability source are a sample
from a bigger population of effects and their estimates are not independent from each other. Gelman and Hill (2007) use the 
terms varying intercepts and varying slopes to overcome the yet unclear and debatable definitions of fixed and random effects.
The model below is a varying-intercept model in the sense that 1) effects are left to vary within a given categorical source of 
variability; 2) the effects within a source of variability are assumed to come from a commom probability distribution; 3) the 
estimate of effects within a source of variability is influenced by the sample size of the source (partial pooling effects).

# Experimental context
Mixed plantation experiment between Eucalyptus urograndis and Acacia Mangium (N-fixing tree). Measurament of nitrous oxide collected from 
chamber over eleveen months for times a day. The stands are: 1) A100 (pure plantion of Acacia Mangium), 2) E100 (Pure plantion of Eucalyptus urograndis), and 3) A50E50 ( mixed plantion between Acacia and Eucalyptus).

# Building the model         
         cat("model{
      
         # General parameters
         n.time = 12
         n.sampling = 4
         n.rain = 50
         n.treat = 5
         n.block = 2
         
         # Likelihood model
         y[i] ~ N( mu[i], tau.res)  # Stochastic (random) part 
         log(mu[i]) <- exp(log.mu[i]) # Link function
         
         # linear predictor
         log.mu[i] <- alpha +
                      beta.time[TIME[i] + 
                      beta.sampling[SAMPLING[i]] + 
                      beta.treat[TREAT[i]] + 
                      beta.block[BLOCK[i]] + 
                      beta.rain[RAINFALL[i]] + 
                      beta.tb[TREAT[i],BLOCK[i]] + 
                      beta.tt[TREAT[i],TIME[i]]  
                             

           # Prior for parameters:

           # intercept fixed effect
            alpha ~ dnorm(0.001,0.001) 

           # random time (month) effects
           for( i in 1:n.time) { beta.time[i] ~ dnorm(0, tau.time)}

           # random sampling effects
           for( i in 1:n.sampling) {beta.sampling[i] ~ dnorm(0, tau.sampling)}

           # random rainfall effects
           for( i in 1:n.rain) {beta.rain[i] ~ dnorm(0, tau.rain)}

           # random treatment effects
           for( i in 1:n.treat) {beta.treat[i] ~ dnorm(0, tau.treat)}
 
           # random block effects 
           for( i in 1:n.block) {beta.block[i] ~ dnorm (sigma.block*sigma.block)}

           # random treatment-time effects
           for( i in 1:n.time) {for( j in j:n.treat) { beta.tt[i,j] ~ dnorm (0, tau.tt)}}

           # random treatment-block effects
           for( i in 1:n.block) {for( j in j:n.treat) {beta.tb[i,j] ~ dnorm (0, tau.tb)}}

           # Prior for variance components
           tau.res <- 1/(sigma.res*sigma.res)           
           sigma.res ~ dunif(0,100)

           tau.time <- 1/(sigma.time*sigma.time)
           sigma.time ~dunif(0,100)

           tau. sampling <- 1/(sigma. sampling *sigma. sampling)
           sigma. sampling ~dunif(0,100)

           tau.treat <- 1/(sigma.treat*sigma.treat)
           sigma.treat ~dunif(0,100)

           tau.block <- 1/(sigma.block*sigma.block)
           sigma.block ~dunif(0,100)

           tau.tt <- 1/(sigma.tt*sigma.tt)
           sigma.tt ~dunif(0,100)

           tau.tb <- 1/(sigma.tb*sigma.tb)
           sigma.tb ~dunif(0,100)

        # finite population standard deviations for all variance sources
          for( i in 1:N) {y.res[i] <- y[i] - mu[i]} # calculates the residuals at individual observation level i
          sd.res <- sd(y.res[]) 
          sd.time <- sd(beta.time[])
          sd.sampling <- sd(beta.sampling[])
          sd.treat <- sd(beta.treat[]) 
          sd.block <- sd(beta.block[])
          sd.rain <- sd(beta.rain[])
          sd.tt <- sd(beta.tt[,])
          sd.tb <- sd(beta.tb[,])

        # Estimating MAIN effects within variability sources (derived quantities)
        # main treatment effects
        for (i in 1:n.treat){g.treat[i] <- beta.treat[i]-mean(beta.treat[])}
        
        # main time effects
         for (i in 1:n.treat){g.time[i] <- beta.time[i]-mean(beta.time[])}
        
        # main block effects
        for (i in 1:n.block){g.block[i] <- beta.block[i]-mean(beta.block[])}}

        # main rainfall effects
        for (i in 1:n.rain){g.rainfall[i] <- beta.rain[i]-mean(beta.rain[])}
        
        # main treatment-time effects
        for (i in 1:n.time){for(j in 1:n.treat){g.tt[i,j] <- beta.tt[i,j]-mean(beta.tt[,])}} 

       # Specific questions:
       # At the what extent are the N2O emissions from pure plantations (A100 vs E100) different from each other?
         diff.AcEuc <- g.treat[1] – g.treat[2]
         prob.diff.AcEuc <- step(diff.AcEuc)

       # At the what extent are the N2O emissions from E100 and A50E50 different from each other?
        diff.EucMix <- g.treat[2] – g.treat[3]
        prob.diff.EucMix <- step(diff.EucMix)
                
       # At the what extent are the N2O emissions from A100 and A50E50 different from each other?
        diff.AcMix <- g.treat[1] – g.treat[3]
        prob.diff.AcMix <- step(diff.AcMix)

        }",file="EUCALEGmod1.txt")
    

    
       # Bundle data
       win.data <- list(y=as.numeric(y), n.time=n.time, n.sampling=n.sampling, n.rain=n.rain, n.treat=n.treat, n.tt=n.tt,n.tb=n.tb, N=N, MONTH = MONTH, SAMPLING = SAMPLING, RAINFALL = RAINFALL,TREAT=TREAT,BLOCK = BLOCK)


      # Inits function for parameters
        inits <- function(){ list(alpha = rnorm(1), 
                          beta.time= rnorm(12,0,1),
                          beta.sampling = rnorm(8,0,1),
                          beta.rain= rlnorm(2,0,1),
                          beta.treat = rnorm(3,0,1),
                          beta.block = rnorm(3,0,1),
                          sigma.time= rlnorm(1),
                          sigma.sampling = rlnorm(1),
                          sigma.rain= rlnorm(1),
                          sigma.treat = rlnorm(1),
                          sigma.block = rlnorm(1),
                          sigma.tt = rlnorm(1),
                          sigma.tb = rlnorm(1))}
           

        # Parameters to monitor (which is going to be displayed)
         params <- c(“sd.time”, “sd.sampling”, “sd.rain”, “sd.treat”, “sd.block”, “sd.tt”, “sd.tb”, “sd.res”, "g.time","g.rain”, “g.treat”, "g.block", “g.tt", “diff.AcMix”, “diff.EucMix”, “diff.AcEuc”, “prob.diff.AcMix”, “prob.diff.EucMix”, “prob.diff.AcEuc”)

        # MCMC settings
        ni <- 50000
        nb <- 20000
        nt <- 25
        nc <- 3

        # Run the MCMC simulation
        data.model <- jags.model(file="EUCALEGmod1.txt", data=win.data ,inits = inits,n.chains=nc, n.adapt=nb)
        # Compile MCMC samples, monitor all parameters,
        # hyperparameters, and replicated values
        data.samples <- coda.samples(data.model,variable.names = params, n.iter=ni, thin=nt, n.adapt=nb)  

        # gelman diagnostics
          gelman.diag(data.samples,multivariate=F)
