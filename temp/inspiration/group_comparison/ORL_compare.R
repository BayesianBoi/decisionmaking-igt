install.packages("pacman")
pacman::p_load(R2jags, parallel, polspline)

set.seed(1983)

setwd('/work/Module4/group_comparison')

# defining a function for calculating the maximum of the posterior density (not exactly the same as the mode)
MPD <- function(x) {
  density(x)$x[which(density(x)$y==max(density(x)$y))]
}

#load control data
ctr_data <- read.table("data/IGTdata_healthy_control.txt",header=TRUE)
opi_data <- read.table("data/IGTdata_heroin.txt",header=TRUE)
amp_data <- read.table("data/IGTdata_amphetamine.txt",header=TRUE)

#----------prepare data for jags models - want trial x subject arrays for choice, gain, and loss ----
# identify and count unique subject IDs
subIDs_ctr <- unique(ctr_data$subjID)
nsubs_ctr <- length(subIDs_ctr)

subIDs_opi <- unique(opi_data$subjID)
nsubs_opi <- length(subIDs_opi)

subIDs_amp <- unique(amp_data$subjID)
nsubs_amp <- length(subIDs_amp)

ntrials_max <- 100

# SCALE THE REWARDS DOWN SO WE DON'T NEED SUCH WIDE PRIORS
# all choices (x) and outcomes (X)
x_raw_ctr <- ctr_data$deck
X_raw_ctr <- (ctr_data$gain + ctr_data$loss)/100 #note the sign!!!!!

x_raw_opi <- opi_data$deck
X_raw_opi <- (opi_data$gain + opi_data$loss)/100 #note the sign!!!!!

x_raw_amp <- amp_data$deck
X_raw_amp <- (amp_data$gain + amp_data$loss)/100 #note the sign!!!!!

#--- assign choices and outcomes in trial x sub matrix

#different number of trials across subjects. We'll need to fix this by padding arrays of < 100
#this is just so we can make the array
#then we'll also need to record number of valid trials for each sub, 
#then run the JAGS model on only valid trials

# empty arrays to fill
ntrials_ctr <- array(0,c(nsubs_ctr))
x_ctr <- array(0,c(nsubs_ctr,ntrials_max))
X_ctr <- array(0,c(nsubs_ctr,ntrials_max))

ntrials_opi <- array(0,c(nsubs_opi))
x_opi <- array(0,c(nsubs_opi,ntrials_max))
X_opi <- array(0,c(nsubs_opi,ntrials_max))

ntrials_amp <- array(0,c(nsubs_amp))
x_amp <- array(0,c(nsubs_amp,ntrials_max))
X_amp <- array(0,c(nsubs_amp,ntrials_max))

# make control data matrices
for (s in 1:nsubs_ctr) {
  
  #record n trials for subject s
  ntrials_ctr[s] <- length(x_raw_ctr[ctr_data$subjID==subIDs_ctr[s]])
  
  #pad trials with NA if n trials < maximum (i.e. 100)
  x_sub_ctr <- x_raw_ctr[ctr_data$subjID==subIDs_ctr[s]] 
  length(x_sub_ctr) <- ntrials_max
  
  X_sub_ctr <- X_raw_ctr[ctr_data$subjID==subIDs_ctr[s]] 
  length(X_sub_ctr) <- ntrials_max
  
  # assign arrays
  x_ctr[s,] <- x_sub_ctr
  X_ctr[s,] <- X_sub_ctr
  
}

# make opioid data matrices
for (s in 1:nsubs_opi) {
  
  #record n trials for subject s
  ntrials_opi[s] <- length(x_raw_opi[opi_data$subjID==subIDs_opi[s]])
  
  #pad trials with NA if n trials < maximum (i.e. 100)
  x_sub_opi <- x_raw_opi[opi_data$subjID==subIDs_opi[s]] 
  length(x_sub_opi) <- ntrials_max
  
  X_sub_opi <- X_raw_opi[opi_data$subjID==subIDs_opi[s]] 
  length(X_sub_opi) <- ntrials_max
  
  # assign arrays
  x_opi[s,] <- x_sub_opi
  X_opi[s,] <- X_sub_opi
  
}

# make amph data matrices
for (s in 1:nsubs_amp) {
  
  #record n trials for subject s
  ntrials_amp[s] <- length(x_raw_amp[amp_data$subjID==subIDs_amp[s]])
  
  #pad trials with NA if n trials < maximum (i.e. 100)
  x_sub_amp <- x_raw_amp[amp_data$subjID==subIDs_amp[s]] 
  length(x_sub_amp) <- ntrials_max
  
  X_sub_amp <- X_raw_amp[amp_data$subjID==subIDs_amp[s]] 
  length(X_sub_amp) <- ntrials_max
  
  # assign arrays
  x_amp[s,] <- x_sub_amp
  X_amp[s,] <- X_sub_amp
  
}


X <- c(-0.5, 0,5)

X <- (X-mean(X))/sd(X)

invSigma <- solve(t(X)%*%X) # required for JZS priors

ngroups=length(X)

# set up jags and run jags model on one subject
data <- list("x_ctr","X_ctr","ntrials_ctr","nsubs_ctr",
             "x_opi","X_opi","ntrials_opi","nsubs_opi", "X", "invSigma", "ngroups") 
params<-c("alpha_a_rew","alpha_a_pun","alpha_K","alpha_omega_f","alpha_omega_p") #,
#"a_rew_ctr","a_pun_ctr","K_ctr","omega_f_ctr","omega_p_ctr",
#"a_rew_opi","a_pun_opi","K_opi","omega_f_opi","omega_p_opi")
start_time = Sys.time()

samples <- jags.parallel(data, inits=NULL, params,
                model.file ="ORL_compare.txt",
                n.chains=3, n.iter=5000, n.burnin=1000, n.thin=1, n.cluster=4)

end_time = Sys.time()
end_time - start_time

BF_effect <- NULL
BF_null <- NULL

# savage dickey plot
plot(density(rnorm(12000,0,1/sqrt(1))),ylim=c(0,1),main=" ")
lines(density(samples$BUGSoutput$sims.list$alpha_a_rew),col="red")

fit.posterior <- logspline(samples$BUGSoutput$sims.list$alpha_a_rew)
null.posterior <- dlogspline(0, fit.posterior)

null.prior     <- dnorm(0,0,(1/sqrt(1)))                   

BF_effect$alpha_a_rew <- null.prior/null.posterior
BF_null$alpha_a_rew <- null.posterior/null.prior

# savage dickey plot
plot(density(rnorm(12000,0,1/sqrt(1))),ylim=c(0,1),main=" ")
lines(density(samples$BUGSoutput$sims.list$alpha_omega_f),col="red")

fit.posterior <- logspline(samples$BUGSoutput$sims.list$alpha_omega_f)
null.posterior <- dlogspline(0, fit.posterior)

null.prior     <- dnorm(0,0,(1/sqrt(1)))                   

BF_effect$alpha_omega_f <- null.prior/null.posterior
BF_null$alpha_omega_f <- null.posterior/null.prior



par(mfrow=c(3,2))
plot(density(rnorm(10000,0,1/sqrt(1))),ylim=c(0,.7),main="a_rew")
lines(density(samples$BUGSoutput$sims.list$alpha_a_rew),col="red")

plot(density(rnorm(10000,0,1/sqrt(1))),ylim=c(0,.7),main="a_pun")
lines(density(samples$BUGSoutput$sims.list$alpha_a_pun),col="red")

plot(density(rnorm(10000,0,1/sqrt(1))),ylim=c(0,.7),main="omega_f")
lines(density(samples$BUGSoutput$sims.list$alpha_omega_f),col="red")

plot(density(rnorm(10000,0,1/sqrt(1))),ylim=c(0,.7),main="omega_p")
lines(density(samples$BUGSoutput$sims.list$alpha_omega_p),col="red")

plot(density(rnorm(10000,0,1/sqrt(1))),ylim=c(0,.7),main="K")
lines(density(samples$BUGSoutput$sims.list$alpha_K),col="red")















############### Look at cumulative balance to justify learning ##############

# calculate cumulative value
xcum_ctr <- array(0,c(nsubs_ctr,100))
for (i in 1:nsubs_ctr) {
  xcum_ctr[i,] <- cumsum(X_ctr[i,]) 
}

xcum_opi <- array(0,c(nsubs_opi,100))
for (i in 1:nsubs_opi) {
  xcum_opi[i,] <- cumsum(X_opi[i,]) 
}

xcum_amp <- array(0,c(nsubs_amp,100))
for (i in 1:nsubs_amp) {
  xcum_amp[i,] <- cumsum(X_amp[i,]) 
}

plot(colMeans(xcum_ctr),ylim=c(-15,8),type='l',lwd=2,  
     xlab = "Trial", ylab = "Cumulative Balance") # should be -1500, 800 if not scaled
lines(colMeans(xcum_opi),col="red",lwd=2)
lines(colMeans(xcum_amp),col="blue",lwd=2)
legend(x='topright', legend=c('healthy controls', 'heroine', 'ampetamine'))

############### compare final means ##############

# set up jags and run jags model on one subject
data <- list("xcum_ctr","nsubs_ctr",
             "xcum_opi","nsubs_opi") 
params<-c("mu","alpha","Smu_ctr","Smu_opi")
temp_samples <- jags.parallel(data, inits=NULL, params,
                model.file ="balance_compare.txt",
                n.chains=3, n.iter=15000, n.burnin=1000, n.thin=3, n.cluster=4)


# savage dickey plot
plot(density(rnorm(12000,0,1/sqrt(1))),main=" ")
lines(density(temp_samples$BUGSoutput$sims.list$alpha),col="red")

fit.posterior <- logspline(temp_samples$BUGSoutput$sims.list$alpha)
null.posterior <- dlogspline(0, fit.posterior)

null.prior     <- dnorm(0,0,(1/sqrt(1)))                   

BF <- null.prior/null.posterior
BF

# savage dickey plot
plot(density(rnorm(12000,0,1/sqrt(1))),main="mu")
lines(density(temp_samples$BUGSoutput$sims.list$mu),col="red")

fit.posterior <- logspline(temp_samples$BUGSoutput$sims.list$mu)
null.posterior <- dlogspline(0, fit.posterior)

null.prior     <- dnorm(0,0,(1/sqrt(1)))                   

BF <- null.prior/null.posterior
BF


