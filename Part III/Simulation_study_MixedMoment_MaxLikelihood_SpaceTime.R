###################################################
#            Complete Simulation Study            #
#               MM vs. ML estimators              #
#   cf. Chapter 6 Section 6.4 Silva Lomba (2023)  #
# Based and adapted from public code by Chen Zhou #
#   https://github.com/zhouchen0527/rainscedasis  #
#             cf. Einmahl et al. (2022)           #
###################################################


################# PREAMBLE #################


# Packages
library(kolmim)
library(reshape)
library(ggplot2)
library(EnvStats)
library(eva)
library(tictoc)
library(foreach)
library(doRNG)
library(doParallel)

# seed for data generation
set.seed(663)


## Parameter specification

# Simulation setup
# The number of stations
m <- 12
# The number of observations over time
n <- 1350
# The number of simulations
sim <- 500

## Model setup

# Dependence type
#dependence<-"indep"; beta<- 0; nBM <- 0
dependence<-"dep"; beta<- 0.05; nBM <- 28
# Number of simulated BM - depends on result of determination below

# Scedasis type
#scedasis<-"stationary"; stationarity<-0
scedasis<- "sin_scedasis"; vos<-0.4; vot<-0.5; shift<-0.5; stationarity<-1
#scedasis<-"linear_scedasis"; stationarity<-2
#scedasis<-"exp_scedasis"; vt<-5; stationarity<-3

# Common EVI 
xi <- 0.05
#xi <- 0.5


################# FUNCTIONS #################

## Simulation

# Simulate Brownian motion with parameter beta at location 1:m
simBM <- function(beta,m,nsim){
  eps <- rnorm(m*nsim, mean=0, sd=sqrt(beta))
  dim(eps) <- c(m,nsim)
  results <- apply(eps, 2, cumsum)
  return(results)
}

# Simulate a Brown-Resnick process with fixed number of BM
simBR<- function(beta,m,nBM){
  BMs <- simBM(beta,m,nBM)
  BMs <- BMs-(1:m)%*%t(rep(1,nBM))*beta/2
  expo <- rexp(nBM)
  cexp <- 1/cumsum(expo)
  results <- exp(BMs)%*%diag(cexp)
  results <- apply(results, 1, max)
  return(results)
}

# Pre-simulate Brown-Resnick processes
# To determine the number of Brownian motion needed
# For any beta and m, always run this functon first
# The optimal nBM should be either the maximum of the results
# or a very high quantile of it
presimBR <- function (beta, m, ntest = 1000, start = 5, end = 30) {
  results <- c()
  for(j in 1:ntest){
    BMs <- simBM(beta,m,end)
    BMs <- BMs-(1:m)%*%t(rep(1,end))*beta/2
    expo <- rexp(end)
    cexp <- 1/cumsum(expo)
    temp <- exp(BMs)%*%diag(cexp)
    tempbr <- apply(temp,1,max)
    temppos <- start:end
    for(l in temppos){
      partial <- apply(temp[,1:l],1,max)
      if(sum(partial!=tempbr)==0) break
    }
    results=c(results,l)
  }
  return(results)
}

summary(presimBR(beta,m))



# Generating scedasis function 

cfun <- function(m,n,stationarity){
  results <- rep(0,m*n)
  dim(results) <- c(n,m)
  if(stationarity==1){
    Cfun <- 1/m*(1+ vos*sin((1:m)/m*(2*pi)))
    for(j in 1:m){
      results[,j]=m*Cfun[j]*(1+vot*sin((1:n)/n*(2*pi)+shift*j))
    }
  }else if(stationarity==2){
    Cfun <- 2/(m*(m+3))*(1+ (1:m)) #not needed
    for(j in 1:m){
      results[,j]=2*(j+2*(1:n)/n)/(3+m)
    }
  }else if(stationarity==3){
    for(j in 1:m){
      results[,j]=2*(vt*j*exp(vt*(1:n)/n))/((exp(vt)-1)*(1+m))
    }
  }  
  else{
    for(j in 1:m){
      results[,j]=1
    }
  }
  return(results)
}



## Data generation
# Data structure, n rows, m columns -> vector

# Dependence (BR process)
gendata <- function(m,n,stationarity, xi, beta, nBM){
  skedasis <- cfun(m,n,stationarity)
  brprocess <- replicate(n, simBR(beta,m,nBM))
  results <- (skedasis*t(brprocess))^xi
  return (as.vector(results))
}

# Independence (GPd)
gendata_indep <- function(m,n,stationarity,xi){
  skedasis <- cfun(m,n,stationarity)
  paretos <- matrix(EnvStats::rpareto(n*m, location=1, shape = 1),n,m)
  results <- ((paretos*skedasis)^xi-1)/xi
  return (as.vector(results))
}



## Estimation

#Mixed Moment estimates for sample fractions from kmin to kmax
Mix_Mom_Est<-function(aggregated_data, kmin,kmax){
  
  X<-sort(aggregated_data)  #ascending order
  N<-length(X)
  K<-kmin:kmax     #Number of upper order statistics
  n_est<-length(K) #Number of estimates 
  
  MjN<-matrix(NA,2,n_est)
  L1N<-numeric(n_est)
  phi_hat<-numeric(n_est)
  xi_hat<-numeric(n_est)
  
  #Computing MjN(k)
  for(j in 1:2){            
    for(k in K){
      MjN[j,(k-kmin+1)]<-(1/k)*sum(log(X[(N-k+1):N]/X[N-k])^j)
    }
  }
  
  #Computing L1N(k), phi_hat(k) and xi_hat(k)
  for(k in K){ 
    L1N[(k-kmin+1)]<-(1/k)*sum(1-rep(X[N-k],k)/X[(N-k+1):N])
    phi_hat[(k-kmin+1)]<-(MjN[1,(k-kmin+1)]-L1N[(k-kmin+1)])/(L1N[(k-kmin+1)]^2)
    xi_hat[(k-kmin+1)]<-(phi_hat[(k-kmin+1)]-1)/(1+2*min(phi_hat[(k-kmin+1)]-1,0))
  }
  
  return(list(K=K,Phi=phi_hat,xi_MM=xi_hat))
}



# Process for estimation either by MM or ML estimators
onesampleEST <- function(data, kmin,kmax, estimator){
  estimates<-NULL
  if(estimator=="MMe"){
    estimates <-  Mix_Mom_Est(data,kmin,kmax)$xi_MM}
  else{
    for (k in kmin:kmax){
      estimates<-c(estimates,eva::gpdFit(data, nextremes=k, method = c("mle"), start=NULL, opt = "Nelder-Mead", maxit = 10000)$par.ests[2])
    }
  }
  return(estimates)
}





################# DATA GENERATION #################

samples<-matrix(NA,sim,m*n)

# independence samples
tic()
for(i in 1:sim){
  samples[i,] = gendata_indep(m,n,stationarity, xi)
}
toc()
beepr::beep(4)


# dependence samples
tic()
for(i in 1:sim){
  samples[i,] = gendata(m,n,stationarity, xi, beta, nBM)
}
toc()
beepr::beep(4)



################# EVI ESTIMATION #################
#----------- with parallel computation ----------#

# sample fraction
kmin=10
kmax=10000


est<-c()

# Choice of estimator
estimator<-"MMe"
#estimator<-"MLe"


#setup parallel backend to use all processors
cores=detectCores()
cl <- makeCluster(cores[1]-1) #not to overload your computer
registerDoParallel(cl)


tic()
est <- foreach(i=1:sim, .combine=rbind) %dopar% {
  est_temp = onesampleEST(samples[i,],kmin,kmax,estimator)
  est_temp
}
toc()
beepr::beep(4)

#stop cluster
stopCluster(cl)

dtf<-data.frame(est)
colnames(dtf)<-c(kmin:kmax)


write.csv(dtf,paste0("Sim_",dependence,"_",scedasis,"_",estimator,"_",xi,"_FULL.csv"),row.names = F)



###################################################
#            Complete Simulation Study            #
#               MM vs. ML estimators              #
#   cf. Chapter 6 Section 6.4 Silva Lomba (2023)  #
# Based and adapted from public code by Chen Zhou #
#   https://github.com/zhouchen0527/rainscedasis  #
#             cf. Einmahl et al. (2022)           #
###################################################
