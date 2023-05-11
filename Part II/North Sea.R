###################################################################
######################### North Sea Data ########################## 
############# 628 storm peak significant wave heights ############# 
########################### ONE LOCATION ##########################
#################### October 1964 - March 1995 ####################
################## Only periods October - March ###################
#################### From Northrop et al 2017 ##################### 
###################################################################


# Packages
library(evir)
library(evd)
library(threshr)
library(eva)
library(spatialEco)
library(SpatialExtremes)

# import L-moments estimators
source('Sample L-Moments.R')
# import source code for STSM - Northrop and Coleman (2014) 
source("NorthropColeman2014.fns")

# load data
data('ns')
x<-ns

# sample size
n<-length(x)


#---------------------------------------------------------  Histogram
# Figure 3.9
#####

hist(x, probability=T, breaks = 60,main="", xlab= 'significant wave height (m)',
     ylim = c(0,0.45),mgp=c(1.8, .5, 0), col='gray', density=30, angle=30,border='black')  #same as reference paper
box(lty="solid",col='black')
axis(side=3, at=quantile(ns,probs = c(.25,.50,.65,.75,.80,.85,.90,.95,.975,.99)), 
     labels=c(25,50,65,75,80,85,90,95,97.5,99)
     ,tick=T, mgp=c(3, .5, 0),cex.axis=0.7)

# estimated density
lines(density(x), col='red')
#####


#---------------------------------------------------------  MRL plot
# Langousis et al. (2016)
# Figure A.6 + Table 3.5
#####


# ascending order statistics (OS)
u<-x[1:(n-10)]
# u[n-10]=8.634

# mean excess function at OS
me_u<-c()
for(i in 1:(n-10)){
  me_u[i]<-mean(x[x>u[i]]-u[i])
}

# candidates
u_star=seq(1.652,8.552,by=0.1)   #length(u_star)=70 

# mean excess function at candidates
me_u_star<-c()
for(i in 1:(length(u_star))){
  me_u_star[i]<-mean(x[x>u_star[i]]-u_star[i])
}

# weighted least squares
w_star<-c()
for(j in 1:(length(u_star))){
  w_star[j]<-(length(x[x>u_star[j]]))/var(x[x>u_star[j]]-u_star[j])
}

wls_fit_star<-vector(mode='list', length = length(u_star))
wls_WMSE_star<-c()

for(j in 1:(length(u_star))){
  wls_fit_star[[j]]<-lm(me_u_star[j:(length(u_star))] ~ u_star[j:(length(u_star))],
                        weights=w_star[j:(length(u_star))])
  
  wls_WMSE_star[j]<-sum(w_star[j:(length(u_star))]*(wls_fit_star[[j]]$residuals)^2
  )/sum(w_star[j:(length(u_star))])
}

# choice of threshold
min_WMSE_star<-local.min.max(wls_WMSE_star, dev = mean, plot = F)$minima[1]
which(round(wls_WMSE_star,9)==round(min_WMSE_star,9))



# plot
plot(u,me_u,col='blue',xlim=c(1,9),ylim=c(0.6,2.2),cex=1.1,lwd=0.5,
     xaxt='n',yaxt='n', xlab=expression(x[i:n]), ylab = 'Mean Excess',mgp=c(1.8,0, 0))
axis(side=2, at=seq(0.6,2.2,by=0.2), labels=seq(0.6,2.2,by=0.2),tick=T,mgp=c(1.8, .5, 0),cex.axis=0.7)
axis(side=1, at=seq(1.652,8.552,by=0.1), labels=NA,tick=T, col='black',
     mgp=c(3, .5, 0),cex.axis=0.7, lty=3)
axis(side=1, at=seq(1,9,by=1), labels=seq(1,9,by=1),tick=T,mgp=c(1.8, .5, 0),cex.axis=0.7)
axis(side=3, at=quantile(x,probs = c(.25,.50,.65,.80,.90,.95)), 
     labels=c(25,50,65,80,90,95),tick=T, mgp=c(1.8, .5, 0),cex.axis=0.6)
title('Mean Excess Plots for North Sea Data')

points(u_star,me_u_star,pch=4, cex=0.9, lwd=2)

y0_star<-wls_fit_star[[4]]$coefficients[1]
m_star<-wls_fit_star[[4]]$coefficients[2]
segments(u_star[4], u_star[4]*m_star+y0_star, x1 = u_star[70], y1 = u_star[70]*m_star+y0_star,
         lwd=3)


# Parameter estimates - 1st local min WMSE - BY PLOT LINE ADJUSTMENT
u_opt_WLS2<-u_star[4]
Nu_WLS2<-length(x[x>u_opt_WLS2])

fit_WLS_xi2<-m_star/(1+m_star)   #xi = A/(1+A) = -0.2468237
fit_WLS_sigma2<-y0_star*(1-fit_WLS_xi2)+fit_WLS_xi2*u_opt_WLS2  #a_u*=B(1-xi)+xi*u =2.706154

#100, 1000, 1000 year return levels (20.26 observations per year)
Q_100_WLS2<- u_opt_WLS2 + fit_WLS_sigma2/fit_WLS_xi2 * (((n/Nu_WLS2)*(1/2026))^(-fit_WLS_xi2)-1)
Q_1000_WLS2<- u_opt_WLS2 + fit_WLS_sigma2/fit_WLS_xi2 * (((n/Nu_WLS)*(1/20260))^(-fit_WLS_xi2)-1)
Q_10000_WLS2<- u_opt_WLS2 + fit_WLS_sigma2/fit_WLS_xi2 * (((n/Nu_WLS)*(1/202600))^(-fit_WLS_xi2)-1)


# Parameter estimates - 1st local min WMSE -- MAXIMUM LIKELIHOOD
u_opt_WLS<-u_star[4]
Nu_WLS<-length(x[x>u_opt_WLS])
fit_WLS <- evir::gpd(x,u_opt_WLS)
fit_WLS_xi<-fit_WLS$par.ests[[1]]
fit_WLS_sigma<-fit_WLS$par.ests[[2]]

#100, 1000, 1000 year return levels (20.26 observations per year)
Q_100_WLS<- u_opt_WLS + fit_WLS_sigma/fit_WLS_xi * (((n/Nu_WLS)*(1/2026))^(-fit_WLS_xi)-1)
Q_1000_WLS<- u_opt_WLS + fit_WLS_sigma/fit_WLS_xi * (((n/Nu_WLS)*(1/20260))^(-fit_WLS_xi)-1)
Q_10000_WLS<- u_opt_WLS + fit_WLS_sigma/fit_WLS_xi * (((n/Nu_WLS)*(1/202600))^(-fit_WLS_xi)-1)
#####



#-------------------------------------------------------- GP QQ plot
# Table 3.5
#####


# to test correlation in the GP qq-plots above each threshold, 
# estimation of xi is needed
xi_fit_k<-c()
qqplot_cor<-c()

# fitting the GPd and computing correlation
for(k in 1:(n-10)){
  fit_gpd<- evir::gpd(x,sort(x)[k])
  xi_fit_k[k]<-fit_gpd$par.ests[[1]]
  
  plotpos<-(1:(n-k))/(n-k+1)
  qqplot_cor[k]<-cor(evir::qgpd(plotpos,xi_fit_k[k]), sort(x)[(k+1):n]-sort(x)[k])
}

plot(sort(x)[1:(n-10)],qqplot_cor, main= 'GP QQ-plot correlations North Sea', 
     xlab=bquote('Threshold '(x[k:n])), ylab='Correlation', xaxt='n')
axis(1,at=sort(x)[-c(n:(n-9))], labels=sort(x)[-c(n:(n-9))],cex.axis=1, mgp=c(3, .5, 0))

# choosing threshold
k_qqplot<-which.max(qqplot_cor)
u_qqplot<-sort(x)[which.max(qqplot_cor)]

abline(h=max(qqplot_cor), lty=2, col='red')
abline(v=u_qqplot, lty=2, col='red')
text(1.6,1,round(max(qqplot_cor),3), col='red', cex=0.8)
text(4.5,0.955,expression('   4.3 \n (k=448)'), col='red', cex=0.8)



# plotting 'best' GP QQ-plot
plotpos_qqplot<-(1:(n-k_qqplot))/(n-k_qqplot+1)
quantis_qqplot<- evir::qgpd(plotpos_qqplot,xi_fit_k[k_qqplot])

plot(quantis_qqplot,sort(x)[(k_qqplot+1):n]-sort(x)[k_qqplot], main='GP QQplot for maximum corelation North Sea - u = 4.3', 
     xlab=bquote(((1-p[i])^-hat(xi)-1)/hat(xi)), ylab=expression(x[i:n]))

fit_qqplot<-lm(sort(x)[(k_qqplot+1):n]-sort(x)[k_qqplot]~quantis_qqplot-1)
abline(fit_qqplot, col='red', lty=2)
text(2.25,6,round(max(qqplot_cor),3), col='red', cex=0.8)
text(2.25,0,bquote(hat(xi) == -0.3274747), col='red', cex=0.8)


# Parameter estimates - max correlation
u_opt_QQ<-u_qqplot
Nu_QQ<-length(x[x>u_opt_QQ])
fit_QQ <- evir::gpd(x,u_opt_QQ)
fit_QQ_xi<-fit_QQ$par.ests[[1]]
fit_QQ_sigma<-fit_QQ$par.ests[[2]]

#100, 1000, 1000 year return levels (20.26 observations per year)
Q_100_QQ<- u_opt_QQ + fit_QQ_sigma/fit_QQ_xi * (((n/Nu_QQ)*(1/2026))^(-fit_QQ_xi)-1)
Q_1000_QQ<- u_opt_QQ + fit_QQ_sigma/fit_QQ_xi * (((n/Nu_QQ)*(1/20260))^(-fit_QQ_xi)-1)
Q_10000_QQ<- u_opt_QQ + fit_QQ_sigma/fit_QQ_xi * (((n/Nu_QQ)*(1/202600))^(-fit_QQ_xi)-1)

#####



#---------------------------------------------------- Surprise Plots
# Lee et al. (2015) - base code made available by Dr. Kate J. E. Lee
# Table 3.5 + Figure 3.12
#####

# Posterior in log
lnPost <- function(y,u,xi,sigma){
  f = sum(dgpd(y[y>u],loc=u,scale=sigma,shape=xi,log=TRUE))-log(sigma)-log(1+xi)-0.5*log(1+2*xi); return(f)
}

# Metropolis-Hastings updates 
MHSim <- function(y,u,T,batch,lnPost){
  xi0=0.02; sig0 =8; # initial parameter values; xi0(shape) and sig0(scale)
  Wsig=2; Wxi=1; # initial proposal walk size for sig and xi
  lnP0=lnPost(y,u,xi0,sig0) # log-posterior for the initial parameter values
  accept=0; # number of accepted proposals 
  Para=matrix(,T,2); # MCMC outputs for parameter values [shape,scale]
  
  ### Metropolis-Hastings update ####
  for (t in 1:T){  
    sig=exp(rnorm(1,log(sig0),Wsig)); xi=rnorm(1,xi0, Wxi); xi=xi*(xi>(-0.5))+xi0*(xi<(-0.5))
    lnP=lnPost(y,u,xi,sig)
    if(runif(1)<exp(lnP-lnP0)){sig0=sig; xi0=xi; lnP0 = lnP; accept=accept+1; }
    Para[t,]=c(sig0,xi0); 
    
    # readjust proposal walk size 
    if((t/batch)==floor(t/batch)){ 
      Wsig=max(0.1,2.38*sd(log(Para[(t-batch+1):t,1])))
      Wxi=max(0.1,2.38*sd(Para[(t-batch+1):t,2]))}}
  
  out=list(Para=Para,accept); return(out)
}

# candidate thresholds
U<- quantile(x,probs = seq(0.25,0.925,by=0.075), names=F)  #10 quantiles
#U<- quantile(x,probs = seq(0.25,0.953,by=0.037), names=F)  #20 quantiles

Ni0=50; # number replicates
N=9e3; # number of posterior samples used for the p-value estimate
MS=matrix(,length(U),Ni0) #posterior predictive p-value for Ni0 replicates

T=1e4; # Total number of iterations (length of Markov chain)
batch=1e3; # the proposal walk size is adjusted at every "batch"-iteration

start.time <- Sys.time()

for(i0 in 1:Ni0){
  
  Ny=length(x)
  y=x
  
  for(i in 1:length(U)){
    # Simulate T - posterior samples 
    out=MHSim(y,U[i],T,batch,lnPost)
    
    # Posterior predictive p-value calculation
    PosP=0; 
    for(j in 1:N){ ind=sample(c(1e3:T),1) 
    # simulated data using a posterior sample
    ysim=rgpd(sum(y>U[i]),loc=U[i],scale=out$Para[ind,1],shape=out$Para[ind,2])
    # Here the test statistic is a likelihood
    Tsim=mean(dgpd(ysim,loc=U[i],scale=out$Para[ind,1],shape=out$Para[ind,2]))
    Tobs=mean(dgpd(y[y>U[i]],loc=U[i],scale=out$Para[ind,1],shape=out$Para[ind,2]))
    PosP= PosP+as.numeric(Tsim<Tobs)
    }
    MS[i,i0]=PosP/N #p-value 
  }
}

boxplot(t(MS),xaxt='n',cex.axis=0.8, mgp=c(0,0,1),cex.lab=0.8)
axis(side=1, at=c(1:10),labels=seq(25,92.5,7.5),cex.axis=0.8, mgp=c(1.8, .5, 0))
title(ylab='posterior predictive p-values', mgp=c(2.2, .5, 0))
title(xlab='candidate thresholds - sample quantiles (%)', mgp=c(2, .5, 0), cex=0.8)
abline(h=0.5, lty=2, col='red')

end.time <- Sys.time()
difftime(end.time, start.time, units='secs')

#####



#------------------------------------------------------------- aSTSM
# Northrop and Coleman (2014)
# Table 3.5 + Figure 3.13
#####

# sort sample
x<-sort(x)

## Candidate thresholds
# u<-x[c(1:(n-10))]
## testing all sample points as thresholds
## error in score.fitrange : thresholds are not in increasing order (equal observations)
# u<-c(u[which(u[1:617]-u[2:618]<0)],u[618])  #580 thresholds

#u<- quantile(x,probs = seq(0.25,0.925,by=0.075), names=F)  #10 quantiles
u<- quantile(x,probs = seq(0.25,0.953,by=0.037), names=F)  #20 quantiles


# Score test
northrop <- tryCatch(score.fitrange(x,u), error= function(e) NA) 

# automatic choice
if(length(northrop)==10){
  i=0
  aux=F
  while(i<(length(u)-1) && aux==F){
    i<-i+1
    if(northrop$e.p.values[i] > 0.05 & all(northrop$e.p.values[i:(length(u)-1)] > 0.05)){aux=TRUE}
  }
  (u_opt_N<-u[i])
  Nu_opt_N<-length(x[x>u_opt_N])
  } else{
  u_opt_N<-NA
  Nu_opt_N<-NA
}

# add info to p-values plot - Figure 3.11
abline(v=u_opt_N,col='red',lty=2)
#text(floor(u_opt_N),1,u_opt_N, col='red')
text(3.2,1,u_opt_N, col='red')

#alternative automatic selection rule
u_opt_N_alt<-u[which(northrop$e.p.values > 0.05)][1]
abline(v=u_opt_N_alt,col='blue',lty=2)
text(ceiling(u_opt_N_alt),1,u_opt_N_alt, col='blue')

# ## Parameter estimates - 10 candidates
# # A
# u_opt_N_10<-u_opt_N
# Nu_N_10<-length(x[x>u_opt_N_10])
# fit_N_10<- evir::gpd(x,u_opt_N_10)
# fit_N_10_xi<-fit_N_10$par.ests[[1]]
# fit_N_10_sigma<-fit_N_10$par.ests[[2]]
# # B
# u_opt_N_10_2<-u_opt_N_alt
# Nu_N_10_2<-length(x[x>u_opt_N_10_2])
# fit_N_10_2<- evir::gpd(x,u_opt_N_10_2)
# fit_N_10_2_xi<-fit_N_10_2$par.ests[[1]]
# fit_N_10_2_sigma<-fit_N_10_2$par.ests[[2]]
# 
# # 100, 1000, 1000 year return levels (20.26 observations per year)
# # A
# Q_100_N_10<- u_opt_N_10 + fit_N_10_sigma/fit_N_10_xi * (((n/Nu_N_10)*(1/2026))^(-fit_N_10_xi)-1)
# Q_1000_N_10<- u_opt_N_10 + fit_N_10_sigma/fit_N_10_xi * (((n/Nu_N_10)*(1/20260))^(-fit_N_10_xi)-1)
# Q_10000_N_10<- u_opt_N_10 + fit_N_10_sigma/fit_N_10_xi * (((n/Nu_N_10)*(1/202600))^(-fit_N_10_xi)-1)
# # B
# Q_100_N_10_2<- u_opt_N_10_2 + fit_N_10_2_sigma/fit_N_10_2_xi * (((n/Nu_N_10_2)*(1/2026))^(-fit_N_10_2_xi)-1)
# Q_1000_N_10_2<- u_opt_N_10_2 + fit_N_10_2_sigma/fit_N_10_2_xi * (((n/Nu_N_10_2)*(1/20260))^(-fit_N_10_2_xi)-1)
# Q_10000_N_10_2<- u_opt_N_10_2 + fit_N_10_2_sigma/fit_N_10_2_xi * (((n/Nu_N_10_2)*(1/202600))^(-fit_N_10_2_xi)-1)

## Parameter estimates - 20 candidates
# A
u_opt_N_20<-u_opt_N
Nu_N_20<-length(x[x>u_opt_N_20])
fit_N_20<- evir::gpd(x,u_opt_N_20)
fit_N_20_xi<-fit_N_20$par.ests[[1]]
fit_N_20_sigma<-fit_N_20$par.ests[[2]]

# B
u_opt_N_20_2<-u_opt_N_alt
Nu_N_20_2<-length(x[x>u_opt_N_20_2])
fit_N_20_2<- evir::gpd(x,u_opt_N_20_2)
fit_N_20_2_xi<-fit_N_20_2$par.ests[[1]]
fit_N_20_2_sigma<-fit_N_20_2$par.ests[[2]]
 
# 100, 1000, 1000 year return levels (20.26 observations per year)
# A
Q_100_N_20<- u_opt_N_20 + fit_N_20_sigma/fit_N_20_xi * (((n/Nu_N_20)*(1/2026))^(-fit_N_20_xi)-1)
Q_1000_N_20<- u_opt_N_20 + fit_N_20_sigma/fit_N_20_xi * (((n/Nu_N_20)*(1/20260))^(-fit_N_20_xi)-1)
Q_10000_N_20<- u_opt_N_20 + fit_N_20_sigma/fit_N_20_xi * (((n/Nu_N_20)*(1/202600))^(-fit_N_20_xi)-1)
# B
Q_100_N_20_2<- u_opt_N_20_2 + fit_N_20_2_sigma/fit_N_20_2_xi * (((n/Nu_N_20_2)*(1/2026))^(-fit_N_20_2_xi)-1)
Q_1000_N_20_2<- u_opt_N_20_2 + fit_N_20_2_sigma/fit_N_20_2_xi * (((n/Nu_N_20_2)*(1/20260))^(-fit_N_20_2_xi)-1)
Q_10000_N_20_2<- u_opt_N_20_2 + fit_N_20_2_sigma/fit_N_20_2_xi * (((n/Nu_N_20_2)*(1/202600))^(-fit_N_20_2_xi)-1)

#####



#------------------------------------------------------------- SGFSM
# Bader et al. (2018)
# Table 3.5
#####

## Candidate thresholds
#u<-sort(x)[c(1:(n-10))]
#u<-c(u[which(u[1:617]-u[2:618]<0)],u[618])  #580 thresholds

u<- quantile(x,probs = seq(0.25,0.925,by=0.075), names=F)  #10 quantis
#u<- quantile(x,probs = seq(0.25,0.953,by=0.037), names=F)  #20 quantis


bader_p_val <- tryCatch(eva::gpdSeqTests(x,u,method='ad'), error= function (e) NA, warning = function(w) NA)
FS <- tryCatch(pSeqStop(bader_p_val$p.values), error= function (e) NA, warning = function(w) NA)

#automatic choice
if(length(FS)>1){
  if(sum(FS$ForwardStop<=0.05)!=0){
    i0_fs<-which.max(u[FS$ForwardStop<=0.05])+1
    u_opt_B<- u[min(i0_fs,length(u))]
    
    rm(i0_fs)
  }
  else{u_opt_B<-u[1]}
  Nu_opt_B<-length(x[x>u_opt_B])}
if(length(FS)==1){
  u_opt_B<- NA
}

## Parameter estimates - 10 candidates
u_opt_B_10<-u_opt_B
Nu_B_10<-length(x[x>u_opt_B_10])
fit_B_10<- evir::gpd(x,u_opt_B_10)
fit_B_10_xi<-fit_B_10$par.ests[[1]]
fit_B_10_sigma<-fit_B_10$par.ests[[2]]

# 100, 1000, 1000 year return levels (20.26 observations per year)
Q_100_B_10<- u_opt_B_10 + fit_B_10_sigma/fit_B_10_xi * (((n/Nu_B_10)*(1/2026))^(-fit_B_10_xi)-1)
Q_1000_B_10<- u_opt_B_10 + fit_B_10_sigma/fit_B_10_xi * (((n/Nu_B_10)*(1/20260))^(-fit_B_10_xi)-1)
Q_10000_B_10<- u_opt_B_10 + fit_B_10_sigma/fit_B_10_xi * (((n/Nu_B_10)*(1/202600))^(-fit_B_10_xi)-1)

## Parameter estimates - 20 candidates
u_opt_B_20<-u_opt_B
Nu_B_20<-length(x[x>u_opt_B_20])
fit_B_20<- evir::gpd(x,u_opt_B_20)
fit_B_20_xi<-fit_B_20$par.ests[[1]]
fit_B_20_sigma<-fit_B_20$par.ests[[2]]

# 100, 1000, 1000 year return levels (20.26 observations per year)
Q_100_B_20<- u_opt_B_20 + fit_B_20_sigma/fit_B_20_xi * (((n/Nu_B_20)*(1/2026))^(-fit_B_20_xi)-1)
Q_1000_B_20<- u_opt_B_20 + fit_B_20_sigma/fit_B_20_xi * (((n/Nu_B_20)*(1/20260))^(-fit_B_20_xi)-1)
Q_10000_B_20<- u_opt_B_20 + fit_B_20_sigma/fit_B_20_xi * (((n/Nu_B_20)*(1/202600))^(-fit_B_20_xi)-1)

#####




#-------------------------------------------- L-Moment Ratio Diagram
#  Table 3.5
#####

#sort sample
x<-sort(x)

#GPd theoretical curve & plot
tau.4<-function(tau.3){tau.3*(1+5*tau.3)/(5+tau.3)}
curve(tau.4(x),0,0.6,xlab = expression(tau[3]), #might need to adjust the window
      ylab = expression(tau[4]))
abline(h=0,v=0,col="grey")


#space for sample l-mom ratios
t3 <- t4 <-c()
#space for corresponding l-mom ratios on GPd curve
t4_line<- t3_line <- numeric()

## Candidate thresholds
#u<-x[c(1:(n-10))] #618
u<- quantile(x,probs = seq(0.25,0.925,by=0.075), names=F)  #10 quantis
#u<- quantile(x,probs = seq(0.25,0.953,by=0.037), names=F)  #20 quantis


I<-1:length(u)

for (i in I) {
  x_exc <- x[x > u[i]]   #exceedances of the current level
  Xu=x_exc-u[i]          #corresponding excesses
  
  #sample l-mom ratios for the current threshold
  t4[i]=l4(Xu)/l2(Xu)
  t3[i]=l3(Xu)/l2(Xu)
  
  #optimize min distance in 2dim to the theoretical line
  ff <- function(x){x*(1+5*x)/(5+x)}    # tau3 vs tau4
  obj <- function(x, xi,yi) {
    sum((xi-x)^2+(yi-ff(x))^2)
  }       #euclidian distance of (xi,yi) to the GP curve
  
  xmin <- optimize(f=obj, c(-1, 1), tol = 0.0000001, xi= t3[i],yi= t4[i])
  
  
  t3_line[i] <-xmin$minimum
  t4_line[i] <- ff(xmin$minimum)
  points(t3_line[i],t4_line[i],col="blue",pch=16)
  
  rm(x_exc)
}

#add sequence of sample l-mom ratios to the plot
lines(t3[I],t4[I], type="p",col="navy") #(t3,t4)_u_i

#find closest pair to the curve
error_line <- (t3-t3_line)^2+(t4-t4_line)^2
index_line=1:length(error_line)
i0_line <- index_line[error_line[index_line]==min(error_line)]

#mark optimal pair on plot
points(t3_line[i0_line],t4_line[i0_line],col="red",pch=20, cex=2)
(u_opt_L<-u[i0_line]) # optimal threshold


# # samplepath of shape estimates
# out <- gpd(x,u_opt_L) # Fits GPD to excess over u[i0_boot]
# n_exc <- length(x[x>u_opt_L])
# shape(x,reverse = FALSE, end = length(x))
# abline(v=n_exc )

title("L-moment Ratio Diagram - 10 Candidate Thresholds")
text(0.8,1,"Selected Threshold = 4.8088")
text(0.8,0.96,"(k~488)")

# title("L-moment Ratio Diagram - 20 Candidate Thresholds")
# text(0.8,1,"Selected Threshold = 5.11294")
# text(0.8,0.96,"(k~506)")

# title("L-moment Ratio Diagram - 618 Candidate Thresholds")
# text(0.8,1,"Selected Threshold = 1.87")
# text(0.8,0.96,"(k=70/71)")



## Parameter estimates - 10 candidates
u_opt_L_10<-u_opt_L
Nu_L_10<-length(x[x>u_opt_L_10])
fit_L_10<- evir::gpd(x,u_opt_L_10)
fit_L_10_xi<-fit_L_10$par.ests[[1]]
fit_L_10_sigma<-fit_L_10$par.ests[[2]]

# 100, 1000, 1000 year return levels (20.26 observations per year)
Q_100_L_10<- u_opt_L_10 + fit_L_10_sigma/fit_L_10_xi * (((n/Nu_L_10)*(1/2026))^(-fit_L_10_xi)-1)
Q_1000_L_10<- u_opt_L_10 + fit_L_10_sigma/fit_L_10_xi * (((n/Nu_L_10)*(1/20260))^(-fit_L_10_xi)-1)
Q_10000_L_10<- u_opt_L_10 + fit_L_10_sigma/fit_L_10_xi * (((n/Nu_L_10)*(1/202600))^(-fit_L_10_xi)-1)


# ## Parameter estimates - 20 candidates
# u_opt_L_20<-u_opt_L
# Nu_L_20<-length(x[x>u_opt_L_20])
# fit_L_20<- evir::gpd(x,u_opt_L_20)
# fit_L_20_xi<-fit_L_20$par.ests[[1]]
# fit_L_20_sigma<-fit_L_20$par.ests[[2]]
# 
# # 100, 1000, 1000 year return levels (20.26 observations per year)
# Q_100_L_20<- u_opt_L_20 + fit_L_20_sigma/fit_L_20_xi * (((n/Nu_L_20)*(1/2026))^(-fit_L_20_xi)-1)
# Q_1000_L_20<- u_opt_L_20 + fit_L_20_sigma/fit_L_20_xi * (((n/Nu_L_20)*(1/20260))^(-fit_L_20_xi)-1)
# Q_10000_L_20<- u_opt_L_20 + fit_L_20_sigma/fit_L_20_xi * (((n/Nu_L_20)*(1/202600))^(-fit_L_20_xi)-1)


# ## Parameter estimates - 618 candidates
# u_opt_L_618<-u_opt_L[1]
# Nu_L_618<-length(x[x>u_opt_L_618])
# fit_L_618<- evir::gpd(x,u_opt_L_618)
# fit_L_618_xi<-fit_L_618$par.ests[[1]]
# fit_L_618_sigma<-fit_L_618$par.ests[[2]]
# 
# #100, 1000, 1000 year return levels (20.26 observations per year)
# Q_100_L_618<- u_opt_L_618 + fit_L_618_sigma/fit_L_618_xi * (((n/Nu_L_618)*(1/2026))^(-fit_L_618_xi)-1)
# Q_1000_L_618<- u_opt_L_618 + fit_L_618_sigma/fit_L_618_xi * (((n/Nu_L_618)*(1/20260))^(-fit_L_618_xi)-1)
# Q_10000_L_618<- u_opt_L_618 + fit_L_618_sigma/fit_L_618_xi * (((n/Nu_L_618)*(1/202600))^(-fit_L_618_xi)-1)


#####



#-------------------------------------------- Additional validation
#  Table 3.5 + Section 3.4.3
#####

u<- quantile(x,probs = seq(0.25,0.925,by=0.075), names=F)  #10 quantis
#u<- quantile(x,probs = seq(0.25,0.953,by=0.037), names=F)  #20 quantis

bader_p_val <- tryCatch(eva::gpdSeqTests(x,u,method='ad'), error= function (e) NA, warning = function(w) NA)
FS <- tryCatch(pSeqStop(bader_p_val$p.values), error= function (e) NA, warning = function(w) NA)

# GOF ADtest all candidates
round(bader_p_val,4)

# GOF ADtest MRL WLS choice
eva::gpdAd(x[x>1.952]-1.952)

# GOF ADtest Max cor QQplot choice
eva::gpdAd(x[x>4.300]-4.300)

# GOF ADtest ALRSM choice I=618
eva::gpdAd(x[x>1.870]-1.870)


# CCDF - Figure 3.14

#Weibull plotting-positions
pp<-(1:n)/(n+1)
exp_quantiles<- log(1-pp)

plot(sort(ns),exp_quantiles, ylab = bquote(log[ ] (1-p[i])), xlab=expression(x[i:n]), 
     mgp=c(1.8, .5, 0), col='red',pch=18)  #same as reference paper
box(lty='1111',col='black')

axis(side=3, at=quantile(ns,probs = c(.25,.50,.65,.80,.90,.95,.975,.99)), 
     labels=c(25,50,65,80,90,95,97.5,99)
     ,tick=T, mgp=c(3, .5, 0),cex.axis=0.7)

abline(v=quantile(ns,probs = c(.25,.50,.65,.80,.90,.95,.975,.99)), col='lightgray',lty=3)
grid(NA,NULL)  



###################################################################
######################### North Sea Data ########################## 
############# 628 storm peak significant wave heights ############# 
########################### ONE LOCATION ##########################
#################### October 1964 - March 1995 ####################
################## Only periods October - March ###################
#################### From Northrop et al 2017 ##################### 
###################################################################
