##########################################
#        Complete Simulation Study       #
#        ALRSM vs aSTSM vs SGFSM         #
# cf. Silva Lomba and Fraga Alves (2020) #
#       Northrop and Coleman (2014)      #
#           Bader et al (2018)           #
##########################################


################# PREAMBLE #################

# Packages
library(evir)
library(eva)

# import L-moments estimators
source('Sample_L-Moments.R')

# import hybrid random generator
source('rhybrid.R')

# import source code for STSM  
source("NorthropColeman2014.fns")


################# SELECTION FUNCTIONS #################
# x - sample
# u - candidate  thresholds

#-----------------------------------------------------
# aSTSM Northrop and Coleman (2014)

aSTSM<-function(x,u){
  northrop <- tryCatch(score.fitrange(x,u), error= function(e) NA)               
  
  ## automatized choice
  if(length(northrop)==10){
    i=0
    aux=F
    while(i<(length(u)-1) && aux==F){
      i<-i+1
      if(northrop$e.p.values[i] > 0.05 & all(northrop$e.p.values[i:(length(u)-1)] > 0.05)){aux=TRUE}
    }
    u_star<-u[i]
    Nu_star<-length(x[x>u_star])}
  else{
    u_star<-NA
    Nu_star<-NA
  }
  fit_xi_star<-tryCatch(northrop$xi.mle[i], error= function (e) NA)
  fit_sigma_star<-tryCatch(northrop$sigma.mle[i], error= function (e) NA)
  
  choice<-data.frame("threshold" = u_star, "Nexceedances" = Nu_star, 
                     "EVI" = fit_xi_star, "scale" = fit_sigma_star)
  return(choice)
} 


#-----------------------------------------------------
# SGFSM Bader et al. (2018)

SGFSM<-function(x,u){
  bader_p_val <- tryCatch(eva::gpdSeqTests(x,u,method='ad'), error= function (e) NA, warning = function(w) NA)
  FS <- tryCatch(pSeqStop(bader_p_val$p.values), error= function (e) NA, warning = function(w) NA)
  
  ## automatized choice
  if(length(FS)>1){
    if(sum(FS$ForwardStop<=0.05)!=0){
      i0_fs<-which.max(u[FS$ForwardStop<=0.05])+1
      u_star<- u[min(i0_fs,length(u))]
      
      rm(i0_fs)
    }
    else{u_star<-u[1]}
    Nu_star<-length(x[x>u_star])}
  if(length(FS)==1){
    u_star<- NA
    Nu_star<-NA
  }
  
  ## GPd fit at chosen level
  fit <- tryCatch(evir::gpd(x,u_star), error= function (e) NA, warning = function(w) NA)
  fit_xi_star<-tryCatch(fit$par.ests[[1]], error= function (e) NA)
  fit_sigma_star<-tryCatch(fit$par.ests[[2]], error= function (e) NA)
  
  choice<-data.frame("threshold" = u_star, "Nexceedances" = Nu_star, 
                     "EVI" = fit_xi_star, "scale" = fit_sigma_star)
  return(choice)
} 


#-----------------------------------------------------
# ALRSM Silva Lomba e Fraga Alves (2020)

ALRSM<-function(x,u){
  ## computing distance to GPd-specific curve of L-statistics
  i0_line<-t3<-t4<-c()
  t4_line<- t3_line <- numeric()
  I_seq<-1:length(u)
  
  for (i in I_seq) {
    x_exc <- x[x > u[i]]   #exceedances of the current level
    Xu=x_exc-u[i]          #corresponding excesses
    
    ## sample l-mom ratios for the current threshold
    t4[i]=l4(Xu)/l2(Xu)
    t3[i]=l3(Xu)/l2(Xu)
    
    ## optimize min distance in 2dim to the theoretical line
    ff <- function(x){x*(1+5*x)/(5+x)}    # tau3 vs tau4
    obj <- function(x, xi,yi) {
      sum((xi-x)^2+(yi-ff(x))^2)
    }       ## euclidian distance of (xi,yi) to the GP curve
    
    xmin <- optimize(f=obj, c(-1, 1), tol = 0.0000001, xi= t3[i],yi= t4[i])
    
    t3_line[i] <-xmin$minimum
    t4_line[i] <- ff(xmin$minimum)
    
    rm(x_exc)
  }
  
  ## automatized choice
  error_line <- (t3-t3_line)^2+(t4-t4_line)^2
  index_line=1:length(error_line)
  i0_line <- index_line[error_line[index_line]==min(error_line)]
  u_star<-u[i0_line]
  Nu_star<-length(x[x>u_star])
  
  ## GPd fit at chosen level
  fit <- tryCatch(evir::gpd(x,u_star), error= function (e) NA, warning = function(w) NA)
  fit_xi_star<-tryCatch(fit$par.ests[[1]], error= function (e) NA)
  fit_sigma_star<-tryCatch(fit$par.ests[[2]], error= function (e) NA)
  
  
  choice<-data.frame("threshold" = u_star, "Nexceedances" = Nu_star, 
                     "EVI" = fit_xi_star, "scale" = fit_sigma_star)
  return(choice)
} 




################# Choose distribution #################

distribution<- "hybrid"        #Northrop and Coleman (2014)
#distribution<- "gev"           #Generalized Extreme Value Distribution
#distribution<- "quant_hybrid"  #Rounded Hybrid to the second decimal 

################# Choose method #################

#method<-"aSTSM" #Northrop and Coleman (2014)
#method<-"SGFSM" #Bader et al. (2018)
method<-"ALRSM" #Silva Lomba and Fraga Alves (2020)


################# SIMULATION FUNCTION #################
# B - number of runs
# n - sample size
# xi - true EVI value
# true_u - true threshold (Hybrid)
# probs - sequence of probabilities for candidate thresholds at sample quantiles
# seed - recreate generated samples 

run_sim_ATS<-function(distribution, method, B, n, xi, true_u, probs, seed) {
  
  if (!is.character(my_string) || !(distribution %in% c("hybrid", "gev", "quant_hybrid"))) {
    stop("distribution must be a string equal to 'hybrid', 'gev' or 'quant_hybrid'.")
  }
  
  if (!is.character(method) || !(method %in% c("aSTSM", "SGFSM", "ALRSM"))) {
    stop("method must be a string equal to 'aSTSM', 'SGFSM' or 'ALRSM'.")
  }
  
  ## true 100/1000/10000 obs return level Hybrid
  if(distribution == "gev"){
    ## scale and location for GEVd
    sigma = (1-true_u)^(xi+1) 
    mu = true_u+(1-true_u)/xi*((1-true_u)^xi-1)

    Q_0.01_true<-evir::qgev(0.99,xi,mu,sigma)
    Q_0.001_true<-evir::qgev(0.999,xi,mu,sigma)
    Q_0.0001_true<-evir::qgev(0.9999,xi,mu,sigma)
    
  } else{
    Q_0.01_true<-qhybrid(0.99,xi2=xi,qq=true_u)
    Q_0.001_true<-qhybrid(0.999,xi2=xi,qq=true_u)
    Q_0.0001_true<-qhybrid(0.9999,xi2=xi,qq=true_u)
    
  }
  
  ## estimated RLs
  Q_0.01_est<-Q_0.001_est<-Q_0.0001_est<-c()
  ## estimated threshold and exceedance number
  u_opt<-Nu_opt<-c()
  ## estimated GPd parameters
  fit_xi<-fit_sigma<-c()
  
  
  set.seed(seed)  
  start.time <- Sys.time()
  
  for(b in 1:B){
    
    ## generate and order data 
    if(distribution == "gev"){
      x<- sort(rgev(n,xi,mu,sigma))
    } else if(distribution == "hybrid"){
      x<- sort(rhybrid(n, xi2=xi,qq=true_u))
    } else if(distribution == "quant_hybrid"){
      x<- sort(round(rhybrid(n, xi2=xi,qq=true_u),2))
    }
    
    ## candidate thresholds 
    u<- quantile(x, probs, names=F)
    
    ## automatized choice
    if(method == "aSTSM"){
      
      choice<-aSTSM(x,u)
      
    } else if(method == "SGFSM"){
      
      choice<-SGFSM(x,u)
      
    }else if(method == "ALRSM"){
      
      choice<-ALRSM(x,u)
      
    }
    
    u_opt[b]<-choice$threshold
    Nu_opt[b]<-choice$Nexceedances
    
    ## GPd fit at chosen level
    fit_xi[b]<-choice$EVI
    fit_sigma[b]<-choice$scale
    
    ## estimated quantiles at chosen level
    Q_0.01_est[b]<- u_opt[b] +  fit_sigma[b]/fit_xi[b] * ((n*0.01/Nu_opt[b])^(-fit_xi[b])-1)
    Q_0.001_est[b]<- u_opt[b] +  fit_sigma[b]/fit_xi[b] * ((n*0.001/Nu_opt[b])^(-fit_xi[b])-1)
    Q_0.0001_est[b]<- u_opt[b] +  fit_sigma[b]/fit_xi[b] * ((n*0.0001/Nu_opt[b])^(-fit_xi[b])-1)
  
    rm(choice)
  }
  
  u_med<-mean(u_opt[is.na(u_opt)==F])
  RMSE_u_opt<-sqrt(mean((u_opt[is.na(u_opt)==F]-true_u)^2))
  
  xi_med<- mean(fit_xi[is.na(fit_xi)==F])
  RMSE_xi<-sqrt(mean((fit_xi[is.na(fit_xi)==F]-xi)^2))
  
  Q_0.01_est_med<-mean(Q_0.01_est[is.na(Q_0.01_est)==F])
  Q_0.001_est_med<-mean(Q_0.001_est[is.na(Q_0.001_est)==F])
  Q_0.0001_est_med<-mean(Q_0.0001_est[is.na(Q_0.0001_est)==F])
  
  Q_0.01_r<- Q_0.01_est_med / Q_0.01_true
  Q_0.001_r<- Q_0.001_est_med / Q_0.001_true
  Q_0.0001_r<- Q_0.0001_est_med / Q_0.0001_true
  
  RMSE_0.01<-sqrt(mean((Q_0.01_est[is.na(Q_0.01_est)==F]-Q_0.01_true)^2))
  RMSE_0.001<-sqrt(mean((Q_0.001_est[is.na(Q_0.001_est)==F]-Q_0.001_true)^2))
  RMSE_0.0001<-sqrt(mean((Q_0.0001_est[is.na(Q_0.0001_est)==F]-Q_0.0001_true)^2))
  
  end.time <- Sys.time()
  
  ## Reporting the results
  if(distribution == "gev"){
    cat('Distribution/ method/ xi/ n/ I:  ', 
      distribution,'/', method,'/',
      xi,'/',n,'/',length(u), '\n',
      'Mean Selected Threshold:  ', u_med , '\n',
      'Mean EVI Bias:  ', xi_med - xi, '\n',
      'RMSE EVI:  ', RMSE_xi, '\n',
      'Mean Est. Quantile / True Quantile \n', 
      Q_0.01_r,'\n', Q_0.001_r, '\n',Q_0.0001_r, '\n',
      'RMSE Quantiles \n', 
      RMSE_0.01, '\n', RMSE_0.001, '\n',RMSE_0.0001, '\n',
      'Run time: ',format(difftime(end.time, start.time, units='secs')), '\n',
      'Failures: ', length(Q_0.001_est[is.na(Q_0.001_est)==T]), '\n')
  } else{
    cat('Distribution/ method/ xi/ n/ I:  ', 
        distribution,'/', method,'/',
        xi,'/',n,'/',length(u), '\n',
        'Mean Threshold Bias:  ', u_med - true_u, '\n',
        'RMSE Threshold:  ', RMSE_u_opt, '\n',
        'Mean EVI Bias:  ', xi_med - xi, '\n',
        'RMSE EVI:  ', RMSE_xi, '\n',
        'Mean Est. Quantile / True Quantile \n', 
        Q_0.01_r,'\n', Q_0.001_r, '\n',Q_0.0001_r, '\n',
        'RMSE Quantiles \n', 
        RMSE_0.01, '\n', RMSE_0.001, '\n',RMSE_0.0001, '\n',
        'Run time: ',format(difftime(end.time, start.time, units='secs')), '\n',
        'Failures: ', length(Q_0.001_est[is.na(Q_0.001_est)==T]), '\n')
  }

}




################# SET PARAMETERS AND RUN SIMULATION #################

# seed for multiple runs
# re-set seed for every run
seed=66

# proposed scenarios 
## Number of runs
B=1500

## sample size
#n=1000
#n=500
n=200

## EVI value
#xi=0.5
xi=0.2
#xi=-0.2

## true threshold for hybrid distribution
#true_u=0.5
true_u=0.75

## number of candidate thresholds
I=10; probs <- seq(0.25,0.925,by=0.075)
#I=20; probs <- seq(0.25,0.953,by=0.037)


### RUN 
run_sim_ATS(distribution, method, B, n, xi, true_u, probs, seed)


##########################################
#        Complete Simulation Study       #
#        ALRSM vs aSTSM vs SGFSM         #
# cf. Silva Lomba and Fraga Alves (2020) #
#       Northrop and Coleman (2014)      #
#           Bader et al (2018)           #
##########################################
