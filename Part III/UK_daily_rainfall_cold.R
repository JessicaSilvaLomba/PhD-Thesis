##############################################################
#           UK RAINFALL DATA ANALYSIS - COLD SEASON          #
#             cf. Chapter 7 Silva Lomba (2023)               #
# Based and adapted from public code by C. Zhou and C. Neves #
#        https://github.com/zhouchen0527/rainscedasis        #
#                  cf. Einmahl et al. (2022)                 #
##############################################################



################# PREAMBLE #################

# Packages

library(beepr)
library(ggplot2)
library(viridis)
library(ggmap)
library(maps)
library(mapdata)
library(devtools)
library(stringr)
library(tidyverse)
library(RColorBrewer)
library(dplyr)
library(patchwork)
library(purrr)
library(grid)
library(gridExtra)
library(eva)
library(stats)
library(reshape)
library(tictoc)
library(exdex)
library(Rcpp)
library(SpatialExtremes)
library(sp)
# Basic maps
library("sf")            
library("rnaturalearth")  
library("rnaturalearthdata")

################# FUNCTIONS #################


# Mix Moment EVI estimator
# computes estimates for a series for thresholds 
# between X[N-kmax : N] and X[N-kmin : N]
Mix_Mom_Est<-function(pooled_data, kmin,kmax){
  
  X<-sort(pooled_data)  #ascending order
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



# asymptotic variance of MM under independence, given an array of estimates of xi
var_mom_asym_0<-function(EVI){ #functions for vectors or matrices
  l<-dim(as.array(EVI))
  sigma_sqr<-array(NA,l)
  sigma_sqr[EVI<0]<-(1-2*EVI[EVI<0])^4*(1-EVI[EVI<0])^2*
    (6*EVI[EVI<0]^2-EVI[EVI<0]+1)/((1-2*EVI[EVI<0])^3*(1-3*EVI[EVI<0])*(1-4*EVI[EVI<0]))
  sigma_sqr[EVI>0]<-(1+EVI[EVI>0])^2
  return(sigma_sqr)
}



# Obtain sum of covariance matrix M(j/k,l/k) cf. pg 153
# Double sum of covariances omega_j1(s)omega_j2(t)
# ==== REQUIRED Cj1 =======
# Total Integrated Scedasis per station Cj1=1/k*(Sum excesses)

M_st <- function(j, l, k, evi){  #k must be the same as used for Cj1
  
  # if(season=="Summer"){Cj1<-Cj1_Summer}  #this must make the algorithm slower.
  # else{if(season=="Winter"){Cj1<-Cj1_Winter}
  #   else{stop('season must be "Summer" or "Winter"')}}
  
  c1 = matrix(1,n,m)
  c2 = matrix(1,n,m)
  c3 = matrix(1,n,m)
  c1[data<X[N-j+1]]<- 0 # Select j exceedances above X_{N-j:N} 
  c2[data<X[N-l+1]]<- 0 # Select l exceedances above X_{N-l:N} 
  c3[data<X[N-k+1]]<- 0 # Select k exceedances above the baseline overall threshold X_{N-k:N}
  ajl=matrix(0,nrow = m, ncol = m)  #this is an mXm matrix - entries of M(s,t)
  
  if(evi>0){ # xi pos
    for (j1 in 1:m){
      for (j2 in 1:m){
        #if(j2!=j1){
        ajl[j1,j2]= (j*l/k^2)^(-1)*sum(as.vector(c1[,j1]*c2[,j2]))/k +
          sum(as.vector(c3[,j1]*c3[,j2]))/k - 
          (j/k)^(-1) * sum(as.vector(c1[,j1]*c3[,j2]))/k - 
          (l/k)^(-1)* sum(as.vector(c3[,j1]*c2[,j2]))/k
      }#else{# Main diagonal
      ajl[j1,j1]= ((j*l/k^2)^(-1 ) * min(j,l)/k - 1)* Cj1[j1]
      #}
      #}
    }
  }
  else{ #xi neg
    for (j1 in 1:m){
      for (j2 in 1:m){
        #if(j2!=j1){
        ajl[j1,j2]= (j*l/k^2)^(-evi-1)*sum(as.vector(c1[,j1]*c2[,j2]))/k +
          sum(as.vector(c3[,j1]*c3[,j2]))/k - 
          (j/k)^(-evi-1) * sum(as.vector(c1[,j1]*c3[,j2]))/k - 
          (l/k)^(-evi-1)* sum(as.vector(c3[,j1]*c2[,j2]))/k
      }#else{# Main diagonal
      ajl[j1,j1]= ((j*l/k^2)^(-1 ) * min(j,l)/k + 1 - 
                     (j/k)^(-evi)-(l/k)^(-evi))* Cj1[j1]
      #}
      #}
    }
  }
  
  #rm(c1,c2,c3)
  aux= sum(ajl) # dependence across space is assumed present
  #rm(ajl)
  
  return(aux) 
}


# Asymptotic std deviation of MM for an EVI estimate
# uses covariance matrix M(j/k,l/k) -- M_st above
f_aux<-function(s,evi){return(1-(1+2*evi)*s^evi)}    #positive xi
g_aux<-function(s,evi){return((1-2*evi)*s^(-evi)-1)} #negative xi

library(compiler)
g_aux<- cmpfun(g_aux)
f_aux<- cmpfun(f_aux)

EVI_Mix_Mom_StDv<- function(k0, evi0){  
  #already sigma/sqrt(k)
  aux1= 0
  if(evi0>0){
    for (j in 1:(k0-1)){
      print(paste("current j is", j))
      for(l in 1:(k0-1)){
        inner= M_st(j, l, k0, evi0)
        aux1= aux1 +
          f_aux((j/k0),evi0) * inner * f_aux((l/k0),evi0) /k0^2
      }   
    }
    evistdv=  (evi0+1)^2/evi0* sqrt(aux1/k0)
  }
  else{ 
    for (j in 1:(k0-1)){
      print(paste("current j is", j))
      for(l in 1:(k0-1)){
        inner= M_st(j, l, k0, evi0)
        aux1= aux1 +
          g_aux((j/k0),evi0) * inner * g_aux((l/k0),evi0) /k0^2
      }   
    }
    evistdv=  sqrt((1-2*evi0)^4)*sqrt((1-evi0)^4/(evi0^2*(1-2*evi0)^2))* sqrt(aux1/k0)
  }
  return(evistdv)
}




# TESTING THE SCEDASIS

# station-wise trend test
# testing Cj(t)=tCj(1) for 0<=t<=1
# cf. Einmahl et al. (2022)

testC_stationwise <- function(data, n, m , k){  #at fixed overall threshold
  N<-n*m
  s <- seq(1/n, 1, 1/n)
  results <- array(dim=c(m,2))
  thres<- sort(data)[N-k]
  for(j in 1:m){
    tempdata <- data[,j]
    Cest <- cumsum((tempdata>thres))/k      
    results[j,1] = sqrt(k*Cj1[j])* max(abs(Cest/Cj1[j]-s))   # KS Test statistic
    results[j,2] <- k*Cj1[j]/(n) * (t(Cest/Cj1[j]-s)%*% (Cest/Cj1[j]-s)) # AD statistic
  }
  return(results)
}


# Testing total integrated scedasis is constant over space
# overall station test for a given k, elapsed time
# testing Cj(1)=1/m for all j=1,...m
# cf. Einmahl et al. (2022)

testC_total_time <- function(data, n, m, k){
  klen <- length(k)
  results <- numeric(klen)
  for(i in 1:klen){
    thresh <- sort(data)[N-k[i]]
    potdata <- data>thresh
    matsigma <- t(potdata)%*%potdata/k[i]  # Hat_Sigma(1,1,1,1)
    B1 <- cbind(diag(1,m-1,m-1),rep(-1,m-1))
    Sigma <- B1%*%matsigma%*%t(B1)
    Cones <- apply(potdata,2,sum)/k[i]
    compare <- B1%*%Cones
    results[i] <- tryCatch(k[i]*t(compare)%*%solve(Sigma)%*%(compare), error = function(e) NA) # solve is inverse
  }
  return(cbind(k,results)) #->d to Chi^2(m-1)
}




################# LOAD DECLUSTERED DATA AND STATION SELECTION #################

# information on the 81 stations -- lat and long
load("locMat_81stations.rda")

# declustered daily records for cold months, 81 stations
# 1346x81
load("Data_cold_EqSpacDeclFALSE.rda")
ALL_WINTER<-Data


season<-"Winter"; initial<-"W"
m<-81


# REGION 1 (41 previously selected stations)
R1<-c(4,5,6,8,9,10,12,13,14,17,24,26,27,35,36,40,43,
      45,49,53,54,55,56,57,58,61,62,63,64,67,68,69,70,72,74,
      75,76,77,79,80,81)
REGION1<-paste0("st.",R1) 

# REGION 2: (14 previously selected stations)
R2<-c(1,3,11,20,21,22,23,30,31,32,41,42,50,51)
REGION2<-paste0("st.",R2)

# add station id.
locMat81$Station<-paste0("st.",c(1:81))

# add region id. per station
locMat81$Region[R1]<-1
locMat81$Region[R2]<-2

# remove unselected stations -- comment for the plot of all stations
locMat81<-subset(locMat81, is.na(Region)==F)

Winter1<-Data[,R1]
colnames(Winter1)<-REGION1

Winter2<-Data[,R2]
colnames(Winter2)<-REGION2

n<-dim(Data)[1]
rm(Data)



# MAP ALL STATIONS / MAP REGIONS

sbbox_all <- make_bbox(lon = c(-4, 2), lat = c(52.5,56.5), f = .1)
# get map
region_all <- get_map(location=sbbox_all, zoom=7, source = "osm",maptype="roadmap")
# create map
map_region_all <- ggmap(region_all)  +
  xlab("Longitude") + ylab("Latitude") 

# Complete map 81 -- before removing unselected stations
# (map_region_all_stations<- map_region_all + ggtitle("Gauging Stations")+
#     geom_point(data = locMat81, mapping = aes(x = lon, y = lat), color = "red", shape=4, stroke=3)+
#     geom_text(data = locMat81,aes(x = lon, y = lat, label = row.names(locMat81), 
#                                   angle = 0, vjust = 1.3, hjust=-0.4), size=4))

# Selcted regions 1 and 2
(map_regions2<- map_region_all + ggtitle(paste0("Gauging Stations - ",season))+
    geom_point(data = locMat81, mapping = aes(x = lon, y = lat, color=as.factor(Region)), shape=4, stroke=3)+
    geom_point(data = locMat81, mapping = aes(x = lon, y = lat, color=as.factor(Region)), shape=4,colour = "black")+
    geom_text(data = locMat81,aes(x = lon, y = lat, label = row.names(locMat81),
                                  angle = 0, vjust = 1.5, hjust=-0.3), size=4,fontface=2)+
    scale_colour_manual(name="Region",values=c("Green","Blue"))+ theme_gray(base_size=20)+ theme(plot.title = element_text(size=22)))

#ggsave("W_REGIONS.pdf", width = 6, height = 6, dpi = 100)





# Box and Violin-plots of station MM EVI for station choice

# stability region
kmin<-230
kmax<-280
kseq<-kmin:kmax


# station-wise EVI estimation through MM
MMe_station<-matrix(nrow=m,ncol=length(kseq))

for(s in 1:m){
  MMe_station[s,]<-Mix_Mom_Est(ALL_WINTER[,s],kmin,kmax)$xi_MM
}

df.MMe.station<-data.frame(kseq,t(MMe_station))
colnames(df.MMe.station)<-c("k_star",paste0("st.",seq(1:81)))
df.MMe.station_81<-melt(df.MMe.station , id.vars="k_star", variable_name = 'Station')
df.MMe.station_81$id<-as.numeric(df.MMe.station_81$Station)

# Region 1
(boxplot_81_reg1<-ggplot(data=df.MMe.station_81, 
                         aes(reorder(id, value, mean),value))+
    geom_boxplot()+geom_boxplot(data=subset(df.MMe.station_81,Station %in% REGION1),fill="Green")+
    #labs(x = "Station", y = bquote(hat(xi)[j]^MM*(k^s))) +
    labs(x = " ", y = bquote(hat(xi)[j]^MM*(k^s))) +
    #ggtitle(paste0(season,", k^s=",kmin,"...",kmax))+
    #ggtitle(paste0("Region, k^s=",kmin,"...",kmax))
    scale_y_continuous(limits=c(-0.04,0.24),breaks = seq(-0.14,0.4,0.02))+
    theme(panel.grid.major.y = element_line(color="white",size=1))+
    theme_gray(base_size=20)+
    theme(axis.text.x = element_text(angle = 60, vjust =1, hjust=1, size=10))+
    geom_hline(yintercept=0.07 , color="Green",size=1,linetype='dashed')+
    geom_hline(yintercept=0.14,color="Green",size=1,linetype='dashed')+
    theme(axis.title.x = element_blank())+
    theme(axis.title.y = element_text(vjust =-1.5, size=17))) 

# Region 2
(boxplot_81_reg2<-ggplot(data=df.MMe.station_81, 
                         aes(reorder(id, value, mean),value))+
    geom_boxplot()+geom_boxplot(data=subset(df.MMe.station_81,Station %in% REGION2),fill="Blue")+
    labs(x = "Station", y = bquote(hat(xi)[j]^MM*(k^s))) +
    #ggtitle(paste0(season,", k^s=",kmin,"...",kmax))+
    #ggtitle(paste0("Region, k^s=",kmin,"...",kmax))
    scale_y_continuous(limits=c(-0.04,0.24),breaks = seq(-0.14,0.4,0.02))+
    theme(panel.grid.major.y = element_line(color="white",size=1))+
    theme_gray(base_size=20)+
    theme(axis.text.x = element_text(angle = 60, vjust =1, hjust=1, size=10))+
    geom_hline(yintercept=0.13 , color="Blue",size=1,linetype='dashed')+
    geom_hline(yintercept=0.20,color="Blue",size=1,linetype='dashed')+
    theme(axis.title.y = element_text(vjust =-1.5, size=17))+
    theme(axis.title.x = element_text(vjust =1.5, size=17)))

grouped<- boxplot_81_reg1 / boxplot_81_reg2 & theme(plot.margin = unit(c(0.1,0.1,0.1,0.1), "cm"))
grouped + plot_annotation(
  title = bquote(.(season)*":"~ k^s == .(kmin)*",...,"*.(kmax)),
  theme = theme(plot.title = element_text(size = 20)))

#pdf(file=paste0(initial,"_BP_xi.pdf"), width=15, height=9.5)
grouped + plot_annotation(
  title = bquote(.(season)*":"~ k^s == .(kmin)*",...,"*.(kmax)),
  theme = theme(plot.title = element_text(size = 20)))
dev.off()

# Region 2
(violinplot_81_reg1<-ggplot(data=df.MMe.station_81, 
                            aes(reorder(id, value, mean),value))+
    geom_violin()+geom_violin(data=subset(df.MMe.station_81,Station %in% REGION1),fill="Green")+
    #labs(x = "Station", y = bquote(hat(xi)[j]^MM*(k^s))) +
    labs(x = " ", y = bquote(hat(xi)[j]^MM*(k^s))) +
    #ggtitle(paste0(season,", k^s=",kmin,"...",kmax))+
    #ggtitle(paste0("Region, k^s=",kmin,"...",kmax))
    scale_y_continuous(limits=c(-0.04,0.24),breaks = seq(-0.14,0.4,0.02))+
    theme(panel.grid.major.y = element_line(color="white",size=1))+
    stat_summary(fun=mean, geom="point", aes(group=Station),color="black")+
    theme_gray(base_size=20)+
    theme(axis.text.x = element_text(angle = 60, vjust =1, hjust=1, size=10))+
    geom_hline(yintercept=0.07 , color="Green",size=1,linetype='dashed')+
    geom_hline(yintercept=0.14,color="Green",size=1,linetype='dashed')+
    theme(axis.title.x = element_blank())+
    theme(axis.title.y = element_text(vjust =-1.5, size=17)) ) 

# Region 2
(violinplot_81_reg2<-ggplot(data=df.MMe.station_81, 
                            aes(reorder(id, value, mean),value))+
    geom_violin()+geom_violin(data=subset(df.MMe.station_81,Station %in% REGION2),fill="Blue")+
    labs(x = "Station", y = bquote(hat(xi)[j]^MM*(k^s))) +
    #ggtitle(paste0(season,", k^s=",kmin,"...",kmax))+
    #ggtitle(paste0("Region, k^s=",kmin,"...",kmax))
    scale_y_continuous(limits=c(-0.04,0.24),breaks = seq(-0.14,0.4,0.02))+
    theme(panel.grid.major.y = element_line(color="white",size=1))+
    stat_summary(fun=mean, geom="point", aes(group=Station),color="black")+
    theme_gray(base_size=20)+
    theme(axis.text.x = element_text(angle = 60, vjust =1, hjust=1, size=10))+
    geom_hline(yintercept=0.13 , color="Blue",size=1,linetype='dashed')+
    geom_hline(yintercept=0.20,color="Blue",size=1,linetype='dashed')+
    theme(axis.title.y = element_text(vjust =-1.5, size=17))+
    theme(axis.title.x = element_text(vjust =1.5, size=17)))

grouped<- violinplot_81_reg1 / violinplot_81_reg2 & theme(plot.margin = unit(c(0.1,0.1,0.1,0.1), "cm"))
grouped + plot_annotation(
  title = bquote(.(season)*":"~ k^s == .(kmin)*",...,"*.(kmax)),
  theme = theme(plot.title = element_text(size = 20)))

#pdf(file=paste0(initial,"_VP_xi.pdf"), width=15, height=9.5)
grouped + plot_annotation(
  title = bquote(.(season)*":"~ k^s == .(kmin)*",...,"*.(kmax)),
  theme = theme(plot.title = element_text(size = 20)))
dev.off()




# Line plots of doubtful stations in each region

kmin<-20
kmax<-round(n*0.75)

kseq<-kmin:kmax

# MM and ML station-wise estimates of EVI
MMe_station<-matrix(nrow=m,ncol=length(kseq))
#MLe_station<-matrix(nrow=m,ncol=length(kseq))


for(s in 1:m){
  MMe_station[s,]<-Mix_Mom_Est(ALL_WINTER[,s],kmin,kmax)$xi_MM
  
  # for(i in (1:(kmax-kmin+1))){
  #   MLe_station[s,i]<-gpdFit(as.vector(ALL_WINTER[,s]), nextremes=kseq[i], method = c("mle"), start=NULL, opt = "Nelder-Mead", maxit = 10000)$par.ests[2]
  # }
}
beep(4)

df.MMe.station<-data.frame(kseq,t(MMe_station))
colnames(df.MMe.station)<-c("k_star",paste0("st.",seq(1:81)))
df.MMe.station_81<-melt(df.MMe.station , id.vars="k_star", variable_name = 'Station')
#write.csv(df.MMe.station_81,file=paste0(initial,"_full_estimates_81stations_MM.csv"),row.names = F)

df.MLe.station<-data.frame(kseq,t(MLe_station))
colnames(df.MLe.station)<-c("k_star",paste0("st.",seq(1:81)))
df.MLe.station_81<-melt(df.MLe.station , id.vars="k_star", variable_name = 'Station')
#write.csv(df.MLe.station_81,file=paste0(initial,"_full_estimates_81stations_ML.csv"),row.names = F)

df.estimators.station<-cbind(df.MMe.station_81,df.MLe.station_81[,-c(1,2)])
colnames(df.estimators.station)<-c("k_star",'Station',"MM","ML")
#write.csv(df.estimators.station,file=paste0(initial,"_full_estimates_81stations.csv"),row.names = F)

###

#df.estimators.station<-read.csv(file=paste0(initial,"_full_estimates_81stations.csv"))
#df.MMe.station_81<-read.csv(file=paste0(initial,"_full_estimates_81stations_MM.csv"))
#df.MLe.station_81<-read.csv(file=paste0(initial,"_full_estimates_81stations_ML.csv"))

# line plots of EVI MM and ML estimates for a specific station
region<-2
(pg4<-ggplot(subset(df.estimators.station, Station %in% c("st.41")))+
    geom_line(aes(k_star,MM,color="MM"),size=1) +
    geom_line(aes(k_star,ML,color="ML"),size=1)+
    scale_y_continuous(limits=c(-0.3,0.6), breaks=seq(-0.3, 0.9, 0.15))+
    scale_color_manual(name="Estimator",values=c("MM"="blue", "ML"="red"))+
    #gtitle(paste0(season,": Region ",region))+
    ggtitle(paste0(season,": Region ",region, ", Station 41"))+
    labs(x = bquote(k^s), y = bquote(hat(xi)[j]^"MM"*(k^s)))+
    theme_gray(base_size = 20)+
    #geom_vline(xintercept = 230, size=0.5,linetype='dashed')+
    #geom_vline(xintercept = 280, size=0.5,linetype='dashed')+
    # geom_vline(xintercept = 245, color="grey", size=0.5,linetype='dashed')+
    # geom_vline(xintercept = 250, color="grey", size=0.5,linetype='dashed')+
    # geom_vline(xintercept = 260, color="grey", size=0.5,linetype='dashed')+
    #annotate("text", x=190, y=-0.3, label=paste("~k^s==~", 230),parse=TRUE, size=5)+
    #annotate("text", x=320, y=-0.3, label=paste("~k^s==~", 280),parse=TRUE, size=5)+
    theme(axis.title.y = element_text(vjust =-1.5, size=17))+
    theme(axis.title.x = element_text(vjust =1.5, size=17))+
    theme(plot.title = element_text(size = 20)))

#png(filename=paste0(initial,"_Reg_",region,"_specific_stations_EVI_FULL.png"), width=5500, height=5000,res=300)
combined <- pg1 + pg2 + pg3 + pg4 & theme(legend.position = "bottom") & theme(plot.margin = unit(c(0.1,0.1,0.1,0.1), "cm"))
combined + plot_layout(guides = "collect")
dev.off()


# sample paths of doubtful stations' MM-EVI over all lines for region 1
region<-1
(pg1<-ggplot(subset(df.MMe.station_81, Station %in% REGION1), aes(k_star,value))+
    geom_line(size=.1) +
    scale_y_continuous(limits=c(-0.3,0.6), breaks=seq(-0.3, 0.9, 0.15))+
    ggtitle(paste0(season,": Region ",region))+
    labs(x = bquote(k^s), y = bquote(hat(xi)[j]^"MM"*(k^s)))+
    theme_gray(base_size = 20)+
    geom_vline(xintercept = 230, size=0.5,linetype='dashed')+
    geom_vline(xintercept = 280, size=0.5,linetype='dashed')+
    # geom_vline(xintercept = 245, color="grey", size=0.5,linetype='dashed')+
    # geom_vline(xintercept = 250, color="grey", size=0.5,linetype='dashed')+
    # geom_vline(xintercept = 260, color="grey", size=0.5,linetype='dashed')+
    annotate("text", x=190, y=-0.3, label=paste("~k^s==~", 230),parse=TRUE, size=5)+
    annotate("text", x=320, y=-0.3, label=paste("~k^s==~", 280),parse=TRUE, size=5)+
    theme(axis.title.y = element_text(vjust =-1.5, size=17))+
    theme(axis.title.x = element_text(vjust =1.5, size=17))+
    theme(plot.title = element_text(size = 20)))

pg2<-pg1+geom_line(data=subset(df.MMe.station_81, Station %in% paste0("st.",c(29,66))), 
                   size=0.5, aes(color=Station)) + 
  scale_color_manual(values=c("st.29"="yellow", "st.66"="cyan"))
pg2

#pdf(file=paste0(initial,"_line_plot_reg",region,"_xi.pdf"), width=12, height=5)
pg2 & theme(plot.margin = unit(c(0.1,0.1,0.1,0.1), "cm"))
dev.off()

# sample paths of doubtful stations' MM-EVI over all lines for region 2
region<-2
(pg1<-ggplot(subset(df.MMe.station_81, Station %in% REGION2), aes(k_star,value))+
    geom_line(size=.1) +
    scale_y_continuous(limits=c(-0.3,0.6), breaks=seq(-0.3, 0.9, 0.15))+
    ggtitle(paste0(season,": Region ",region))+
    labs(x = bquote(k^s), y = bquote(hat(xi)[j]^"MM"*(k^s)))+
    theme_gray(base_size = 20)+
    geom_vline(xintercept = 230, size=0.5,linetype='dashed')+
    geom_vline(xintercept = 280, size=0.5,linetype='dashed')+
    # geom_vline(xintercept = 245, color="grey", size=0.5,linetype='dashed')+
    # geom_vline(xintercept = 250, color="grey", size=0.5,linetype='dashed')+
    # geom_vline(xintercept = 260, color="grey", size=0.5,linetype='dashed')+
    annotate("text", x=190, y=-0.3, label=paste("~k^s==~", 230),parse=TRUE, size=5)+
    annotate("text", x=320, y=-0.3, label=paste("~k^s==~", 280),parse=TRUE, size=5)+
    theme(axis.title.y = element_text(vjust =-1.5, size=17))+
    theme(axis.title.x = element_text(vjust =1.5, size=17))+
    theme(plot.title = element_text(size = 20)))

pg2<-pg1+geom_line(data=subset(df.MMe.station_81, Station %in% paste0("st.",c(52,65,73))), 
                   size=0.5, aes(color=Station)) + 
  scale_color_manual(values=c("st.52"="yellow", "st.65"="cyan", "st.73"="deeppink"))
pg2

#pdf(file=paste0(initial,"_line_plot_reg",region,"_xi.pdf"), width=12, height=5)
pg2 & theme(plot.margin = unit(c(0.1,0.1,0.1,0.1), "cm"))
dev.off()








################# REGIONAL ANALYSIS #################

# Toggle between regions

Data<-Winter1; region<-1; R<-R1; REGION<-REGION1
#Data<-Winter2; region<-2; R<-R2; REGION<-REGION2

m<-dim(Data)[2]
n<-dim(Data)[1] #1346
N<-m*n



# Map stationwise MM-EVI estimates

kmin<-200
kmax<-300
kseq<-kmin:kmax
MMe_station<-matrix(nrow=m,ncol=length(kseq))
for(s in 1:m){
  MMe_station[s,]<-Mix_Mom_Est(Data[,s],kmin,kmax)$xi_MM
}
df.MMe.station<-data.frame(kseq,t(MMe_station))
colnames(df.MMe.station)<-c("k_star",REGION)

#locMat81$MMe230<-NA

k<-230
locMat81$MMe230[which(locMat81$Region==region)]<-t(df.MMe.station[(k-kmin+1),-1])


theme_set(theme_bw())    #dark-on-light theme
uk<-ne_countries(country='United Kingdom' ,scale = "medium", returnclass = "sf")
map_region<- ggplot(data = uk) + geom_sf() +
  coord_sf(xlim = c(-4, 2), ylim = c(52,57), expand = FALSE) + 
  xlab("Longitude") + ylab("Latitude") + 
  ggtitle(season) 

(map1_bw<- map_region + 
    geom_point(data = subset(locMat81,Station %in% REGION) , mapping = aes(x = lon, y = lat, color = MMe230), size=9)+
    geom_text(data = subset(locMat81,Station %in% REGION),aes(x = lon, y = lat, 
                                                              label = row.names(subset(locMat81,Station %in% REGION)), angle = 0, vjust = -0.5, hjust=1.5), size=4, fontface=2)+
    geom_text(data = subset(locMat81,Station %in% REGION),aes(x = lon, y = lat, 
                                                              label = round(MMe230,3), angle = 0, vjust = 2.5, hjust=0), size=3)+ 
    ggtitle(paste0(season,": Region ",region))+ 
    #scale_colour_viridis_c(name= paste0("MMe k^s=",k), option = "plasma"))
    scale_colour_distiller(name= bquote(k^s==.(k)), limits=c(-0.2,0.4),palette = "Spectral")+
    theme_bw(base_size=20))
theme_set(theme_gray())

#pdf(file=paste0(initial,"_reg_",region,"_map_k_",k,".pdf"))
#map2_bw
#dev.off()





# TESTING RESIDUAL SERIAL DEPENDENCE (TIME)



# Testing the Series of daily maxima for residual DEPENDENCE over TIME
# Daily maxima series regardless of station
maxdata<-as.vector(apply(Data,1,FUN=max))
maxdata2<-data.frame("day"=rownames(Data),"DailyMax"=maxdata)

# Plot series
reg1<-ggplot(maxdata2,aes(day,DailyMax))+geom_point()+
  geom_linerange(aes(x=day, ymax=DailyMax,ymin=0))+
  scale_x_discrete(breaks=maxdata2$day[seq(1,n,100)])+
  theme_gray(base_size = 20)+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.7, hjust=0.7, size=10))+#element_blank())+
  labs(x=' ', y="Daily Maximum Rainfall")+
  geom_hline(yintercept=0)+
  ggtitle(paste0(season," - Region ",region))

#(combine<-reg1 / reg2 & theme(plot.margin = unit(c(0.1,0.1,0.1,0.1), "cm"),
                              axis.title.x = element_blank(), axis.title.y = element_blank()))

#pdf(file=paste0(initial,"_DailyMax.pdf"), width=15, height=9.5)
#combine
#dev.off()




# Estimating theta for daily maxima series

# theta <- spm(maxdata,b=50,bias_adjust = "BB3", varN=FALSE, which_dj = "first")
# # Estimates: BB2018b is BB2018 - 1/b
# theta
# summary(theta)
# plot(confint(theta, interval_type = 'both'))
# 
# seqb<-seq(from=10, to=100, by=1)  # block size
# lb<-length(seqb)                  # #tested blocks sizes
# 
# theta_block<- choose_b(maxdata, seqb,interval_type = 'norm', bias_adjust="BB3",
#                        varN = FALSE, )
# plot(theta_block)

seqb<-seq(from=5, to=85, by=10)  # block size
lb<-length(seqb)

N2015<-data.frame("b"=NA,"theta"=NA, "se"=NA, "UBn"=NA, "LBn"=NA,"UBl"=NA, "LBl"=NA)
BB2018b<-data.frame("b"=NA,"theta"=NA, "se"=NA, "UBn"=NA, "LBn"=NA,"UBl"=NA, "LBl"=NA)

for(i in 1:lb){
  b <- seqb[i]
  theta<-spm(data = maxdata, b, bias_adjust = "BB3", varN = T)
  confs<-confint(theta, interval_type = 'both', maxima=c("sliding"))
  N2015[i,]<-c(b,theta$theta_sl[1],theta$se_sl[1],confs$cis[1,2],confs$cis[1,1],confs$cis[4,2],confs$cis[4,1])
  BB2018b[i,]<-c(b,theta$theta_sl[3],theta$se_sl[3],confs$cis[3,2],confs$cis[3,1],confs$cis[6,2],confs$cis[6,1])
}

round(N2015,2)
round(BB2018b,2)



u <- quantile(maxdata, probs = c(seq(0.5, 0.98, by = 0.025), 0.99), na.rm = TRUE)
imt_theta <- choose_uk(maxdata, u = u, k = 1:5)
plot(imt_theta, uprob = F, alpha = 0.05)
summary(kgaps(maxdata, 19, k = 1)) #Reg1 19 #Reg2 14
summary(exdex::dgaps(maxdata, u=19, D = 1))



# Station-wise Extremal Index of declustered series

current.station<-Data[,1] # see for all stations
u <- quantile(current.station, probs = c(seq(0.5, 0.98, by = 0.025), 0.99))
imt_theta <- choose_uk(current.station, u = u, k = 1:5)
plot(imt_theta, uprob = F, alpha = 0.05/m)
summary(kgaps(current.station, 11, k = 1)) #~90%


theta_station_kgaps<-matrix(nrow=ncol(Data), ncol=2)
theta_station_Dgaps<-matrix(nrow=ncol(Data), ncol=2)
colnames(theta_station_kgaps)<-c("theta", "se")
colnames(theta_station_Dgaps)<-c("theta", "se")
for (i in 1:ncol(Data)) {
  current.station<-Data[,i]
  aux<-kgaps(current.station, 11, k = 1) #Reg1 11 Reg2 9
  aux2<-dgaps(current.station, 11, D = 1)
  theta_station_kgaps[i,] <- c(aux$theta,aux$se)
  theta_station_Dgaps[i,] <- c(aux2$theta,aux2$se)
}
theta_station_kgaps; theta_station_Dgaps



# STUDYING SPATIAL DEPENDENCE - non-parametric approach




#k<-1150  # number of order statistics for pointwise estimation R1

kmin<-100; kmax<-4000  # range of thresholds considered
N<-n*m

# Pairwise dependence for all pairs of stations

X<-sort(as.vector(Data))
pairs_OG<-combn(R,2) # station combinations (original numbers)
pairs<-t(combn(1:m,2))    # stations combination sequential numbers, as in col(Summer)
colnames(pairs)<-c('j1','j2')

# msig_j1j2<- matrix(0,nrow=dim(pairs)[1], ncol=(kmax-kmin+1)) # pairwise dep. for all pairs
# # varying k between kmin and kmax
# colnames(msig_j1j2)<-seq(kmin,kmax,1)
# 
# tic()
# for(i in 1:dim(pairs)[1]){
#   j1<-pairs[i,1] ; j2<-pairs[i,2]  # current pair
#   for(k0 in kmin:kmax){
#     ck<- matrix(1,n,m)
#     ck[Data<X[N-k0+1]]<- 0
#     msig_j1j2[i,k0-kmin+1]= sum(as.vector(ck[,j1]*ck[,j2]))/k0 * m # REMOVE * m for sigma
#     rm(ck)
#   }
# }
# toc()
# beep(4)
# 
# write.csv(msig_j1j2,paste0(initial,"_Reg_",region,"_msigma_matrix.csv"),row.names = F)
# write.csv(msig_j1j2,paste0(initial,"_Reg_",region,"_Sigma_matrix.csv"),row.names = F)

#Sig_j1j2<-read.csv(paste0(initial,"_Reg_",region,"_Sigma_matrix.csv"))
msig_j1j2<-read.csv(paste0(initial,"_Reg_",region,"_msigma_matrix.csv"))


pair_choice<-1:dim(pairs)[1] # work with all pairs, change here to plot specific pairs

#pair_choice<-sort(rdunif(100,b=dim(pairs)[1]))  
pair_choice_original<-rbind(REGION[pairs[pair_choice,1]],
                            REGION[pairs[pair_choice,2]]) # retrieves original id of stations


aux_plot_data<-data.frame((kmin:kmax), t(msig_j1j2)[,pair_choice])
colnames(aux_plot_data)<-c("k",paste("p",as.character(pair_choice)))

aux_plot_data2<-melt(aux_plot_data, id="k",variable='pairs')

(pair_dep_plot<-ggplot(aux_plot_data2, aes(k, value)) +
    ggtitle(season)+labs(x=bquote(k[N]),y=bquote(m * hat(sigma)[.("j1,j2")](k[N])))+
    geom_line(aes(group=pairs, colour = pairs), size=1)+
    #scale_y_continuous(limits=c(0,1.6), breaks=seq(0, 1.6, 0.2)) + 
    scale_x_continuous(limits=c(kmin,kmax), breaks=seq(0,kmax,200)) + 
    scale_colour_discrete(name="Stations\nj1-j2", type=plasma(4),
                          labels=c(paste(pair_choice_original[1,],"-",pair_choice_original[2,]))))& 
  theme(legend.position = "none")


#png(filename=paste0(initial,"_Reg_",region,"_mSigma_FULL.png"), width=2300, height=1800,res=300)
#pair_dep_plot & theme(legend.position = "none")
#dev.off()


# Calculate euclidean distance between stations and plot dependence vs. distance

subset_locMat<-as.matrix(locMat81[locMat81$Region==region,c(1,2)])
# vector of distances in degrees
# distance(subset_locMat)
# matrix of distances in km - with 0 diagonal
dist_mat_km<-spDists(subset_locMat, subset_locMat, longlat = T, segments = FALSE, diagonal = FALSE)
colnames(dist_mat_km)<-REGION
rownames(dist_mat_km)<-REGION


dist_pairs<-data.frame('j1'=pair_choice_original[1,], 'j2'=pair_choice_original[2,], "dist"=NA)
for(i in 1:(dim(pairs)[1])){
  dist_pairs$dist[i]<-dist_mat_km[dist_pairs$j1[i],dist_pairs$j2[i]]
}

#ksigma<-1150   #R1 - threshold = 18.5
ksigma<-500    #R1 - threshold = 22.3
#ksigma<-1400   #R2 - threshold = 11.6
#ksigma<-250    #R2 - threshold = 20.7


col_ksigma<-paste0("X",ksigma)
dist_pairs$msig<-msig_j1j2[,col_ksigma]

(psig<-ggplot(dist_pairs, aes(dist,msig))+geom_point(size=2,shape=4,col="blue",stroke=2)+
    ggtitle(paste0(season,": Region ",region))+
    scale_y_continuous(limits=c(0,1), breaks=seq(0, 1, 0.2))+
    scale_x_continuous(limits=c(min(dist_pairs$dist),max(dist_pairs$dist)), breaks=seq(0,kmax,50))+
    labs(x="d(j1,j2) (km)",y=bquote(m * hat(sigma)[.("j1,j2")](.(ksigma))))+
    theme_gray(base_size = 20) & theme(plot.margin = unit(c(0.1,0.1,0.1,0.1), "cm")))

#pdf(file=paste0(initial,"_Reg_",region,"_msig_",ksigma,".pdf"))
psig
#dev.off()


# Line-plots of dependence for some pairs

# REG 1 stations further away (589.9 km) 9-74 (184); 
# REG 1 all side by side pairs are approx. at the same distance (55.3 approx.), 
# REG 1 choose 3: 8-9 (118); 67-68 (755); 80-81 (820)

# REG 2 stations further away (392.2 km) 1-51 (13); 
# REG 2 all side by side pairs are at the same distance (55.3 approx.), 
# REG 2 choose 2: 21-22 (47); 50-51 (91)


 pair_choice<-c(184,118,755,820)
# pair_choice<-c(13,47,91) 
 

# comparison between seasons
# WR2 vs SR1
pair_choice<-c(13,56)  #13: 1-51 56:22-23 

# WR1 vs SR2
pair_choice<-c(668,577)  #668: 57-58 577: 49-64


pair_choice_original<-rbind(REGION[pairs[pair_choice,1]],
                            REGION[pairs[pair_choice,2]])#retrieves original id of stations


aux_plot_data<-data.frame((kmin:kmax), t(msig_j1j2)[,pair_choice])
# aux_plot_data<-data.frame((kmin:kmax), t(Sig_j1j2)[,pair_choice])
colnames(aux_plot_data)<-c("k",paste("p",as.character(pair_choice)))

aux_plot_data2<-melt(aux_plot_data, id="k",variable='pairs')

(pair_dep_plot<-ggplot(aux_plot_data2, aes(k, value)) +
    ggtitle(paste0(season,": Region ",region))+
    labs(x=bquote(k[N]),y=bquote(m * hat(sigma)[.("j1,j2")](k[N])))+
    #labs(x=bquote(k[N]),y=bquote(hat(sigma)[.("j1,j2")](k[N])))+
    geom_line(aes(group=pairs, colour = pairs), size=1)+
    scale_y_continuous(limits=c(0,1), breaks=seq(0, 1.6, 0.1)) +
    #scale_y_continuous(limits=c(0,0.06), breaks=seq(0, 0.06, 0.01)) +
    scale_x_continuous(limits=c(kmin,1400), breaks=seq(0,kmax,200)) + 
    scale_colour_discrete(name="Stations\n j1-j2", type=plasma(4),
                          labels=c(paste0(pair_choice_original[1,],"\n",pair_choice_original[2,])))+
    theme_gray(base_size = 20)+
    plot_layout(guides = "collect",heights = c(4,0.2)) &
    theme(plot.margin = unit(c(0.1,0.1,-0.4,0.1), "cm"), 
          legend.position = "bottom",legend.box = "vertical",legend.margin=margin()))

#pdf(file=paste0(initial,"_Reg_",region,"_Sig_lines_COMP.pdf"))
pair_dep_plot
dev.off()




# TESTING SCEDASIS (THROUGH TIME)

Pooled_Data<-as.vector(Data)
X<-sort(Pooled_Data)
k<-950 # threshold level R1
# k<-810 # threshold level R2
threshold <- X[N-k]
print(paste("overall thershold for", season, "is",threshold,"mm"))


# Total Integrated Scedasis per station Cj(1)=1/k*(Sum excesses)
Cj1 <- matrix(0,m,0)
# count total excesses of threshold per station and average by k
Cj1<- cbind(Cj1,apply(Data>threshold,2,sum)/k) 
maxC=max(Cj1)
print(paste("max frequency for", season, "is",maxC*k))


# Test constant scedasis at each station; simulates critical values from Brownian bridge
# cf. Zhou and Neves (github)
set.seed(123467)

grid <- 5000  
simbm <- 2000
increments <- rnorm(grid*simbm)
dim(increments) <- c(grid,simbm)
simbms <- apply(increments, 2, cumsum)/sqrt(grid) # W(t), t in (1,5000)

gstep <- (1:grid)/grid   # t/T in(0,1)

simbms <- simbms - gstep %*% t(simbms[grid,]) # 2000 replicas of B(t)= W(t)-(t/T)*W(T)

t1c <- apply(abs(simbms),2,max)
t2c <- apply(simbms^2, 2, sum)/grid  # for the A-D goodness-of-fit test


testing <- testC_stationwise(Data,n,m,k)
pvalue_C_j <- testing

for(j in 1:m){
  pvalue_C_j[j,1] <- sum(t1c>testing[j,1])/simbm
  pvalue_C_j[j,2] <- sum(t2c>testing[j,2])/simbm
}

locMat81$KStest<-NA
locMat81$ADtest<-NA
locMat81$KStest[which(locMat81$Region==region)]= pvalue_C_j[,1] 
locMat81$ADtest[which(locMat81$Region==region)]= pvalue_C_j[,2]


# plotting KS p-values on Map

b<- c( 0.05, 0.15, 0.25, 0.5, 0.75) #p-value scale
(mapKS<- map_region_all + 
    geom_point(data = subset(locMat81, Region %in% region) , mapping = aes(x = lon, y = lat, color = KStest), size=9)+
    geom_text(data = subset(locMat81, Region %in% region),aes(x = lon, y = lat, label = R, angle = 0, vjust = -0.5, hjust=1.5), size=4,fontface=2)+
    geom_text(data = subset(locMat81, Region %in% region),aes(x = lon, y = lat, label = round(KStest,2), angle = 0, vjust = 2.5, hjust=0), size=3.5)+ 
    #scale_colour_viridis_c(name= "K-S\np-value", breaks= b, option = "plasma"))
    scale_color_gradient2(name= "K-S\np-value", mid = "red",  high = "green", breaks= b)+
    ggtitle(paste0(season,": Region ",region))+ 
    theme_bw(base_size=20))


#pdf(file=paste0(initial,"_reg_",region,"_map_KS_timescedasis.pdf"))
mapKS
dev.off()





# TESTING SCEDASIS (THROUGH SPACE)


# Testing Cj1=1/m across stations for different values of k & plot
kminTS<-300 #R1
# kminTS<-100 #R2
kmaxTS<-1500
lseq <- seq(from=kminTS, to = kmaxTS, by=1)
tic()
Cfun <- testC_total_time(Data, n, m, lseq) 
toc()
pvalue_C_k <- 1-pchisq(Cfun[,2],df=m-1)
dt1<- data.frame(lseq, pvalue_C_k)
colnames(dt1) <- c("K", "pvalue")
#rm(lseq)

#png(filename=paste0("Scedasis_Space_",season,".png"), width=2300, height=1500,res=300)

(p1 <- ggplot(dt1, aes(K, pvalue)) + geom_line(color="blue", size=1) + 
    labs(x = "k", y = " ")+ 
    scale_y_continuous(limits=c(0,0.4), breaks=seq(0, 1.0, 0.1))+
    geom_hline(yintercept = 0.05, color="grey", size=0.5,linetype='dashed') +  
    scale_x_continuous(limits=c(kminTS,kmaxTS), breaks=seq(100,kmaxTS, 200))+
    labs(x = bquote(k[N]), y =bquote("p-values for "~T[n](k[N])) )+
    ggtitle(paste0(season,": Region ",region))+ 
    theme_grey(base_size=20))


#pdf(file=paste0(initial,"_reg_",region,"scedasis_space.pdf"))
p1
dev.off()






# SEASONAL, REGIONAL EVI ESRIMATION THROUGH MM AND ML ESTIMATORS

N<-m*n 

kmin<- 100
kmax<- round(N*0.75)

# CN suggested  
# kmin<- 700
# kmax<- 4000

 
#sorted sample, ascending order
X<-sort(Pooled_Data)


k<-950 # threshold level R1
# k<-810 # threshold level R2
threshold <- X[N-k]
print(paste("overall thershold for", season, "is",threshold,"mm"))
print(paste("thersholds considered for", season,"in the range between",X[N-kmax],"mm and", X[N-kmin], "mm"))



 
# Computing MM and ML estimates
# tic()
# Mix_Mom<-Mix_Mom_Est(Pooled_Data,kmin,kmax)  # return k, hat(phi), hat(xi)
# beep(4)
# toc()
# 
# 
# tic()
# ML<-numeric()
# for (i in kmin:kmax){
#   ML<-c(ML,gpdFit(Rank_Pooled_Data, nextremes=i, method = c("mle"), start=NULL, opt = "Nelder-Mead", maxit = 10000)$par.ests[2])
# }
# beep(4)
# toc()
# 
# 
# dtf<-data.frame(kmin:kmax,Mix_Mom$xi_MM, ML)
# colnames(dtf) <- c("k", "MM", "ML")

# write.csv(dtf,paste0(initial,"_Reg_",region,"_full_estimates.csv"),row.names = F)


dtf_evi<-read.csv(paste0(initial,"_Reg_",region,"_full_estimates.csv"))
dtf.2<-melt(dtf_evi,id.vars = "k",variable_name = "Estimator")

# Point estimates
(overall_EVI_k<-dtf_evi$MM[dtf_evi$k==k])
(overall_EVI_k_ML<-dtf_evi$ML[dtf_evi$k==k])

# Estimators sample path - full range no vertical line
(pg1 <- ggplot(dtf.2, aes(k, value))+
    scale_y_continuous(limits=c(-0.15,0.75), breaks=seq(-0.3, 2, 0.15))+ 
    scale_x_continuous(limits=c(0,kmax), breaks=c(seq(0,10000,2500),seq(10000,kmax,5000)))+ #REG1
    #scale_x_continuous(limits=c(0,kmax), breaks=c(seq(0,5000,1000),seq(5000,kmax,2500)))+ #REG2
    geom_line(aes(colour = Estimator), size=1)+ #scale_color_discrete(type=plasma(3))+
    scale_color_manual(values=c("MM"="blue", "ML"="red"))+
    #geom_vline(xintercept = 385, color="grey") +  
    #geom_line(aes(y = EVI_MM), color="#00AFBB", size=0.5) +
    #geom_line(aes(y = EVI_ML), color="red", size=0.5) +
    ggtitle(paste0(season, ": Region ",region))+labs(x = bquote(k[N]), y = bquote(hat(xi)^E*(k[N])))+
    #ggtitle(paste0(season, ": Region ",region))+labs(x = bquote(k[N]), y = "")+
    #geom_vline(xintercept = 1150, color="black", size=0.8,linetype='dashed')+
    #geom_vline(xintercept = 1700, color="black", size=0.5,linetype='dashed'))
    geom_hline(yintercept = 0, color="black", size=1,linetype='dashed')+
    #annotate("text",x=1250,y=0.145,label=bquote("k"[N]*"=1150"), color="black",size=7)+ 
    theme_gray(base_size = 20))


# pdf(file=paste0(initial,"_Reg_",region,"_overall_EVI_FULL.pdf"))
# pg1
# dev.off()
pg1_w_full<-pg1 
pg2_w_full<-pg1


# (grouped<- (pg1_w_full + pg2_w_full) /guide_area() + plot_layout(guides = "collect",heights = c(3.3,0.3))&
#    theme(plot.margin =unit(c(0.1,0.1,0,0.1), "cm"), legend.position = "bottom", 
#          axis.text.x = element_text(angle = 60, vjust =1, hjust=1, size=15)))
# axis.title.x = element_text(vjust=2)))

#pdf(file=paste0(initial,"_overall_EVI_FULL_xi.pdf"), width=14, height=6)
#grouped
#dev.off()



#Estimators sample path - smaller range no horizontal line
(pg <- ggplot(dtf.2, aes(k, value))+
    scale_y_continuous(limits=c(-0.15,0.15), breaks=seq(-2, 2, 0.02))+ 
    #scale_y_continuous(limits=c(-0.1,0.15), breaks=seq(-2, 2, 0.02))+ 
    scale_x_continuous(limits=c(100,5000), breaks=seq(0,5000,500))+
    #scale_x_continuous(limits=c(100,1800), breaks=seq(0,5000,400))+
    #scale_x_continuous(limits=c(500,5000), breaks=seq(0,5000,500))+
    geom_line(aes(colour = Estimator), size=1)+ #scale_color_discrete(type=plasma(3))+
    scale_color_manual(values=c("MM"="blue", "ML"="red"))+
    #geom_vline(xintercept = 385, color="grey") +  
    #geom_line(aes(y = EVI_MM), color="#00AFBB", size=0.5) +
    #geom_liplot.mplot.mplot.m, color="red", size=0.5) +
    ggtitle(paste0(season, ": Region ",region))+labs(x = bquote(k[N]), y = bquote(hat(xi)^E*(k[N])))+
    #geom_vline(xintercept = k, color="black", size=0.8,linetype='dashed')+
    #geom_vline(xintercept = 830, color="black", size=0.5,linetype='dashed')+
    geom_hline(yintercept = 0, color="black", size=1,linetype='dashed')+
    #annotate("text",x=1250,y=0.145,label=bquote("k"[N]==.(k)), color="black",size=7)+ 
    theme_gray(base_size = 20)+
    theme(axis.text.x = element_text(angle = 60, vjust =1, hjust=1, size=15)))

pg1_w_short<-pg & theme(legend.position = "none")
pg2_w_short<-pg & theme(legend.position = "none")

# pdf(file=paste0(initial,"_overall_EVI_ZOOM_xi.pdf"),width=14, height=6)
# ((pg1_w_short + pg2_w_short)/guide_area())  + plot_layout(guides = "collect",heights = c(3.3,0.3)) &
#  theme(plot.margin = unit(c(0.1,0.1,0.1,0.1), "cm"), legend.position = "bottom") 
# dev.off()



# MM ESTIMATION VALIDATION LEAVE ONE STATION OUT 
N=n*(m-1)
kmin<- 100
kmax<- round(N*0.75)

# Estimation with one station removed at a time
# dtf<-data.frame(kmin:kmax)
# colnames(dtf) <- c("k")
# 
# tic()
# for( out in 1:m){
#   Pooled_Data<-as.vector(Data[,-out])    
#   
#   Mix_Mom<-Mix_Mom_Est(Pooled_Data,kmin,kmax)  #return k, hat(phi), hat(xi)
#   dtf<-cbind(dtf,Mix_Mom$xi_MM)
#   
# }
# toc()
# beep(4)
# #colnames(dtf) <- c("k",REGION1)
# colnames(dtf) <- c("k",REGION2)
#write.csv(dtf,paste0(initial,"_Reg_",region,"_full_estimates_CROSS_VAL.csv"),row.names = F)

dtf2<-read.csv(paste0(initial,"_Reg_",region,"_full_estimates_CROSS_VAL.csv"))

dtf.3<-melt(dtf2,id.vars = "k",variable_name = "Removed_Station")

(pg1 <- ggplot(dtf.3, aes(k, value))+
    #scale_y_continuous(limits=c(-0.1,0.15), breaks=seq(-2, 2, 0.02))+
    scale_x_continuous(limits=c(0,kmax), breaks=c(seq(0,10000,2500),seq(10000,kmax,5000)))+ #REG1
    #scale_x_continuous(limits=c(0,kmax), breaks=c(seq(0,5000,1000),seq(5000,kmax,2500)))+ #REG2
    geom_line(aes(colour = Removed_Station), size=.1)+
    scale_color_discrete(name="Removed\nStation",type=plasma(14))+
    #scale_color_manual(values=c("MM"="blue", "ML"="red"))+
    #geom_vline(xintercept = 385, color="grey") +
    #geom_line(aes(y = EVI_MM), color="#00AFBB", size=0.5) +
    #geom_line(aes(y = EVI_ML), color="red", size=0.5) +
    ggtitle(paste0(season, ": Region ",region))+labs(x = bquote(k[N]), y = bquote(hat(xi)^(MM)*(k[N])))+
    #geom_vline(xintercept = 1150, color="black", size=0.8,linetype='dashed')+
    #geom_vline(xintercept = 1700, color="black", size=0.5,linetype='dashed'))
    geom_hline(yintercept = 0, color="black", size=1,linetype='dashed')+
    #annotate("text",x=1250,y=0.145,label=bquote("k"[N]*"=1150"), color="black",size=7)+
    theme_gray(base_size = 20))

k_select<-c(750,950,1150)
# k_select<-c(780,810,1100)
ests_full<-dtf_evi$MM[dtf_evi$k %in% k_select]

(pg2<-ggplot(data=subset(dtf.3,k %in% k_select), aes(as.factor(k), value))+
    geom_boxplot()+
    geom_point(data = data.frame(k = factor(k_select), ests = ests_full),
               aes(x=k, y=ests), color = 'blue',stroke=4,shape=20)+
    labs(x = bquote(k[N]), y = bquote(hat(xi)[m-1]^(MM)*(k[N]))) +
    ggtitle(paste0(season, ": Region ",region))+
    #scale_y_continuous(limits=c(0.04,0.1),breaks = seq(-0.14,0.6,0.02))+
    scale_y_continuous(limits=c(-0.02,0.08),breaks = seq(-0.14,0.6,0.02))+
    theme(panel.grid.major.y = element_line(color="white",size=1))+
    theme_gray(base_size=20))

#pdf(file=paste0(initial,"_Reg_",region,"_L1O_val_xi.pdf"))
#pg2
#dev.off()

# One-by-one box-plots
(pg2<-ggplot(data=subset(dtf.3,k %in% 900:1700), aes(Removed_Station, value))+
    geom_boxplot()+
    #labs(x = "Station", y = bquote(hat(xi)[j]^MM*(k^s))) +
    #ggtitle(paste0(season,", k^s=",kmin,"...",kmax))+
    #scale_y_continuous(limits=c(-0.02,0.35),breaks = seq(-0.14,0.4,0.02))+
    theme(panel.grid.major.y = element_line(color="white",size=1))+
    theme_gray(base_size=20)+
    theme(axis.text.x = element_text(angle = 60, vjust =1, hjust=1, size=10)))





# MM-BASED CONFIDENCE BANDS 


kmin<-700
kmax<-1200
#kmin<-300
#kmax<-1100
Est<- numeric(kmin-1)
lv<- numeric(kmin-1)
UB <-numeric(kmax)
LB <-numeric(kmax)
UBInd <- numeric(kmax)
LBInd <- numeric(kmax)



ESTIMATES<-Mix_Mom_Est(Pooled_Data,kmin,kmax)
Est<-c(Est,ESTIMATES$xi_MM)
lv<-c(lv,ESTIMATES$K)


X<-sort(Pooled_Data)
k<-950 #threshold level R1
#k<-810 #threshold level R2
threshold <- X[N-k]


# PARALLEL COMPUTING
# library(foreach)
# library(doParallel)
# 
# #setup parallel backend to use many processors
# cores=detectCores()
# cl <- makeCluster(cores[1]-1) #not to overload your computer
# registerDoParallel(cl)
# 
# tic()
# se <- foreach(i=kmin:kmax, .combine=cbind) %dopar% {
#   se_temp = EVI_Mix_Mom_StDv(i, Est[i])
#   se_temp #Equivalent to se = cbind(se,se_temp)
# }
# toc()
# 
# #stop cluster
# stopCluster(cl)
# 
# 
# 
# for (i in kmin:kmax){
#   se_aux= se[i-kmin+1]
#   se_ind= sqrt(var_mom_asym_0(Est[i])/lv[i])
#   UB[i] = Est[i] + qnorm(0.975)* se_aux
#   LB[i] = Est[i] - qnorm(0.975)* se_aux
#   UBInd[i] <- Est[i] + qnorm(0.975)*se_ind  # Lower bound for the CI in the independence setting
#   LBInd[i] <- Est[i] - qnorm(0.975)*se_ind # Upper bound for the CI in the independence setting
#   
# }
# 
# 
# ML<-numeric(kmin-1)
# for (i in kmin:kmax){
#   ML<-c(ML,gpdFit(Pooled_Data, nextremes=i, method = c("mle"), start=NULL, opt = "Nelder-Mead", maxit = 10000)$par.ests[2])
# }
# 
# 
# 
# dt2<- data.frame(lv,Est, LB, UB, LBInd, UBInd, Summer_ML)
# #print(paste("overall EVI estimate:", Est[k], LB[k], UB[k], LBInd[k], UBInd[k]))
# colnames(dt2) <- c("K", "EVI", "LB", "UB", "LBInd", "UBInd", "EVI_ML")
# write.csv(dt2,paste0(season,"_EVI_BANDS.csv"), row.names=F)
# rm(lv)


bands_full<-read.csv(paste0(initial,"_Reg_",region,"_EVI_BANDS.csv"),header=T,sep=",")
kmin<-700
kmax<-1200
#kmin<-300 
#kmax<-1100


(p21 <- ggplot(bands_full, aes(K, EVI)) + 
    labs(x = bquote(k[N]), y = bquote(hat(xi)^E*(k[N])))+
    #scale_y_continuous(limits=c(-0.35,0.3), breaks=seq(-1.5, 1.5, 0.1))+ 
    scale_y_continuous(limits=c(-0.17,0.17), breaks=seq(-2, 2, 0.02))+ 
    scale_x_continuous(limits=c(kmin,kmax), breaks=seq(0,kmax,100))+
    geom_ribbon(aes(ymin = LB, ymax = UB), fill = "blue4", alpha= 0.8) +
    geom_ribbon(aes(ymin = LBInd, ymax = UBInd), fill = "darkslategray1", alpha= 0.6) +
    geom_line(aes(y = EVI, color="MM"), size=1) +
    geom_line(aes(y = EVI_ML, color="ML"), size=0.5) +
    #annotate(geom="text", x=1270, y=-0.095, label="MLe",  color="orange",size=8) +
    #annotate(geom="text", x=1270, y=0.035, label="MMe",  color="red",size=8)+
    #geom_vline(xintercept = 810, color="grey", size=0.5,linetype='dashed')+
    geom_hline(yintercept = 0, color="grey",linetype='dashed')+
    scale_color_manual(values = c('MM' = 'blue','ML' = 'red')) +
    ggtitle(paste0(season, ": Region ",region))+
    labs(color = 'Estimator')+
    theme_gray(base_size = 20))


# pdf(file=paste0(initial,"_Reg_",region,"_BANDS.pdf"))
# p21 + plot_layout(guides = "collect",heights = c(3.3,0.3)) & 
#  theme(plot.margin =unit(c(0.1,0.1,0,0.1), "cm"), legend.position = "bottom")
# dev.off()



# pdf(file=paste0(initial,"_Bands_xi.pdf"),width=14, height=6)
# ((p21 + p22)/guide_area())  + plot_layout(guides = "collect",heights = c(3.3,0.3)) &
#  theme(plot.margin = unit(c(0.1,0.1,0.1,0.1), "cm"), legend.position = "bottom") 
# dev.off()





##############################################################
#           UK RAINFALL DATA ANALYSIS - COLD SEASON          #
#             cf. Chapter 7 Silva Lomba (2023)               #
# Based and adapted from public code by C. Zhou and C. Neves #
#        https://github.com/zhouchen0527/rainscedasis        #
#                  cf. Einmahl et al. (2022)                 #
##############################################################
