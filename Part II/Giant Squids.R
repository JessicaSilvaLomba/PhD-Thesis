###################################################################
######################## GIANT SQUID DATA #########################
############# Total Length and Mantle Length Analysis ############# 
#######################  From Paxton 2016  ######################## 
###################################################################

library(evd)
library(pracma)
library(ggplot2)
library(spatialEco)
library(memisc)
library(evir)
library(lmom)


source("NorthropColeman2014.fns")
source('Sample_L-Moments.R')
source('L-moment Ratios Asymptotic Conditional Variances.R')



#--------------- PRELIMINARY ANALYSIS 



#--- FUNCTIONS
# Intergral of E_n statistic - Dietrich et al 2002, by Husler and Li 2006

integranda<-function(t,X,k){
        X<-sort(X);n<-length(X)
        hill_est<-max(1/k*sum(log(X[(n-k+1):n])-log(X[(n-k)])),0)
        A<-(log(X[(n-round(k*t))])-log(X[(n-k)]))/hill_est
        neg_mom_est<-min(1-1/2*(1-(hill_est)^2/(1/k*sum((log(X[(n-k+1):n])-log(X[(n-k)]))^2)))^(-1),0)
        B<-((t^(-neg_mom_est)-1)/neg_mom_est)*(1-neg_mom_est)
        return((A-B)^2*t^2)
}
gama<-function(X,k){  #Moment estimator of EVI
        X<-sort(X);n<-length(X)
        hill_est<-max(1/k*sum(log(X[(n-k+1):n])-log(X[(n-k)])),0)
        neg_mom_est<-min(1-1/2*(1-(hill_est)^2/(1/k*sum((log(X[(n-k+1):n])-log(X[(n-k)]))^2)))^(-1),0)
        return(hill_est+neg_mom_est)
}
Critical<-function(Gama,values){  #Interpolated Critical values of E_gamma depending on the Moment estimate
        cases(
                is.na(Gama)==T -> 0,
                (Gama>=0) -> values[1],
                (Gama<0 & Gama>= -0.1) -> values[2]+(Gama+0.1)/0.1*(values[2]-values[1]),
                (Gama< -.1 & Gama>= -0.2) -> values[3]+(Gama+0.2)/0.1*(values[3]-values[2]),
                (Gama< -.2 & Gama>= -0.3) -> values[4]+(Gama+0.3)/0.1*(values[4]-values[3]),
                (Gama< -.3 & Gama>= -0.4) -> values[5]+(Gama+0.4)/0.1*(values[5]-values[4]),
                (Gama< -.4 & Gama>= -0.5) -> values[6]+(Gama+0.5)/0.1*(values[6]-values[5]),
                (Gama< -.5 & Gama>= -0.6) -> values[7]+(Gama+0.6)/0.1*(values[7]-values[6]),
                (Gama< -.6 & Gama>= -0.7) -> values[8]+(Gama+0.7)/0.1*(values[8]-values[7]),
                (Gama< -0.7) -> values[8]
        )}
values95<-c(.150,.144,.141,.140,.141,.141,.144,.147) #critical values for E_gamma
                                                        #for varying gamma at 95%


# MOMENT without logs order r  (k o.s.)   
# as Eq. 21 Fraga Alves and Neves 2017

# main functional of Location Invariant Moment Estimator (EVI<0)
Nk_r <- function(x,r,k) { # x=vector# 1<= k <= n-1; r integer 
        n=length(x)
        X<-sort(x,decreasing = F)
        Nr <-mean((X[(n-k+1):n]-X[n-k])^r)
        return(Nr)
} 

# EVI<0 Moment Estimator (fixed order k)    
# Location Invariant Moment Estimator (EVI<0)
# as Eq. 19 Fraga Alves and Neves 2017
Mom.est.neg<-function(x,k) {# 2<=k k <= n-1; 
        N1=Nk_r(x,1,k);N2=Nk_r(x,2,k)
        Mom.k= 1-0.5*(1-N1^2/N2)^(-1)  
        return(Mom.k)
}

#  EVI<0 Moment  EstimatoR -  (2<k<n-1  )
# Location Invariant Moment Estimator 
MOM.neg <- function(data){
        g_MOM.neg<-c()
        n=length(data)
        for(k in 2:(n-1)) {
                g_MOM.neg[k] <- Mom.est.neg(data,k)
        }
        return(g_MOM.neg)
}

# EVI<0  Moment Endpoint Estimator 
# as x^* in Fraga Alves and Neves 2017, all k
# Eq. 18

endpoint.mom.neg<- function(data){# 2<= k <= n-1;
        X<-sort(data,decreasing = F)
        n=length(data)
        xF_mom.neg <- c()
        for(k in 2:(n-1)) {
                g.k <- Mom.est.neg(X,k)
                N1=Nk_r(X,1,k)
                xF_mom.neg[k] <- X[n-k] - N1*(1-g.k)/g.k
        }
        return(xF_mom.neg)
}



# xF_hat:= endpoint of data   (k = 1,2,...,n0:=floor(n/2))  
# X_FAN estimator at fixed k -- Fraga Alves and Neves 2017

xF_FAN_k<- function(data,k){ #k0 <= n/2
        k=floor(k/2)
        X<-sort(data)
        n=length(data)
        xF_hat <- c()
        aik<-c()
        i<-1:k 
        aik=-log((k+i-1)/(k+i))/log(2)
        xF_hat <- X[n]+sum(aik*(X[n-k]-X[n-k-i+1]))
        return(xF_hat)
}

#  all sample 

xF_FAN <- function(data){
        xF_hat <-c()
        n<-length(data)
        kpar<-seq(from = 2, to = n-2,by=2)
        for(k in kpar) {
                xF_hat[k]<-xF_FAN_k(data,k)
        } 
        return(xF_hat)
}


# h(gamma)   
# function for bias reduced version 
h <- function(g){((2^(-g)-1)/(g*log(2))+1)/g}


# a0(n/k) Estimator  at k    (2<k<n-1  ) #
# as in Eq. 20 of Fraga Alves and Neves 2017

a0_k<-function(x,k) {# 2<=k k <= n-1; 
        X <- sort(x)
        # M1 = Mk_r(X,1,k);
        N1=Nk_r(X,1,k)
        # g= Mom_est(X,k)
        g <- Mom.est.neg(X,k)
        # a0k=X[n-k]*M1*(1-g)
        a0k=N1*(1-g)
        return(a0k)
}

# Red-Bias 1 
xF_RB1 <- function(x){
        X<-sort(x)
        n=length(X)
        RB1<-c()
        for(k in 2:(n-1)) {
                # g = Mom_est(X,k) 
                g = Mom.est.neg(X,k)
                xF = xF_FAN_k(X,k)
                a0 = a0_k(X,k)
                RB1[k] <- xF - a0*h(g)
        }
        return(RB1)
}

# Red-Bias 2  
# as Eq. 13 of Fraga Alves and Neves 2017

xF_RB2 <- function(x){
        X<-sort(x)
        n=length(X)
        RB2<-c()
        for(k in 2:(n-1)) {
                #g = Mom_est(X,k) 
                g = Mom.est.neg(X,k)
                xF = xF_FAN_k(X,k)
                a0 = a0_k(X,k)
                if(g>-1/2){
                        RB2[k] <- xF - a0*h(g) - a0 * gamma(1-g)/g * (k^g)
                }
                else{
                        RB2[k] <- xF- a0*h(g)
                }
                RB2[k] <-max(RB2[k],X[n])  # RB2admissible
        }
        return(RB2)
}

# Approximated 100(1-alpha)%-Confidence Upper Bound for endpoint
# as Eq. 14 Fraga Alves and Neves 2017

UB_xF_FAN<-function(x,alpha){
        X<-sort(x)
        n=length(X)
        UB<-c()
        for(k in 2:(n-1)) {
                g = Mom.est.neg(X,k)
                xF = xF_FAN_k(X,k)
                a0 = a0_k(X,k)
                q_alpha = (-log(alpha))^(-g)/g
                UB[k] <- xF - a0*(h(g)+(k^g)*q_alpha)
        }
        return(UB)
}



#----------------------------------#
#------------- TL -----------------#
#----------------------------------#



##### Data available in the Supplementary Material of Paxton (2016)
squid <- read.table ("data squid original.csv",header=T, sep=",") 
squidtemp<-squid[is.na(squid$EL)==F,]

TL<-squidtemp$EL #Full 
TL_72<-squidtemp$EL[squidtemp$EL<16.81] #Only Reliable measurements

#hist(TL, probability=T,xlab = "Total Length (m)",main='')
#hist(TL_72, probability=T,xlab = "Total Length (m)",main='')

n<-length(TL)
n_72<-length(TL_72)


# CDF plot on log-y scale
# Weibull plotting-positions

#pdf(file="CCDF_TL.pdf")
#par(mfrow=c(1,3), mai = c(0.4, 0.4, 0.2, 0.1))
pp<-(1:n)/(n+1)
exp_quantiles<- log(1-pp)

plot(sort(TL),exp_quantiles, ylab = bquote(log[] (1-p[i])), xlab=expression(TL[i:74]), 
     mgp=c(1.8, .5, 0), col='red',pch=18, cex.lab=1.2)  #same as reference paper
box(lty='1111',col='black')

axis(side=3, at=quantile(TL,probs = c(.20,.50,.65,.80,.85,.90,.95,.99)), 
     labels=c(20,50,65,80,85,90,95,99)
     ,tick=T, mgp=c(3, .5, 0),cex.axis=1)

abline(v=quantile(TL,probs = c(.20,.50,.65,.80,.85,.90,.95,.99)), col='lightgray',lty=3)
grid(NA,NULL)  
legend("bottomleft","max(TL) = 19 m",bg="white",cex=1.5)

#dev.off()


#pdf(file="CCDF_TL_72.pdf")
pp<-(1:n_72)/(n_72+1)
exp_quantiles<- log(1-pp)

plot(sort(TL_72),exp_quantiles, ylab = bquote(log[] (1-p[i])), xlab=expression(TL[i:72]), 
     mgp=c(1.8, .5, 0), col='red',pch=18, cex.lab=1.2)  #same as reference paper
box(lty='1111',col='black')

axis(side=3, at=quantile(TL_72,probs = c(.20,.50,.65,.80,.85,.90,.95,.99)), 
     labels=c(20,50,65,80,85,90,95,99)
     ,tick=T, mgp=c(3, .5, 0),cex.axis=1)

abline(v=quantile(TL_72,probs = c(.20,.50,.65,.80,.85,.90,.95,.99)), col='lightgray',lty=3)
grid(NA,NULL)
legend("bottomleft","max(TL) = 16.764 m",bg="white", cex=1.5)
#dev.off()


# MRLP
#pdf(file="MRLP_TL.pdf")
#x<-sort(TL);n<-length(x)

x<-sort(TL_72);n<-length(x)

u<-x[1:(n-10)] # thresholds = ascending order statistics (o.s)


me_u<-c()
for(i in 1:(n-10)){
        me_u[i]<-mean(x[x>u[i]]-u[i])
}

plot(u,me_u,col='blue',xlim=c(0.5,12), xaxt='n',yaxt='n', xlab=expression(paste("threshold ",TL[i:74])), ylab = 'Mean Excess TL')
axis(side=2, at=seq(2.8,7.8,by=0.4), labels=seq(2.8,7.8,by=0.4),tick=T)
axis(side=1, at=seq(0,12,by=1), labels=seq(0,12,by=1),tick=T)
axis(side=3, at=quantile(TL,probs = c(.20,.50,.65,.80,.85,.90,.95,.99)), 
     labels=c(20,50,65,80,85,90,95,99)
     ,tick=T, mgp=c(3, .5, 0),cex.axis=1)

abline(v=quantile(TL,probs = c(.20,.50,.65,.80,.85,.90,.95,.99)), col='lightgray',lty=3)
grid(NA,NULL)
#grid(NULL,NULL)

me_u_3<-c()
for(i in 1:(n-10)){
        me_u_3[i]<-mean(x[x>u[i]]-u[i])
}

w_3<-c()
for(j in 1:(length(u))){
        w_3[j]<-(length(x[x>u[j]]))/var(x[x>u[j]]-u[j])
}


wls_fit_3<-vector(mode='list', length = (length(u)-10) )
wls_WMSE_3<-c()


for(j in (1:(n-20))){   
        wls_fit_3[[j]]<-lm(me_u_3[j:(n-10)] ~ u[j:(n-10)],
                           weights=w_3[j:(n-10)])
        
        wls_WMSE_3[j]<-sum(w_3[j:(n-10)]*(wls_fit_3[[j]]$residuals)^2
        )/sum(w_3[j:(n-10)])
}


#points(u[1:(n-20)],wls_WMSE_3, col='blue')
min_WMSE_3<-local.min.max(wls_WMSE_3, dev = mean, plot = F)$minima[1]
which(round(wls_WMSE_3,10)==round(min_WMSE_3,10))
#abline(v=u[33], col='blue')  #1.062


y0_3<-wls_fit_3[[14]]$coefficients[1]
m_3<-wls_fit_3[[14]]$coefficients[2]
segments(u[14], u[14]*m_3+y0_3, x1 = x[n-10], y1 = x[n-10]*m_3+y0_3,
         lwd=2, lty=1, col='black')



abline(v=u[14],col='black',lty=3,lwd=0.5)
text(6.3,2.9,"u*=5.51",col='black',cex=1.2)

legend("topright", legend="Best WMSE-LS fit to Mean Excess above sample points\n (Langousis et al. 2016)",
       cex=1, bg="white",
       col=c('black'),lty=1,
       lwd=2,seg.len = 1.5)
legend("bottomleft","max(TL) = 19 m",bg="white")

#dev.off()

#pdf(file="MRLP_TL_72.pdf")

x<-sort(TL_72);n<-length(x)

u<-x[1:(n-10)] # thresholds = ascending order statistics (o.s)


me_u<-c()
for(i in 1:(n-10)){
        me_u[i]<-mean(x[x>u[i]]-u[i])
}

plot(u,me_u,col='blue',xlim=c(0.5,12), xaxt='n',yaxt='n', xlab=expression(paste("threshold ",TL[i:72])), ylab = 'Mean Excess TL')
axis(side=2, at=seq(2.2,7.8,by=0.4), labels=seq(2.2,7.8,by=0.4),tick=T)
axis(side=1, at=seq(0,12,by=1), labels=seq(0,12,by=1),tick=T)
axis(side=3, at=quantile(TL_72,probs = c(.20,.50,.65,.80,.85,.90,.95,.99)), 
     labels=c(20,50,65,80,85,90,95,99)
     ,tick=T, mgp=c(3, .5, 0),cex.axis=1)

abline(v=quantile(TL_72,probs = c(.20,.50,.65,.80,.85,.90,.95,.99)), col='lightgray',lty=3)
grid(NA,NULL)

me_u_3<-c()
for(i in 1:(n-10)){
        me_u_3[i]<-mean(x[x>u[i]]-u[i])
}

w_3<-c()
for(j in 1:(length(u))){
        w_3[j]<-(length(x[x>u[j]]))/var(x[x>u[j]]-u[j])
}


wls_fit_3<-vector(mode='list', length = (length(u)-10) )
wls_WMSE_3<-c()


for(j in (1:(n-20))){   
        wls_fit_3[[j]]<-lm(me_u_3[j:(n-10)] ~ u[j:(n-10)],
                           weights=w_3[j:(n-10)])
        
        wls_WMSE_3[j]<-sum(w_3[j:(n-10)]*(wls_fit_3[[j]]$residuals)^2
        )/sum(w_3[j:(n-10)])
}


#points(u[1:(n-20)],wls_WMSE_3, col='blue')
min_WMSE_3<-local.min.max(wls_WMSE_3, dev = mean, plot = F)$minima[1]
which(round(wls_WMSE_3,10)==round(min_WMSE_3,10))
#abline(v=u[33], col='blue')  #1.062


y0_3<-wls_fit_3[[14]]$coefficients[1]
m_3<-wls_fit_3[[14]]$coefficients[2]
segments(u[14], u[14]*m_3+y0_3, x1 = x[n-10], y1 = x[n-10]*m_3+y0_3,
         lwd=2, lty=1, col='black')



abline(v=u[14],col='black',lty=3,lwd=0.5)
text(6.3,2.3,"u*=5.51",col='black',cex=1.2)

legend("topright", legend="Best WMSE-LS fit to Mean Excess above sample points\n (Langousis et al. 2016)",
       cex=1, bg="white",
       col=c('black'),lty=1,
       lwd=2,seg.len = 1.5)
legend("bottomleft","max(TL) = 16.764 m",bg="white")


#dev.off()




# SEMI-PARAMETRIC
# Testing extreme value condition E_n(k)

#pdf(file="En_k_TL.pdf", width=6*1.1,height=5*1.1)

n<-length(TL)
Gama<-c()
Enk<-c();
#Enk1<-Enk2<-Enk3<-c()
for(k in 2:(n-1)){
        Enk[k]<-integrate(integranda,0,1,X=TL,k,subdivisions = 100000)$value*k
        # Enk1[k]<-integral(integranda,0,1,method="Kronrod",X=EL,k=k)*k
        # Enk2[k]<-integral(integranda,0,1,method="Clenshaw",X=EL,k=k)*k
        # Enk3[k]<-integral(integranda,0,1,method="Simpson",X=EL,k=k)*k
        Gama[k]<-gama(TL,k)}

critical95<-Critical(Gama,values95)

plot(Enk,type = 'l',ylab=NA,xlab=NA, col="blue",lwd=2)
title(ylab=expression(E[n](k)), xlab='k',line = 2.2)
polygon(c(2,73,seq(73,2)),c(0.72,0.72,critical95[seq(73,2)]),lty = 1, 
        col="burlywood",angle=45,density=20, border=NA)

lines(Enk,type = 'l',col="blue",lwd=2)
lines(c(NaN,critical95[-1]),col='black',lty=2)

n<-length(TL_72)
Gama_72<-c()
#Enk<-c();
Enk1<-Enk2<-Enk3<-c()
for(k in 2:(n-1)){
        # Enk[k]<-integrate(integranda,0,1,X=TL_72,k,subdivisions = 100000)$value*k
        # Enk1[k]<-integral(integranda,0,1,method="Kronrod",X=TL_72,k=k)*k
        # Enk2[k]<-integral(integranda,0,1,method="Clenshaw",X=EL,k=k)*k
        Enk3[k]<-integral(integranda,0,1,method="Simpson",X=TL_72,k=k)*k
        Gama_72[k]<-gama(TL_72,k)}

critical95_72<-Critical(Gama_72,values95) #since the critical values are all very similar,
                                        # it is not significant to draw both lines

#plot(Enk,type = 'l',ylab=expression(paste("Test Statistic ",E[n](k))), xlab='k', col="blue")
#polygon(c(2,73,73,2),c(0.72,0.72,critical95[73],critical95[2]),lty = 1, 
#        col="lightsteelblue1",angle=45,density=30)
lines(Enk3, col='orange',lwd=2)
#lines(c(NaN,critical95[-1]),col='red',lty=2)
legend("topleft", legend=c(expression(paste(E[n](k)," with max(TL) = 16.764 m")),expression(paste(E[n](k)," with max(TL) = 19 m")),expression(paste("95% critical values of ",E[hat(xi)]," "))),
       col=c('orange','blue','black'),cex=1,lty=c(1,1,2),lwd=c(2,2,1),bg="white", bty='o',title.adj=0.1)

#dev.off()




# Choice of Domain of attraction
# n<-length(TL)
n<-length(TL_72)
K<-1:(n-1)

#XX<-sort(TL)   #ordered sample
XX<-sort(TL_72)

# main functional of Location Invariant Moment Estimator
Njk<-matrix(seq(1:(2*(n-1))),2,n-1)
for(j in 1:2){
        for(k in K){
                Njk[j,k]<-Nk_r(XX,j,k)
        }
}



#pdf(file="Choice_Domain_TL_72.pdf", width=9*1.1,height=3*1.1)
par(mfrow=c(1,3),mai = c(0.5, 0.55, 0.1, 0.05))
#Greenwood Statistic Gr=N2/(N1^2) ;  Gr*=sqrt(k/4)*(Gr-2)
Gr<-c()
Gr_ast<-c()
for(k in K){
        Gr[k]<-Njk[2,k]/(Njk[1,k]^2)
        Gr_ast[k]<-sqrt(k/4)*(Gr[k]-2)
}
plot(K,Gr_ast,type='n', ylim=c(-3.5,0.5), xlab=NA,
     ylab=NA)
title(ylab=expression(paste(Gr[n]^"*"*(k))), line=2.2, cex.lab=1.2, xlab="k")
axis(2,at=seq(-3.5,0.5,0.5), labels = seq(-3.5,0.5,0.5) )
polygon(c(0,72,72,0),c(qnorm(0.05),qnorm(0.05),-3.6,-3.6),lty = 1, 
        col="lightsteelblue1",angle=45,density=30)
points(K,Gr_ast, type='l', col='blue',lwd=2)
abline(h=qnorm(0.05), lty=2, lwd=2)


text(67,-1.45,expression(z[0.05]),cex=1.5)
legend("bottomleft","max(TL) = 16.764 m",bg="white",cex=1.2)

#Hasofer-Wang Statistic W=(1/k)*(1/(Gr-1)) ; W*=sqrt(k/4)*(k*W-1)
W<-c()
W_ast<-c()
for(k in K){
        W[k]<-(1/k)*(1/(Gr[k]-1))
        W_ast[k]<-sqrt(k/4)*(k*W[k]-1)
}
plot(K,W_ast,type='n', ylim=c(-2,15), xlab=NA,
     ylab=NA, yaxt="n" )
title(ylab=expression(paste(W[n]^"*"*(k))), line=2.2, cex.lab=1.2,xlab="k")
axis(2,at=seq(-2,15,2), labels = seq(-2,15,2) )
polygon(c(0,72,72,0),c(qnorm(0.95),qnorm(0.95),15.5,15.5),lty = 1, 
        col="peachpuff",angle=45,density=30)
points(K,W_ast, type='l', col='orange',lwd=2)
abline(h=qnorm(0.95), lty=2, lwd=2)
text(67,2.5,expression(z[0.95]),cex=1.5)
legend("bottomleft","max(TL) = 16.764 m",bg="white",cex=1.2)






#Ratio Statistic R=(XX[n]-XX[n-k])/N1 ;  R*=R-log(k)
R<-c()
R_ast<-c()
for(k in K){
        R[k]<-(XX[n]-XX[n-k])/Njk[1,k]
        R_ast[k]<-R[k]-log(k)
}

plot(K,R_ast,type='n', ylim=c(-2,1.5), xlab=NA,
     ylab=NA, yaxt="n" )
title(ylab=expression(paste(R[n]^"*"*(k))), line=2.2, cex.lab=1.2,xlab="k")
axis(2,at=seq(-2,1.5,0.5), labels = seq(-2,1.5,0.5) )
polygon(c(0,72,72,0),c(qgumbel(0.05),qgumbel(0.05),-2.1,-2.1),lty = 1, 
        col="palegreen1",angle=45,density=30)
points(K,R_ast, type='l', col='forestgreen',lwd=2)
abline(h=qgumbel(0.05), lty=2, lwd=2)
text(67,-0.95,expression(g[0.05]),cex=1.5)
legend("bottomleft","max(TL) = 16.764 m",bg="white",cex=1.2)

#dev.off()


#bilateral
#pdf(file="Choice_Domain_Bilateral_TL_72.pdf", width=9*1.1,height=3*1.1)
par(mfrow=c(1,3),mai = c(0.5, 0.55, 0.1, 0.05))
#Greenwood Statistic Gr=N2/(N1^2) ;  Gr*=sqrt(k/4)*(Gr-2)
plot(K,Gr_ast,type='n', ylim=c(-3.5,2), xlab=NA,
     ylab=NA)
title(ylab=expression(paste(Gr[n]^"*"*(k))), line=2.2, cex.lab=1.2, xlab="k")
axis(2,at=seq(-3.5,2,0.5), labels = seq(-3.5,2,0.5) )
polygon(c(0,72,72,0),c(qnorm(1-0.05/2),qnorm(1-0.05/2),2.2,2.2),lty = 1, 
        col="lightsteelblue1",angle=45,density=30)
polygon(c(0,72,72,0),c(-qnorm(1-0.05/2),-qnorm(1-0.05/2),-3.6,-3.6),lty = 1, 
        col="lightsteelblue1",angle=45,density=30)
points(K,Gr_ast, type='l', col='blue',lwd=2)
abline(h=qnorm(1-0.05/2), lty=2, lwd=2)
abline(h=-qnorm(1-0.05/2), lty=2, lwd=2)
text(67,-1.7,expression(z[0.025]),cex=1.5)
text(67,1.7,expression(z[0.975]),cex=1.5)
legend("topleft","max(TL) = 16.764 m",bg="white",cex=1.2)

#Hasofer-Wang Statistic W=(1/k)*(1/(Gr-1)) ; W*=sqrt(k/4)*(k*W-1)
plot(K,W_ast,type='n', ylim=c(-3,15), xlab=NA,
     ylab=NA, yaxt="n" )
title(ylab=expression(paste(W[n]^"*"*(k))), line=2.2, cex.lab=1.2,xlab="k")
axis(2,at=seq(-2,15,2), labels = seq(-2,15,2) )
polygon(c(0,72,72,0),c(qnorm(0.975),qnorm(0.975),15.5,15.5),lty = 1, 
        col="peachpuff",angle=45,density=30)
polygon(c(0,72,72,0),c(qnorm(0.025),qnorm(0.025),-3.5,-3.5),lty = 1, 
        col="peachpuff",angle=45,density=30)
points(K,W_ast, type='l', col='orange',lwd=2)
abline(h=qnorm(1-0.05/2), lty=2, lwd=2)
abline(h=-qnorm(1-0.05/2), lty=2, lwd=2)
text(67,-1.5,expression(z[0.025]),cex=1.5)
text(67,1.5,expression(z[0.975]),cex=1.5)
legend("bottomleft","max(TL) = 16.764 m",bg="white",cex=1.2)


#Ratio Statistic R=(XX[n]-XX[n-k])/N1 ;  R*=R-log(k)
plot(K,R_ast,type='n', ylim=c(-2,4), xlab=NA,
     ylab=NA, yaxt="n" )
title(ylab=expression(paste(R[n]^"*"*(k))), line=2.2, cex.lab=1.2,xlab="k")
axis(2,at=seq(-2,4,1), labels = seq(-2,4,1) )
polygon(c(0,72,72,0),c(qgumbel(0.025),qgumbel(0.025),-2.1,-2.1),lty = 1, 
        col="palegreen1",angle=45,density=30)
polygon(c(0,72,72,0),c(qgumbel(0.975),qgumbel(0.975),4.1,4.1),lty = 1, 
        col="palegreen1",angle=45,density=30)
points(K,R_ast, type='l', col='forestgreen',lwd=2)
abline(h=qgumbel(0.025), lty=2, lwd=2)
abline(h=qgumbel(0.975), lty=2, lwd=2)
text(67,-1.1,expression(g[0.025]),cex=1.5)
text(67,3.6,expression(g[0.975]),cex=1.5)
legend("topleft","max(TL) = 16.764 m",bg="white",cex=1.2)

#dev.off()











n<-length(TL)
K<-1:(n-1)

XX<-sort(TL)   #ordered sample


# main functional of Location Invariant Moment Estimator
Njk<-matrix(seq(1:(2*(n-1))),2,n-1)
for(j in 1:2){
        for(k in K){
                Njk[j,k]<-Nk_r(XX,j,k)
        }
}


#pdf(file="Choice_Domain_TL.pdf", width=9*1.1,height=3*1.1)
par(mfrow=c(1,3),mai = c(0.5, 0.55, 0.1, 0.05))

# Greenwood Statistic Gr=N2/(N1^2) ;  Gr*=sqrt(k/4)*(Gr-2)
Gr<-c()
Gr_ast<-c()
for(k in K){
        Gr[k]<-Njk[2,k]/(Njk[1,k]^2)
        Gr_ast[k]<-sqrt(k/4)*(Gr[k]-2)
}
plot(K,Gr_ast,type='n', ylim=c(-3.5,0.5), xlab=NA,
     ylab=NA)
title(ylab=expression(paste(Gr[n]^"*"*(k))), line=2.2, cex.lab=1.2, xlab="k")
axis(2,at=seq(-3.5,0.5,0.5), labels = seq(-3.5,0.5,0.5) )
polygon(c(0,74,74,0),c(qnorm(0.05),qnorm(0.05),-3.6,-3.6),lty = 1, 
        col="lightsteelblue1",angle=45,density=30)
points(K,Gr_ast, type='l', col='blue',lwd=2)
abline(h=qnorm(0.05), lty=2, lwd=2)


text(67,-1.45,expression(z[0.05]),cex=1.5)
legend("bottomleft","max(TL) = 19 m",bg="white",cex=1.2)

# Hasofer-Wang Statistic W=(1/k)*(1/(Gr-1)) ; W*=sqrt(k/4)*(k*W-1)
W<-c()
W_ast<-c()
for(k in K){
        W[k]<-(1/k)*(1/(Gr[k]-1))
        W_ast[k]<-sqrt(k/4)*(k*W[k]-1)
}
plot(K,W_ast,type='n', ylim=c(-2,15), xlab=NA,
     ylab=NA, yaxt="n" )
title(ylab=expression(paste(W[n]^"*"*(k))), line=2.2, cex.lab=1.2,xlab="k")
axis(2,at=seq(-2,15,2), labels = seq(-2,15,2) )
polygon(c(0,74,74,0),c(qnorm(0.95),qnorm(0.95),15.5,15.5),lty = 1, 
        col="peachpuff",angle=45,density=30)
points(K,W_ast, type='l', col='orange',lwd=2)
abline(h=qnorm(0.95), lty=2, lwd=2)
text(67,2.5,expression(z[0.95]),cex=1.5)
legend("bottomleft","max(TL) = 19 m",bg="white",cex=1.2)


 
# Ratio Statistic R=(XX[n]-XX[n-k])/N1 ;  R*=R-log(k)
R<-c()
R_ast<-c()
for(k in K){
        R[k]<-(XX[n]-XX[n-k])/Njk[1,k]
        R_ast[k]<-R[k]-log(k)
}

plot(K,R_ast,type='n', ylim=c(-2,1.5), xlab=NA,
     ylab=NA, yaxt="n" )
title(ylab=expression(paste(R[n]^"*"*(k))), line=2.2, cex.lab=1.2,xlab="k")
axis(2,at=seq(-2,1.5,0.5), labels = seq(-2,1.5,0.5) )
polygon(c(0,74,74,0),c(qgumbel(0.05),qgumbel(0.05),-2.1,-2.1),lty = 1, 
        col="palegreen1",angle=45,density=30)
points(K,R_ast, type='l', col='forestgreen',lwd=2)
abline(h=qgumbel(0.05), lty=2, lwd=2)
text(67,-0.95,expression(g[0.05]),cex=1.5)
legend("bottomleft","max(TL) = 19 m",bg="white",cex=1.2)

#dev.off()




#bilateral
#pdf(file="Choice_Domain_Bilateral_TL.pdf", width=9*1.1,height=3*1.1)
par(mfrow=c(1,3),mai = c(0.5, 0.55, 0.1, 0.05))

# Greenwood Statistic Gr=N2/(N1^2) ;  Gr*=sqrt(k/4)*(Gr-2)
plot(K,Gr_ast,type='n', ylim=c(-3.5,2), xlab=NA,
     ylab=NA)
title(ylab=expression(paste(Gr[n]^"*"*(k))), line=2.2, cex.lab=1.2, xlab="k")
axis(2,at=seq(-3.5,2,0.5), labels = seq(-3.5,2,0.5) )
polygon(c(0,74,74,0),c(qnorm(1-0.05/2),qnorm(1-0.05/2),2.2,2.2),lty = 1, 
        col="lightsteelblue1",angle=45,density=30)
polygon(c(0,74,74,0),c(-qnorm(1-0.05/2),-qnorm(1-0.05/2),-3.6,-3.6),lty = 1, 
        col="lightsteelblue1",angle=45,density=30)
points(K,Gr_ast, type='l', col='blue',lwd=2)
abline(h=qnorm(1-0.05/2), lty=2, lwd=2)
abline(h=-qnorm(1-0.05/2), lty=2, lwd=2)
text(67,-1.7,expression(z[0.025]),cex=1.5)
text(67,1.7,expression(z[0.975]),cex=1.5)
legend("topleft","max(TL) = 19 m",bg="white",cex=1.2)

# Hasofer-Wang Statistic W=(1/k)*(1/(Gr-1)) ; W*=sqrt(k/4)*(k*W-1)
plot(K,W_ast,type='n', ylim=c(-3,15), xlab=NA,
     ylab=NA, yaxt="n" )
title(ylab=expression(paste(W[n]^"*"*(k))), line=2.2, cex.lab=1.2,xlab="k")
axis(2,at=seq(-2,15,2), labels = seq(-2,15,2) )
polygon(c(0,74,74,0),c(qnorm(0.975),qnorm(0.975),15.5,15.5),lty = 1, 
        col="peachpuff",angle=45,density=30)
polygon(c(0,74,74,0),c(qnorm(0.025),qnorm(0.025),-3.5,-3.5),lty = 1, 
        col="peachpuff",angle=45,density=30)
points(K,W_ast, type='l', col='orange',lwd=2)
abline(h=qnorm(1-0.05/2), lty=2, lwd=2)
abline(h=-qnorm(1-0.05/2), lty=2, lwd=2)
text(67,-1.5,expression(z[0.025]),cex=1.5)
text(67,1.5,expression(z[0.975]),cex=1.5)
legend("topleft","max(TL) = 19 m",bg="white",cex=1.2)


# Ratio Statistic R=(XX[n]-XX[n-k])/N1 ;  R*=R-log(k)
plot(K,R_ast,type='n', ylim=c(-2,4), xlab=NA,
     ylab=NA, yaxt="n" )
title(ylab=expression(paste(R[n]^"*"*(k))), line=2.2, cex.lab=1.2,xlab="k")
axis(2,at=seq(-2,4,1), labels = seq(-2,4,1) )
polygon(c(0,74,74,0),c(qgumbel(0.025),qgumbel(0.025),-2.1,-2.1),lty = 1, 
        col="palegreen1",angle=45,density=30)
polygon(c(0,74,74,0),c(qgumbel(0.975),qgumbel(0.975),4.1,4.1),lty = 1, 
        col="palegreen1",angle=45,density=30)
points(K,R_ast, type='l', col='forestgreen',lwd=2)
abline(h=qgumbel(0.025), lty=2, lwd=2)
abline(h=qgumbel(0.975), lty=2, lwd=2)
text(67,-1.1,expression(g[0.025]),cex=1.5)
text(67,3.6,expression(g[0.975]),cex=1.5)
legend("topleft","max(TL) = 19 m",bg="white",cex=1.2)

#dev.off()








# Testing Finiteness of endpoint sample path
#n<-length(TL)
n<-length(TL_72)
K<-1:(n-1)

#XX<-sort(TL)   #ordered sample
XX<-sort(TL_72)

# Equivalent to Moment statistics M_nk^(j)
logNjk<-matrix(seq(1:(2*(n-1))),2,n-1)
for(j in 1:2){ for(k in K){
        logNjk[j,k]<-mean((log(XX[(n-k+1):n])-log(XX[n-k]))^j)
}}

a_n<-c();TT1<-c();TT1_ast<-c()

for(k in K){
        a_n[k]<-XX[n-k]*logNjk[1,k]*0.5*(1-logNjk[1,k]^2/logNjk[2,k])^(-1)
        TT1[k]<-(1/k)*(sum((XX[(n-k):(n-1)]-XX[n-k]-a_n[k])/(XX[n]-XX[n-k])))
        TT1_ast[k]<-sqrt(k)*log(k)*TT1[k]
}


TT2<-c();TT2_ast<-c()
for(k in (2:(n-1))){
        i<-1:(k-1)
        TT2[k]<-(1/k)*(sum(i*(XX[n-i+1]-XX[n-i])/XX[n-k]))
        TT2_ast[k]<-sqrt(k)*(log(n/k)*TT2[k]-1)
}

#pdf(file="Test_Finiteness_Endpoint_TL_72.pdf", width=6*1.1,height=5*1.1)
plot(K, TT1_ast, ylim=c(-24,4), xlab=NA, ylab=NA,yaxt='n',
     type='n')
title(ylab=expression(T[1]*"* and "*T[2]*"* statistics"), line=2.2, cex.lab=1.2,xlab="k")
axis(2,at=seq(-24,4,4), labels = seq(-24,4,4) )

polygon(c(0,72,72,0),c(qnorm(0.025),qnorm(0.025),-24.1,-24.1),lty = 1, 
        col="lightsteelblue1",angle=45,density=30)
polygon(c(0,72,72,0),c(qnorm(0.975),qnorm(0.975),4.1,4.1),lty = 1, 
        col="lightsteelblue1",angle=45,density=30)
lines(TT1_ast,  col='navy',lwd=2,lty=1)
lines(TT2_ast, col='blue',lwd=2)
(krej1<-which(TT1_ast<=qnorm(0.05)))
(krej2<-which(TT2_ast>=qnorm(1-0.05)))
abline(h=qnorm(0.025), lty=2, lwd=2)
abline(h=qnorm(1-0.025), lty=2, lwd=2)

text(67,-3,expression(z[0.025]),cex=1.5)
text(67,1,expression(z[0.975]),cex=1.5)
legend("bottomleft",title="max(TL) = 16.764 m",bg="white",cex=1.2, title.adj = 0.5,
       legend=c(expression(T[1]*"*"),expression(T[2]*"*")),col=c("navy","blue"),
       lwd=c(2,2),lty=c(1,1))
#dev.off()




n<-length(TL)
K<-1:(n-1)

XX<-sort(TL)   #ordered sample

# Equivalent to Moment statistics M_nk^(j)
logNjk<-matrix(seq(1:(2*(n-1))),2,n-1)
for(j in 1:2){ for(k in K){
        logNjk[j,k]<-mean((log(XX[(n-k+1):n])-log(XX[n-k]))^j)
}}

a_n<-c();TT1<-c();TT1_ast<-c()

for(k in K){
        a_n[k]<-XX[n-k]*logNjk[1,k]*0.5*(1-logNjk[1,k]^2/logNjk[2,k])^(-1)
        TT1[k]<-(1/k)*(sum((XX[(n-k):(n-1)]-XX[n-k]-a_n[k])/(XX[n]-XX[n-k])))
        TT1_ast[k]<-sqrt(k)*log(k)*TT1[k]
}


TT2<-c();TT2_ast<-c()
for(k in (2:(n-1))){
        i<-1:(k-1)
        TT2[k]<-(1/k)*(sum(i*(XX[n-i+1]-XX[n-i])/XX[n-k]))
        TT2_ast[k]<-sqrt(k)*(log(n/k)*TT2[k]-1)
}

#pdf(file="Test_Finiteness_Endpoint_TL.pdf", width=6*1.1,height=5*1.1)
plot(K, TT1_ast, ylim=c(-24,4), xlab=NA, ylab=NA,yaxt='n',
     type='n')
title(ylab=expression(T[1]*"* and "*T[2]*"* statistics"), line=2.2, cex.lab=1.2,xlab="k")
axis(2,at=seq(-24,4,4), labels = seq(-24,4,4) )

polygon(c(0,74,74,0),c(qnorm(0.025),qnorm(0.025),-24.1,-24.1),lty = 1, 
        col="lightsteelblue1",angle=45,density=30)
polygon(c(0,74,74,0),c(qnorm(0.975),qnorm(0.975),4.1,4.1),lty = 1, 
        col="lightsteelblue1",angle=45,density=30)
lines(TT1_ast,  col='navy',lwd=2,lty=1)
lines(TT2_ast, col='blue',lwd=2)
(krej1<-which(TT1_ast<=qnorm(0.05)))
(krej2<-which(TT2_ast>=qnorm(1-0.05)))
abline(h=qnorm(0.025), lty=2, lwd=2)
abline(h=qnorm(1-0.025), lty=2, lwd=2)

text(67,-3,expression(z[0.025]),cex=1.5)
text(67,1,expression(z[0.975]),cex=1.5)
legend("bottomleft",title="max(TL) = 19 m",bg="white",cex=1.2, title.adj = 0.5,
       legend=c(expression(T[1]*"*"),expression(T[2]*"*")),col=c("navy","blue"),
       lwd=c(2,2),lty=c(1,1))
#dev.off()




# Estimation of right endpoint sample path
#n<-length(TL)
n<-length(TL_72)

#XX<-sort(TL)   #ordered sample
XX<-sort(TL_72)


pdf(file="Endpoint_Estimate_TL_72.pdf", width=6*1.1,height=5*1.1)
endpoint_TL_72<-xF_FAN(XX)
BR_endpoint_TL_72<-xF_RB2(XX)
#UB_TL_72<-UB_xF_FAN(XX,0.05)
plot(BR_endpoint_TL_72,ylim=c(16,45),type='l', col='blue',lwd=2, ylab=NA, xlab=NA)
title(ylab=expression(hat(x)^F),xlab="k",line=2.2, cex.lab=1.2)
kpar<-seq(from = 2, to = n-2,by=2)
lines(kpar,endpoint_TL_72[is.na(endpoint_TL_72)==F],col="orange",lwd=2)
#lines(kpar,UB_TL_72[kpar],col="red")
abline(h=max(TL_72),lty=2,lwd=2)
legend("topleft",bg="white",cex=1.2,
       legend=c("max(TL) = 16.764 m",expression(hat(x)[FAN]^F),expression(hat(x)[RB2]^F)),
       col=c("black","orange","blue"),
       lwd=c(2,2,2),lty=c(2,1,1))

dev.off()

pdf(file="Endpoint_Estimate_ZOOM_TL_72.pdf", width=6*1.1,height=5*1.1)
endpoint_TL_72<-xF_FAN(XX)
BR_endpoint_TL_72<-xF_RB2(XX)
#UB_TL_72<-UB_xF_FAN(XX,0.05)
plot(BR_endpoint_TL_72,ylim=c(15,25),type='l', col='blue',lwd=2, ylab=NA, xlab=NA)
title(ylab=expression(hat(x)^F),xlab="k",line=2.2, cex.lab=1.2)
kpar<-seq(from = 2, to = n-2,by=2)
lines(kpar,endpoint_TL_72[is.na(endpoint_TL_72)==F],col="orange",lwd=2)
#lines(kpar,UB_TL_72[kpar],col="red")
abline(h=max(TL_72),lty=2,lwd=2)
abline(h=max(TL),lty=3,lwd=2)
legend("topleft",bg="white",cex=1.2,
       legend=c(expression(hat(x)[FAN]^F),expression(hat(x)[RB2]^F)),
       col=c("orange","blue"),
       lwd=c(2,2),lty=c(1,1))
legend("bottomleft", legend = "max(TL) = 16.764 m",col="black",lwd=2,lty=2)
text(69,19.3,"19 m",col="black",cex=1.2)

dev.off()




n<-length(TL)
XX<-sort(TL)

pdf(file="Endpoint_Estimate_ZOOM_TL.pdf", width=6*1.1,height=5*1.1)
endpoint_TL<-xF_FAN(XX)
BR_endpoint_TL<-xF_RB2(XX)
#UB_TL_72<-UB_xF_FAN(XX,0.05)
plot(BR_endpoint_TL,ylim=c(15,30),type='l', col='blue',lwd=2, ylab=NA, xlab=NA)
title(ylab=expression(hat(x)^F),xlab="k",line=2.2, cex.lab=1.2)
kpar<-seq(from = 2, to = n-2,by=2)
lines(kpar,endpoint_TL[is.na(endpoint_TL)==F],col="orange",lwd=2)
#lines(kpar,UB_TL_72[kpar],col="red")
abline(h=max(TL),lty=2,lwd=2)
legend("bottomright",bg="white",cex=1.2,
       legend=c(expression(hat(x)[FAN]^F),expression(hat(x)[RB2]^F)),
       col=c("orange","blue"),
       lwd=c(2,2),lty=c(1,1))
legend("bottomleft", legend = "max(TL) = 19 m",col="black",lwd=2,lty=2)


dev.off()




# General Right Endpoint -based test of EVI
# General Right Endpoint Estimator xF
n<-length(TL_72)
XX<-sort(TL_72)

K2<-1:(n/2-1)
gen_xF<-c()
for(k in K2){
        gen_xF[k]<-XX[n]+XX[n-k]-(1/log(2))*sum(log(1+1/(k+(0:(k-1))))*XX[(n-k):(n-k-k+1)])
}

# Gnk Test Statistics
Gnk<-c()
Gnk0_ast<-c()
for(k in K2){
        Gnk[k]<-(gen_xF[k]-XX[n-k])/(XX[n-k]-XX[n-2*k])
        Gnk0_ast[k]<-log(2)*Gnk[k]-(log(k)+log(2)/2)
        
}

pdf(file="Test_Domain_Endpoint_TL_72.pdf", width=6*1.1,height=5*1.1)
plot(2*K2,Gnk0_ast,type='n',col='red3', lwd=2, ylab=NA, xlab=NA)
title(ylab=expression(G[list(n,k)]^paste('*')*paste('(0)')),xlab="k",line=2.2, cex.lab=1.2)
polygon(c(0,72,72,0),c(qgumbel(0.05),qgumbel(0.05),-2.6,-2.6),lty = 1, 
        col="lightsalmon",angle=45,density=30)
lines(2*K2,Gnk0_ast,type='l',col='red3', lwd=2)
legend("bottomleft","max(TL) = 16.764 m",bg="white",cex=1.2)
abline(h=qgumbel(0.05), lty=2,lwd=2)
text(68,-0.95,expression(g[0.05]),cex=1.5)
dev.off()












#----------------------------------#
#------------- MtL -----------------#
#----------------------------------#



##### Reload data
squid <- read.table ("data squid original.csv",header=T, sep=",")
squidtemp<-squid[is.na(squid$FinalML)==F,]

ML<-squidtemp$FinalML[squidtemp$FinalML<3] #Only Reliable measurements

#hist(ML, probability=T,xlab = "Total Length (m)",main='')

n<-length(ML)



# CDF plot on log-y scale
# Weibull plotting-positions

#pdf(file="CCDF_ML.pdf")
pp<-(1:n)/(n+1)
exp_quantiles<- log(1-pp)

plot(sort(ML),exp_quantiles, ylab = bquote(log[] (1-p[i])), xlab=expression(MtL[i:163]), 
     mgp=c(1.8, .5, 0), col='red',pch=18,cex.lab=1.2)  #same as reference paper
box(lty='1111',col='black')

axis(side=3, at=quantile(ML,probs = c(.20,.50,.65,.80,.85,.90,.95,.99)), 
     labels=c(20,50,65,80,85,90,95,99)
     ,tick=T, mgp=c(3, .5, 0),cex.axis=1)

abline(v=quantile(ML,probs = c(.20,.50,.65,.80,.85,.90,.95,.99)), col='lightgray',lty=3)
grid(NA,NULL)  
legend("bottomleft","max(MtL) = 2.794 m",bg="white", cex=1.5)


#dev.off()


# MRLP

#pdf(file="MRLP_ML.pdf")
x<-sort(ML);n<-length(x)


u<-x[1:(n-10)] # tresholds = ascending order statistics (o.s)


me_u<-c()
for(i in 1:(n-10)){
        me_u[i]<-mean(x[x>u[i]]-u[i])
}

plot(u,me_u,col='blue',xlim=c(0.01,2.15), xaxt='n',yaxt='n', xlab=expression(paste("threshold ",MtL[i:163])), ylab = 'Mean Excess MtL')
axis(side=2, at=seq(0,1.4,by=0.2), labels=seq(0,1.4,by=0.2),tick=T)
axis(side=1, at=seq(0,2.15,by=0.15), labels=seq(0,2.15,by=0.15),tick=T)
axis(side=3, at=quantile(ML,probs = c(.20,.50,.65,.80,.85,.90,.95,.99)), 
     labels=c(20,50,65,80,85,90,95,99)
     ,tick=T, mgp=c(3, .5, 0),cex.axis=1)

abline(v=quantile(ML,probs = c(.20,.50,.65,.80,.85,.90,.95,.99)), col='lightgray',lty=3)
grid(NA,NULL)

me_u_3<-c()
for(i in 1:(n-10)){
        me_u_3[i]<-mean(x[x>u[i]]-u[i])
}

w_3<-c()
for(j in 1:(length(u))){
        w_3[j]<-(length(x[x>u[j]]))/var(x[x>u[j]]-u[j])
}


wls_fit_3<-vector(mode='list', length = (length(u)-10) )
wls_WMSE_3<-c()


for(j in (1:(n-20))){   
        wls_fit_3[[j]]<-lm(me_u_3[j:(n-10)] ~ u[j:(n-10)],
                           weights=w_3[j:(n-10)])
        
        wls_WMSE_3[j]<-sum(w_3[j:(n-10)]*(wls_fit_3[[j]]$residuals)^2
        )/sum(w_3[j:(n-10)])
}


#points(u[1:(n-20)],wls_WMSE_3, col='blue')
min_WMSE_3<-local.min.max(wls_WMSE_3, dev = mean, plot = F)$minima[1]
which(round(wls_WMSE_3,10)==round(min_WMSE_3,10))
#abline(v=u[33], col='blue')  #1.062


y0_3<-wls_fit_3[[42]]$coefficients[1]
m_3<-wls_fit_3[[42]]$coefficients[2]
segments(u[42], u[42]*m_3+y0_3, x1 = x[n-10], y1 = x[n-10]*m_3+y0_3,
         lwd=2, lty=1, col='black')



abline(v=u[42],col='black',lty=3,lwd=0.5)
text(1.25,0.22,"u*=1.07",col='black',cex=1.2)

legend("topright", legend="Best WMSE-LS fit to Mean Excess above sample points\n (Langousis et al. 2016)",
       cex=1, bg="white",
       col=c('black'),lty=1,
       lwd=2,seg.len = 1.5)
legend("bottomleft","max(MtL) = 2.794 m",bg="white")

dev.off()



# SEMI-PARAMETRIC
# Testing extreme value condition E_n(k)

#pdf(file="En_k_ML.pdf", width=6*1.1,height=5*1.1)

n<-length(ML)
Gama<-c()
#Enk<-c();
Enk1<-Enk2<-Enk3<-c()
for(k in 2:(n-1)){
        # Enk[k]<-integrate(integranda,0,1,X=ML,k,subdivisions = 100000)$value*k
        # Enk1[k]<-integral(integranda,0,1,method="Kronrod",X=EL,k=k)*k
        # Enk2[k]<-integral(integranda,0,1,method="Clenshaw",X=EL,k=k)*k
         Enk3[k]<-integral(integranda,0,1,method="Simpson",X=ML,k=k)*k
        Gama[k]<-gama(ML,k)}


critical95<-Critical(Gama,values95)

plot(Enk3,type = 'l',ylab=NA,xlab=NA, col="forestgreen",lwd=2, xaxt='n')
axis(1,at=seq(0,170,20), labels = seq(0,170,20))
title(ylab=expression(E[n](k)), xlab='k',line = 2.2)
polygon(c(2,162,seq(162,2)),c(0.595,0.595,critical95[seq(162,2)]),lty = 1, 
        col="darkolivegreen1",angle=45,density=20, border=NA)
lines(c(NaN,critical95[-1]),col='black',lty=2)
lines(Enk3,type = 'l', col="forestgreen",lwd=2)


legend("topleft", legend=c(expression(paste(E[n](k)," with max(MtL) = 2.794 m")),expression(paste("95% critical values of ",E[hat(xi)]," "))),
       col=c('forestgreen','black'),cex=1,lty=c(1,2),lwd=c(2,1), bty='o',bg="white",title.adj=0.1)

#dev.off()


# Choice of Domain of attraction
n<-length(ML)
K<-1:(n-1)

XX<-sort(ML)   #ordered sample


# main functional of Location Invariant Moment Estimator
Njk<-matrix(seq(1:(2*(n-1))),2,n-1)
for(j in 1:2){
        for(k in K){
                Njk[j,k]<-Nk_r(XX,j,k)
        }
}


#pdf(file="Choice_Domain_ML.pdf", width=9*1.1,height=3*1.1)
par(mfrow=c(1,3),mai = c(0.5, 0.55, 0.1, 0.05))

# Greenwood Statistic Gr=N2/(N1^2) ;  Gr*=sqrt(k/4)*(Gr-2)
Gr<-c()
Gr_ast<-c()
for(k in K){
        Gr[k]<-Njk[2,k]/(Njk[1,k]^2)
        Gr_ast[k]<-sqrt(k/4)*(Gr[k]-2)
}
plot(K,Gr_ast,type='n', ylim=c(-5.5,0.5), xlab=NA,
     ylab=NA, yaxt="n", xaxt='n')
title(ylab=expression(paste(Gr[n]^"*"*(k))), line=2.2, cex.lab=1.2, xlab="k")
axis(2,at=seq(-5.5,0.5,1), labels = seq(-5.5,0.5,1) )
axis(1,at=seq(0,170,20), labels = seq(0,170,20))
polygon(c(0,163,163,0),c(qnorm(0.05),qnorm(0.05),-5.6,-5.6),lty = 1, 
        col="lightsteelblue1",angle=45,density=30)
points(K,Gr_ast, type='l', col='blue',lwd=2)
abline(h=qnorm(0.05), lty=2, lwd=2)
text(155,-1.4,expression(z[0.05]),cex=1.5)
legend("bottomleft","max(MtL) = 2.794 m",bg="white",cex=1.2)

# Hasofer-Wang Statistic W=(1/k)*(1/(Gr-1)) ; W*=sqrt(k/4)*(k*W-1)
W<-c()
W_ast<-c()
for(k in K){
        W[k]<-(1/k)*(1/(Gr[k]-1))
        W_ast[k]<-sqrt(k/4)*(k*W[k]-1)
}
plot(K,W_ast,type='n', ylim=c(-2,38), xlab=NA,
     ylab=NA, yaxt="n", xaxt='n')
title(ylab=expression(paste(W[n]^"*"*(k))), line=2.2, cex.lab=1.2,xlab="k")
axis(2,at=seq(-5,38,5), labels = seq(-5,38,5) )
axis(1,at=seq(0,170,20), labels = seq(0,170,20))
polygon(c(0,163,163,0),c(qnorm(0.95),qnorm(0.95),38.5,38.5),lty = 1, 
        col="peachpuff",angle=45,density=30)
points(K,W_ast, type='l', col='orange',lwd=2)
abline(h=qnorm(0.95), lty=2, lwd=2)
text(155,3,expression(z[0.95]),cex=1.5)
legend("topleft","max(MtL) = 2.794 m",bg="white",cex=1.2)


# Ratio Statistic R=(XX[n]-XX[n-k])/N1 ;  R*=R-log(k)
R<-c()
R_ast<-c()
for(k in K){
        R[k]<-(XX[n]-XX[n-k])/Njk[1,k]
        R_ast[k]<-R[k]-log(k)
}

plot(K,R_ast,type='n', ylim=c(-3,1.5), xlab=NA,
     ylab=NA, yaxt="n" , xaxt='n')
axis(1,at=seq(0,170,20), labels = seq(0,170,20))
title(ylab=expression(paste(R[n]^"*"*(k))), line=2.2, cex.lab=1.2,xlab="k")
axis(2,at=seq(-3.5,1.5,0.5), labels = seq(-3.5,1.5,0.5) )
polygon(c(0,163,163,0),c(qgumbel(0.05),qgumbel(0.05),-3.1,-3.1),lty = 1, 
        col="palegreen1",angle=45,density=30)
points(K,R_ast, type='l', col='forestgreen',lwd=2)
abline(h=qgumbel(0.05), lty=2, lwd=2)
text(155,-0.95,expression(g[0.05]),cex=1.5)
legend("bottomleft","max(MtL) = 2.794 m",bg="white",cex=1.2)

#dev.off()




# bilateral
#pdf(file="Choice_Domain_Bilateral_ML.pdf", width=9*1.1,height=3*1.1)
par(mfrow=c(1,3),mai = c(0.5, 0.55, 0.1, 0.05))

# Greenwood Statistic Gr=N2/(N1^2) ;  Gr*=sqrt(k/4)*(Gr-2)
plot(K,Gr_ast,type='n', ylim=c(-5.5,2.5), xlab=NA,
     ylab=NA, yaxt="n", xaxt='n')
title(ylab=expression(paste(Gr[n]^"*"*(k))), line=2.2, cex.lab=1.2, xlab="k")
axis(2,at=seq(-5.5,2.5,1), labels = seq(-5.5,2.5,1) )
axis(1,at=seq(0,170,20), labels = seq(0,170,20))
polygon(c(0,163,163,0),c(qnorm(1-0.05/2),qnorm(1-0.05/2),2.7,2.7),lty = 1, 
        col="lightsteelblue1",angle=45,density=30)
polygon(c(0,163,163,0),c(-qnorm(1-0.05/2),-qnorm(1-0.05/2),-5.6,-5.6),lty = 1, 
        col="lightsteelblue1",angle=45,density=30)
points(K,Gr_ast, type='l', col='blue',lwd=2)
abline(h=qnorm(1-0.05/2), lty=2, lwd=2)
abline(h=-qnorm(1-0.05/2), lty=2, lwd=2)
text(155,-1.45,expression(z[0.025]),cex=1.5)
text(155,1.45,expression(z[0.975]),cex=1.5)
legend("bottomleft","max(MtL) = 2.794 m",bg="white",cex=1.2)




# Hasofer-Wang Statistic W=(1/k)*(1/(Gr-1)) ; W*=sqrt(k/4)*(k*W-1)
plot(K,W_ast,type='n', ylim=c(-2,38), xlab=NA,
     ylab=NA, yaxt="n", xaxt='n')
title(ylab=expression(paste(W[n]^"*"*(k))), line=2.2, cex.lab=1.2,xlab="k")
axis(2,at=seq(-5,38,5), labels = seq(-5,38,5) )
axis(1,at=seq(0,170,20), labels = seq(0,170,20))
polygon(c(0,163,163,0),c(qnorm(0.975),qnorm(0.975),38.5,38.5),lty = 1, 
        col="peachpuff",angle=45,density=30)
polygon(c(0,163,163,0),c(qnorm(0.025),qnorm(0.025),-3,-3),lty = 1, 
        col="peachpuff",angle=45,density=30)
points(K,W_ast, type='l', col='orange',lwd=2)
abline(h=qnorm(1-0.05/2), lty=2, lwd=2)
abline(h=-qnorm(1-0.05/2), lty=2, lwd=2)
text(155,0,expression(z[0.025]),cex=1.5)
text(155,4,expression(z[0.975]),cex=1.5)
legend("topleft","max(MtL) = 2.794 m",bg="white",cex=1.2)


# Ratio Statistic R=(XX[n]-XX[n-k])/N1 ;  R*=R-log(k)
plot(K,R_ast,type='n', ylim=c(-3,4), xlab=NA,
     ylab=NA, yaxt="n" , xaxt='n')
axis(1,at=seq(0,170,20), labels = seq(0,170,20))
title(ylab=expression(paste(R[n]^"*"*(k))), line=2.2, cex.lab=1.2,xlab="k")
axis(2,at=seq(-3.5,4,1), labels = seq(-3.5,4,1) )
polygon(c(0,163,163,0),c(qgumbel(0.025),qgumbel(0.025),-3.1,-3.1),lty = 1, 
        col="palegreen1",angle=45,density=30)
polygon(c(0,163,163,0),c(qgumbel(0.975),qgumbel(0.975),4.1,4.1),lty = 1, 
        col="palegreen1",angle=45,density=30)
points(K,R_ast, type='l', col='forestgreen',lwd=2)
abline(h=qgumbel(0.025), lty=2, lwd=2)
abline(h=qgumbel(0.975), lty=2, lwd=2)
text(155,-1.1,expression(g[0.025]),cex=1.5)
text(155,3.55,expression(g[0.975]),cex=1.5)
legend("bottomleft","max(MtL) = 2.794 m",bg="white",cex=1.2)

#dev.off()


# Testing Finiteness of Endpoint
n<-length(ML)
K<-1:(n-1)

XX<-sort(ML)   #ordered sample

# Equivalent to Moment statistics M_nk^(j)
logNjk<-matrix(seq(1:(2*(n-1))),2,n-1)
for(j in 1:2){ for(k in K){
        logNjk[j,k]<-mean((log(XX[(n-k+1):n])-log(XX[n-k]))^j)
}}

a_n<-c();TT1<-c();TT1_ast<-c()

for(k in K){
        a_n[k]<-XX[n-k]*logNjk[1,k]*0.5*(1-logNjk[1,k]^2/logNjk[2,k])^(-1)
        TT1[k]<-(1/k)*(sum((XX[(n-k):(n-1)]-XX[n-k]-a_n[k])/(XX[n]-XX[n-k])))
        TT1_ast[k]<-sqrt(k)*log(k)*TT1[k]
}


TT2<-c();TT2_ast<-c()
for(k in (2:(n-1))){
        i<-1:(k-1)
        TT2[k]<-(1/k)*(sum(i*(XX[n-i+1]-XX[n-i])/XX[n-k]))
        TT2_ast[k]<-sqrt(k)*(log(n/k)*TT2[k]-1)
}

#pdf(file="Test_Finiteness_Endpoint_ML.pdf", width=6*1.1,height=5*1.1)
plot(K, TT1_ast, ylim=c(-24,4), xlab=NA, ylab=NA,yaxt='n',
     type='n', xaxt='n')
axis(1,at=seq(0,170,20), labels = seq(0,170,20))
title(ylab=expression(T[1]*"* and "*T[2]*"* statistics"), line=2.2, cex.lab=1.2,xlab="k")
axis(2,at=seq(-24,4,4), labels = seq(-24,4,4) )

polygon(c(0,163,163,0),c(qnorm(0.025),qnorm(0.025),-24.1,-24.1),lty = 1, 
        col="lightsteelblue1",angle=45,density=30)
polygon(c(0,163,163,0),c(qnorm(0.975),qnorm(0.975),4.1,4.1),lty = 1, 
        col="lightsteelblue1",angle=45,density=30)
lines(TT1_ast,  col='navy',lwd=2,lty=1)
lines(TT2_ast, col='blue',lwd=2)
(krej1<-which(TT1_ast<=qnorm(0.05)))
(krej2<-which(TT2_ast>=qnorm(1-0.05)))
abline(h=qnorm(0.025), lty=2, lwd=2)
abline(h=qnorm(1-0.025), lty=2, lwd=2)

text(155,-3,expression(z[0.025]),cex=1.5)
text(155,1,expression(z[0.975]),cex=1.5)
legend("bottomleft",title="max(MtL) = 2.794 m",bg="white",cex=1.2, title.adj = 0.5,
       legend=c(expression(T[1]*"*"),expression(T[2]*"*")),col=c("navy","blue"),
       lwd=c(2,2),lty=c(1,1))
#dev.off()



# Estimation of right endpoint sample path
n<-length(ML)
XX<-sort(ML)   #ordered sample


pdf(file="Endpoint_Estimate_ML.pdf", width=6*1.1,height=5*1.1)
endpoint_ML<-xF_FAN(XX)
BR_endpoint_ML<-xF_RB2(XX)
#UB_TL_72<-UB_xF_FAN(XX,0.05)
plot(BR_endpoint_ML,type='l', col='blue',lwd=2, ylab=NA, xlab=NA, xaxt='n')
axis(1,at=seq(0,170,20), labels = seq(0,170,20))
title(ylab=expression(hat(x)^F),xlab="k",line=2.2, cex.lab=1.2)
kpar<-seq(from = 2, to = n-2,by=2)
lines(kpar,endpoint_ML[is.na(endpoint_ML)==F],col="orange",lwd=2)
#lines(kpar,UB_TL_72[kpar],col="red")
abline(h=max(ML),lty=2,lwd=2)
legend("topright",bg="white",cex=1.2,
       legend=c("max(MtL) = 2.794 m",expression(hat(x)[FAN]^F),expression(hat(x)[RB2]^F)),
       col=c("black","orange","blue"),
       lwd=c(2,2,2),lty=c(2,1,1))

dev.off()

pdf(file="Endpoint_Estimate_ZOOM_ML.pdf", width=6*1.1,height=5*1.1)
plot(BR_endpoint_ML,ylim=c(2.5,5),type='l', col='blue',lwd=2, ylab=NA, xlab=NA,xaxt='n')
axis(1,at=seq(0,170,20), labels = seq(0,170,20))
title(ylab=expression(hat(x)^F),xlab="k",line=2.2, cex.lab=1.2)
kpar<-seq(from = 2, to = n-2,by=2)
lines(kpar,endpoint_ML[is.na(endpoint_ML)==F],col="orange",lwd=2)
#lines(kpar,UB_TL_72[kpar],col="red")
abline(h=max(ML),lty=2,lwd=2)
abline(h=3.35,lty=3,lwd=2)
legend("topright",bg="white",cex=1.2,
       legend=c(expression(hat(x)[FAN]^F),expression(hat(x)[RB2]^F)),
       col=c("orange","blue"),
       lwd=c(2,2),lty=c(1,1))
legend("bottomleft", legend = "max(MtL) = 2.794 m",col="black",lwd=2,lty=2)
text(155,3.45,"3.35 m",col="black",cex=1.2)

dev.off()




n<-length(ML)
XX<-sort(ML)

K2<-1:(n/2)
gen_xF<-c()
for(k in K2){
        gen_xF[k]<-XX[n]+XX[n-k]-(1/log(2))*sum(log(1+1/(k+(0:(k-1))))*XX[(n-k):(n-k-k+1)])
}

# Gnk Test Statistics
Gnk<-c()
Gnk0_ast<-c()
for(k in K2){
        Gnk[k]<-(gen_xF[k]-XX[n-k])/(XX[n-k]-XX[n-2*k])
        Gnk0_ast[k]<-log(2)*Gnk[k]-(log(k)+log(2)/2)
        
}

pdf(file="Test_Domain_Endpoint_ML.pdf", width=6*1.1,height=5*1.1)
plot(2*K2,Gnk0_ast,type='n',col='red3', lwd=2, ylab=NA, xlab=NA,xaxt='n')
axis(1,at=seq(0,170,20), labels = seq(0,170,20))
title(ylab=expression(G[list(n,k)]^paste('*')*paste('(0)')),xlab="k",line=2.2, cex.lab=1.2)
polygon(c(0,163,163,0),c(qgumbel(0.05),qgumbel(0.05),-5,-5),lty = 1, 
        col="lightsalmon",angle=45,density=30)
lines(2*K2,Gnk0_ast,type='l',col='red3', lwd=2)
legend("topright","max(MtL) = 2.794 m",bg="white",cex=1.2)
abline(h=qgumbel(0.05), lty=2,lwd=2)
text(155,-0,expression(g[0.05]),cex=1.5)
dev.off()









#--------------- THRESHOLD AND ESTIMATION










#----------------------------------#
#------------- TL -----------------#
#----------------------------------#


##### Reload data
squid <- read.table ("data squid original.csv",header=T, sep=",")
squidtemp<-squid[is.na(squid$EL)==F,]

TL<-squidtemp$EL #Full 
TL_72<-squidtemp$EL[squidtemp$EL<16.81] #Only Reliable measurements

#hist(TL, probability=T,xlab = "Total Length (m)",main='')
#hist(TL_72, probability=T,xlab = "Total Length (m)",main='')

n<-length(TL)
n_72<-length(TL_72)


# ALRSM
#pdf(file="Lmom_TL_72_62cand.pdf")

# LMRD
tau.4<-function(tau.3){tau.3*(1+5*tau.3)/(5+tau.3)}
curve(tau.4(x),-0.2,0.4,xlab = expression(tau[3]), 
      ylab = expression(tau[4]),cex.lab=1.5)
abline(h=0,v=0,col="grey")



x<-sort(TL_72)
#x<-sort(TL)
n<-length(x)

t4_line<- t3_line <- numeric()

tt3<-c();tt4<-c() 
var_t4dadott3<-c()
var_t3dadott4<-c()

#Candidate thresholds
#u<-x[c(1:(n-10))]
#u<- quantile(x,probs = seq(0.25,0.925,by=0.075), names=F)  #10 quantis
#u<- quantile(x,probs = seq(0.25,0.953,by=0.037), names=F)   #20 quantis
#u<-quantile(x, probs = seq(0.20,0.85,by=0.05), names=F)   #14 quantis
u<-x[1:62]   #62 candidates (Leave 10 out)

I<-1:length(u)

for (i in I) {
        x_exc <- x[x > u[i]] 
        Xu=x_exc-u[i]          
        nu<-length(Xu)
        
        #sample l-mom ratios for the current threshold
        tt4[i]=l4(Xu)/l2(Xu)
        tt3[i]=l3(Xu)/l2(Xu)
       
        
        #optimize min distance in 2dim to the theoretical line
        ff <- function(x){x*(1+5*x)/(5+x)}    #relação tau3 tau4
        #ff <- function(x){x^2+1}
        obj <- function(x, xi,yi) {
                sum((xi-x)^2+(yi-ff(x))^2)
        }       
        
        xmin <- optimize(f=obj, c(-1, 1), tol = 0.0000001, xi= tt3[i],yi= tt4[i])

        
        t3_line[i] <-xmin$minimum
        t4_line[i] <- ff(xmin$minimum)
        points(t3_line[i],t4_line[i],col="blue",pch=16)
        
        #### Estimated Root Variance with t4=tau.4(t3)
        TT1<- T_mat_1(Xu)
        #EVal_t4dadott3[i] <- tt4[i]  #this way the points would always fall inside
        var_t4dadott3[i] <- sqrt(TT1[2,2] /nu ) * sqrt(1-(TT1[1,2]^2/(TT1[1,1]*TT1[2,2])))
        
        #### Estimated Root Variance with  t3= tau.3(t4)
        TT2<- T_mat_2(Xu)
        var_t3dadott4[i] <- sqrt(TT2[1,1] /nu ) * sqrt(1-(TT2[1,2]^2/(TT2[1,1]*TT2[2,2])))
        
        #points(tau_3,tau_4,col="red",pch=20)
        rm(x_exc)
        rm(Xu)
}

# length(t3)
# length(t4)
# cat(t3,"\n",t4)

#pontos originais
lines(tt3[I],tt4[I], type="p",col="navy") #, pch=20)

error_line <- (tt3-t3_line)^2+(tt4-t4_line)^2
index_line=1:length(error_line)
i0_line <- index_line[error_line[index_line]==min(error_line[which(is.na(error_line)==F)])]
#points(t3_line[i0_line],t4_line[i0_line],col="red",pch=20, cex=2)
points(tt3[i0_line],tt4[i0_line],col="red",pch=20, cex=1)
(u_opt_L<-u[i0_line][1]) # optimal threshold



#title("L-moment Ratio Diagram - 10 Candidate Thresholds")
#legend(0.25,0.5,expression(paste("(", t[3], "," ,t[4],") for u* = ", 10.44691 %~~% TL[57:74])),cex = 0.8, pt.cex = 1,
#       bty='n',pch=21, col=c("navy"), pt.bg=c("red"))
lines(tt3,tau.4(tt3)+qnorm(0.975)*var_t4dadott3, type="l", lty= 2, lwd= 2 ,col="orange") 
lines(tt3,tau.4(tt3)+qnorm(0.025)*var_t4dadott3, type="l", lty= 2, lwd= 2 ,col="orange")


#ESTIMATED Variance with t4 and t3=tau.3(t4)
#no mesmo plot t3 dado t4
lines(tau.3(tt4[I])+qnorm(0.975)*var_t3dadott4, tt4[I], type="l", lty= 2, lwd= 2 ,col="green") 
lines(tau.3(tt4[I])+qnorm(0.025)*var_t3dadott4, tt4[I], type="l", lty= 2, lwd= 2 ,col="green")


legend("topleft",legend=c(expression(paste("95% Conf. Band for ",tau[4]," given ",t[3])),
                       expression(paste("95% Conf. Band for ",tau[3]," given ",t[4])), 
                       expression(paste("(", t[3], "," ,t[4],") for u* = ", 6.2))),
       title= "max(TL) = 16.764 m; 62 candidates",title.adj = 0.15,
       col=c('orange','green','navy'), lty =c(2,2,0),pch = c(NA,NA,21), bty='o', pt.cex = c(1,1,1.2), 
       pt.bg=c(NA,NA,'red'),cex=1.2,bg="white",lwd=c(2,2,NA))

dev.off()


############################ GPD fit above u_10
############################
u_10<-u_opt_L
(fit_10<-evd::fpot(TL_72,u_10,model="gpd",npp=1,std.err = T, mper= 1000))
(fit_10_scale<-evd::fpot(TL_72,u_10,model="gpd",std.err = T))
(xF_10<-u_10-fit_10_scale$scale/fit_10_scale$param[2])
(sd_xF_10<-sqrt(fit_10_scale$var.cov[1,1]*(-1/fit_10_scale$param[2])^2 + 
        fit_10_scale$var.cov[2,2]*(fit_10_scale$param[1]/fit_10_scale$param[2]^2)^2+
        2* fit_10_scale$var.cov[1,2]*(-1/fit_10_scale$param[2])*(fit_10_scale$param[1]/fit_10_scale$param[2]^2)))

par(mfrow=c(2,2))
plot(fit_10)
par(mfrow=c(1,1))


############################ GPD fit above u_20
############################
u_20<-u_opt_L
(fit_20<-evd::fpot(TL_72,u_20,model="gpd",npp=1,std.err = T, mper= 1000))
(fit_20_scale<-evd::fpot(TL_72,u_20,model="gpd",std.err = T))
(xF_20<-u_20-fit_20_scale$scale/fit_20_scale$param[2])
(sd_xF_20<-sqrt(fit_20_scale$var.cov[1,1]*(-1/fit_20_scale$param[2])^2 + 
                fit_20_scale$var.cov[2,2]*(fit_20_scale$param[1]/fit_20_scale$param[2]^2)^2+
                2* fit_20_scale$var.cov[1,2]*(-1/fit_20_scale$param[2])*(fit_20_scale$param[1]/fit_20_scale$param[2]^2)))

par(mfrow=c(2,2))
plot(fit_20)
par(mfrow=c(1,1))



############################ GPD fit above u_14
############################
u_14<-u_opt_L
(fit_14<-evd::fpot(TL_72,u_14,model="gpd",npp=1,std.err = T, mper= 1000))
(fit_14_scale<-evd::fpot(TL_72,u_14,model="gpd",std.err = T))
(xF_14<-u_14-fit_14_scale$scale/fit_14_scale$param[2])
(sd_xF_14<-sqrt(fit_14_scale$var.cov[1,1]*(-1/fit_14_scale$param[2])^2 + 
                fit_14_scale$var.cov[2,2]*(fit_14_scale$param[1]/fit_14_scale$param[2]^2)^2+
                2* fit_14_scale$var.cov[1,2]*(-1/fit_14_scale$param[2])*(fit_14_scale$param[1]/fit_14_scale$param[2]^2)))

par(mfrow=c(2,2))
plot(fit_14)
par(mfrow=c(1,1))





# ALRSM
pdf(file="Lmom_TL.pdf")
# LMDR

tau.4<-function(tau.3){tau.3*(1+5*tau.3)/(5+tau.3)}
curve(tau.4(x),-0.1,0.55,xlab = expression(tau[3]), 
      ylab = expression(tau[4]),cex.lab=1.5)
abline(h=0,v=0,col="grey")


x<-sort(TL)
n<-length(TL)

t4_line<- t3_line <- numeric()

tt3<-c();tt4<-c() 
var_t4dadott3<-c()
var_t3dadott4<-c()

#Candidate thresholds
#u<-x[c(1:(n-10))]
#u<- quantile(x,probs = seq(0.25,0.925,by=0.075), names=F)  #10 quantis
#u<- quantile(x,probs = seq(0.25,0.953,by=0.037), names=F)   #20 quantis
u<-quantile(x, probs = seq(0.20,0.85,by=0.05), names=F)   #14 quantis


I<-1:length(u)

for (i in I) {
        x_exc <- x[x > u[i]] 
        Xu=x_exc-u[i]          
        nu<-length(Xu)
        
        #sample l-mom ratios for the current threshold
        tt4[i]=l4(Xu)/l2(Xu)
        tt3[i]=l3(Xu)/l2(Xu)
        
        
        #optimize min distance in 2dim to the theoretical line
        ff <- function(x){x*(1+5*x)/(5+x)}    
        #ff <- function(x){x^2+1}
        obj <- function(x, xi,yi) {
                sum((xi-x)^2+(yi-ff(x))^2)
        }       
        
        xmin <- optimize(f=obj, c(-1, 1), tol = 0.0000001, xi= tt3[i],yi= tt4[i])
        
        
        t3_line[i] <-xmin$minimum
        t4_line[i] <- ff(xmin$minimum)
        points(t3_line[i],t4_line[i],col="blue",pch=16)
        
        #### Estimated Root Variance with  t4=tau.4(t3)
        TT1<- T_mat_1(Xu)
        #EVal_t4dadott3[i] <- tt4[i]  #this way the points would always fall inside
        var_t4dadott3[i] <- sqrt(TT1[2,2] /nu ) * sqrt(1-(TT1[1,2]^2/(TT1[1,1]*TT1[2,2])))
        
        #### Estimated Root Variance with  t3= tau.3(t4)
        TT2<- T_mat_2(Xu)
        var_t3dadott4[i] <- sqrt(TT2[1,1] /nu ) * sqrt(1-(TT2[1,2]^2/(TT2[1,1]*TT2[2,2])))
        
        #points(tau_3,tau_4,col="red",pch=20)
        rm(x_exc)
        rm(Xu)
}

# length(t3)
# length(t4)
# cat(t3,"\n",t4)

#pontos originais
lines(tt3[I],tt4[I], type="p",col="navy") #, pch=20)

error_line <- (tt3-t3_line)^2+(tt4-t4_line)^2
index_line=1:length(error_line)
i0_line <- index_line[error_line[index_line]==min(error_line[which(is.na(error_line)==F)])]
#points(t3_line[i0_line],t4_line[i0_line],col="red",pch=20, cex=2)
points(tt3[i0_line],tt4[i0_line],col="red",pch=20, cex=1)
(u_opt_L<-u[i0_line]) # optimal threshold



#title("L-moment Ratio Diagram - 10 Candidate Thresholds")
#legend(0.25,0.5,expression(paste("(", t[3], "," ,t[4],") for u* = ", 10.44691 %~~% TL[57:74])),cex = 0.8, pt.cex = 1,
#       bty='n',pch=21, col=c("navy"), pt.bg=c("red"))
lines(tt3,tau.4(tt3)+qnorm(0.975)*var_t4dadott3, type="l", lty= 2, lwd= 2 ,col="orange") 
lines(tt3,tau.4(tt3)+qnorm(0.025)*var_t4dadott3, type="l", lty= 2, lwd= 2 ,col="orange")


#ESTIMATED Variance with t4 and t3=tau.3(t4)
#no mesmo plot t3 dado t4
lines(tau.3(tt4[I])+qnorm(0.975)*var_t3dadott4, tt4[I], type="l", lty= 2, lwd= 2 ,col="green") 
lines(tau.3(tt4[I])+qnorm(0.025)*var_t3dadott4, tt4[I], type="l", lty= 2, lwd= 2 ,col="green")


legend("topleft",legend=c(expression(paste("95% Conf. Band for ",tau[4]," given ",t[3])),
                          expression(paste("95% Conf. Band for ",tau[3]," given ",t[4])), 
                          expression(paste("(", t[3], "," ,t[4],") for u* = ", 6.102))),
       title= "max(TL) = 19 m; 14 candidates",title.adj = 0.15,
       col=c('orange','green','navy'), lty =c(2,2,0),pch = c(NA,NA,21), bty='o', pt.cex = c(1,1,1.2), 
       pt.bg=c(NA,NA,'red'),cex=1.2,bg="white",lwd=c(2,2,NA))

dev.off()


############################ GPD fit above u_10
############################
u_10<-u_opt_L
(fit_10<-evd::fpot(TL,u_10,model="gpd",npp=1,std.err = T, mper= 1000))
(fit_10_scale<-evd::fpot(TL,u_10,model="gpd",std.err = T))
(xF_10<-u_10-fit_10_scale$scale/fit_10_scale$param[2])
(sd_xF_10<-sqrt(fit_10_scale$var.cov[1,1]*(-1/fit_10_scale$param[2])^2 + 
                fit_10_scale$var.cov[2,2]*(fit_10_scale$param[1]/fit_10_scale$param[2]^2)^2+
                2* fit_10_scale$var.cov[1,2]*(-1/fit_10_scale$param[2])*(fit_10_scale$param[1]/fit_10_scale$param[2]^2)))

par(mfrow=c(2,2))
plot(fit_10)
par(mfrow=c(1,1))


############################ GPD fit above u_20
############################
u_20<-u_opt_L[1]
(fit_20<-evd::fpot(TL,u_20,model="gpd",npp=1,std.err = T, mper= 1000))
(fit_20_scale<-evd::fpot(TL,u_20,model="gpd",std.err = T))
(xF_20<-u_20-fit_20_scale$scale/fit_20_scale$param[2])
(sd_xF_20<-sqrt(fit_20_scale$var.cov[1,1]*(-1/fit_20_scale$param[2])^2 + 
                fit_20_scale$var.cov[2,2]*(fit_20_scale$param[1]/fit_20_scale$param[2]^2)^2+
                2* fit_20_scale$var.cov[1,2]*(-1/fit_20_scale$param[2])*(fit_20_scale$param[1]/fit_20_scale$param[2]^2)))

par(mfrow=c(2,2))
plot(fit_20)
par(mfrow=c(1,1))



############################ GPD fit above u_14
############################
u_14<-u_opt_L
(fit_14<-evd::fpot(TL,u_14,model="gpd",npp=1,std.err = T, mper= 1000))
(fit_14_scale<-evd::fpot(TL,u_14,model="gpd",std.err = T))
(xF_14<-u_14-fit_14_scale$scale/fit_14_scale$param[2])
(sd_xF_14<-sqrt(fit_14_scale$var.cov[1,1]*(-1/fit_14_scale$param[2])^2 + 
                fit_14_scale$var.cov[2,2]*(fit_14_scale$param[1]/fit_14_scale$param[2]^2)^2+
                2* fit_14_scale$var.cov[1,2]*(-1/fit_14_scale$param[2])*(fit_14_scale$param[1]/fit_14_scale$param[2]^2)))

par(mfrow=c(2,2))
plot(fit_14)
par(mfrow=c(1,1))




####################################################################
############################################  Bader et al. (2018)
####################################################################
x<-sort(TL_72);n<-length(TL_72)
#u<- quantile(x,probs = seq(0.25,0.925,by=0.075), names=F)  #10 quantis
#u<- quantile(x,probs = seq(0.25,0.953,by=0.037), names=F)   #20 quantis
u<-quantile(x, probs = seq(0.20,0.85,by=0.05), names=F)   #14 quantis
bader_p_val <- tryCatch(eva::gpdSeqTests(x,u,method='ad', nsim=100), error= function (e) NA, warning = function(w) NA)
(FS <- tryCatch(eva::pSeqStop(bader_p_val$p.values), error= function (e) NA, warning = function(w) NA))



####################################################################
############################################  Northrop e Coleman
####################################################################
x<-sort(TL_72);n<-length(TL_72)
u1<- quantile(x,probs = seq(0.25,0.925,by=0.075), names=F)  #10 quantis
u2<- quantile(x,probs = seq(0.25,0.953,by=0.037), names=F)   #20 quantis
u3<-quantile(x, probs = seq(0.20,0.85,by=0.05), names=F)   #14 quantis


pdf(file="Northrop_TL_72.pdf",width=9*.9,height=3*.9)
par(mfrow=c(1,3),mai=c(0.55,0.5,0.25,0.05))
northrop <- score.fitrange(x,u1)
abline(h=0.05,lty=3,col='grey')
text(11.5,0.1,"0.05", col="grey",cex=1.2)
abline(v=u_10,lty=2,col='red')
northrop <- score.fitrange(x,u2)
abline(h=0.05,lty=3,col='grey')
text(12.5,0.1,"0.05", col="grey",cex=1.2)
abline(v=u_20,lty=2,col='red')
northrop <- score.fitrange(x,u3)
abline(h=0.05,lty=3,col='grey')
text(10,0.1,"0.05", col="grey",cex=1.2)
abline(v=u_14,lty=2,col='red')
dev.off()

pdf(file="Northrop_TL_72_62.pdf",width=6,height=6)
u<-x[1:62]   #62 candidates (Leave 10 out)
northrop <- score.fitrange(x,unique(u)) #56 unique threshold
abline(h=0.05,lty=3,col='grey')
text(11,0.1,"0.05", col="grey",cex=1.2)
abline(v=unique(u)[25],lty=2,col='red')
dev.off()

############################ GPD fit above u_10
############################
u_10<-u1[2]
(fit_10<-evd::fpot(x,u_10,model="gpd",npp=1,std.err = T, mper= 1000))
(fit_10_scale<-evd::fpot(x,u_10,model="gpd",std.err = T))
(xF_10<-u_10-fit_10_scale$scale/fit_10_scale$param[2])
(sd_xF_10<-sqrt(fit_10_scale$var.cov[1,1]*(-1/fit_10_scale$param[2])^2 + 
                        fit_10_scale$var.cov[2,2]*(fit_10_scale$param[1]/fit_10_scale$param[2]^2)^2+
                        2* fit_10_scale$var.cov[1,2]*(-1/fit_10_scale$param[2])*(fit_10_scale$param[1]/fit_10_scale$param[2]^2)))
# (dRL_dsigma_10<-(1/fit_10_scale$param[2])*(((72/16)*(1/1000))^(-fit_10_scale$param[2])-1))
# (dRL_dxi_10<-(-fit_10_scale$scale/fit_10_scale$param[2]^2)*(((72/16)*(1/1000))^(-fit_10_scale$param[2])-1)-
#         (fit_10_scale$scale/fit_10_scale$param[2])*(((72/16)*(1/1000))^(-fit_10_scale$param[2]))*log((72/16)*(1/1000)))
# (sd_RL_10<-sqrt(fit_10_scale$var.cov[1,1]*dRL_dsigma_10^2 + fit_10_scale$var.cov[2,2]*dRL_dxi_10^2 +
#                         2*dRL_dsigma_10*dRL_dxi_10*fit_10_scale$var.cov[2,1]))


############################ GPD fit above u_20
############################
u_20<-u2[3]
(fit_20<-evd::fpot(x,u_20,model="gpd",npp=1,std.err = T, mper= 1000))
(fit_20_scale<-evd::fpot(x,u_20,model="gpd",std.err = T))
(xF_20<-u_20-fit_20_scale$scale/fit_20_scale$param[2])
(sd_xF_20<-sqrt(fit_20_scale$var.cov[1,1]*(-1/fit_20_scale$param[2])^2 + 
                        fit_20_scale$var.cov[2,2]*(fit_20_scale$param[1]/fit_20_scale$param[2]^2)^2+
                        2* fit_20_scale$var.cov[1,2]*(-1/fit_20_scale$param[2])*(fit_20_scale$param[1]/fit_20_scale$param[2]^2)))

############################ GPD fit above u_14
############################
u_14<-u3[4]
(fit_14<-evd::fpot(x,u_14,model="gpd",npp=1,std.err = T, mper= 1000))
(fit_14_scale<-evd::fpot(x,u_14,model="gpd",std.err = T))
(xF_14<-u_14-fit_14_scale$scale/fit_14_scale$param[2])
(sd_xF_14<-sqrt(fit_14_scale$var.cov[1,1]*(-1/fit_14_scale$param[2])^2 + 
                        fit_14_scale$var.cov[2,2]*(fit_14_scale$param[1]/fit_14_scale$param[2]^2)^2+
                        2* fit_14_scale$var.cov[1,2]*(-1/fit_14_scale$param[2])*(fit_14_scale$param[1]/fit_14_scale$param[2]^2)))





####################################################################
############################################  Bader et al. (2018)
####################################################################
x<-sort(TL);n<-length(TL)
#u<- quantile(x,probs = seq(0.25,0.925,by=0.075), names=F)  #10 quantis
#u<- quantile(x,probs = seq(0.25,0.953,by=0.037), names=F)   #20 quantis
u<-quantile(x, probs = seq(0.20,0.85,by=0.05), names=F)   #14 quantis
bader_p_val <- tryCatch(eva::gpdSeqTests(x,u,method='ad', nsim=100), error= function (e) NA, warning = function(w) NA)
(FS <- tryCatch(eva::pSeqStop(bader_p_val$p.values), error= function (e) NA, warning = function(w) NA))
############################ GPD fit above u_14
############################
u_14<-u[4]
(fit_14<-evd::fpot(x,u_14,model="gpd",npp=1,std.err = T, mper= 1000))
(fit_14_scale<-evd::fpot(x,u_14,model="gpd",std.err = T))
(xF_14<-u_14-fit_14_scale$scale/fit_14_scale$param[2])
(sd_xF_14<-sqrt(fit_14_scale$var.cov[1,1]*(-1/fit_14_scale$param[2])^2 + 
                        fit_14_scale$var.cov[2,2]*(fit_14_scale$param[1]/fit_14_scale$param[2]^2)^2+
                        2* fit_14_scale$var.cov[1,2]*(-1/fit_14_scale$param[2])*(fit_14_scale$param[1]/fit_14_scale$param[2]^2)))



####################################################################
############################################  Northrop e Coleman
####################################################################
x<-sort(TL);n<-length(TL)
u<- quantile(x,probs = seq(0.25,0.925,by=0.075), names=F)  #10 quantis
#u<- quantile(x,probs = seq(0.25,0.953,by=0.037), names=F)   #20 quantis
#u<-quantile(x, probs = seq(0.20,0.85,by=0.05), names=F)   #14 quantis

pdf(file="Northrop_TL.pdf",width=6,height=6)
northrop <- score.fitrange(x,u)
abline(h=0.05,lty=3,col='grey')
text(10,0.1,"0.05", col="grey",cex=1.2)
abline(v=u_14,lty=2,col='red')
dev.off()

northrop <- score.fitrange(x,u)
abline(h=0.05,lty=3,col='grey')

############################ GPD fit above u_10
############################
u_10<-u[2]
(fit_10<-evd::fpot(x,u_10,model="gpd",npp=1,std.err = T, mper= 1000))
(fit_10_scale<-evd::fpot(x,u_10,model="gpd",std.err = T))
(xF_10<-u_10-fit_10_scale$scale/fit_10_scale$param[2])
(sd_xF_10<-sqrt(fit_10_scale$var.cov[1,1]*(-1/fit_10_scale$param[2])^2 + 
                        fit_10_scale$var.cov[2,2]*(fit_10_scale$param[1]/fit_10_scale$param[2]^2)^2+
                        2* fit_10_scale$var.cov[1,2]*(-1/fit_10_scale$param[2])*(fit_10_scale$param[1]/fit_10_scale$param[2]^2)))
(dRL_dsigma_10<-(1/fit_10_scale$param[2])*(((72/16)*(1/1000))^(-fit_10_scale$param[2])-1))
(dRL_dxi_10<-(-fit_10_scale$scale/fit_10_scale$param[2]^2)*(((72/16)*(1/1000))^(-fit_10_scale$param[2])-1)-
                (fit_10_scale$scale/fit_10_scale$param[2])*(((72/16)*(1/1000))^(-fit_10_scale$param[2]))*log((72/16)*(1/1000)))
(sd_RL_10<-sqrt(fit_10_scale$var.cov[1,1]*dRL_dsigma_10^2 + fit_10_scale$var.cov[2,2]*dRL_dxi_10^2 +
                        2*dRL_dsigma_10*dRL_dxi_10*fit_10_scale$var.cov[2,1]))


############################ GPD fit above u_20
############################
u_20<-u[3]
(fit_20<-evd::fpot(x,u_20,model="gpd",npp=1,std.err = T, mper= 1000))
(fit_20_scale<-evd::fpot(x,u_20,model="gpd",std.err = T))
(xF_20<-u_20-fit_20_scale$scale/fit_20_scale$param[2])
(sd_xF_20<-sqrt(fit_20_scale$var.cov[1,1]*(-1/fit_20_scale$param[2])^2 + 
                        fit_20_scale$var.cov[2,2]*(fit_20_scale$param[1]/fit_20_scale$param[2]^2)^2+
                        2* fit_20_scale$var.cov[1,2]*(-1/fit_20_scale$param[2])*(fit_20_scale$param[1]/fit_20_scale$param[2]^2)))

############################ GPD fit above u_14
############################
u_14<-u[4]
(fit_14<-evd::fpot(x,u_14,model="gpd",npp=1,std.err = T, mper= 1000))
(fit_14_scale<-evd::fpot(x,u_14,model="gpd",std.err = T))
(xF_14<-u_14-fit_14_scale$scale/fit_14_scale$param[2])
(sd_xF_14<-sqrt(fit_14_scale$var.cov[1,1]*(-1/fit_14_scale$param[2])^2 + 
                        fit_14_scale$var.cov[2,2]*(fit_14_scale$param[1]/fit_14_scale$param[2]^2)^2+
                        2* fit_14_scale$var.cov[1,2]*(-1/fit_14_scale$param[2])*(fit_14_scale$param[1]/fit_14_scale$param[2]^2)))









#----------------------------------#
#------------- MtL ----------------#
#----------------------------------#



##### Reload data
squid <- read.table ("data squid original.csv",header=T, sep=",")
squidtemp<-squid[is.na(squid$FinalML)==F,]

ML<-squidtemp$FinalML[squidtemp$FinalML<3] #Only Reliable measurements

#hist(ML, probability=T,xlab = "Total Length (m)",main='')

n<-length(ML)

# ALRSM
pdf(file="Lmom_ML_153.pdf")
# LMRD
tau.4<-function(tau.3){tau.3*(1+5*tau.3)/(5+tau.3)}
curve(tau.4(x),-0.1,0.55,xlab = expression(tau[3]), 
      ylab = expression(tau[4]),cex.lab=1.5)
abline(h=0,v=0,col="grey")


x<-sort(ML)
n<-length(ML)

t4_line<- t3_line <- numeric()

tt3<-c();tt4<-c() 
var_t4dadott3<-c()
var_t3dadott4<-c()

#Candidate thresholds
#u<-x[c(1:(n-10))]
#u<- quantile(x,probs = seq(0.25,0.925,by=0.075), names=F)  #10 quantis
u<- quantile(x,probs = seq(0.25,0.953,by=0.037), names=F)   #20 quantis
#u<-x[1:153]

I<-1:length(u)

for (i in I) {
        x_exc <- x[x > u[i]] 
        Xu=x_exc-u[i]         
        nu<-length(Xu)
        
        #sample l-mom ratios for the current threshold
        tt4[i]=l4(Xu)/l2(Xu)
        tt3[i]=l3(Xu)/l2(Xu)
        
        
        #optimize min distance in 2dim to the theoretical line
        ff <- function(x){x*(1+5*x)/(5+x)}   
        #ff <- function(x){x^2+1}
        obj <- function(x, xi,yi) {
                sum((xi-x)^2+(yi-ff(x))^2)
        }       
        
        xmin <- optimize(f=obj, c(-1, 1), tol = 0.0000001, xi= tt3[i],yi= tt4[i])
        
        
        t3_line[i] <-xmin$minimum
        t4_line[i] <- ff(xmin$minimum)
        points(t3_line[i],t4_line[i],col="blue",pch=16)
        
        #### Estimated Root Variance with t4=tau.4(t3)
        TT1<- T_mat_1(Xu)
        #EVal_t4dadott3[i] <- tt4[i]  #this way the points would always fall inside
        var_t4dadott3[i] <- sqrt(TT1[2,2] /nu ) * sqrt(1-(TT1[1,2]^2/(TT1[1,1]*TT1[2,2])))
        
        #### Estimated Root Variance with t3= tau.3(t4)
        TT2<- T_mat_2(Xu)
        var_t3dadott4[i] <- sqrt(TT2[1,1] /nu ) * sqrt(1-(TT2[1,2]^2/(TT2[1,1]*TT2[2,2])))
        
        #points(tau_3,tau_4,col="red",pch=20)
        rm(x_exc)
        rm(Xu)
}

# length(t3)
# length(t4)
# cat(t3,"\n",t4)

#pontos originais
lines(tt3[I],tt4[I], type="p",col="navy") #, pch=20)

error_line <- (tt3-t3_line)^2+(tt4-t4_line)^2
index_line=1:length(error_line)
i0_line <- index_line[error_line[index_line]==min(error_line[which(is.na(error_line)==F)])]
#points(t3_line[i0_line],t4_line[i0_line],col="red",pch=20, cex=2)
points(tt3[i0_line],tt4[i0_line],col="red",pch=20, cex=1)
(u_opt_L<-u[i0_line]) # optimal threshold



#title("L-moment Ratio Diagram - 10 Candidate Thresholds")
#legend(0.25,0.5,expression(paste("(", t[3], "," ,t[4],") for u* = ", 10.44691 %~~% TL[57:74])),cex = 0.8, pt.cex = 1,
#       bty='n',pch=21, col=c("navy"), pt.bg=c("red"))
lines(tt3,tau.4(tt3)+qnorm(0.975)*var_t4dadott3, type="l", lty= 2, lwd= 2 ,col="orange") 
lines(tt3,tau.4(tt3)+qnorm(0.025)*var_t4dadott3, type="l", lty= 2, lwd= 2 ,col="orange")


#ESTIMATED Variance with t4 and t3=tau.3(t4)
#no mesmo plot t3 dado t4
lines(tau.3(tt4[I])+qnorm(0.975)*var_t3dadott4, tt4[I], type="l", lty= 2, lwd= 2 ,col="green") 
lines(tau.3(tt4[I])+qnorm(0.025)*var_t3dadott4, tt4[I], type="l", lty= 2, lwd= 2 ,col="green")

#10 candidates
legend("topleft",legend=c(expression(paste("95% Conf. Band for ",tau[4]," given ",t[3])),
                          expression(paste("95% Conf. Band for ",tau[3]," given ",t[4])), 
                          expression(paste("(", t[3], "," ,t[4],") for u* = ", 1.829))),
       title= "max(MtL) = 2.794 m; 20 candidates",title.adj = 0.15,
       col=c('orange','green','navy'), lty =c(2,2,0),pch = c(NA,NA,21), bty='o', pt.cex = c(1,1,1.2), 
       pt.bg=c(NA,NA,'red'),cex=1.2,bg="white",lwd=c(2,2,NA))


dev.off()


############################ GPD fit above u_10
############################
u_10<-u_opt_L
(fit_10<-evd::fpot(ML,u_10,model="gpd",npp=1,std.err = T, mper= 1000))
(fit_10_scale<-evd::fpot(ML,u_10,model="gpd",std.err = T))
(xF_10<-u_10-fit_10_scale$scale/fit_10_scale$param[2])
(sd_xF_10<-sqrt(fit_10_scale$var.cov[1,1]*(-1/fit_10_scale$param[2])^2 + 
                        fit_10_scale$var.cov[2,2]*(fit_10_scale$param[1]/fit_10_scale$param[2]^2)^2+
                        2* fit_10_scale$var.cov[1,2]*(-1/fit_10_scale$param[2])*(fit_10_scale$param[1]/fit_10_scale$param[2]^2)))

par(mfrow=c(2,2))
plot(fit_10)
par(mfrow=c(1,1))


############################ GPD fit above u_20
############################
u_20<-u_opt_L[1]
(fit_20<-evd::fpot(ML,u_20,model="gpd",npp=1,std.err = T, mper= 1000))
(fit_20_scale<-evd::fpot(ML,u_20,model="gpd",std.err = T))
(xF_20<-u_20-fit_20_scale$scale/fit_20_scale$param[2])
(sd_xF_20<-sqrt(fit_20_scale$var.cov[1,1]*(-1/fit_20_scale$param[2])^2 + 
                        fit_20_scale$var.cov[2,2]*(fit_20_scale$param[1]/fit_20_scale$param[2]^2)^2+
                        2*fit_20_scale$var.cov[1,2]*(-1/fit_20_scale$param[2])*(fit_20_scale$param[1]/fit_20_scale$param[2]^2)))

par(mfrow=c(2,2))
plot(fit_20)
par(mfrow=c(1,1))



####################################################################
############################################  Bader et al. (2018)
####################################################################
x<-sort(ML);n<-length(ML)
#u<- quantile(x,probs = seq(0.25,0.925,by=0.075), names=F)  #10 quantis
u<- quantile(x,probs = seq(0.25,0.953,by=0.037), names=F)   #20 quantis
bader_p_val <- tryCatch(eva::gpdSeqTests(x,u,method='ad', nsim=100), error= function (e) NA, warning = function(w) NA)
(FS <- tryCatch(eva::pSeqStop(bader_p_val$p.values), error= function (e) NA, warning = function(w) NA))


####################################################################
############################################  Northrop e Coleman
####################################################################
x<-sort(ML);n<-length(ML)
u1<- quantile(x,probs = seq(0.25,0.925,by=0.075), names=F)  #10 quantis
u2<- quantile(x,probs = seq(0.25,0.953,by=0.037), names=F)   #20 quantis
#u<-quantile(x, probs = seq(0.20,0.85,by=0.05), names=F)   #14 quantis
northrop <- score.fitrange(x,u)
abline(h=0.05,lty=3,col='grey')


pdf(file="Northrop_ML.pdf",width=8,height=4)
par(mfrow=c(1,2),mai=c(0.8,0.78,0.35,0.05))
northrop <- score.fitrange(x,u1)
abline(h=0.05,lty=3,col='grey')
text(1.15,0.1,"0.05", col="grey",cex=1.2)
abline(v=u_10,lty=2,col='red')
northrop <- score.fitrange(x,u2)
abline(h=0.05,lty=3,col='grey')
text(1.95,0.1,"0.05", col="grey",cex=1.2)
abline(v=u_20,lty=2,col='red')
dev.off()

############################ GPD fit above u_10
############################
u_10<-u1[8]
(fit_10<-evd::fpot(x,u_10,model="gpd",npp=1,std.err = T, mper= 1000))
(fit_10_scale<-evd::fpot(x,u_10,model="gpd",std.err = T))
(xF_10<-u_10-fit_10_scale$scale/fit_10_scale$param[2])
(sd_xF_10<-sqrt(fit_10_scale$var.cov[1,1]*(-1/fit_10_scale$param[2])^2 + 
                        fit_10_scale$var.cov[2,2]*(fit_10_scale$param[1]/fit_10_scale$param[2]^2)^2+
                        2* fit_10_scale$var.cov[1,2]*(-1/fit_10_scale$param[2])*(fit_10_scale$param[1]/fit_10_scale$param[2]^2)))


############################ GPD fit above u_20
############################
u_20<-u2[14]
(fit_20<-evd::fpot(x,u_20,model="gpd",npp=1,std.err = T, mper= 1000))
(fit_20_scale<-evd::fpot(x,u_20,model="gpd",std.err = T))
(xF_20<-u_20-fit_20_scale$scale/fit_20_scale$param[2])
(sd_xF_20<-sqrt(fit_20_scale$var.cov[1,1]*(-1/fit_20_scale$param[2])^2 + 
                        fit_20_scale$var.cov[2,2]*(fit_20_scale$param[1]/fit_20_scale$param[2]^2)^2+
                        2* fit_20_scale$var.cov[1,2]*(-1/fit_20_scale$param[2])*(fit_20_scale$param[1]/fit_20_scale$param[2]^2)))





###################################################################
######################## GIANT SQUID DATA #########################
############# Total Length and Mantle Length Analysis ############# 
#######################  From Paxton 2016  ######################## 
###################################################################
