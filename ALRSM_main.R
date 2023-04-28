####################################################
#    Automatic L-moment Ration Selection Method    #
####################################################

#import L-moments estimators
source('ALRSM_Sample L-Moments.R')

#'data' - your array of iid observations
x<-sort(data)
n<-length(x)

#GPd theoretical curve & plot
tau.4<-function(tau.3){tau.3*(1+5*tau.3)/(5+tau.3)}
curve(tau.4(x),0,0.6,xlab = expression(tau[3]), #might need to adjust the window
      ylab = expression(tau[4]))
abline(h=0,v=0,col="grey")


#space for sample l-mom ratios
t3 <- t4 <-c()
#space for corresponding l-mom ratios on GPd curve
t4_line<- t3_line <- numeric()

#Candidate thresholds
#u<- quantile(x,probs = seq(0.25,0.925,by=0.075), names=F)  #10 quantiles
u<- quantile(x,probs = seq(0.25,0.953,by=0.037), names=F)  #20 quantiles

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


####################################################
#    Automatic L-moment Ration Selection Method    #
####################################################