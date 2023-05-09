###########################################
#           Hybrid Distribution           #
#          Northrop et al. (2017)         #
###########################################


### from sim_2.R of Northrop et al. (2017)
qhybrid <- function(x,sigma1=1,xi1=-1,xi2=0.2,qq=0.75){
  uu <- sigma1*((1-qq)^(-xi1)-1)/xi1  #quantile 0.75 de GP(-1,1)
  sigma2 <- sigma1+xi1*uu             
  pp <- (1+xi1*uu/sigma1)^(-1/xi1)    #exceedance prob of uu for GP(-1,1)
                                      #pp = 1 - qq = 0.25
  xd <- (x-1+pp)/pp                   #xd = (x - 0.75)/0.25
  ifelse(x<qq,sigma1*((1-x)^(-xi1)-1)/xi1,uu+sigma2*((1-xd)^(-xi2)-1)/xi2)
}

### random generator from the Uniform-GP hybrid distribution
### based in Northrop et al. (2017)
### U(0,1) up to its qq% quantile = qq
### GP(xi2, mu=0, sigma=1-qq ) from qq forward  
### . n - sample size
### . qq - probability of the threshold quantile
### . xi2 - shape of the GP distribution
### . default values - qq=0.75, xi=0.1 
### . GP with set location and scale

rhybrid <- function(n, xi2=0.1,qq=0.75){
  aux<-runif(n)
  qhybrid(aux,xi2=xi2,qq=qq)
  
}

# ###### test trial ########
#n=1000
#x<-rhybrid(n,xi2=0.2,qq=0.75)
#hist(x,probability=T, breaks=30, xlim=c(0,3))
#lines(density(x), col='red')
# ###### empirical CDF check #########
# yy<-seq(1/n,1,by=1/(n))
# plot(sort(x),yy)
# 
# h <- function (x) {
#   hx<-c()
#   for(i in 1:length(x)){
#   if(x[i]>=0 & x[i]<=0.75) {
#     hx[i]=1
#   }
#   else {
#     hx[i]=(0.4+0.8*x[i])^(-6)
#   }}
#   return(hx)
# }
# 
# 
# #curve(h(x),0,3,add=T, col='red')
# mu_gev<-0.75-(0.25/0.2)*((-log(0.75))^(-0.2)-1)
# #curve(dgev(x,xi=0.2,mu=mu_gev,sigma=0.25),add=T, col='blue')
# pgev(0.75,xi=0.2,mu=mu_gev,sigma=0.25)
# 
# 
# 
# xx=seq(0,4,length.out = 10000)
# hxx <- h(xx)
# hxxgev<-dgev(xx,xi=0.2,mu=mu_gev,sigma=0.25)
# hxxgev2<-dgev(xx,xi=0.2,mu=0.75,sigma=0.25)
# hxxgev3<-dgev(xx,xi=0.2,mu=0,sigma=1)

#plot(xx,hxx,type="l",xlim=c(0,4),col="red",xlab=bquote(x),ylab="Hybrid(x|0.75,0.2)",cex.lab=1.2 )
#plot(xx,hxx,type="l",xlim=c(u0,2.5),ylim=c(0,5.5),col="red",xlab=bquote(x),ylab="Hybrid(x|0.75,0.2)",cex.lab=1.2 )
#abline(h=0, col="grey",lty=3)
#segments(0.75,1,
#         0.75,0, lty=2)

#lines(xx,hxxgev,col='blue')
#lines(xx,hxxgev2,col='green')
#lines(xx,hxxgev3,col='orange')



###########################################
#           Hybrid Distribution           #
#          Northrop et al. (2017)         #
###########################################