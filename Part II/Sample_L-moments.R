#############################
#      Sample L-Moments     #         
#     cf. Hosking (1986)    #       
# Hosking and Wallis (1997) #
#############################


#estimator of PWM beta0 := M[1,0,0]
b0<-function(x){
  return(mean(x))
}
#estimator of PWM beta1 := M[1,1,0]
b1<-function(x){
  n=length(x)
  x<-sort(x)
  return(sum(((2:n)-1)*x[2:n])/(n*(n-1)))
}
#estimator of PWM beta2 := M[1,2,0]
b2<-function(x){
  n=length(x)
  x<-sort(x)
  return(sum(((3:n)-1)*((3:n)-2)*x[3:n])/(n*(n-1)*(n-2)))
}
#estimator of PWM beta3 := M[1,3,0]
b3<-function(x){
  n=length(x)
  x<-sort(x)
  return(sum(((4:n)-1)*((4:n)-2)*((4:n)-3)*x[4:n])/(n*(n-1)*(n-2)*(n-3)))
}



#estimator of PWM alpha0 := M[1,0,0]
a0<-function(x){
  return(b0(x))
}
#estimator of PWM alpha1 := M[1,0,1]
a1<-function(x){
  return(b0(x)-b1(x))
}
#estimator of PWM alpha2 := M[1,0,2]
a2<-function(x){
  return(b0(x)-2*b1(x)+b2(x))
}
#estimator of PWM alpha3 := M[1,0,3]
a3<-function(x){
  return(b0(x)-3*b1(x)+3*b2(x)-b3(x))
}




#estimator of 1st L-moment lambda1
l1<-function(x){
  return(a0(x))
}
#estimator of 2nd L-moment lambda2
l2<-function(x){
  return(a0(x)-2*a1(x))
}
#estimator of 3rd L-moment lambda3
l3<-function(x){
  return(a0(x)-6*a1(x)+6*a2(x))
}
#estimator of 4th L-moment lambda4
l4<-function(x){
  return(a0(x)-12*a1(x)+30*a2(x)-20*a3(x))
}



#estimator of L-skewness
t3<-function(x){
  return(l3(x)/l2(x))
}

#estimator of L-kurtosis
t4<-function(x){
  return(l4(x)/l2(x))
}


##########################
#    Sample L-Moments    #
##########################
