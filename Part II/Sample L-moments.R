##########################
#    Sample L-Moments    #
##########################

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


#estimator of 1st L-moment lambda1
l1<-function(x){
  return(b0(x))
}

#estimator of 2nd L-moment lambda2
l2<-function(x){
  return(2*b1(x)-b0(x))
}

#estimator of 3rd L-moment lambda3
l3<-function(x){
  return(6*b2(x)-6*b1(x)+b0(x))
}

#estimator of 4th L-moment lambda4
l4<-function(x){
  return(20*b3(x)-30*b2(x)+12*b1(x)-b0(x))
}

##########################
#    Sample L-Moments    #
##########################
