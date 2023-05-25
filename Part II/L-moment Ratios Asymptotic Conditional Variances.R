#############################################
#     L-statistics Asymptotic Variances     #
#            cf. Hosking (1986)             #
#         Hosking and Wallis (1997)         #
#############################################

# THEORETICAL L-moments and L-mom ratios for GP(sigma, xi)
# lambda1<-sigma/(1-xi)
# lambda2<-sigma/((1-xi)*(2-xi))
# tau3<-(1+xi)/(3-xi)
# tau4<-(1+xi)*(2+xi)/((3-xi)*(4-xi))

# THEORETICAL asymptotic co-variances of (a_r,a_s) - unbiased estimators of PWMs: (alpha_r,alpha_s)
# A_asym<-matrix(ncol=4,nrow = 4)
# for(r in 0:3){
#   for(s in 0:3){
#     A_asym[r+1,s+1]<-sigma^2/((r+1-xi)*(s+1-xi)*(r+s+1-2*xi))
#   }
# }

M<-matrix(c(1,1,1,1,0,-2,-6,-12,0,0,6,30,0,0,0,-20),4,4)
Mt<-t(M)

# THEORETICAL asymptotic co-variances of (l1,l2,l3,l4) - unbiased estimators of 4 l-moments
# Lambda_asym<-M%*%A_asym%*%Mt


# T_asym - Aymptotic var/cov of t3, t4  
#   with THEORETICAL parameters        

# T_asym<-matrix(ncol=2,nrow=2)
# T_asym[1,1]<-(Lambda_asym[3,3] - 2*tau3*Lambda_asym[2,3] + tau3^2*Lambda_asym[2,2])/(lambda2^2)
# T_asym[1,2]<-(Lambda_asym[3,4] - tau3*Lambda_asym[2,3] - tau4*Lambda_asym[2,4] + 
#                tau3*tau4*Lambda_asym[2,2])/(lambda2^2)
# T_asym[2,1]<-T_asym[1,2]
# T_asym[2,2]<-(Lambda_asym[4,4] - 2*tau4*Lambda_asym[2,4] + tau4^2*Lambda_asym[2,2])/(lambda2^2)


# Curve of t4 as function of t3 (for GP model)
# tau4=g(tau3)
tau.4<-function(tau.3){tau.3*(1+5*tau.3)/(5+tau.3)}

# Curve of t3 as function of t4 (for GP model with t4>=0)
tau.3<-function(tau.4){(tau.4-1)/10 + (1/10)*sqrt(tau.4^2+98*tau.4+1)}


#  GPd Parameter Estimates based on PWM  

sigma_pwm<-function(x){
  return(2*a0(x)*a1(x)/l2(x))
}

xi_pwm<-function(x){
  return(2 - a0(x)/l2(x))
}


# A_mat(x) - matrix cov(ar,a_s) with estimated parameters    
#    Estimator of A_asym        

A_mat<- function(x){
  A<-matrix(ncol=4,nrow = 4)
  for(r in 0:3){
    for(s in 0:3){
      A[r+1,s+1]<-sigma_pwm(x)^2/((r+1-xi_pwm(x))*(s+1-xi_pwm(x))*(r+s+1-2*xi_pwm(x)))
    }
  }
  return(A)
}


# Lambda(x) - matrix cov(l1,l2,l3,l4) with estimated parameters         
#    Estimator of Lambda_asym        

Lambda<-function(x){
  return(M%*%A_mat(x)%*%Mt)
}


# T_mat(x) - Aymptotic var/cov of (t3,t4) with estimated parameters            
#    Estimator of T_asym            

# t4 estimated as known function of t3 for the GP model
# t4= tau.4(t3)

T_mat_1<-function(x){
  TT<-matrix(ncol=2,nrow=2)
  TT[1,1]<-(Lambda(x)[3,3] - 2*t3(x)*Lambda(x)[2,3] + t3(x)^2*Lambda(x)[2,2])/(l2(x)^2)
  TT[1,2]<-(Lambda(x)[3,4] - t3(x)*Lambda(x)[2,3] - tau.4(t3(x))*Lambda(x)[2,4] +
              t3(x)*tau.4(t3(x))*Lambda(x)[2,2])/(l2(x)^2)
  TT[2,1]<-TT[1,2]
  TT[2,2]<-(Lambda(x)[4,4] - 2*tau.4(t3(x))*Lambda(x)[2,4] + tau.4(t3(x))^2*Lambda(x)[2,2])/(l2(x)^2)
  return(TT)
}


# t3 estimated as known function of t4 for the GP model
# t3= tau.3(t4)

T_mat_2<-function(x){
  TT<-matrix(ncol=2,nrow=2)
  TT[1,1]<-(Lambda(x)[3,3] - 2*tau.3(t4(x))*Lambda(x)[2,3] + tau.3(t4(x))^2*Lambda(x)[2,2])/(l2(x)^2)
  TT[1,2]<-(Lambda(x)[3,4] - tau.3(t4(x))*Lambda(x)[2,3] - t4(x)*Lambda(x)[2,4] +
              t4(x)*tau.3(t4(x))*Lambda(x)[2,2])/(l2(x)^2)
  TT[2,1]<-TT[1,2]
  TT[2,2]<-(Lambda(x)[4,4] - 2*t4(x)*Lambda(x)[2,4] + t4(x)^2*Lambda(x)[2,2])/(l2(x)^2)
  return(TT)
}





# Conditional expected values (theoretical)
EVal_t4giventt3<-function(tau3,tau4,tt3){
  return(tau4 + (T_asym[1,2]/T_asym[1,1])*(tt3-tau3))
}

EVal_t3giventt4<-function(tau3,tau4,tt4){
  return(tau3 + (T_asym[1,2]/T_asym[2,2])*(tt4-tau4))
}


# Conditional asymptotic variances (theoretical)

var_t4giventt3_teorica<-function(nu){
  return(sqrt(T_asym[2,2] /nu ) * sqrt(1-(T_asym[1,2]^2/(T_asym[1,1]*T_asym[2,2]))))
}


var_t3giventt4_teorica<-function(nu){
  return(sqrt(T_asym[1,1] /nu ) * sqrt(1-(T_asym[1,2]^2/(T_asym[1,1]*T_asym[2,2]))))
} 

# Conditional asymptotic variances (estimated)

var_t4giventt3<- function(Xu,nu){
  TT<- T_mat_1(Xu)
  return(sqrt(TT[2,2] /nu ) * sqrt(1-(TT[1,2]^2/(TT[1,1]*TT[2,2]))))
} 


var_t3giventt4<- function(Xu,nu){ 
  TT<- T_mat_2(Xu)
  return(sqrt(TT[1,1] /nu ) * sqrt(1-(TT[1,2]^2/(TT[1,1]*TT[2,2]))))
}



#############################################
#     L-statistics Asymptotic Variances     #
#            cf. Hosking (1986)             #
#         Hosking and Wallis (1997)         #
#############################################
