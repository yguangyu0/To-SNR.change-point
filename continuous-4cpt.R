############# Try 4-Knot Our Proposed Method Simulation Result ###################
set.seed(20190101)

##########################################################################
#############   Functions  ########################
NR_newfunction<-function(X, Y, tau00)
{
  erro = FALSE
  eta = rep(1,4); j = 1
  tau0 = tau00
  while (any(abs(eta) > 1e-06) & j<=1000 & erro==FALSE)
  {
    tau0_1 = tau0[1];tau0_2 = tau0[2];tau0_3 = tau0[3];tau0_4 = tau0[4]
    ### generate inital values ###
    X_tau_1 =  as.numeric(X > tau0_1) * ( X - tau0_1 )
    X_tau_2 =  as.numeric(X > tau0_2) * ( X - tau0_2 )
    X_tau_3 =  as.numeric(X > tau0_3) * ( X - tau0_3 )
    X_tau_4 =  as.numeric(X > tau0_4) * ( X - tau0_4 )    
    lin.inital <- lm(Y~ X+ X_tau_1+X_tau_2+ X_tau_3+X_tau_4)
    
    theta0 = c(as.numeric(coef(lin.inital)),tau0_1,tau0_2,tau0_3,tau0_4)
    y_pred0 = predict(lin.inital)
    
    I_tau_1 = as.numeric(X > tau0_1);I_tau_2 = as.numeric(X > tau0_2);
    I_tau_3 = as.numeric(X > tau0_3);I_tau_4 = as.numeric(X > tau0_4);
    sig_error = Y-y_pred0
    ### Calculate U and J ###
    Un=c(theta0[3]*mean(sig_error*I_tau_1),
         theta0[4]*mean(sig_error*I_tau_2),
         theta0[5]*mean(sig_error*I_tau_3),
         theta0[6]*mean(sig_error*I_tau_4))
    Un = matrix(Un,nrow=4)
    
    Jn =-c(theta0[3]^2*mean(I_tau_1),                  theta0[3]*theta0[4]*mean(I_tau_1*I_tau_2),  theta0[3]*theta0[5]*mean(I_tau_1*I_tau_3),   theta0[3]*theta0[6]*mean(I_tau_1*I_tau_4),
           theta0[4]*theta0[3]*mean(I_tau_2*I_tau_1),  theta0[4]^2*mean(I_tau_2)  ,                theta0[4]*theta0[5]*mean(I_tau_2*I_tau_3),   theta0[4]*theta0[6]*mean(I_tau_2*I_tau_4),                                                      
           theta0[5]*theta0[3]*mean(I_tau_3*I_tau_1),  theta0[5]*theta0[4]*mean(I_tau_3*I_tau_2),  theta0[5]^2*mean(I_tau_3),                   theta0[5]*theta0[6]*mean(I_tau_3*I_tau_4),
           theta0[6]*theta0[3]*mean(I_tau_4*I_tau_1),  theta0[6]*theta0[4]*mean(I_tau_4*I_tau_2),  theta0[6]*theta0[5]*mean(I_tau_4*I_tau_3),   theta0[6]^2*mean(I_tau_4)
    )
    
    Jn = matrix(Jn,ncol=4,nrow=4)
    
    tryCatch( {eta = solve(Jn)%*%Un},error=function(cond) erro = TRUE)
    
    if (erro == FALSE)
    {
      ### Update Tau ###
      tau_up = tau0 + eta
      tau0 = tau_up
      
      #tau1_c = c(tau1_c, tau0[1])
      #tau2_c = c(tau2_c,tau0[2])
      #tau3_c = c(tau3_c, tau0[3])
      #tau4_c = c(tau4_c,tau0[4])
    }
    j=j+1}
  
  tau.out = c(NA,NA,NA,NA)
  if (erro == FALSE & max(tau0,na.rm=T)<max(X) & min(tau0,na.rm=T)>min(X))
  { tau.out = tau0}
  
  return(tau.out)
}

sd_function<- function(tau.out)
{
  tau0_1 =  tau.out[1];tau0_2 =  tau.out[2];
  tau0_3 =  tau.out[3];tau0_4 =  tau.out[4];
  # tau0_1 =  0.1;tau0_2 = 0.3;
  # tau0_3 =  0.5;tau0_4 = 0.7;
  erro = FALSE
  X_tau_1 =  as.numeric(X > tau0_1) * ( X - tau0_1 )
  X_tau_2 =  as.numeric(X > tau0_2) * ( X - tau0_2 )
  X_tau_3 =  as.numeric(X > tau0_3) * ( X - tau0_3 )
  X_tau_4 =  as.numeric(X > tau0_4) * ( X - tau0_4 )    
  lin.inital <- lm(Y~ X+ X_tau_1+X_tau_2+ X_tau_3+X_tau_4)
  sig.hat = sigma(lin.inital)
  theta0 = c(as.numeric(coef(lin.inital)),tau0_1,tau0_2,tau0_3,tau0_4)
  y_pred0 = predict(lin.inital)
  sig_error = Y-y_pred0
  
  I_tau_1 = as.numeric(X > tau0_1);I_tau_2 = as.numeric(X > tau0_2);
  I_tau_3 = as.numeric(X > tau0_3);I_tau_4 = as.numeric(X > tau0_4);
  
  Jn1 =2*c(1,             mean(X),         mean(X_tau_1),          mean(X_tau_2),         mean(X_tau_3),         mean(X_tau_4),
           mean(X),       mean(X^2),       mean(X*X_tau_1),        mean(X*X_tau_2),       mean(X*X_tau_3),       mean(X*X_tau_4),
           mean(X_tau_1), mean(X*X_tau_1), mean(X_tau_1^2),        mean(X_tau_1*X_tau_2), mean(X_tau_1*X_tau_3), mean(X_tau_1*X_tau_4),
           mean(X_tau_2), mean(X*X_tau_2), mean(X_tau_1*X_tau_2),  mean(X_tau_2^2),       mean(X_tau_2*X_tau_3), mean(X_tau_2*X_tau_4),
           mean(X_tau_3), mean(X*X_tau_3), mean(X_tau_1*X_tau_3),  mean(X_tau_2*X_tau_3), mean(X_tau_3^2),       mean(X_tau_3*X_tau_4),
           mean(X_tau_4), mean(X*X_tau_4), mean(X_tau_1*X_tau_4),  mean(X_tau_2*X_tau_4), mean(X_tau_3*X_tau_4), mean(X_tau_4^2)
  )
  Jn1 = matrix( Jn1 , ncol=6,nrow=6)
  
  Jn2=-2*c(theta0[3]*mean(I_tau_1),                             theta0[4]*mean(I_tau_2),                          theta0[5]*mean(I_tau_3),                           theta0[6]*mean(I_tau_4),
           theta0[3]*mean(X*I_tau_1),                           theta0[4]*mean(X*I_tau_2),                        theta0[5]*mean(X*I_tau_3),                         theta0[6]*mean(X*I_tau_4),
           theta0[3]*mean(X_tau_1)-mean(I_tau_1*sig_error),     theta0[4]*mean(X_tau_1*I_tau_2),                  theta0[5]*mean(X_tau_1*I_tau_3),                   theta0[6]*mean(X_tau_1*I_tau_4),
           theta0[3]*mean(X_tau_2*I_tau_1),                     theta0[4]*mean(X_tau_2)-mean(I_tau_2*sig_error),  theta0[5]*mean(X_tau_2*I_tau_3),                   theta0[6]*mean(X_tau_2*I_tau_4),
           theta0[3]*mean(X_tau_3*I_tau_1),                     theta0[4]*mean(X_tau_3*I_tau_2),                  theta0[5]*mean(X_tau_3)-mean(I_tau_3*sig_error),   theta0[6]*mean(X_tau_3*I_tau_4),
           theta0[3]*mean(X_tau_4*I_tau_1),                     theta0[4]*mean(X_tau_4*I_tau_2),                  theta0[5]*mean(X_tau_4*I_tau_3),                   theta0[6]*mean(X_tau_4)-mean(I_tau_4*sig_error)
  )
  Jn2 = matrix( Jn2 , ncol=4,nrow=6,byrow = T)
  
  Jn3 = t(Jn2)
  
  Jn4 =2*c(theta0[3]^2*mean(I_tau_1),                  theta0[3]*theta0[4]*mean(I_tau_1*I_tau_2),  theta0[3]*theta0[5]*mean(I_tau_1*I_tau_3),   theta0[3]*theta0[6]*mean(I_tau_1*I_tau_4),
           theta0[4]*theta0[3]*mean(I_tau_2*I_tau_1),  theta0[4]^2*mean(I_tau_2)  ,                theta0[4]*theta0[5]*mean(I_tau_2*I_tau_3),   theta0[4]*theta0[6]*mean(I_tau_2*I_tau_4),                                                      
           theta0[5]*theta0[3]*mean(I_tau_3*I_tau_1),  theta0[5]*theta0[4]*mean(I_tau_3*I_tau_2),  theta0[5]^2*mean(I_tau_3),                   theta0[5]*theta0[6]*mean(I_tau_3*I_tau_4),
           theta0[6]*theta0[3]*mean(I_tau_4*I_tau_1),  theta0[6]*theta0[4]*mean(I_tau_4*I_tau_2),  theta0[6]*theta0[5]*mean(I_tau_4*I_tau_3),   theta0[6]^2*mean(I_tau_4)
  )
  
  Jn4 = matrix(Jn4,ncol=4,nrow=4)
  J = solve(Jn4 - Jn3%*%solve(Jn1)%*%Jn2)
  tryCatch( {V = 2*sig.hat^2*J/n},error=function(cond) erro = TRUE)
  sd.out = c(NA,NA,NA,NA)
  if (erro == FALSE )
  {
    sd_1 = sqrt(V[1,1])
    sd_2 = sqrt(V[2,2])
    sd_3 = sqrt(V[3,3])
    sd_4 = sqrt(V[4,4])
  }
  sd.out = c(sd_1, sd_2, sd_3, sd_4)
  return(sd.out)
}


theor_sd_function<- function()
{
  tau0_1 =  theta[7];tau0_2 = theta[8];
  tau0_3 =  theta[9];tau0_4 =  theta[10];
  
  EX_tau_1 =  1/2*(1-tau0_1^2)-tau0_1*(1-tau0_1)
  EX_tau_2 =  1/2*(1-tau0_2^2)-tau0_2*(1-tau0_2)
  EX_tau_3 =  1/2*(1-tau0_3^2)-tau0_3*(1-tau0_3)
  EX_tau_4 =  1/2*(1-tau0_4^2)-tau0_4*(1-tau0_4)  
  
  theta0 = theta
  
  I_tau_1 = as.numeric(X > tau0_1);I_tau_2 = as.numeric(X > tau0_2);
  I_tau_3 = as.numeric(X > tau0_3);I_tau_4 = as.numeric(X > tau0_4);
  
  E_r1r2<-function(r1,r2)
  {
    r = max(r1,r2)
    out = 1/3-(r1+r2)/2+r1*r2-(1/3*r^3-1/2*(r1+r2)*r^2+r1*r2*r)
    return(out)
  }
  Jn10=2*c(1,       1/2,                     EX_tau_1,                        EX_tau_2,                         EX_tau_3,                        EX_tau_4,
           1/2,      1/3,                     1/3-tau0_1/2+tau0_1^3/6,         1/3-tau0_2/2+tau0_2^3/6,          1/3-tau0_3/2+tau0_3^3/6,         1/3-tau0_4/2+tau0_4^3/6,
           EX_tau_1, 1/3-tau0_1/2+tau0_1^3/6, 1/3-tau0_1+tau0_1^2-1/3*tau0_1^3,E_r1r2(tau0_1,tau0_2),            E_r1r2(tau0_1,tau0_3),           E_r1r2(tau0_1,tau0_4),
           EX_tau_2, 1/3-tau0_2/2+tau0_2^3/6, E_r1r2(tau0_1,tau0_2),           1/3-tau0_2+tau0_2^2-1/3*tau0_2^3, E_r1r2(tau0_2,tau0_3),           E_r1r2(tau0_2,tau0_4),
           EX_tau_3, 1/3-tau0_3/2+tau0_3^3/6, E_r1r2(tau0_1,tau0_3),           E_r1r2(tau0_2,tau0_3),            1/3-tau0_3+tau0_3^2-1/3*tau0_3^3,E_r1r2(tau0_3,tau0_4),
           EX_tau_4, 1/3-tau0_4/2+tau0_4^3/6, E_r1r2(tau0_1,tau0_4),           E_r1r2(tau0_2,tau0_4),            E_r1r2(tau0_3,tau0_4),           1/3-tau0_4+tau0_4^2-1/3*tau0_4^3
  )
  Jn10 = matrix( Jn10 , ncol=6,nrow=6)
  
  E_Xr1_Ir2<-function(r1,r2)
  {
    r = max(r1,r2)
    out = 1/2-r1-(r^2/2-r1*r)
    return(out)
  }
  
  Jn20=-2*c(theta0[3]*(1-tau0_1),                theta0[4]*(1-tau0_2),               theta0[5]*(1-tau0_3),               theta0[6]*(1-tau0_4),
            theta0[3]*(1-tau0_1^2)/2,            theta0[4]*(1-tau0_2^2)/2,           theta0[5]*(1-tau0_3^2)/2,           theta0[6]*(1-tau0_4^2)/2,
            theta0[3]*EX_tau_1,                  theta0[4]*E_Xr1_Ir2(tau0_1,tau0_2), theta0[5]*E_Xr1_Ir2(tau0_1,tau0_3), theta0[6]*E_Xr1_Ir2(tau0_1,tau0_4),
            theta0[3]*E_Xr1_Ir2(tau0_2,tau0_1),  theta0[4]*EX_tau_2,                 theta0[5]*E_Xr1_Ir2(tau0_2,tau0_3), theta0[6]*E_Xr1_Ir2(tau0_2,tau0_4),
            theta0[3]*E_Xr1_Ir2(tau0_3,tau0_1),  theta0[4]*E_Xr1_Ir2(tau0_3,tau0_2), theta0[5]*EX_tau_3,                 theta0[6]*E_Xr1_Ir2(tau0_3,tau0_4),
            theta0[3]*E_Xr1_Ir2(tau0_4,tau0_1),  theta0[4]*E_Xr1_Ir2(tau0_4,tau0_2), theta0[5]*E_Xr1_Ir2(tau0_4,tau0_3), theta0[6]*EX_tau_4
  )
  Jn20 = matrix( Jn20 , ncol=4,nrow=6,byrow = T)
  
  Jn30= t(Jn20)
  
  E_Ir1_Ir2<-function(r1,r2)
  {
    r = max(r1,r2)
    out = 1-r
    return(out)
  }
  
  Jn40 = 2*c (theta0[3]^2*(1-tau0_1),                        theta0[3]*theta0[4]*E_Ir1_Ir2(tau0_1,tau0_2),  theta0[3]*theta0[5]*E_Ir1_Ir2(tau0_1,tau0_3),   theta0[3]*theta0[6]*E_Ir1_Ir2(tau0_1,tau0_4),
              theta0[4]*theta0[3]*E_Ir1_Ir2(tau0_1,tau0_2),  theta0[4]^2*(1-tau0_2)  ,                      theta0[4]*theta0[5]*E_Ir1_Ir2(tau0_2,tau0_3),   theta0[4]*theta0[6]*E_Ir1_Ir2(tau0_2,tau0_4),                                                      
              theta0[5]*theta0[3]*E_Ir1_Ir2(tau0_1,tau0_3),  theta0[5]*theta0[4]*E_Ir1_Ir2(tau0_2,tau0_3),  theta0[5]^2*(1-tau0_3),                         theta0[5]*theta0[6]*E_Ir1_Ir2(tau0_3,tau0_4),
              theta0[6]*theta0[3]*E_Ir1_Ir2(tau0_1,tau0_4),  theta0[6]*theta0[4]*E_Ir1_Ir2(tau0_2,tau0_4),  theta0[6]*theta0[5]*E_Ir1_Ir2(tau0_3,tau0_4),   theta0[6]^2*(1-tau0_4)
  )
  
  Jn40 = matrix(Jn40,ncol=4,nrow=4)
  J = solve(Jn40 - Jn30%*%solve(Jn10)%*%Jn20)
  V = 2*sig^2*J/n
  
  sd_1 = sqrt(V[1,1])
  sd_2 = sqrt(V[2,2])
  sd_3 = sqrt(V[3,3])
  sd_4 = sqrt(V[4,4])
  
  sd.out = c(sd_1, sd_2, sd_3, sd_4)
  return(sd.out)
}
##########################################################################
##########################################################################
##########################################################################
##########################################################################




#################################### Settings ############################
######################## true values beta0, betas,tau ####
sig = 0.05
##################### Setup 4.1 ########################################
theta = c(0.3, 1, -2, 4,-5, 3 , 0.2, 0.4, 0.6, 0.8)
tau00=c(0.1,0.3,0.5,0.9)
tau0_1=0.1; tau0_2 = 0.3; tau0_3=0.5; tau0_4 = 0.9;


# ##################### Setup 4.2 ########################################
# theta = c(5, -1, 3, 5,-9,6 , 0.3, 0.6, 0.8, 0.9)
# tau0_1=0.2; tau0_2 = 0.5; tau0_3=0.85; tau0_4 = 0.95;

# ###################### Setup 4.3 ########################################
# theta = c(9, -5, 6, 7,-8,-3 , 0.1, 0.3, 0.5, 0.7)
# tau0_1=0.05; tau0_2 = 0.4; tau0_3=0.55; tau0_4 = 0.8;


Nsim = 1000;  n = 5000
i = 1 ; o = 0
sd_est1 = NULL; sd_est2 = NULL; sd_est3 = NULL; sd_est4 = NULL;
tau_est1 = NULL;tau_est2 = NULL;tau_est3 = NULL;tau_est4 = NULL;
cov1 = 0; cov2 = 0; cov3 = 0; cov4 = 0;
####### Simulation Procedure ###########

while (i <= Nsim)
{
  ### generate dataset ###
  X = runif (n,0,1)
  e = rnorm (n,0,sd=sig)
  Y = theta[1] + theta[2]*X + theta[3] * as.numeric(X > theta[7]) * ( X - theta[7] ) + theta[4] * as.numeric(X > theta[8]) * ( X - theta[8] ) + theta[5] * as.numeric(X > theta[9]) * ( X - theta[9] )  + theta[6] * as.numeric(X > theta[10]) * ( X - theta[10] )+e
  
  #ptm <- proc.time()
  tau.out = NR_newfunction(X,Y,tau00)#tau0_1,tau0_2,tau0_3,tau0_4)
  #proc.time() - ptm
  
  ### Result ###
  if (any(!is.na(tau.out)))
  {  tau.out = sort(tau.out)
  tau_est1[i] = tau.out[1]
  tau_est2[i] = tau.out[2]
  tau_est3[i] = tau.out[3]
  tau_est4[i] = tau.out[4]
  
  sd.out=sd_function(tau.out)
  sd_est1[i] = sd.out[1]
  sd_est2[i] = sd.out[2]
  sd_est3[i] = sd.out[3]
  sd_est4[i] = sd.out[4]
  
  cov1 = cov1 + as.numeric(tau.out[1]-1.96*sd.out[1] <= theta[7] & tau.out[1]+1.96*sd.out[1]>=theta[7])
  cov2 = cov2 + as.numeric(tau.out[2]-1.96*sd.out[2] <= theta[8] & tau.out[2]+1.96*sd.out[2]>=theta[8])
  cov3 = cov3 + as.numeric(tau.out[3]-1.96*sd.out[3] <= theta[9] & tau.out[3]+1.96*sd.out[3]>=theta[9])
  cov4 = cov4 + as.numeric(tau.out[4]-1.96*sd.out[4] <= theta[10] & tau.out[4]+1.96*sd.out[4]>=theta[10])
  
  
  
  i=i+1
  if (i%% 100 == 0) {print (i)}
  }
  
  o=o+1
}

#### Coverage  Probability ####
c(cov1,cov2,cov3,cov4)/Nsim*100
##### Mean & SD (MTCL) ####
bias_tau1 =  mean(tau_est1) - theta[7] 
sd_tau1 = sd(tau_est1)
mean_sd1 = mean(sd_est1)
round(c(bias_tau1,sd_tau1,mean_sd1 )*1000,2)


bias_tau2 = mean(tau_est2)- theta[8] 
sd_tau2 = sd(tau_est2)
mean_sd2 = mean(sd_est2)
round(c(bias_tau2,sd_tau2,mean_sd2 )*1000,2)

bias_tau3 = mean(tau_est3) -theta[9] 
sd_tau3 = sd(tau_est3)
mean_sd3 = mean(sd_est3)
round(c(bias_tau3,sd_tau3,mean_sd3 )*1000,2)

bias_tau4 = mean(tau_est4) - theta[10] 
sd_tau4 = sd(tau_est4)
mean_sd4 = mean(sd_est4)
round(c(bias_tau4,sd_tau4,mean_sd4 )*1000,2)

###Theoretical S.D. #####

round(theor_sd_function()*1000,2)

### The probability of having the result ###
Nsim/o

par(mfrow=c(2,2))
### Density Plot ####
den<-density (tau_est1)
plot(den,col='blue',lwd=2,main = 'Density Plot for Proposed Method Tau1',
     xlab='X')
abline ( v= theta[5], lty=2)
abline (v=mean(tau_est1), lty=2, col='red')

den<-density (tau_est2)
plot(den,col='blue',lwd=2,main = 'Density Plot for Proposed Method Tau2',
     xlab='X')
abline ( v= theta[6], lty=2)
abline (v=mean(tau_est2), lty=2, col='red')

den<-density (tau_est3)
plot(den,col='blue',lwd=2,main = 'Density Plot for Proposed Method Tau3',
     xlab='X')
abline ( v= theta[5], lty=2)
abline (v=mean(tau_est3), lty=2, col='red')

den<-density (tau_est4)
plot(den,col='blue',lwd=2,main = 'Density Plot for Proposed Method Tau4',
     xlab='X')
abline ( v= theta[6], lty=2)
abline (v=mean(tau_est4), lty=2, col='red')

##### Sample Plot ###########
par(mfrow=c(1,1))
plot (X,Y,col='blue',main='Example Question')
curve(theta[1] + theta[2]*x + theta[3] * as.numeric(x > theta[7]) * ( x - theta[7] ) + theta[4] * as.numeric(x > theta[8]) * ( x - theta[8] ) + theta[5] * as.numeric(x > theta[9]) * ( x - theta[9] )  + theta[6] * as.numeric(x > theta[10]) * ( x - theta[10] ),
      col='blue', add = T,lwd=2,lty=2)
abline(v=c(0.2,0.4,0.6,0.8),lty=2)

curve(theta0[1] + theta0[2]*x + theta0[3] * as.numeric(x > theta0[7]) * ( x - theta0[7] )+ theta0[4] * as.numeric(x > theta0[8]) * ( x - theta0[8] ) + theta0[5] * as.numeric(x > theta0[9]) * ( x - theta0[9] ) + theta0[6] * as.numeric(x > theta0[10]) * ( x - theta0[10] ),
      col='blue', add = T,lwd=1,lty=2)

##### Sample Plot ###########
par(mfrow=c(1,3))
plot (X,Y,col='orange',main="Setup 4.3")
curve(theta[1] + theta[2]*x + theta[3] * as.numeric(x > theta[7]) * ( x - theta[7] ) + theta[4] * as.numeric(x > theta[8]) * ( x - theta[8] ) + theta[5] * as.numeric(x > theta[9]) * ( x - theta[9] )  + theta[6] * as.numeric(x > theta[10]) * ( x - theta[10] ),
      add = T,lwd=2,lty=2)

abline(v=c(theta[7:10]),lty=2,col='red')


