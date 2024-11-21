###################################################################################
########### To-SNR Algorithm: Continuous Outcome With Four Change Points ##########
######            Each Covariate Has Two Change Points                       ###### 
###################################################################################

#set.seed(20231113)

############# Functions For To-SNR Algorithm ####################
ToSNR_continuous<-function(tau0)
{
  erro = FALSE
  eta = rep(1,4); j = 1
  
  while (any(abs(eta) > 1e-7) & j<=100 & erro==FALSE)
  {
    tau0_1 = tau0[1];tau0_2 = tau0[2];
    tau0_3 = tau0[3];tau0_4 = tau0[4]
    X1_tau1 =  as.numeric(X1 > tau0_1) * ( X1 - tau0_1 )
    X1_tau2 =  as.numeric(X1 > tau0_2) * ( X1 - tau0_2 )
    X2_tau1 =  as.numeric(X2 > tau0_3) * ( X2 - tau0_3 )
    X2_tau2 =  as.numeric(X2 > tau0_4) * ( X2 - tau0_4 )
    
    I1_tau1 = as.numeric(X1 > tau0_1)
    I1_tau2 = as.numeric(X1 > tau0_2)
    I2_tau1 = as.numeric(X2 > tau0_3)
    I2_tau2 = as.numeric(X2 > tau0_4)
    
    lm.inital <- lm(Y~X1+ X1_tau1+X1_tau2+X2+X2_tau1+X2_tau2+Z)
    theta0 =coef( lm.inital) 
    be.tau = c(theta0['X1_tau1'],theta0['X1_tau2'],theta0['X2_tau1'],theta0['X2_tau2'])
    mu0 = predict(lm.inital)
    ### Calculate U and J ###
    Un=c(be.tau[1]*mean((mu0-Y)*I1_tau1),be.tau[2]*mean((mu0-Y)*I1_tau2),be.tau[3]*mean((mu0-Y)*I2_tau1),be.tau[4]*mean((mu0-Y)*I2_tau2))
    Un = matrix(Un,nrow=4)
    
    Jn = c(be.tau[1]^2*mean(I1_tau1),                  be.tau[1]*be.tau[2]*mean(I1_tau1*I1_tau2),  be.tau[1]*be.tau[3]*mean(I1_tau1*I2_tau1),   be.tau[1]*be.tau[4]*mean(I1_tau1*I2_tau2),
           be.tau[2]*be.tau[1]*mean(I1_tau2*I1_tau1),  be.tau[2]^2*mean(I1_tau2)  ,                be.tau[2]*be.tau[3]*mean(I1_tau2*I2_tau1),   be.tau[2]*be.tau[4]*mean(I1_tau2*I2_tau2),                                                      
           be.tau[3]*be.tau[1]*mean(I2_tau1*I1_tau1),  be.tau[3]*be.tau[2]*mean(I2_tau1*I1_tau2),  be.tau[3]^2*mean(I2_tau1),                   be.tau[3]*be.tau[4]*mean(I2_tau1*I2_tau2),
           be.tau[4]*be.tau[1]*mean(I2_tau2*I1_tau1),  be.tau[4]*be.tau[2]*mean(I2_tau2*I1_tau2),  be.tau[4]*be.tau[3]*mean(I2_tau2*I2_tau1),   be.tau[4]^2*mean(I2_tau2))
    Jn = matrix(Jn,ncol=4,nrow=4)
    
    tryCatch( {eta = solve(Jn)%*%Un},error=function(cond) {erro <<-  TRUE})
    if(erro ==TRUE| any(is.na(eta))) {eta=rep(0,4)}
    if (erro == FALSE)
    {
      ### Update Tau ###
      eta=as.numeric(eta)
      tau0 = tau0 + eta
    }
    
    j=j+1}
  
  tau.out = c(NA,NA,NA,NA)
  Test = (erro == FALSE)
  if ( Test){ tau.out = tau0}
  
  return(tau.out)
}

sd_continuous<- function(tau.out)
{
  erro = FALSE
  tau0_1 =  tau.out[1];tau0_2 =  tau.out[2];
  tau0_3 =  tau.out[3];tau0_4 =  tau.out[4];
  
  X1_tau1 =  as.numeric(X1 > tau0_1) * ( X1 - tau0_1 )
  X1_tau2 =  as.numeric(X1 > tau0_2) * ( X1 - tau0_2 )
  X2_tau1 =  as.numeric(X2 > tau0_3) * ( X2 - tau0_3 )
  X2_tau2 =  as.numeric(X2 > tau0_4) * ( X2 - tau0_4 )
  
  I1_tau1 = as.numeric(X1 > tau0_1)
  I1_tau2 = as.numeric(X1 > tau0_2)
  I2_tau1 = as.numeric(X2 > tau0_3)
  I2_tau2 = as.numeric(X2 > tau0_4)
  
  lm.inital <- lm(Y~X1+ X1_tau1+X1_tau2+X2+X2_tau1+X2_tau2+Z,family="poisson",
                    control = list(maxit = 25))
  theta0 =coef( lm.inital) 
  be.tau = c(theta0['X1_tau1'],theta0['X1_tau2'],theta0['X2_tau1'],theta0['X2_tau2'])
  sig0 = sigma(lm.inital)
  
  Jn1 =   c(     1,         mean(X1),         mean(X1_tau1),          mean(X1_tau2),        mean(X2),            mean(X2_tau1),             mean(X2_tau2),               mean(Z),
            mean(X1),       mean(X1*X1),      mean(X1_tau1*X1),       mean(X1_tau2*X1),     mean(X2*X1),         mean(X2_tau1*X1),          mean(X2_tau2*X1),            mean(Z*X1),
            mean(X1_tau1),  mean(X1*X1_tau1), mean(X1_tau1*X1_tau1),  mean(X1_tau2*X1_tau1), mean(X2*X1_tau1),   mean(X2_tau1*X1_tau1),     mean(X2_tau2*X1_tau1),       mean(Z*X1_tau1),
            mean(X1_tau2),  mean(X1*X1_tau2), mean(X1_tau1*X1_tau2),  mean(X1_tau2*X1_tau2), mean(X2*X1_tau2),   mean(X2_tau1*X1_tau2),     mean(X2_tau2*X1_tau2),       mean(Z*X1_tau2),
            mean(X2),       mean(X1*X2),      mean(X1_tau1*X2),       mean(X1_tau2*X2),      mean(X2*X2),        mean(X2_tau1*X2),          mean(X2_tau2*X2),            mean(Z*X2),
            mean(X2_tau1),  mean(X1*X2_tau1), mean(X1_tau1*X2_tau1),  mean(X1_tau2*X2_tau1), mean(X2*X2_tau1),   mean(X2_tau1*X2_tau1),     mean(X2_tau2*X2_tau1),       mean(Z*X2_tau1),
            mean(X2_tau2),  mean(X1*X2_tau2), mean(X1_tau1*X2_tau2),  mean(X1_tau2*X2_tau2), mean(X2*X2_tau2),   mean(X2_tau1*X2_tau2),     mean(X2_tau2*X2_tau2),       mean(Z*X2_tau2),
            mean(Z),        mean(X1*Z),       mean(X1_tau1*Z),        mean(X1_tau2*Z),       mean(X2*Z),         mean(X2_tau1*Z),           mean(X2_tau2*Z),             mean(Z*Z))
  Jn1 = matrix( Jn1 , ncol=8,nrow=8)
  
  Jn2=-c(be.tau[1]*mean(I1_tau1),                              be.tau[2]*mean(I1_tau2),                          be.tau[3]*mean(I2_tau1),                           be.tau[4]*mean(I2_tau2),
         be.tau[1]*mean(X1*I1_tau1),                           be.tau[2]*mean(X1*I1_tau2),                       be.tau[3]*mean(X1*I2_tau1),                        be.tau[4]*mean(X1*I2_tau2),
         be.tau[1]*mean(X1_tau1),                              be.tau[2]*mean(X1_tau1*I1_tau2),                  be.tau[3]*mean(X1_tau1*I2_tau1),                   be.tau[4]*mean(X1_tau1*I2_tau2),
         be.tau[1]*mean(X1_tau2*I1_tau1),                      be.tau[2]*mean(X1_tau2),                          be.tau[3]*mean(X1_tau2*I2_tau1),                   be.tau[4]*mean(X1_tau2*I2_tau2),
         be.tau[1]*mean(X2*I1_tau1),                           be.tau[2]*mean(X2*I1_tau2),                       be.tau[3]*mean(X2*I2_tau1),                        be.tau[4]*mean(X2*I2_tau2),
         be.tau[1]*mean(X2_tau1*I1_tau1),                      be.tau[2]*mean(X2_tau1*I1_tau2),                  be.tau[3]*mean(X2_tau1),                           be.tau[4]*mean(X2_tau1*I2_tau2),
         be.tau[1]*mean(X2_tau2*I1_tau1),                      be.tau[2]*mean(X2_tau2*I1_tau2),                  be.tau[3]*mean(X2_tau2*I2_tau1),                   be.tau[4]*mean(X2_tau2),
         be.tau[1]*mean(Z*I1_tau1),                            be.tau[2]*mean(Z*I1_tau2),                        be.tau[3]*mean(Z*I2_tau1),                         be.tau[4]*mean(Z*I2_tau2))
  Jn2 = matrix( Jn2 , ncol=4,nrow=8,byrow = T)
  Jn3 = t(Jn2)
  
  Jn4 =c(be.tau[1]^2*mean(I1_tau1),                  be.tau[1]* be.tau[2]*mean(I1_tau1*I1_tau2),  be.tau[1]*be.tau[3]*mean(I1_tau1*I2_tau1),   be.tau[1]*be.tau[4]*mean(I1_tau1*I2_tau2),
         be.tau[2]*be.tau[1]*mean(I1_tau2*I1_tau1),  be.tau[2]^2*mean(I1_tau2)  ,                 be.tau[2]*be.tau[3]*mean(I1_tau2*I2_tau1),   be.tau[2]*be.tau[4]*mean(I1_tau2*I2_tau2),                                                      
         be.tau[3]*be.tau[1]*mean(I2_tau1*I1_tau1),  be.tau[3]* be.tau[2]*mean(I2_tau1*I1_tau2),  be.tau[3]^2*mean(I2_tau1),                   be.tau[3]*be.tau[4]*mean(I2_tau1*I2_tau2),
         be.tau[4]*be.tau[1]*mean(I2_tau2*I1_tau1),  be.tau[4]* be.tau[2]*mean(I2_tau2*I1_tau2),  be.tau[4]*be.tau[3]*mean(I2_tau2*I2_tau1),   be.tau[4]^2*mean(I2_tau2))
  Jn4 = matrix(Jn4,ncol=4,nrow=4)
  
  sd.out =c(NA,NA,NA,NA)
  tryCatch( { V = solve(Jn4 - Jn3%*%solve(Jn1)%*%Jn2)/(n-12)},
            error=function(cond) {erro <<- TRUE})
  
  if (erro == FALSE)
  {
    sd_1 = sig0*sqrt(V[1,1])
    sd_2 = sig0*sqrt(V[2,2])
    sd_3 = sig0*sqrt(V[3,3])
    sd_4 = sig0*sqrt(V[4,4])
    sd.out = c(sd_1, sd_2, sd_3, sd_4)}
  return(sd.out)
}
###################################################
#########        END OF FUNCTIONS    ##############
###################################################



##########################################################################
#####                    Data Generating Procedure               #########
##########################################################################
#true parameters and true change points
theta = c(1.5, 2, -5, 4,-2,6,-8)
tau=c(0.4,0.7,0.8,1.4)
eta = 0.2

#number of MC replicates
Nsim = 1000

#initial value of change points
tau00=c(0.45,0.8,0.9,1.5) 

#sample size: 250, 500 and 1000
nlst=c(250,500,1000)
length.nlst=length(nlst)
for (nn in 1:length.nlst)
{
  n=nlst[nn]
  i = 1 ; o = 0
  tau_est1 = NULL;tau_est2 = NULL;tau_est3 = NULL;tau_est4 = NULL
  sd_est1 = NULL; sd_est2 = NULL;sd_est3 = NULL; sd_est4 = NULL
  cov1 = 0; cov2 = 0; cov3 = 0; cov4 = 0
   ####### Simulation Procedure ###########
  while (o < Nsim)
{
  ### generate dataset ###
  Z = rnorm(n,0,1)
  a1=  rnorm(n,0,1);u = (a1+Z)/sqrt(2)
  a2=  rnorm(n,0,2);v = (a2+Z)/sqrt(5)
  X1= pnorm(u,mean=0,sd=1)     #X1 and Z are correlated and X1~uniform(0,1)
  X2= pnorm(v,mean=0,sd=1)*2   #X2 and Z are correlated and X2~uniform(0,2)
  e=rnorm(n,0,sd=0.2)
  Y=theta[1] + theta[2]*X1 + theta[3] * as.numeric(X1 > tau[1]) * ( X1 - tau[1] ) + theta[4] * as.numeric(X1 > tau[2]) * ( X1 - tau[2] ) + 
    theta[5]*X2 + theta[6]*as.numeric(X2 > tau[3]) * (X2- tau[3])+theta[7] * as.numeric(X2 > tau[4]) * ( X2 - tau[4] )+ eta*Z+e
  
  ### estimated change points via To-SNR algorithm ###
  tau.out = ToSNR_continuous(tau00)
  
  ### calculate standard error of estimated change points ###
  if (!any(is.na(tau.out)))
  {  
    tau.out = sort(tau.out) 
    sd.out= sd_continuous(tau.out)
    if (!any(is.na( sd.out))){  
      tau_est1[i] = tau.out[1]
      tau_est2[i] = tau.out[2]
      tau_est3[i] = tau.out[3]
      tau_est4[i] = tau.out[4]
      
      sd_est1[i] = sd.out[1]
      sd_est2[i] = sd.out[2]
      sd_est3[i] = sd.out[3]
      sd_est4[i] = sd.out[4]
      
      ## whether true change points within the 95% CI ##
      cov1 = cov1 + as.numeric(tau.out[1]-1.96*sd.out[1] <= tau[1] & tau.out[1]+1.96*sd.out[1]>=tau[1])
      cov2 = cov2 + as.numeric(tau.out[2]-1.96*sd.out[2] <= tau[2] & tau.out[2]+1.96*sd.out[2]>=tau[2])
      cov3 = cov3 + as.numeric(tau.out[3]-1.96*sd.out[3] <= tau[3] & tau.out[3]+1.96*sd.out[3]>=tau[3])
      cov4 = cov4 + as.numeric(tau.out[4]-1.96*sd.out[4] <= tau[4] & tau.out[4]+1.96*sd.out[4]>=tau[4])
      
      i=i+1}}
  o=o+1
  
}

print(paste0("sample size: ", n))
print("To-SNR Method")
#first change point results
bias_tau1 = mean(tau_est1) -  tau[1] 
sd_tau1 = sd(tau_est1)
rmse_tau1=sqrt(bias_tau1^2+sd_tau1 ^2)
mean_sd1 = mean(sd_est1)
#Output:  Bias; MSCD; MSE; AVESE; CP%
print("change point 1")
print(c(round(c(bias_tau1,sd_tau1,rmse_tau1,mean_sd1 )*10,3),round(cov1/(i-1)*100,1)))

#second change point results 
bias_tau2 = mean(tau_est2) -  tau[2] 
sd_tau2 = sd(tau_est2)
rmse_tau2=sqrt(bias_tau2^2+sd_tau2 ^2)
mean_sd2 = mean(sd_est2)
print("change point 2")
print(c(round(c(bias_tau2,sd_tau2,rmse_tau2,mean_sd2 )*10,3),round(cov2/(i-1)*100,1))) 

#third change point results
bias_tau3 = mean(tau_est3) -  tau[3]
sd_tau3 = sd(tau_est3)
rmse_tau3=sqrt(bias_tau3^2+sd_tau3 ^2)
mean_sd3 = mean(sd_est3)
print("change point 3")
print(c(round(c(bias_tau3,sd_tau3,rmse_tau3,mean_sd3 )*10,3),round(cov3/(i-1)*100,1)))

#fourth change point results
bias_tau4 = mean(tau_est4) - tau[4]
sd_tau4 = sd(tau_est4)
rmse_tau4=sqrt(bias_tau4^2+sd_tau4 ^2)
mean_sd4 = mean(sd_est4)
#Output:  Bias; MSCD; MSE; AVESE; CP%
print("change point 4")
print(c(round(c(bias_tau4,sd_tau4,rmse_tau4,mean_sd4 )*10,3),round(cov4/(i-1)*100,1)))


cat("\n\n")

}