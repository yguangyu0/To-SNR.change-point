#################################################################################
########### Simulation Codes For Logistic Model With Two Change Points ##########
#################################################################################

#install.packages("segmented")
require("segmented")
#set.seed(20231113)
set.seed(20231006)
NR_semismooth<-function(tau0)
{
  ### generate inital values ###
  X1_tau1 =  as.numeric(X1 > tau0[1]) * ( X1 - tau0[1] )
  X2_tau1 =  as.numeric(X2 > tau0[2]) * ( X2 - tau0[2] )
  glm.inital <- glm(Y~X1+X1_tau1+X2+X2_tau1+Z,family="binomial",control = list(maxit = 25))
  tau0=c(as.numeric(coef( glm.inital)),tau0)
  
  erro = FALSE;eta = rep(1,2); j = 1
  while (any(abs(eta) > 1e-7) & j<=100 & erro==FALSE)
  {
    X1_tau1 =  as.numeric(X1 > tau0[7]) * ( X1 - tau0[7] )
    X2_tau1 =  as.numeric(X2 > tau0[8]) * ( X2 - tau0[8] )
    mui=tau0[1] + tau0[2]*X1 + tau0[3] * X1_tau1+tau0[4]*X2 +tau0[5] * X2_tau1+tau0[6]*Z
    mu=exp(mui)/(1+exp(mui))
    mu0 = mu*(1-mu);A = Y-mu
    
    I1_tau1 = as.numeric(X1 > tau0[7]);I2_tau1 = as.numeric(X2 > tau0[8])
    
    ### Calculate U and J ###
    Un =  c(mean(A),mean(A*X1),mean(A*X1_tau1),mean(A*X2),mean(A*X2_tau1),mean(A*Z),-tau0[3]*mean(A*I1_tau1),-tau0[5]*mean(A*I2_tau1))
    Un = matrix(Un,nrow=8)
    Jn =  c(mean(mu0),         mean(mu0*X1),          mean(mu0*X1_tau1),        mean(mu0*X2),        mean(mu0*X2_tau1),         mean(mu0*Z),         -tau0[3]*mean(mu0*I1_tau1),     -tau0[5]*mean(mu0*I2_tau1), 
            mean(mu0*X1),      mean(mu0*X1^2),        mean(mu0*X1*X1_tau1),     mean(mu0*X1*X2),     mean(mu0*X1*X2_tau1),      mean(mu0*X1*Z),     -tau0[3]*mean(mu0*X1*I1_tau1),  -tau0[5]*mean(mu0*X1*I2_tau1), 
            mean(mu0*X1_tau1), mean(mu0*X1*X1_tau1),  mean(mu0*X1_tau1^2),      mean(mu0*X1_tau1*X2),mean(mu0*X1_tau1*X2_tau1), mean(mu0*X1_tau1*Z),-tau0[3]*mean(mu0*X1_tau1),     -tau0[5]*mean(mu0*X1_tau1*I2_tau1),
            mean(mu0*X2),      mean(mu0*X1*X2),       mean(mu0*X2*X1_tau1),     mean(mu0*X2^2),      mean(mu0*X2*X2_tau1),      mean(mu0*X2*Z),     -tau0[3]*mean(mu0*X2*I1_tau1),  -tau0[5]*mean(mu0*X2*I2_tau1),
            mean(mu0*X2_tau1), mean(mu0*X1*X2_tau1),  mean(mu0*X1_tau1*X2_tau1),mean(mu0*X2*X2_tau1),mean(mu0*X2_tau1^2),       mean(mu0*X2_tau1*Z),-tau0[3]*mean(mu0*X2_tau1*I1_tau1),-tau0[5]*mean(mu0*X2_tau1),
            mean(mu0*Z),       mean(mu0*X1*Z),        mean(mu0*X1_tau1*Z),      mean(mu0*X2*Z),      mean(mu0*X2_tau1*Z),       mean(mu0*Z^2),      -tau0[3]*mean(mu0*Z*I1_tau1),   -tau0[5]*mean(mu0*Z*I2_tau1),
            -tau0[3]*mean(mu0*I1_tau1),-tau0[3]*mean(mu0*X1*I1_tau1),-tau0[3]*mean(mu0*X1_tau1),  -tau0[3]*mean(mu0*X2*I1_tau1),-tau0[3]*mean(mu0*X2_tau1*I1_tau1),-tau0[3]*mean(mu0*Z*I1_tau1),tau0[3]^2*mean(mu0*I1_tau1*I1_tau1),tau0[3]*tau0[5]*mean(mu0*I1_tau1*I2_tau1),
            -tau0[5]*mean(mu0*I2_tau1),-tau0[5]*mean(mu0*X1*I2_tau1),-tau0[5]*mean(mu0*X1_tau1*I2_tau1),-tau0[5]*mean(mu0*X2*I2_tau1),-tau0[5]*mean(mu0*X2_tau1), -tau0[5]*mean(mu0*Z*I2_tau1),tau0[3]*tau0[5]*mean(mu0*I1_tau1*I2_tau1),tau0[5]^2*mean(mu0*I2_tau1*I2_tau1))
    Jn = matrix( Jn , ncol=8,nrow=8)
    
    tryCatch( {eta = solve(Jn)%*%Un},error=function(cond) {erro <<-  TRUE})
    if(erro ==TRUE| any(is.na(eta))) {eta=rep(0,8)}
    if (erro == FALSE)
    {### Update Tau ###
      eta =as.numeric(eta)
      tau0 = tau0 + eta}
    j=j+1}
  
  tau.out = c(NA,NA)
  Test = (erro == FALSE)
  if ( Test  )
  { tau.out = tau0[7:8]}
  
  return(tau.out)
}

ToSNR_algorithm<-function(tau0)
{
  #tau1_c=tau0[1];
  erro = FALSE;eta = rep(1,2); j = 1
  while (any(abs(eta) > 1e-7) & j<=100 & erro==FALSE)
  {
    tau0_1 = tau0[1]; tau0_2 = tau0[2]
    ### generate inital values ###
    X1_tau1 =  as.numeric(X1 > tau0_1) * ( X1 - tau0_1 )
    X2_tau1 =  as.numeric(X2 > tau0_2) * ( X2 - tau0_2 )
    
    I1_tau1 = as.numeric(X1 > tau0_1)
    I2_tau1 = as.numeric(X2 > tau0_2)
    
    glm.inital <- glm(Y~X1+X1_tau1+X2+X2_tau1+Z,family="binomial",
                      control = list(maxit = 25))
    
    theta0 =coef( glm.inital) 
    be.tau = c(theta0['X1_tau1'],theta0['X2_tau1'])
    
    mu = predict(glm.inital, type="response")
    mu0 = mu*(1-mu)
    
    ### Calculate U and J ###
    Un =  c(be.tau[1]*mean((mu-Y)*I1_tau1),
            be.tau[2]*mean((mu-Y)*I2_tau1))
    Un = matrix(Un,nrow=2)
    
    Jn =c(be.tau[1]^2*mean(mu0*I1_tau1*I2_tau1),         be.tau[1]*be.tau[2]*mean(mu0*I1_tau1*I2_tau1),
          be.tau[1]*be.tau[2]*mean(mu0*I1_tau1*I2_tau1), be.tau[2]^2*mean(mu0*I2_tau1*I2_tau1))
    
    Jn = matrix(Jn,ncol=2,nrow=2)
    
    tryCatch( {eta = solve(Jn)%*%Un},error=function(cond) {erro <<-  TRUE})
    if(erro ==TRUE| any(is.na(eta))) {eta=rep(0,2)}
    if (erro == FALSE)
    {### Update Tau ###
      tau_up = tau0 + eta
      tau0 = tau_up}
    j=j+1}
  
  tau.out = c(NA,NA)
  Test = (erro == FALSE && !any(is.na(eta)))
  if ( Test  )
  { tau.out = tau0}
  
  return(tau.out)
}

sd_logistic<- function(tau.out)
{
  erro=FALSE
  
  tau0_1 = tau.out[1]; tau0_2 = tau.out[2]
  X1_tau1 =  as.numeric(X1 > tau0_1) * ( X1 - tau0_1 )
  X2_tau1 =  as.numeric(X2 > tau0_2) * ( X2 - tau0_2 )
  I1_tau1 = as.numeric(X1 > tau0_1)
  I2_tau1 = as.numeric(X2 > tau0_2)
  
  glm.inital <- glm(Y~X1+X1_tau1+X2+X2_tau1+Z,family="binomial",
                    control = list(maxit = 25))
  
  theta0 =coef( glm.inital) 
  be.tau = c(theta0['X1_tau1'],theta0['X2_tau1'])
  mu = predict(glm.inital, type="response")
  
  mu0 = mu*(1-mu)
  Jn1 =  c(mean(mu0),         mean(mu0*X1),          mean(mu0*X1_tau1),        mean(mu0*X2),        mean(mu0*X2_tau1),         mean(mu0*Z),       
           mean(mu0*X1),      mean(mu0*X1^2),        mean(mu0*X1*X1_tau1),     mean(mu0*X1*X2),     mean(mu0*X1*X2_tau1),      mean(mu0*X1*Z), 
           mean(mu0*X1_tau1), mean(mu0*X1*X1_tau1),  mean(mu0*X1_tau1^2),      mean(mu0*X1_tau1*X2),mean(mu0*X1_tau1*X2_tau1), mean(mu0*X1_tau1*Z), 
           mean(mu0*X2),      mean(mu0*X1*X2),       mean(mu0*X2*X1_tau1),     mean(mu0*X2^2),      mean(mu0*X2*X2_tau1),      mean(mu0*X2*Z), 
           mean(mu0*X2_tau1), mean(mu0*X1*X2_tau1),  mean(mu0*X1_tau1*X2_tau1),mean(mu0*X2*X2_tau1),mean(mu0*X2_tau1^2),       mean(mu0*X2_tau1*Z), 
           mean(mu0*Z),       mean(mu0*X1*Z),        mean(mu0*X1_tau1*Z),      mean(mu0*X2*Z),      mean(mu0*X2_tau1*Z),       mean(mu0*Z^2))
  Jn1 = matrix( Jn1 , ncol=6,nrow=6)
  
  Jn2 = -  c(be.tau[1]*mean(mu0*I1_tau1),             be.tau[2]*mean(mu0*I2_tau1),
             be.tau[1]*mean(mu0*X1*I1_tau1),          be.tau[2]*mean(mu0*X1*I2_tau1),
             be.tau[1]*mean(mu0*X1_tau1),             be.tau[2]*mean(mu0*X1_tau1*I2_tau1),
             be.tau[1]*mean(mu0*X2*I1_tau1),          be.tau[2]*mean(mu0*X2*I2_tau1),
             be.tau[1]*mean(mu0*X2_tau1*I1_tau1),     be.tau[2]*mean(mu0*X2_tau1),
             be.tau[1]*mean(mu0*Z*I1_tau1),           be.tau[2]*mean(mu0*Z*I2_tau1))
  Jn2 = matrix( Jn2 , ncol=2,nrow=6,byrow = T)
  Jn3 = t(Jn2)
  
  Jn4 =   c(be.tau[1]^2*mean(mu0*I1_tau1*I1_tau1),              be.tau[1]*be.tau[2]*mean(mu0*I1_tau1*I2_tau1),
            be.tau[2]*be.tau[1]*mean(mu0*I1_tau1*I2_tau1),      be.tau[2]^2*mean(mu0*I2_tau1*I2_tau1))
  Jn4 = matrix(Jn4,ncol=2,nrow=2)
  
  tryCatch( {J_inverse = solve(Jn4 - Jn3%*%solve(Jn1)%*%Jn2)},error=function(cond) {erro <<-  TRUE})
  
  sd.out = c(NA,NA)
  if (erro==FALSE)
  {V = J_inverse/(n)
  
  sd_1 = sqrt(V[1,1])
  sd_2 = sqrt(V[2,2])
  
  sd.out = c(sd_1, sd_2)}
  
  return(sd.out)
}



##### Data Generating Procedure #########
#true parameters and true change points
theta = c(-1.5, 6, -9.5,3,-8);tau=c(0.45,1.0);eta = 0.2

#number of MC replicates
Nsim = 1000

#initial value of change points
tau00=c(0.5,0.9)


nlst=c(250,500,1000)
length.nlst=length(nlst)
for (nn in 1:length.nlst)
{
  #sample size
  n=nlst[nn]
  i=1;i.2=1;i.seg=1; o = 0
  tau_segest1=NULL;tau_segest2=NULL;sd_segest1 =NULL;sd_segest2 =NULL;cov_seg1=0;cov_seg2=0
  tau_est1 = NULL;tau_est2 = NULL;sd_est1 = NULL; sd_est2 = NULL;cov1 = 0; cov2 = 0
  tau_est.1 = NULL;tau_est.2 = NULL;sd_est.1 = NULL; sd_est.2 = NULL;cov.1 = 0; cov.2 = 0
  ####### Simulation Procedure ###########
  while (o < Nsim)
  {
    erro = FALSE
    ### generate dataset ###
    Z = rnorm(n,0,1)
    a1=  rnorm(n,0,1);u = (a1+Z)/sqrt(2)
    a2=  rnorm(n,0,2);v = (a2+Z)/sqrt(5)
    X1= pnorm(u,mean=0,sd=1)     #X1 and Z are correlated and X1~uniform(0,1)
    X2= pnorm(v,mean=0,sd=1)*2   #X2 and Z are correlated and X2~uniform(0,2)
  
    t= theta[1] + theta[2]*X1 + theta[3] * as.numeric(X1 > tau[1]) * ( X1 - tau[1] )+
      theta[4]*X2 +theta[5] * as.numeric(X2 > tau[2]) * ( X2 - tau[2] )+eta*Z
    pr = 1/(1+exp(-t)) 
    Y = rbinom(n,1,pr)
    
    
    ### estimated change points via segmented package ###
    o1<-glm(Y~X1+X2+Z,family=binomial)
    tryCatch( {o.seg=segmented(o1,seg.Z =~X1+X2,psi=list(X1=tau00[1],X2=tau00[2]),
                               seg.control(toll = 1e-7, it.max = 100,nboot=0,random=F))},
              error=function(cond) {erro <<- TRUE})
    if (erro == FALSE & !is.null(o.seg$psi))
    {seg.out = o.seg$psi[,2]
    seg.sd = o.seg$psi[,3]
    
    ### estimated change points ###
    tau_segest1[i.seg]=seg.out[1]
    tau_segest2[i.seg]=seg.out[2]
    
    ### standard error of estimated change points ###
    sd_segest1[i.seg]=seg.sd [1]
    sd_segest2[i.seg]=seg.sd [2]
    
    ## whether true change points within the 95% CI ##
    ci1=confint(o.seg,"X1");ci2=confint(o.seg,"X2")
    cov_seg1 = cov_seg1  + as.numeric(ci1[2] <= tau[1] & ci1[3] >=tau[1])
    cov_seg2 = cov_seg2  + as.numeric(ci2[2] <= tau[2] & ci2[3] >=tau[2])
    
    i.seg=i.seg+1
    }
    
    ### estimated change points via To-SNR algorithm ###
    tau.out=ToSNR_algorithm(tau00)
    ### calculate standard error of estimated change points ###
    if (!any(is.na(tau.out)))
    {  
      sd.out= sd_logistic(tau.out)
      if (!any(is.na( sd.out)))
      {  
        tau_est1[i] = tau.out[1]
        tau_est2[i] = tau.out[2]
        
        sd_est1[i] = sd.out[1]
        sd_est2[i] = sd.out[2]
        
        ## whether true change points within the 95% CI ##
        cov1 = cov1 + as.numeric(tau.out[1]-1.96*sd.out[1] <= tau[1] & tau.out[1]+1.96*sd.out[1]>=tau[1])
        cov2 = cov2 + as.numeric(tau.out[2]-1.96*sd.out[2] <= tau[2] & tau.out[2]+1.96*sd.out[2]>=tau[2])
        
        i=i+1}}
    
    ### estimated change points via To-SNR algorithm ###
    tau.out2=NR_semismooth(tau00)
    ### calculate standard error of estimated change points ###
    if (!any(is.na(tau.out2)))
    {  
      sd.out2= sd_logistic(tau.out2)
      if (!any(is.na( sd.out2)))
      {  
        tau_est.1[i.2] = tau.out2[1]
        tau_est.2[i.2] = tau.out2[2]
        
        sd_est.1[i.2] = sd.out2[1]
        sd_est.2[i.2] = sd.out2[2]
        
        ## whether true change points within the 95% CI ##
        cov.1 = cov.1 + as.numeric(tau.out2[1]-1.96*sd.out2[1] <= tau[1] & tau.out2[1]+1.96*sd.out2[1]>=tau[1])
        cov.2 = cov.2 + as.numeric(tau.out2[2]-1.96*sd.out2[2] <= tau[2] & tau.out2[2]+1.96*sd.out2[2]>=tau[2])
        
        i.2=i.2+1}}
    
    
    
    o=o+1
  }
  print(paste0("sample size: ", n))
  #segmented method
  #first change point results
  bias_tau1_seg = mean(tau_segest1) -  tau[1]
  sd_tau1_seg = sd(tau_segest1)
  rmse_tau1_seg=sqrt(bias_tau1_seg^2+sd_tau1_seg ^2)
  mean_sd1_seg = mean(sd_segest1)
  #second change point results  
  bias_tau2_seg = mean(tau_segest2) -  tau[2]
  sd_tau2_seg = sd(tau_segest2)
  rmse_tau2_seg=sqrt(bias_tau2_seg^2+sd_tau2_seg ^2)
  mean_sd2_seg = mean(sd_segest2)
  #Output:  Bias; MSCD; MSE; AVESE; CP%
  print("Muggeo Method")
  print("change point 1")
  print(c(round(c(bias_tau1_seg,sd_tau1_seg,rmse_tau1_seg,mean_sd1_seg )*10,3),round(cov_seg1/(i.seg-1)*100,1)))
  print("change point 2")
  print(c(round(c(bias_tau2_seg,sd_tau2_seg,rmse_tau2_seg,mean_sd2_seg )*10,3),round(cov_seg2/(i.seg-1)*100,1)))

  
  #To-SNR method
  bias_tau1 = mean(tau_est1) - tau[1] 
  sd_tau1 = sd(tau_est1)
  rmse_tau1=sqrt(bias_tau1^2+sd_tau1 ^2)
  mean_sd1 = mean(sd_est1)
  #second change point results  
  bias_tau2 = mean(tau_est2) - tau[2] 
  sd_tau2 = sd(tau_est2)
  rmse_tau2=sqrt(bias_tau2^2+sd_tau2 ^2)
  mean_sd2 = mean(sd_est2)
  #Output:  Bias; MSCD; MSE; AVESE; CP%
  print("To-SNR Method")
  print("change point 1")
  print(c(round(c(bias_tau1,sd_tau1,rmse_tau1,mean_sd1 )*10,3),round(cov1/(i-1)*100,1)))
  print("change point 2")
  print(c(round(c(bias_tau2,sd_tau2,rmse_tau2,mean_sd2 )*10,3),round(cov2/(i-1)*100,1))) 
  
  #NR-semismooth method
  bias_tau.1 = mean(tau_est.1) - tau[1] 
  sd_tau.1 = sd(tau_est.1)
  rmse_tau.1=sqrt(bias_tau.1^2+sd_tau.1 ^2)
  mean_sd.1 = mean(sd_est.1)
  #second change point results  
  bias_tau.2 = mean(tau_est.2) - tau[2] 
  sd_tau.2 = sd(tau_est.2)
  rmse_tau.2=sqrt(bias_tau.2^2+sd_tau.2 ^2)
  mean_sd.2 = mean(sd_est.2)
  #Output:  Bias; MSCD; MSE; AVESE; CP%
  print("NR-Semi Method")
  print("change point 1")
  print(c(round(c(bias_tau.1,sd_tau.1,rmse_tau.1,mean_sd.1 )*10,3),round(cov.1/(i.2-1)*100,1)))
  print("change point 2")
  print(c(round(c(bias_tau.2,sd_tau.2,rmse_tau.2,mean_sd.2 )*10,3),round(cov.2/(i.2-1)*100,1))) 
  
  cat("\n\n")
}
