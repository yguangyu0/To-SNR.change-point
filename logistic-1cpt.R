#################################################################################
########### Simulation Codes For Logistic Model With One Change Point  ##########
#################################################################################
#install.packages("segmented")
#install.packages("chngpt")

set.seed(20231006)
library(segmented)
library(chngpt)

###################################################
#############   Functions  ########################
###################################################
NR_semismooth<-function(tau0_1)
{
  #generate initial value of betas
  X_tau_1 =  as.numeric(X > tau0_1) * ( X - tau0_1 )
  glm.inital <- glm(Y~0+X+X_tau_1+Z,family="binomial",control = list(maxit = 25))
  tau0=c(as.numeric(coef( glm.inital)),tau0_1)
  
  erro = FALSE;eta = rep(1,4); j = 1
  
  while (any(abs(eta) > 1e-7) & j<=100 & erro==FALSE)
  {
    ### generate initial values ###
    X_tau_1 =  as.numeric(X > tau0[4]) * ( X - tau0[4] )
    I_tau_1 = as.numeric(X > tau0[4])
    
    mui=tau0[1]*X+tau0[2]*X_tau_1+tau0[3]*Z
    mu=exp(mui)/(1+exp(mui))
    mu0 = mu*(1-mu);A = Y-mu
    
    ### Calculate U and J ###
    U =c(mean(A*X),mean(A*X_tau_1),mean(A*Z),-tau0[2]*mean(A*I_tau_1))
    U = matrix(U,ncol=1,nrow=4)
    J =c(mean(mu0*X^2),        mean(mu0*X*X_tau_1),       mean(mu0*X*Z),         - tau0[2]*mean(mu0*X*I_tau_1), 
         mean(mu0*X*X_tau_1),  mean(mu0*X_tau_1^2),       mean(mu0*X_tau_1*Z),   - tau0[2]*mean(mu0*X_tau_1),
         mean(mu0*X*Z),        mean(mu0*X_tau_1*Z),       mean(mu0*Z^2),         - tau0[2]*mean(mu0*Z*I_tau_1),
         - tau0[2]*mean(mu0*X*I_tau_1),- tau0[2]*mean(mu0*X_tau_1),- tau0[2]*mean(mu0*Z*I_tau_1),tau0[2]^2*mean(mu0*I_tau_1))
    J=matrix(J,ncol=4,nrow=4)
    
    tryCatch( {eta =  solve(J)%*%U},error=function(cond) {erro <<-  TRUE})
    if(erro ==TRUE| any(is.na(eta))) {eta=rep(0,4)}
    
    if (erro == FALSE)
    {
      eta =as.numeric(eta)
      ### Update Tau ###
      tau0 = tau0 + eta}
    j=j+1
  }
  Test = (erro == FALSE &  tau0[4]<max(X) &tau0[4]>min(X))
  tau.out =NA
  
  if ( Test )  { tau.out = tau0}
  return(tau.out[4])
}

NR_logistic<-function(tau0)
{
  #tau1_c=tau0[1];
  erro = FALSE
  eta = 1; j = 1
  
  while (abs(eta) > 1e-7 & j<=100 & erro==FALSE)
  {
    tau0_1 = tau0
    ### generate inital values ###
    X_tau_1 =  as.numeric(X > tau0_1) * ( X - tau0_1 )
    
    glm.inital <- glm(Y~0+X+ X_tau_1+Z,family="binomial",
                      control = list(maxit = 25))
    
    theta0 = c(as.numeric(coef( glm.inital)),tau0_1)
    mu0 = predict(glm.inital, type="response")
    
    I_tau_1 = as.numeric(X > tau0_1)
    
    ### Calculate U and J ###
    Un =  -mean((Y-mu0)*I_tau_1)
    Jn =theta0[2]*mean(mu0*(1-mu0)*I_tau_1)
    
    tryCatch( {eta = Un/Jn},error=function(cond) {erro <<-  TRUE})
    if(erro ==TRUE| any(is.na(eta))) {eta=0}
    
    if (erro == FALSE)
    {
      ### Update Tau ###
      tau0 = tau0 + eta
    }
    j=j+1
  }
  tau.out = NA
  Test = (erro == FALSE && !is.na(eta) & max(tau0,na.rm=T)<max(X) & min(tau0,na.rm=T)>min(X))
  if ( Test )
  { tau.out = tau0}
  return(tau.out)
}

sd_logistic<- function(tau.out)
{
  erro=FALSE
  tau0_1 =  tau.out
  
  X_tau_1 =  as.numeric(X > tau0_1) * ( X - tau0_1 )
  I_tau_1 = as.numeric(X > tau0_1);
  
  glm.inital <- glm(Y~0+X+ X_tau_1+Z,family="binomial",
                    control = list(maxit = 25))
  theta0 = c(as.numeric(coef(glm.inital)),tau0_1)
  mu = predict(glm.inital, type="response")
  
  mu0 = mu*(1-mu)
  Jn1 =  c(mean(mu0*X^2),        mean(mu0*X*X_tau_1),       mean(mu0*X*Z), 
           mean(mu0*X*X_tau_1),  mean(mu0*X_tau_1^2),       mean(mu0*X_tau_1*Z),
           mean(mu0*X*Z),        mean(mu0*X_tau_1*Z),       mean(mu0*Z^2)
  )
  Jn1 = matrix( Jn1 , ncol=3,nrow=3)
  
  Jn2 = - theta0[2]*c(mean(mu0*X*I_tau_1),                            
                      mean(mu0*X_tau_1),
                      mean(mu0*Z*I_tau_1)
  )
  Jn2 = matrix( Jn2 , ncol=1,nrow=3,byrow = T)
  Jn3 = t(Jn2)
  
  Jn4 = theta0[2]^2*mean(mu0*I_tau_1)
  Jn4 = matrix(Jn4,ncol=1,nrow=1)
  
  sd.out = NA
  tryCatch( { V = solve(Jn4 - Jn3%*%solve(Jn1)%*%Jn2)/(n)},
            error=function(cond) {erro <<- TRUE})
  
  if(erro==FALSE){  sd.out = sqrt(V)}
  return(sd.out)
}



##########################################################################
##########################################################################
##########################################################################
##########################################################################







############## Settings ##################
#### true values beta0, beta1, beta2, tau ####
Nsim = 1000 
rho=0.3;c=rho*1.6
nlst=c(250,500,1000)
length.nlst=length(nlst)
for (nn in 1:length.nlst)
{
  n=nlst[nn]
  
  i=1;i.2=1
  tau_est1 = NULL;tau_est2 = NULL;tau_seg=NULL;tau_chngpt=NULL
  sd_est1 = NULL; sd_est2 = NULL; sd_seg=NULL;sd_chngpt=NULL
  cov1 = 0; cov2 = 0; cov.seg=0;cov.chngpt=0
  ####### Simulation Procedure ###########
  for(i in 1:Nsim)
  {

    ### generate dataset ###
    XZ=mvrnorm(n,mu=c(0,4.7),Sigma=matrix(c(1,c,c,1.6^2),ncol=2))
    Z=XZ[,1];X=XZ[,2]
    t=0.4*X-0.92*(X-4.5)*as.numeric(X>4.5)+0.34*Z
    pr = 1/(1+exp(-t)) 
    Y = rbinom(n,1,pr)
    dat=data.frame(Y,X,Z)
    ### segmented method ###
    o1<-glm(Y~0+X+Z,family = "binomial")
    o.seg=segmented(o1,seg.Z=~X,psi=5,control=seg.control(n.boot=0,tol=1e-7,it.max=100))
    seg.out=o.seg$psi[2];seg.sd=o.seg$psi[3]
    tau_seg[i]=seg.out;sd_seg[i]=seg.sd
    ci=confint(o.seg)
    cov.seg = cov.seg + as.numeric(ci[2] <= 4.5 & ci[3]>=4.5)

    
    ### chngpt method ###
    fit.chngpt=chngptm(formula.1=Y~0+Z, formula.2=~X, family="binomial", 
                       data=dat, type="segmented",var.type="model",#bootstrap
                       chngpt.init =5,tol = 1e-7,maxit = 100,
                       est.method = 'smoothapprox')#fastgrid2
    summary.chngpt=summary(fit.chngpt)$chngpt 
    tau_chngpt[i]=summary.chngpt[1]
    sd_chngpt[i]=summary.chngpt[2]
    cov.chngpt = cov.chngpt + as.numeric(summary.chngpt[3] <= 4.5 & summary.chngpt[4]>=4.5)
    
    ### proposed To-SNR method ###
    tau.out=NR_logistic(5)
    sd.out= sd_logistic(tau.out)
    sd_est1[i] = sd.out
    tau_est1[i] = tau.out
    cov1 = cov1 + as.numeric(tau.out-1.96*sd.out <= 4.5 & tau.out+1.96*sd.out>=4.5)
    
    ### NR_semismooth method ###
    tau.out2=NR_semismooth(5)
    if (!is.na( tau.out2))
    {
      sd.out2= sd_logistic(tau.out2)
      sd_est2[i.2] = sd.out2
      tau_est2[i.2] = tau.out2
      cov2 = cov2 + as.numeric(tau.out2-1.96*sd.out2 <= 4.5 & tau.out2+1.96*sd.out2>=4.5)
      i.2=i.2+1
    }
    
  }
  print(paste0("sample size: ", n))
  #chngpt method
  bias_chngpt1 = mean(tau_chngpt) - 4.5
  sd_chngpt1 = sd(tau_chngpt)
  mse_chngpt1=sqrt(bias_chngpt1^2+sd_chngpt1^2)
  mean_sd_chngpt = mean(sd_chngpt)
  print("Fong Method")
  print(c(round(c(bias_chngpt1,sd_chngpt1,mse_chngpt1,mean_sd_chngpt)*10,3),round(cov.chngpt/(i-1)*100,1)))
  
  
  #segmented method
  bias_seg1 = mean(tau_seg) - 4.5
  sd_seg1 = sd(tau_seg)
  mse_seg1=sqrt(bias_seg1^2+sd_seg1^2)
  mean_sd_seg = mean(sd_seg)
  print("Muggeo Method")
  print(c(round(c(bias_seg1,sd_seg1,mse_seg1,mean_sd_seg)*10,3),round(cov.seg/(i-1)*100,1)))
  

  #To-SNR method
  bias_tau1 = mean(tau_est1) - 4.5
  sd_tau1 = sd(tau_est1)
  mse_tau1=sqrt(bias_tau1^2+sd_tau1 ^2)
  mean_sd1 = mean(sd_est1)
  print("To-SNR Method")
  print(c(round(c(bias_tau1,sd_tau1,mse_tau1,mean_sd1)*10,3),round(cov1/(i-1)*100,1)))
  
  #NR-semismooth method
  bias_tau2 = mean(tau_est2) - 4.5
  sd_tau2 = sd(tau_est2)
  mse_tau2=sqrt(bias_tau2^2+sd_tau2 ^2)
  mean_sd2 = mean(sd_est2)
  print("NR-Semi Method")
  print(c(round(c(bias_tau2,sd_tau2,mse_tau2,mean_sd2)*10,3),round(cov2/(i.2-1)*100,1)))
  
  cat("\n\n")
}
