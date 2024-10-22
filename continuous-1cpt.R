#################################################################################
###  Simulation Codes For Linear Segmented Model With One Change Point  #########
#################################################################################
#install.packages("segmented")
#install.packages("chngpt")

set.seed(2023106)
library(segmented)
library(chngpt)
###################################################
#############   Functions  ########################
###################################################
NR_semismooth<-function(tau0_1)
{
  #generate initial value of betas
  X_tau_1 =  as.numeric(X > tau0_1) * ( X - tau0_1 )
  lm.inital <- lm(Y~X+ X_tau_1)
  tau0=c(as.numeric(coef(lm.inital)),tau0_1)
  
  erro = FALSE;eta = rep(1,4); j = 1
  while (any(abs(eta) > 1e-7) & j<=100 & erro==FALSE)
  {
    ### generate initial values ###
    X_tau_1 =  as.numeric(X > tau0[4]) * ( X - tau0[4] )
    I_tau_1 = as.numeric(X > tau0[4])
    
    mu0=tau0[1]+tau0[2]*X+tau0[3]*X_tau_1;A = Y-mu0
    
    ### Calculate U and J ###
    U =c(mean(A),mean(A*X),mean(A*X_tau_1),-tau0[3]*mean(A*I_tau_1))
    U = matrix(U,ncol=1,nrow=4)
    J =c(1,             mean(X),        mean(X_tau_1),   - tau0[3]*mean(I_tau_1),
         mean(X),       mean(X*X),      mean(X_tau_1*X), - tau0[3]*mean(X*I_tau_1),        
         mean(X_tau_1), mean(X*X_tau_1), mean(X_tau_1^2),- tau0[3]*mean(X_tau_1),
         -tau0[3]*mean(I_tau_1),-tau0[3]*mean(X*I_tau_1),- tau0[3]*mean(X_tau_1),   tau0[3]^2*mean(I_tau_1))
    J=matrix(J,ncol=4,nrow=4)
    
    tryCatch( {eta =  solve(J)%*%U},error=function(cond) {erro <<-  TRUE})
    if(erro ==TRUE| any(is.na(eta))) {eta=rep(0,4)}
    if (erro == FALSE)
    { ### Update estimates ###
      eta =as.numeric(eta)
      tau0 = tau0 + eta}
    j=j+1
  }
  Test = (erro == FALSE &  tau0[4]<max(X) &tau0[4]>min(X))
  tau.out =NA
  if ( Test )  { tau.out = tau0}
  return(tau.out[4])
}

ToSNR_continuous<-function(tau0)
{
  erro = FALSE;eta = 1; j = 1
  while (any(abs(eta) > 1e-7) & j<=100 & erro==FALSE)
  {
    ### Step 1. Only update coefficients  ###
    X_tau1 =  as.numeric(X > tau0) * ( X - tau0 )
    lm.inital <- lm(Y~X+ X_tau1)
    theta0 =coef(lm.inital) 
    be.tau = theta0['X_tau1']
    mu0 = predict(lm.inital)
    I_tau1 = as.numeric(X > tau0)
    
    ### Step 2. Only update change point  ###
    Un= mean((mu0-Y)*I_tau1)
    Jn = be.tau[1]*mean(I_tau1)
    tryCatch( {eta = Un/Jn},error=function(cond) {erro <<-  TRUE})
    if(erro ==TRUE| any(is.na(eta))) {eta=0}
    if (erro == FALSE){tau0 = tau0 + eta} ### Update change point ###
    j=j+1}
  
  tau.out =NA
  Test = (erro == FALSE && !is.na(eta))
  if ( Test) { tau.out = tau0}
  
  return(tau.out)
}

sd_continuous<- function(tau0_1)
{
  erro=FALSE
  X_tau_1 =  as.numeric(X > tau0_1) * ( X - tau0_1 )
  I_tau_1 = as.numeric(X > tau0_1)
  
  lm.inital <- lm(Y~X+ X_tau_1)
  theta0 =coef(lm.inital) 
  be.tau = theta0['X_tau_1']
  
  J =c(1,             mean(X),        mean(X_tau_1),   - be.tau*mean(I_tau_1),
       mean(X),       mean(X*X),      mean(X_tau_1*X), - be.tau*mean(X*I_tau_1),        
       mean(X_tau_1), mean(X*X_tau_1), mean(X_tau_1^2),- be.tau*mean(X_tau_1),
       -be.tau*mean(I_tau_1),-be.tau*mean(X*I_tau_1),- be.tau*mean(X_tau_1),   be.tau^2*mean(I_tau_1))
  J=matrix(J,ncol=4,nrow=4)
  sd.out =NA
  tryCatch( { V = solve(J)/(n-4)},error=function(cond) {erro <<- TRUE})
  if (erro == FALSE){ sd.out = sqrt(V[4,4])}
  
  return(sd.out)
}

sd_continuous<- function(tau0_1)
{
  erro=FALSE
  X_tau_1 =  as.numeric(X > tau0_1) * ( X - tau0_1 )
  I_tau_1 = as.numeric(X > tau0_1)
  
  lm.inital <- lm(Y~X+ X_tau_1)
  theta0 =coef(lm.inital) 
  be.tau = theta0['X_tau_1']
  sig=sigma(lm.inital)
  J =c(1,             mean(X),        mean(X_tau_1),   - be.tau*mean(I_tau_1),
       mean(X),       mean(X*X),      mean(X_tau_1*X), - be.tau*mean(X*I_tau_1),        
       mean(X_tau_1), mean(X*X_tau_1), mean(X_tau_1^2),- be.tau*mean(X_tau_1),
       -be.tau*mean(I_tau_1),-be.tau*mean(X*I_tau_1),- be.tau*mean(X_tau_1),   be.tau^2*mean(I_tau_1))
  J=matrix(J,ncol=4,nrow=4)
  sd.out =NA
  tryCatch( { V = solve(J)/(n-4)},error=function(cond) {erro <<- TRUE})
  if (erro == FALSE){ sd.out =  sig*sqrt(V[4,4])}
  
  return(sd.out)
}

###################################################
#########        END OF FUNCTIONS    ##############
###################################################




############## Settings ##################
#### true values beta0, beta1, beta2, tau ####

Nsim = 1000
beta=2.5
#beta=1.8
tau=0.75
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
    X = runif(n,0,1)
    e=rnorm(n,0,sd=0.2)
    Y=3.5- 1.5*X + beta* as.numeric(X > tau) * ( X - tau ) + e
    dat=data.frame(X,Y)
    
    ### chngpt method ###
    fit.chngpt=chngptm(formula.1=Y~1, formula.2=~X, data=dat, family="gaussian", 
                       type="segmented",var.type="bootstrap",#bootstrap
                       chngpt.init =0.7,tol = 1e-7,maxit = 100,
                       est.method = 'fastgrid2')#fastgrid2
    summary.chngpt=summary(fit.chngpt)$chngpt 
    tau_chngpt[i]=summary.chngpt[1]
    sd_chngpt[i]=summary.chngpt[2]
    cov.chngpt = cov.chngpt + as.numeric(summary.chngpt[3] <= tau & summary.chngpt[4]>=tau)
    
    ### segmented method ###
    o1<-lm(Y~X)
    o.seg=segmented(o1,seg.Z =~X, psi=0.7,control=seg.control(n.boot=0, tol=1e-7, it.max = 100, display=F))
    seg.out=o.seg$psi[2];seg.sd=o.seg$psi[3]
    tau_seg[i]=seg.out;sd_seg[i]=seg.sd
    ci=confint(o.seg)
    cov.seg = cov.seg + as.numeric(ci[2] <= tau & ci[3]>=tau)
    
    
    ## Proposed Method ##
    tau.out = ToSNR_continuous(0.7)
    sd.out= sd_continuous(tau.out)
    sd_est1[i] = sd.out
    tau_est1[i] = tau.out
    cov1 = cov1 + as.numeric(tau.out-1.96*sd.out <= tau & tau.out+1.96*sd.out>=tau)
    
    ### NR_semismooth method ###
    tau.out2=NR_semismooth(0.7)
    if (!is.na( tau.out2))
    {
      sd.out2= sd_continuous(tau.out2)
      sd_est2[i.2] = sd.out2
      tau_est2[i.2] = tau.out2
      cov2 = cov2 + as.numeric(tau.out2-1.96*sd.out2 <= tau & tau.out2+1.96*sd.out2>=tau)
      i.2=i.2+1
    }
    
    
    i=i+1
    
  }
  
  print(paste0("sample size: ", n))
  #chngpt method
  bias_chngpt1 = mean(tau_chngpt) - tau
  sd_chngpt1 = sd(tau_chngpt)
  mse_chngpt1=sqrt(bias_chngpt1^2+sd_chngpt1^2)
  mean_sd_chngpt = mean(sd_chngpt)
  print("Fong Method")
  print(c(round(c(bias_chngpt1,sd_chngpt1,mse_chngpt1,mean_sd_chngpt)*10,3),round(cov.chngpt/(i-1)*100,1)))
  
  
  #segmented method
  bias_seg1 = mean(tau_seg) - tau
  sd_seg1 = sd(tau_seg)
  mse_seg1=sqrt(bias_seg1^2+sd_seg1^2)
  mean_sd_seg = mean(sd_seg)
  print("Muggeo Method")
  print(c(round(c(bias_seg1,sd_seg1,mse_seg1,mean_sd_seg)*10,3),round(cov.seg/(i-1)*100,1)))
  
  
  #To-SNR method
  bias_tau1 = mean(tau_est1) - tau
  sd_tau1 = sd(tau_est1)
  mse_tau1=sqrt(bias_tau1^2+sd_tau1 ^2)
  mean_sd1 = mean(sd_est1)
  print("To-SNR Method")
  print(c(round(c(bias_tau1,sd_tau1,mse_tau1,mean_sd1)*10,3),round(cov1/(i-1)*100,1)))
  
  #NR-semismooth method
  bias_tau2 = mean(tau_est2) - tau
  sd_tau2 = sd(tau_est2)
  mse_tau2=sqrt(bias_tau2^2+sd_tau2 ^2)
  mean_sd2 = mean(sd_est2)
  print("NR-Semi Method")
  print(c(round(c(bias_tau2,sd_tau2,mse_tau2,mean_sd2)*10,3),round(cov2/(i-1)*100,1)))
  
  cat("\n\n")
}
