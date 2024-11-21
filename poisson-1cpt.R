############################################################################################
###########  Simulation Codes For Poisson Segmented Model With One Change Point  ###########
############################################################################################
#install.packages("segmented")

#set.seed(2023106)
library(segmented)

###################################################
#############   Functions  ########################
###################################################
ToSNR_poisson<-function(tau0)
{
  erro = FALSE;eta = 1; j = 1
  while (any(abs(eta) > 1e-7) & j<=100 & erro==FALSE)
  {
    ### Step 1. Only update coefficients  ###
    X_tau1 =  as.numeric(X > tau0) * ( X - tau0 )
    glm.inital <- glm(Y~X+ X_tau1,family="poisson",control = list(maxit = 25))
    theta0 =coef( glm.inital) 
    be.tau = theta0['X_tau1']
    mu0 = predict(glm.inital, type="response")
    I_tau1 = as.numeric(X > tau0)
    
    ### Step 2. Only update change point  ###
    Un= mean((mu0-Y)*I_tau1)
    Jn = be.tau[1]*mean(mu0*I_tau1)
    tryCatch( {eta = Un/Jn},error=function(cond) {erro <<-  TRUE})
    if(erro ==TRUE| any(is.na(eta))) {eta=0}
    if (erro == FALSE){tau0 = tau0 + eta} ### Update change point ###
    j=j+1}
  
  tau.out =NA
  Test = (erro == FALSE && !is.na(eta))
  if ( Test) { tau.out = tau0}
  
  return(tau.out)
}

sd_poisson<- function(tau0_1)
{
  erro=FALSE
  ### generate inital values ###
  X1_tau1 =  as.numeric(X > tau0_1) * ( X - tau0_1 )
  I1_tau1 = as.numeric(X > tau0_1)
  
  glm.inital <- glm(Y~X+ X1_tau1,family="poisson",control = list(maxit = 25))
  theta0 =coef( glm.inital) 
  be.tau = theta0['X1_tau1']
  mu0 = predict(glm.inital, type="response")
  
  Jn1 =   c(mean(mu0),         mean(mu0*X),         mean(mu0*X1_tau1),        
            mean(mu0*X),       mean(mu0*X*X),      mean(mu0*X1_tau1*X),       
            mean(mu0*X1_tau1),  mean(mu0*X*X1_tau1), mean(mu0*X1_tau1*X1_tau1))
  Jn1 = matrix( Jn1 , ncol=3,nrow=3)
  Jn2=-c(be.tau*mean(mu0*I1_tau1),                              
         be.tau*mean(mu0*X*I1_tau1),                         
         be.tau*mean(mu0*X1_tau1))                      
  Jn2 = matrix( Jn2 , ncol=1,nrow=3,byrow = T)
  Jn3 = t(Jn2)
  Jn4 =be.tau^2*mean(mu0*I1_tau1);Jn4 = matrix(Jn4,ncol=1,nrow=1)
  sd.out =NA
  tryCatch( { V = solve(Jn4 - Jn3%*%solve(Jn1)%*%Jn2)/(n)},error=function(cond) {erro <<- TRUE})
  if (erro == FALSE){ sd.out = sqrt(V[1,1])}
  
  return(sd.out)
}
###################################################
#########        END OF FUNCTIONS    ##############
###################################################




##########################################################################
#####                    Data Generating Procedure               #########
##########################################################################
#number of MC replicates
Nsim = 1000
#beta=2.5
beta=1.8
tau=0.75

#sample size: 250, 500 and 1000
nlst=c(250,500,1000)
length.nlst=length(nlst)
for (nn in 1:length.nlst)
{
  #sample size: n
  n=nlst[nn]
  
  i=1;tau_est1 = NULL;tau_seg=NULL;sd_est1 = NULL; sd_seg=NULL;
  cov1 = 0; cov.seg=0;
  ####### Simulation Procedure ###########
  for(i in 1:Nsim)
  {
    ### generate dataset ###
    X = runif(n,0,1)
    t=3.5- 1.5*X + beta* as.numeric(X > tau) * ( X - tau ) 
    mu = exp(t)
    Y = rpois(n, lambda=mu)
    dat=data.frame(X,Y)
    
    ### segmented method ###
    o1<-glm(Y~X,family = "poisson")
    o.seg=segmented(o1,seg.Z =~X, psi=0.7,control=seg.control(n.boot=0, tol=1e-7, it.max = 100, display=F))
    seg.out=o.seg$psi[2];seg.sd=o.seg$psi[3]
    tau_seg[i]=seg.out;sd_seg[i]=seg.sd
    ci=confint(o.seg)
    cov.seg = cov.seg + as.numeric(ci[2] <= tau & ci[3]>=tau)
    
    
    ## Proposed Method ##
    tau.out = ToSNR_poisson(0.7)
    sd.out= sd_poisson(tau.out)
    sd_est1[i] = sd.out
    tau_est1[i] = tau.out
    cov1 = cov1 + as.numeric(tau.out-1.96*sd.out <= tau & tau.out+1.96*sd.out>=tau)
  
    
  }
  
  print(paste0("sample size: ", n))
  #segmented method
  bias_seg1 = mean(tau_seg) - tau
  sd_seg1 = sd(tau_seg)
  mse_seg1=sqrt(bias_seg1^2+sd_seg1^2)
  mean_sd_seg = mean(sd_seg)
  print("Muggeo Method")
  print(c(round(c(bias_seg1,sd_seg1,mse_seg1,mean_sd_seg)*10,3),round(cov.seg/i*100,1)))
  
  
  #To-SNR method
  bias_tau1 = mean(tau_est1) - tau
  sd_tau1 = sd(tau_est1)
  mse_tau1=sqrt(bias_tau1^2+sd_tau1 ^2)
  mean_sd1 = mean(sd_est1)
  print("To-SNR Method")
  print(c(round(c(bias_tau1,sd_tau1,mse_tau1,mean_sd1)*10,3),round(cov1/i*100,1)))

  
  cat("\n\n")
}
