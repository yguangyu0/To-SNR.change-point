#################################################################################
########### To-SNR Algorithm: Continuous Outcome With One Change Point ##########
#################################################################################


#set.seed(2023106)
############# Functions For To-SNR Algorithm ####################
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
  sig=sigma(lm.inital)
  J =c(1,             mean(X),        mean(X_tau_1),   - be.tau*mean(I_tau_1),
       mean(X),       mean(X*X),      mean(X_tau_1*X), - be.tau*mean(X*I_tau_1),        
       mean(X_tau_1), mean(X*X_tau_1), mean(X_tau_1^2),- be.tau*mean(X_tau_1),
       -be.tau*mean(I_tau_1),-be.tau*mean(X*I_tau_1),- be.tau*mean(X_tau_1),   be.tau^2*mean(I_tau_1))
  J=matrix(J,ncol=4,nrow=4)
  sd.out =NA
  tryCatch( { V = solve(J)/(n)},error=function(cond) {erro <<- TRUE})
  if (erro == FALSE){ sd.out =  sig*sqrt(V[4,4])}
  
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
beta=1.8 #beta=2.5
tau=0.75

#sample size: 250, 500 and 1000
nlst=c(250,500,1000)
length.nlst=length(nlst)
for (nn in 1:length.nlst)
{
  n=nlst[nn]
  
  i=1; tau_est1 = NULL;sd_est1 = NULL; cov1 = 0;
  ####### Simulation Procedure ###########
  for(i in 1:Nsim)
  {
    ### generate dataset ###
    X = runif(n,0,1)
    e=rnorm(n,0,sd=0.2)
    Y=3.5- 1.5*X + beta* as.numeric(X > tau) * ( X - tau ) + e
    dat=data.frame(X,Y)
    
    ## Proposed Method ##
    tau.out = ToSNR_continuous(0.7)
    sd.out= sd_continuous(tau.out)
    sd_est1[i] = sd.out
    tau_est1[i] = tau.out
    cov1 = cov1 + as.numeric(tau.out-1.96*sd.out <= tau & tau.out+1.96*sd.out>=tau)


}
  
  print(paste0("sample size: ", n))
  #To-SNR method
  bias_tau1 = mean(tau_est1) - tau
  sd_tau1 = sd(tau_est1)
  mse_tau1=sqrt(bias_tau1^2+sd_tau1 ^2)
  mean_sd1 = mean(sd_est1)
  print("To-SNR Method")
  print(c(round(c(bias_tau1,sd_tau1,mse_tau1,mean_sd1)*10,3),round(cov1/i*100,1)))
  

  cat("\n\n")
}
