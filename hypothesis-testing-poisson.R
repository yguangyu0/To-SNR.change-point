####################################################################################################
########### Simulation Codes For Hypothesis Testing of Change Points in a Logistic Model  ##########
####################################################################################################

library(segmented)
set.seed(2023106)

#number of MC replicates
Nsim=10000
#sample size
n=500

o=0
D.pvalue=NULL;p.01=0;p.05=0;p.001=0
rho=0.3;c=rho*1.6
M=10
while(o<Nsim)
{
  X = runif(n,0,1)
  t=3.5- 1.5*X 
  mu = exp(t)
  Y = rpois(n, lambda=mu)
  
  o1<-glm(Y~X,family="poisson")
  
  # Davies' test (segmented package)
  davies_test=davies.test(o1,seg.Z =~X, alternative = c("two.sided"))
  D.pvalue=c(D.pvalue,davies_test$p.value)
  
  #proposed method
  mu0 = predict(o1, type="response")
  X.matrix=matrix(c(rep(1,n),X),ncol=2)
  phi_value=quantile(X,probs = c(1:M)/(M+1))
  phi_function<-function(tau){return(as.numeric(X > tau) * ( X - tau))}
  phi_result=as.numeric(unlist(lapply(phi_value, phi_function)))
  aa=matrix(phi_result,ncol=M,nrow=n,byrow = F)
  phi_mean=matrix(rowMeans(aa),ncol=1)
  D=diag(mu0)
  l0=sum((Y-mu0)*phi_mean)
  V1=t(phi_mean)%*%D%*%X.matrix%*%solve(t(X.matrix)%*%D%*%X.matrix)%*%t(X.matrix)%*%D%*%phi_mean
  var_l0 =  t(phi_mean)%*%D%*%phi_mean-V1
  s0=l0/sqrt(var_l0)
  if(abs(s0)>qnorm(0.95,0,1)){p.01=p.01+1}
  if(abs(s0)>qnorm(0.975,0,1)){p.05=p.05+1}
  if(abs(s0)>qnorm(0.995,0,1)){p.001=p.001+1}
  
  o=o+1
  #if(o%%1000==0){print(o)} #count simulation time
}


#Level 0.01
#Davies test; LR test; Score test; proposed test
c(length(which(D.pvalue<0.01)),p.001)/o

#Level 0.05
#Davies test; LR test; Score test; proposed test
c(length(which(D.pvalue<0.05)),p.05)/o

#Level 0.10
#Davies test; LR test; Score test; proposed test
c(length(which(D.pvalue<0.10)),p.01)/o