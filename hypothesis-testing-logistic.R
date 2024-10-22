####################################################################################################
########### Simulation Codes For Hypothesis Testing of Change Points in a Logistic Model  ##########
####################################################################################################
#install.packages("segmented")
#install.packages("chngpt")

library(segmented)
library(chngpt)
set.seed(20231006)

#number of MC replicates
Nsim=10000;
#sample size
n=1000

M = 10;o=0
p.01=0;p.05=0;p.001=0
rho=0.3;c=rho*1.6
LR.pvalue=NULL;S.pvalue=NULL;D.pvalue=NULL;P.pvalue=NULL

while(o<Nsim)
{
  XZ=mvrnorm(n,mu=c(0,4.7),Sigma=matrix(c(1,c,c,1.6^2),ncol=2))
  Z=XZ[,1];X=XZ[,2]
  t=0.34*Z+0.4*X
  pr = 1/(1+exp(-t)) 
  Y = rbinom(n,1,pr)
  dat=data.frame(Y,X,Z)
  o1<-glm(Y~0+X+Z,family="binomial")
  
  # Davies' test (segmented package)
  davies_test=davies.test(o1,seg.Z =~X, alternative = c("two.sided"))
  D.pvalue=c(D.pvalue,davies_test$p.value)
  #Score & LR test (chngpt package)
  chngpt_test.S=chngpt.test(formula.null=Y~Z, formula.chngpt=~X,family="binomial", 
                          data=dat, type="segmented", test.statistic="score",p.val.method="MC")
  S.pvalue=c(S.pvalue,chngpt_test.S$p.value)
  
  chngpt_test.LR=chngpt.test(formula.null=Y~Z, formula.chngpt=~X,family="binomial", 
                          data=dat, type="segmented", test.statistic="lr",p.val.method="MC")
  LR.pvalue=c(LR.pvalue,chngpt_test.LR$p.value)
  
  #proposed average score test
  mu1 = predict(o1, type="response")
  mu0 = mu1*(1-mu1)
  X.matrix=matrix(c(X,Z),ncol=2)
  phi_value=quantile(X,probs = c(1:M)/(M+1))
  phi_function<-function(tau){return(as.numeric(X > tau) * ( X - tau))}
  phi_result=as.numeric(unlist(lapply(phi_value, phi_function)))
  aa=matrix(phi_result,ncol=M,nrow=n,byrow = F)
  phi_mean=matrix(rowMeans(aa),ncol=1)
  D=diag(mu0)
  l0=sum((Y-mu1)*phi_mean)
  V1=t(phi_mean)%*%D%*%X.matrix%*%solve(t(X.matrix)%*%D%*%X.matrix)%*%t(X.matrix)%*%D%*%phi_mean
  var_l0 =  t(phi_mean)%*%D%*%phi_mean-V1
  s0=l0/sqrt(var_l0)
  if(abs(s0)>qnorm(0.95,0,1)){p.01=p.01+1}
  if(abs(s0)>qnorm(0.975,0,1)){p.05=p.05+1}
  if(abs(s0)>qnorm(0.995,0,1)){p.001=p.001+1}
  
  o=o+1
  #if(o%%100==0){print(o)} #count simulation time
}

#Level 0.01
#Davies test; LR test; Score test; proposed test
c(length(which(D.pvalue<0.01)),length(which(LR.pvalue<0.01)),length(which(S.pvalue<0.01)),p.001)/o

#Level 0.05
#Davies test; LR test; Score test; proposed test
c(length(which(D.pvalue<0.05)),length(which(LR.pvalue<0.05)),length(which(S.pvalue<0.05)),p.05)/o

#Level 0.10
#Davies test; LR test; Score test; proposed test
c(length(which(D.pvalue<0.10)),length(which(LR.pvalue<0.10)),length(which(D.pvalue<0.10)),p.01)/o


