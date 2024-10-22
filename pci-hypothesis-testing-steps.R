########## HYPOTHESIS TESTING OF EXISTENCE CHANGE POINTS ############
####################    REAL APPLICATION.   #########################
#install.packages(c("haven", "sas7bdat"))
#### Application PCI Data #####
library(haven)
library(sas7bdat)


###################################################################
#################### Read In PCI Data #############################
###################################################################
dat<- read_sas("/Users/yguangyu/Desktop/Theoretic\ Work/GLSM/GLSM2024/2024_application/pcidata.sas7bdat")

## select variables in PCI data ##
variable = c('IH_GI','IH_TX','GFR','Age','gender','HX_HYP','combine_dm',
             'HX_SMK','hx_ptca','hx_cabg','race_white','race_black',
             'hx_mi','CI_PRIOR','bmi','wt_kg','ht_cm')
dat<-dat[,variable]
dat<-dat[complete.cases(dat),]

#only include subjects with height >100cm & height <= 250cm and CI_PRIOR!=0
wt=dat$wt_kg
ht=dat$ht_cm
index.ht = which(ht<=250 & ht>=100 & dat$CI_PRIOR!=0)#63156 observation
ht.2 = ht[index.ht]
wt.2 = wt[index.ht]
dat<-dat[index.ht,]

#define outcome Y, factor of interest with change point effects X1, X2
Y = ifelse(dat$IH_GI==0 & dat$IH_TX==0,0,1)#GI Bleed or Tranfustion (binary outcome)
X1 = dat$bmi  #BMI
X2 = dat$GFR  #Glomerular Filtration Rate
#other covariates
Z1 = dat$Age
Z2 = dat$gender
Z3 = dat$race_white
Z4 = dat$race_black
Z5 = dat$HX_SMK #Smoking Status
Z6 = dat$HX_HYP #History of Hypertension
Z7 = dat$combine_dm #History of Diabetes
Z8 = dat$hx_ptca #Previous PCI(percutaneous coronary intervention)
Z9 = dat$hx_cabg #Previous CABG
Z10 = dat$CI_PRIOR #priority: emergent/urgent/nonurgent
Z10 = as.factor(Z10)

n=length(Y) #sample size


### STEP 1. TEST WHETHER ONE CHANGE POINT IN BMI ###
#M: default at 10 
pvalue.step1<-function(M)
{
  o.0<-glm(Y~X1+X2+Z1+Z2+Z3+Z4+Z5+Z6+Z7+Z8+Z9+Z10,family="binomial")
  #test statistics
  mu1 = predict(o.0, type="response")
  mu0 = mu1*(1-mu1)
  X.matrix=matrix(c(rep(1,n),X1,X2,Z1,Z2,Z3,Z4,Z5,Z6,Z7,Z8,Z9,Z10),ncol=13)
  phi_value=quantile(X1,probs = c(1:M)/(M+1))
  phi_function<-function(tau){return(as.numeric(X1 > tau) * ( X1 - tau))}
  phi_result=as.numeric(unlist(lapply(phi_value, phi_function)))
  aa=matrix(phi_result,ncol=M,nrow=length(Y),byrow = F)
  phi_mean=matrix(rowMeans(aa),ncol=1)
  D=diag(mu0)
  l0=sum((Y-mu1)*phi_mean)
  V1=t(phi_mean)%*%D%*%X.matrix%*%solve(t(X.matrix)%*%D%*%X.matrix)%*%t(X.matrix)%*%D%*%phi_mean
  var_l0 =  t(phi_mean)%*%D%*%phi_mean-V1
  s0=l0/sqrt(var_l0)
  
  return((1-pnorm(abs(s0),0,1))*2)
}

pvalue.step1(10)
#p-value<0.0001
  



### STEP 2. TEST WHETHER A SECOND CHANGE POINT IN BMI ###
pvalue.step2<-function(M)
{
  X1_tau1 =  as.numeric(X1 > 32.3) * ( X1 -32.3)
  o.1<-glm(Y~X1+X1_tau1+X2+Z1+Z2+Z3+Z4+Z5+Z6+Z7+Z8+Z9+Z10,family="binomial")
  #test statistics
  mu1 = predict(o.1, type="response")
  mu0 = mu1*(1-mu1)
  X.matrix=matrix(c(rep(1,n),X1,X1_tau1,X2,Z1,Z2,Z3,Z4,Z5,Z6,Z7,Z8,Z9,Z10),ncol=14)
  phi_value=quantile(X1,probs = c(1:M)/(M+1))
  phi_function<-function(tau){return(as.numeric(X1 > tau) * ( X1 - tau))}
  phi_result=as.numeric(unlist(lapply(phi_value, phi_function)))
  aa=matrix(phi_result,ncol=M,nrow=length(Y),byrow = F)
  phi_mean=matrix(rowMeans(aa),ncol=1)
  D=diag(mu0)
  l0=sum((Y-mu1)*phi_mean)
  V1=t(phi_mean)%*%D%*%X.matrix%*%solve(t(X.matrix)%*%D%*%X.matrix)%*%t(X.matrix)%*%D%*%phi_mean
  var_l0 =  t(phi_mean)%*%D%*%phi_mean-V1
  s0=l0/sqrt(var_l0)
  
  return((1-pnorm(abs(s0),0,1))*2)
}
pvalue.step2(10)
#p-value: 0.398146

  
  
### STEP 3. TEST WHETHER ONE CHANGE POINT IN GFR, GIVEN ONE CHANGE POINT IN BMI ###
pvalue.step3<-function(M)
{
  X1_tau1 =  as.numeric(X1 > 32.3) * ( X1 -32.3)
  o.1<-glm(Y~X1+X1_tau1+X2+Z1+Z2+Z3+Z4+Z5+Z6+Z7+Z8+Z9+Z10,family="binomial")
  #test statistics
  mu1 = predict(o.1, type="response")
  mu0 = mu1*(1-mu1)
  X.matrix=matrix(c(rep(1,n),X1,X1_tau1,X2,Z1,Z2,Z3,Z4,Z5,Z6,Z7,Z8,Z9,Z10),ncol=14)
  phi_value=quantile(X2,probs = c(1:M)/(M+1))
  phi_function<-function(tau){return(as.numeric(X2 > tau) * ( X2 - tau))}
  phi_result=as.numeric(unlist(lapply(phi_value, phi_function)))
  aa=matrix(phi_result,ncol=M,nrow=length(Y),byrow = F)
  phi_mean=matrix(rowMeans(aa),ncol=1)
  D=diag(mu0)
  l0=sum((Y-mu1)*phi_mean)
  V1=t(phi_mean)%*%D%*%X.matrix%*%solve(t(X.matrix)%*%D%*%X.matrix)%*%t(X.matrix)%*%D%*%phi_mean
  var_l0 =  t(phi_mean)%*%D%*%phi_mean-V1
  s0=l0/sqrt(var_l0)
  
  return((1-pnorm(abs(s0),0,1))*2)
}
pvalue.step3(10)
#p-value: 4.8107e-06


  



### STEP 4. TEST A SECOND CHANGE POINT IN GFR, GIVEN ONE CHANGE POINT IN BMI ###
pvalue.step4<-function(M)
{
  X1_tau1 =  as.numeric(X1 > 32.3) * ( X1 -32.3)
  X2_tau1 = as.numeric(X2 > 54.2) * (X2 - 54.2)
  o.4<-glm(Y~X1+X1_tau1+X2+X2_tau1+Z1+Z2+Z3+Z4+Z5+Z6+Z7+Z8+Z9+Z10,family="binomial")
  #test statistics
  mu1 = predict(o.4, type="response")
  mu0 = mu1*(1-mu1)
  X.matrix=matrix(c(rep(1,n),X1,X1_tau1,X2,X2_tau1,Z1,Z2,Z3,Z4,Z5,Z6,Z7,Z8,Z9,Z10),ncol=15)
  phi_value=quantile(X2,probs = c(1:M)/(M+1))
  phi_function<-function(tau){return(as.numeric(X2 > tau) * ( X2 - tau))}
  phi_result=as.numeric(unlist(lapply(phi_value, phi_function)))
  aa=matrix(phi_result,ncol=M,nrow=length(Y),byrow = F)
  phi_mean=matrix(rowMeans(aa),ncol=1)
  D=diag(mu0)
  l0=sum((Y-mu1)*phi_mean)
  V1=t(phi_mean)%*%D%*%X.matrix%*%solve(t(X.matrix)%*%D%*%X.matrix)%*%t(X.matrix)%*%D%*%phi_mean
  var_l0 =  t(phi_mean)%*%D%*%phi_mean-V1
  s0=l0/sqrt(var_l0)
  
  return((1-pnorm(abs(s0),0,1))*2)
}
pvalue.step4(10)
#p-value: 0.01338537



### STEP 5. TEST WHETHER A THIRD CHANGE POINT IN GFR, GIVEN ONE CHANGE POINT IN BMI###
pvalue.step5<-function(M)
{
  X1_tau1 =  as.numeric(X1 > 32.3) * ( X1 -32.3)
  X2_tau1 = as.numeric(X2 > 54.2) * (X2 - 54.2)
  X2_tau2 = as.numeric(X2 > 114.7) * (X2 - 114.7)
  o.5<-glm(Y~X1+X1_tau1+X2+X2_tau1+X2_tau2+Z1+Z2+Z3+Z4+Z5+Z6+Z7+Z8+Z9+Z10,family="binomial")
  #test statistics
  mu1 = predict(o.5, type="response")
  mu0 = mu1*(1-mu1)
  X.matrix=matrix(c(rep(1,n),X1,X1_tau1,X2,X2_tau1,X2_tau2,Z1,Z2,Z3,Z4,Z5,Z6,Z7,Z8,Z9,Z10),ncol=16)
  phi_value=quantile(X2,probs = c(1:M)/(M+1))
  phi_function<-function(tau){return(as.numeric(X2 > tau) * ( X2 - tau))}
  phi_result=as.numeric(unlist(lapply(phi_value, phi_function)))
  aa=matrix(phi_result,ncol=M,nrow=length(Y),byrow = F)
  phi_mean=matrix(rowMeans(aa),ncol=1)
  D=diag(mu0)
  l0=sum((Y-mu1)*phi_mean)
  V1=t(phi_mean)%*%D%*%X.matrix%*%solve(t(X.matrix)%*%D%*%X.matrix)%*%t(X.matrix)%*%D%*%phi_mean
  var_l0 =  t(phi_mean)%*%D%*%phi_mean-V1
  s0=l0/sqrt(var_l0)
  
  return((1-pnorm(abs(s0),0,1))*2)
}
pvalue.step5(10)
#p-value: 0.4070974
