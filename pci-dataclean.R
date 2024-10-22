##########################################
#### Code For Application (PCI Data) #####
##########################################
library(haven)
library(sas7bdat)
library(boot)
library(arm)
library(pROC)


###################################################################
###  To-SNR algorithm: Estimated Change Points For BMI & GFR  #####
###################################################################
#To-SNR algorithm:
NR_logistic_Z<-function(tau0)
{
  erro = FALSE;eta = rep(1,3); j = 1
  while (any(abs(eta) > 1e-4) & j<=100 & erro==FALSE)
  {
    tau0_1 = tau0[1];tau0_2 = tau0[2];tau0_3 = tau0[3]
    X1_tau1 =  as.numeric(X1 > tau0_1) * ( X1 - tau0_1 )
    X2_tau1 =  as.numeric(X2 > tau0_2) * ( X2 - tau0_2 )
    X2_tau2 =  as.numeric(X2 > tau0_3) * ( X2 - tau0_3 )
    
    I1_tau1 = as.numeric(X1 > tau0_1)
    I2_tau1 = as.numeric(X2 > tau0_2)
    I2_tau2 = as.numeric(X2 > tau0_3)
    
    glm.inital <- glm(Y~X1+X1_tau1+X2+X2_tau1+X2_tau2+Z1+Z2+Z3+Z4+Z5+Z6+Z7+Z8+Z9+Z10,
                      family="binomial", control = list(maxit = 25))
    theta0 =coef( glm.inital) 
    be.tau = c(theta0['X1_tau1'],theta0['X2_tau1'],theta0['X2_tau2'])
    mu = predict(glm.inital, type="response")
    mu0 = mu*(1-mu)
    
    ### Calculate U and J ###
    Un =  c(be.tau[1]*mean((mu-Y)*I1_tau1),
            be.tau[2]*mean((mu-Y)*I2_tau1),
            be.tau[3]*mean((mu-Y)*I2_tau2))
    Un = matrix(Un,nrow=3)
    
    Jn =c(be.tau[1]^2*mean(mu0*I1_tau1*I2_tau1),         be.tau[1]*be.tau[2]*mean(mu0*I1_tau1*I2_tau1), be.tau[1]*be.tau[3]*mean(mu0*I1_tau1*I2_tau2),
          be.tau[1]*be.tau[2]*mean(mu0*I1_tau1*I2_tau1), be.tau[2]^2*mean(mu0*I2_tau1*I2_tau1),         be.tau[2]*be.tau[3]*mean(mu0*I2_tau2*I2_tau1),
          be.tau[1]*be.tau[3]*mean(mu0*I1_tau1*I2_tau2), be.tau[2]*be.tau[3]*mean(mu0*I2_tau2*I2_tau1), be.tau[3]*be.tau[3]*mean(mu0*I2_tau2*I2_tau2))
    
    Jn = matrix(Jn,ncol=3,nrow=3)
    
    tryCatch( {eta = solve(Jn)%*%Un},error=function(cond) {erro <<-  TRUE})
    if(erro ==TRUE| any(is.na(eta))) {eta=rep(0,3)}
    if (erro == FALSE)
    {### Update Tau ###
      tau_up = tau0 + eta
      tau0 = tau_up}
    j=j+1
  }
  
  tau.out = c(NA,NA,NA)
  if ( erro == FALSE )
  { tau.out = tau0}
  
  return(tau.out)
}

#standard error of estimated change points:
sd_logistic_Z<- function(tau.out)
{
  erro=FALSE
  
  tau0_1 = tau.out[1];tau0_2 = tau.out[2];tau0_3 = tau.out[3]
  X1_tau1 =  as.numeric(X1 > tau0_1) * ( X1 - tau0_1 )
  X2_tau1 =  as.numeric(X2 > tau0_2) * ( X2 - tau0_2 )
  X2_tau2 =  as.numeric(X2 > tau0_3) * ( X2 - tau0_3 )
  
  I1_tau1 = as.numeric(X1 > tau0_1)
  I2_tau1 = as.numeric(X2 > tau0_2)
  I2_tau2 = as.numeric(X2 > tau0_3)
  
  glm.inital <- glm(Y~X1+X1_tau1+X2+X2_tau1+X2_tau2+Z1+Z2+Z3+Z4+Z5+Z6+Z7+Z8+Z9+Z10,
                    family="binomial", control = list(maxit = 25))
  theta0 =coef( glm.inital) 
  be.tau = c(theta0['X1_tau1'],theta0['X2_tau1'],theta0['X2_tau2'])
  mu = predict(glm.inital, type="response")
  mu0 = mu*(1-mu)
  
  n.tau = length(tau.out)
  
  Z102 = ifelse(Z10==2, 1, 0)
  Z103 = ifelse(Z10==3, 1, 0)
  Jn1 =  c(mean(mu0),         mean(mu0*X1),          mean(mu0*X1_tau1),         mean(mu0*X2),         mean(mu0*X2_tau1),             mean(mu0*X2_tau2),             mean(mu0*Z1),         mean(mu0*Z2),         mean(mu0*Z3),         mean(mu0*Z4),        mean(mu0*Z5),        mean(mu0*Z6),         mean(mu0*Z7),        mean(mu0*Z8),         mean(mu0*Z9),        mean(mu0*Z102),        mean(mu0*Z103),             
           mean(mu0*X1),      mean(mu0*X1^2),        mean(mu0*X1*X1_tau1),      mean(mu0*X1*X2),      mean(mu0*X1*X2_tau1),          mean(mu0*X1*X2_tau2),          mean(mu0*X1*Z1),      mean(mu0*X1*Z2),      mean(mu0*X1*Z3),      mean(mu0*X1*Z4),     mean(mu0*X1*Z5),     mean(mu0*X1*Z6),      mean(mu0*X1*Z7),     mean(mu0*X1*Z8),      mean(mu0*X1*Z9),     mean(mu0*X1*Z102),     mean(mu0*X1*Z103),
           mean(mu0*X1_tau1), mean(mu0*X1*X1_tau1),  mean(mu0*X1_tau1^2),       mean(mu0*X1_tau1*X2), mean(mu0*X1_tau1*X2_tau1),     mean(mu0*X1_tau1*X2_tau2),     mean(mu0*X1_tau1*Z1), mean(mu0*X1_tau1*Z2), mean(mu0*X1_tau1*Z3), mean(mu0*X1_tau1*Z4),mean(mu0*X1_tau1*Z5),mean(mu0*X1_tau1*Z6), mean(mu0*X1_tau1*Z7),mean(mu0*X1_tau1*Z8), mean(mu0*X1_tau1*Z9),mean(mu0*X1_tau1*Z102),mean(mu0*X1_tau1*Z103),
           mean(mu0*X2),      mean(mu0*X1*X2),       mean(mu0*X2*X1_tau1),      mean(mu0*X2*X2),      mean(mu0*X2*X2_tau1),          mean(mu0*X2*X2_tau2),          mean(mu0*X2*Z1),      mean(mu0*X2*Z2),      mean(mu0*X2*Z3),      mean(mu0*X2*Z4),     mean(mu0*X2*Z5),     mean(mu0*X2*Z6),      mean(mu0*X2*Z7),     mean(mu0*X2*Z8),      mean(mu0*X2*Z9),     mean(mu0*X2*Z102),     mean(mu0*X2*Z103),
           mean(mu0*X2_tau1), mean(mu0*X1*X2_tau1),  mean(mu0*X1_tau1*X2_tau1), mean(mu0*X2_tau1*X2), mean(mu0*X2_tau1*X2_tau1),     mean(mu0*X2_tau1*X2_tau2),     mean(mu0*X2_tau1*Z1), mean(mu0*X2_tau1*Z2), mean(mu0*X2_tau1*Z3), mean(mu0*X2_tau1*Z4),mean(mu0*X2_tau1*Z5),mean(mu0*X2_tau1*Z6), mean(mu0*X2_tau1*Z7),mean(mu0*X2_tau1*Z8), mean(mu0*X2_tau1*Z9),mean(mu0*X2_tau1*Z102),mean(mu0*X2_tau1*Z103),
           mean(mu0*X2_tau2), mean(mu0*X1*X2_tau2),  mean(mu0*X1_tau1*X2_tau2), mean(mu0*X2*X2_tau2), mean(mu0*X2_tau1*X2_tau2),     mean(mu0*X2_tau2*X2_tau2),     mean(mu0*X2_tau2*Z1), mean(mu0*X2_tau2*Z2), mean(mu0*X2_tau2*Z3), mean(mu0*X2_tau2*Z4),mean(mu0*X2_tau2*Z5),mean(mu0*X2_tau2*Z6), mean(mu0*X2_tau2*Z7),mean(mu0*X2_tau2*Z8), mean(mu0*X2_tau2*Z9),mean(mu0*X2_tau2*Z102),mean(mu0*X2_tau2*Z103),
           mean(mu0*Z1),     mean(mu0*X1*Z1),        mean(mu0*X1_tau1*Z1),      mean(mu0*X2*Z1),      mean(mu0*X2_tau1*Z1),          mean(mu0*X2_tau2*Z1),          mean(mu0*Z1^2),       mean(mu0*Z1*Z2),      mean(mu0*Z1*Z3),      mean(mu0*Z1*Z4),     mean(mu0*Z1*Z5),     mean(mu0*Z1*Z6),      mean(mu0*Z1*Z7),     mean(mu0*Z1*Z8),      mean(mu0*Z1*Z9),     mean(mu0*Z1*Z102),     mean(mu0*Z1*Z103),
           mean(mu0*Z2),     mean(mu0*X1*Z2),        mean(mu0*X1_tau1*Z2),      mean(mu0*X2*Z2),      mean(mu0*X2_tau1*Z2),          mean(mu0*X2_tau2*Z2),          mean(mu0*Z1*Z2),      mean(mu0*Z2^2),       mean(mu0*Z2*Z3),      mean(mu0*Z2*Z4),     mean(mu0*Z2*Z5),     mean(mu0*Z2*Z6),      mean(mu0*Z2*Z7),     mean(mu0*Z2*Z8),      mean(mu0*Z2*Z9),     mean(mu0*Z2*Z102),     mean(mu0*Z2*Z103),
           mean(mu0*Z3),     mean(mu0*X1*Z3),        mean(mu0*X1_tau1*Z3),      mean(mu0*X2*Z3),      mean(mu0*X2_tau1*Z3),          mean(mu0*X2_tau2*Z3),          mean(mu0*Z1*Z3),      mean(mu0*Z2*Z3),      mean(mu0*Z3^2),       mean(mu0*Z3*Z4),     mean(mu0*Z3*Z5),     mean(mu0*Z3*Z6),      mean(mu0*Z3*Z7),     mean(mu0*Z3*Z8),      mean(mu0*Z3*Z9),     mean(mu0*Z3*Z102),     mean(mu0*Z3*Z103),
           mean(mu0*Z4),     mean(mu0*X1*Z4),        mean(mu0*X1_tau1*Z4),      mean(mu0*X2*Z4),      mean(mu0*X2_tau1*Z4),          mean(mu0*X2_tau2*Z4),          mean(mu0*Z1*Z4),      mean(mu0*Z2*Z4),      mean(mu0*Z3*Z4),      mean(mu0*Z4^2),      mean(mu0*Z4*Z5),     mean(mu0*Z4*Z6),      mean(mu0*Z4*Z7),     mean(mu0*Z4*Z8),      mean(mu0*Z4*Z9),     mean(mu0*Z4*Z102),     mean(mu0*Z4*Z103),
           mean(mu0*Z5),     mean(mu0*X1*Z5),        mean(mu0*X1_tau1*Z5),      mean(mu0*X2*Z5),      mean(mu0*X2_tau1*Z5),          mean(mu0*X2_tau2*Z5),          mean(mu0*Z1*Z5),      mean(mu0*Z2*Z5),      mean(mu0*Z3*Z5),      mean(mu0*Z4*Z5),     mean(mu0*Z5^2),      mean(mu0*Z5*Z6),      mean(mu0*Z5*Z7),     mean(mu0*Z5*Z8),      mean(mu0*Z5*Z9),     mean(mu0*Z5*Z102),     mean(mu0*Z5*Z103),
           mean(mu0*Z6),     mean(mu0*X1*Z6),        mean(mu0*X1_tau1*Z6),      mean(mu0*X2*Z6),      mean(mu0*X2_tau1*Z6),          mean(mu0*X2_tau2*Z6),          mean(mu0*Z1*Z6),      mean(mu0*Z2*Z6),      mean(mu0*Z3*Z6),      mean(mu0*Z4*Z6),     mean(mu0*Z5*Z6),     mean(mu0*Z6*Z6),      mean(mu0*Z6*Z7),     mean(mu0*Z6*Z8),      mean(mu0*Z6*Z9),     mean(mu0*Z6*Z102),     mean(mu0*Z6*Z103),
           mean(mu0*Z7),     mean(mu0*X1*Z7),        mean(mu0*X1_tau1*Z7),      mean(mu0*X2*Z7),      mean(mu0*X2_tau1*Z7),          mean(mu0*X2_tau2*Z7),          mean(mu0*Z1*Z7),      mean(mu0*Z2*Z7),      mean(mu0*Z3*Z7),      mean(mu0*Z4*Z7),     mean(mu0*Z5*Z7),     mean(mu0*Z7*Z6),      mean(mu0*Z7*Z7),     mean(mu0*Z7*Z8),      mean(mu0*Z7*Z9),     mean(mu0*Z7*Z102),     mean(mu0*Z7*Z103),
           mean(mu0*Z8),     mean(mu0*X1*Z8),        mean(mu0*X1_tau1*Z8),      mean(mu0*X2*Z8),      mean(mu0*X2_tau1*Z8),          mean(mu0*X2_tau2*Z8),          mean(mu0*Z1*Z8),      mean(mu0*Z2*Z8),      mean(mu0*Z3*Z8),      mean(mu0*Z4*Z8),     mean(mu0*Z5*Z8),     mean(mu0*Z8*Z6),      mean(mu0*Z8*Z7),     mean(mu0*Z8*Z8),      mean(mu0*Z8*Z9),     mean(mu0*Z8*Z102),     mean(mu0*Z8*Z103),
           mean(mu0*Z9),     mean(mu0*X1*Z9),        mean(mu0*X1_tau1*Z9),      mean(mu0*X2*Z9),      mean(mu0*X2_tau1*Z9),          mean(mu0*X2_tau2*Z9),          mean(mu0*Z1*Z9),      mean(mu0*Z2*Z9),      mean(mu0*Z3*Z9),      mean(mu0*Z4*Z9),     mean(mu0*Z5*Z9),     mean(mu0*Z9*Z6),      mean(mu0*Z9*Z7),     mean(mu0*Z9*Z8),      mean(mu0*Z9*Z9),     mean(mu0*Z9*Z102),     mean(mu0*Z9*Z103),
           mean(mu0*Z102),   mean(mu0*X1*Z102),      mean(mu0*X1_tau1*Z102),    mean(mu0*X2*Z102),    mean(mu0*X2_tau1*Z102),        mean(mu0*X2_tau2*Z102),        mean(mu0*Z1*Z102),    mean(mu0*Z2*Z102),    mean(mu0*Z3*Z102),    mean(mu0*Z4*Z102),   mean(mu0*Z5*Z102),   mean(mu0*Z102*Z6),    mean(mu0*Z102*Z7),   mean(mu0*Z102*Z8),    mean(mu0*Z102*Z9),   mean(mu0*Z102*Z102),   mean(mu0*Z102*Z103),
           mean(mu0*Z103),   mean(mu0*X1*Z103),      mean(mu0*X1_tau1*Z103),    mean(mu0*X2*Z103),    mean(mu0*X2_tau1*Z103),        mean(mu0*X2_tau2*Z103),        mean(mu0*Z1*Z103),    mean(mu0*Z2*Z103),    mean(mu0*Z3*Z103),    mean(mu0*Z4*Z103),   mean(mu0*Z5*Z103),   mean(mu0*Z103*Z6),    mean(mu0*Z103*Z7),   mean(mu0*Z103*Z8),    mean(mu0*Z103*Z9),   mean(mu0*Z103*Z102),   mean(mu0*Z103*Z103))
  n.Jn1 = sqrt(length(Jn1))
  Jn1 = matrix( Jn1 , ncol=n.Jn1,nrow=n.Jn1)
  
  Jn2 = - c(be.tau[1]*mean(mu0*I1_tau1),             be.tau[2]*mean(mu0*I2_tau1),             be.tau[3]*mean(mu0*I2_tau2),
            be.tau[1]*mean(mu0*X1*I1_tau1),          be.tau[2]*mean(mu0*X1*I2_tau1),          be.tau[3]*mean(mu0*X1*I2_tau2),
            be.tau[1]*mean(mu0*X1_tau1),             be.tau[2]*mean(mu0*X1_tau1*I2_tau1),     be.tau[3]*mean(mu0*X1_tau1*I2_tau2),
            be.tau[1]*mean(mu0*X2*I1_tau1),          be.tau[2]*mean(mu0*X2*I2_tau1),          be.tau[3]*mean(mu0*X2*I2_tau2),
            be.tau[1]*mean(mu0*X2_tau1*I1_tau1),     be.tau[2]*mean(mu0*X2_tau1),             be.tau[3]*mean(mu0*X2_tau1*I2_tau2),
            be.tau[1]*mean(mu0*X2_tau2*I1_tau1),     be.tau[2]*mean(mu0*X2_tau2*I2_tau1),     be.tau[3]*mean(mu0*X2_tau2),
            be.tau[1]*mean(mu0*Z1*I1_tau1),          be.tau[2]*mean(mu0*Z1*I2_tau1),          be.tau[3]*mean(mu0*Z1*I2_tau2),
            be.tau[1]*mean(mu0*Z2*I1_tau1),          be.tau[2]*mean(mu0*Z2*I2_tau1),          be.tau[3]*mean(mu0*Z2*I2_tau2),
            be.tau[1]*mean(mu0*Z3*I1_tau1),          be.tau[2]*mean(mu0*Z3*I2_tau1),          be.tau[3]*mean(mu0*Z3*I2_tau2),
            be.tau[1]*mean(mu0*Z4*I1_tau1),          be.tau[2]*mean(mu0*Z4*I2_tau1),          be.tau[3]*mean(mu0*Z4*I2_tau2),
            be.tau[1]*mean(mu0*Z5*I1_tau1),          be.tau[2]*mean(mu0*Z5*I2_tau1),          be.tau[3]*mean(mu0*Z5*I2_tau2),
            be.tau[1]*mean(mu0*Z6*I1_tau1),          be.tau[2]*mean(mu0*Z6*I2_tau1),          be.tau[3]*mean(mu0*Z6*I2_tau2),
            be.tau[1]*mean(mu0*Z7*I1_tau1),          be.tau[2]*mean(mu0*Z7*I2_tau1),          be.tau[3]*mean(mu0*Z7*I2_tau2),
            be.tau[1]*mean(mu0*Z8*I1_tau1),          be.tau[2]*mean(mu0*Z8*I2_tau1),          be.tau[3]*mean(mu0*Z8*I2_tau2),
            be.tau[1]*mean(mu0*Z9*I1_tau1),          be.tau[2]*mean(mu0*Z9*I2_tau1),          be.tau[3]*mean(mu0*Z9*I2_tau2),
            be.tau[1]*mean(mu0*Z102*I1_tau1),        be.tau[2]*mean(mu0*Z102*I2_tau1),        be.tau[3]*mean(mu0*Z102*I2_tau2),
            be.tau[1]*mean(mu0*Z103*I1_tau1),        be.tau[2]*mean(mu0*Z103*I2_tau1),        be.tau[3]*mean(mu0*Z103*I2_tau2)
  )
  Jn2 = matrix( Jn2 , ncol=n.tau,nrow=n.Jn1,byrow = T)
  Jn3 = t(Jn2)
  
  Jn4 =  c(be.tau[1]^2*mean(mu0*I1_tau1*I2_tau1),         be.tau[1]*be.tau[2]*mean(mu0*I1_tau1*I2_tau1), be.tau[1]*be.tau[3]*mean(mu0*I1_tau1*I2_tau2),
           be.tau[1]*be.tau[2]*mean(mu0*I1_tau1*I2_tau1), be.tau[2]^2*mean(mu0*I2_tau1*I2_tau1),         be.tau[2]*be.tau[3]*mean(mu0*I2_tau2*I2_tau1),
           be.tau[1]*be.tau[3]*mean(mu0*I1_tau1*I2_tau2), be.tau[2]*be.tau[3]*mean(mu0*I2_tau2*I2_tau1), be.tau[3]*be.tau[3]*mean(mu0*I2_tau2*I2_tau2))
  Jn4 = matrix(Jn4,ncol=n.tau,nrow=n.tau)
  
  tryCatch( {J_inverse = solve(Jn4 - Jn3%*%solve(Jn1)%*%Jn2)},error=function(cond) {erro <<-  TRUE})
  
  sd.out = c(NA,NA)
  if (erro==FALSE)
  {V = J_inverse/(n)
  
  sd_1 = sqrt(V[1,1])
  sd_2 = sqrt(V[2,2])
  sd_3 = sqrt(V[3,3])
  
  sd.out = c(sd_1, sd_2, sd_3)}
  
  return(sd.out)
}




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



###############################################################
######## Estimated Change Point via To-SNR Algorithm ##########
###############################################################
#initial values of change points 
tau0=c(30,60,105)


out = NR_logistic_Z(tau0)#estimated change points
sd = sd_logistic_Z(out) #standard error 

#calculate 95% CIs
low.95 = out -1.96*sd
upp.95 = out +1.96*sd

#estimated change point of BMI
#change point, standard error, 95%CI (low,upper)
round(c(out[1],sd[1],low.95[1], upp.95[1]),1)

#estimated change points of GFR
#change point, standard error, 95%CI (low,upper)
round(c(out[2],sd[2],low.95[2], upp.95[2]),1)
round(c(out[3],sd[3],low.95[3], upp.95[3]),1)

##Final Fitted Model With Estimated Change Points
X1_tau1 = as.numeric(X1 > out[1]) * (X1 -out[1])
X2_tau1 = as.numeric(X2 > out[2]) * (X2 - out[2])
X2_tau2 = as.numeric(X2 > out[3]) * (X2 - out[3])
lin.function<-glm(Y~X1+X1_tau1+X2+X2_tau1+X2_tau2+Z1+Z2+Z3+Z4+Z5+Z6+Z7+Z8+Z9+Z10,family=binomial)

#measure of goodness of fit: c-statistics
Y.pred = predict(lin.function,type="response")
round(auc(Y,Y.pred),2)#0.78

#summary of final fitted model
sum.lin<-summary(lin.function)


#slopes of BMI/GFR before & after change points
coef.lin = sum.lin$coefficients
vcov.lin = vcov(lin.function)
#############################
###### For X1: BMI ##########
#############################
#before change point: odds ratio & p-value
be1 = coef.lin['X1','Estimate']
se_b1 = coef.lin['X1','Std. Error']
p.be1 = coef.lin['X1','Pr(>|z|)']
print("odds ratio before change point (BMI)")
round(exp(be1),2)
print("p-value before change point (BMI)")
round(p.be1,3)
#after change point: odds ratio & p-value
be2 = coef.lin['X1','Estimate']+coef.lin['X1_tau1','Estimate']
se_b2 = sqrt(vcov.lin['X1','X1']+vcov.lin['X1_tau1','X1_tau1']+2*vcov.lin['X1','X1_tau1'])
Z_b2 = be2/se_b2
p.be2= 2*(1-pnorm(abs(Z_b2)))
print("odds ratio after change point (BMI)")
round(exp(be2),2)
print("p-value after change point (BMI)")
round(p.be2,3)
#############################
###### For X2: GFR ##########
#############################
#before change point 1: odds ratio & p-value
be1 = coef.lin['X2','Estimate']
se_b1 = coef.lin['X2','Std. Error']
p.be1 = coef.lin['X2','Pr(>|z|)']
print("odds ratio before change point 1 (GFR)")
round(exp(be1),2)
print("p-value before change point 1 (GFR)")
round(p.be1,3)
#between change point 1 and change point 2: odds ratio & p-value
be2 = coef.lin['X2','Estimate']+coef.lin['X2_tau1','Estimate']
se_b2 = sqrt(vcov.lin['X2','X2']+vcov.lin['X2_tau1','X2_tau1']+2*vcov.lin['X2','X2_tau1'])
Z_b2 = be2/se_b2
p.be2= 2*(1-pnorm(abs(Z_b2)))
round(c(be2,exp(be2), se_b2,p.be2),4)
print("odds ratio between two change points (GFR)")
round(exp(be2),3)
print("p-value between two change points (GFR)")
round(p.be2,3)
#after change point 2: odds ratio & p-value
be3 = coef.lin['X2','Estimate']+coef.lin['X2_tau1','Estimate']+coef.lin['X2_tau2','Estimate']
se_b3 = sqrt(vcov.lin['X2','X2']+vcov.lin['X2_tau1','X2_tau1']+vcov.lin['X2_tau2','X2_tau2']+
               2*vcov.lin['X2','X2_tau1']+2*vcov.lin['X2','X2_tau2']+2*vcov.lin['X2_tau2','X2_tau1'])
Z_b3 = be3/se_b3
p.be3= 2*(1-pnorm(abs(Z_b3)))
print("odds ratio after change point 2 (GFR)")
round(exp(be3),2)
print("p-value after change point 2 (GFR)")
round(p.be3,3)
