rm(list=ls())

load("../Updated_resulttable_rstvar_gcm_empirical_IIV_final_0703_GCM-RS2.Rdata")

N = 108
O = 224
parameter <-c("beta00","sd_beta0","beta10",#"sd_beta1",
                  "AR0","sd_AR",
                  "gamma0_2","sd_gamma_2",#"shift",
                  "alpha","IIV_0",
                  #"betaIIV00",#"sd_betaIIV0",#"betaIIV10","sd_betaIIV1",
                  "midx",
                  "deviance")

npar=length(parameter)

bias=function(x,truex){
  mean((x-truex), na.rm = T)
}

relbias = function(x,truex){
  mean((x-truex)/truex, na.rm = T)
}

rmse = function(x,truex){
  sqrt(mean((x-truex)^2, na.rm = T))
}

power = function(lower, upper,isRound=0,digit){
  
  nrep = length(lower)
  
  if (isRound==0){
    includeFlag = ifelse(lower<= 0 & upper >= 0,1,0) #Does the interval include 0? If so, return a 1, else return a 0
  } else{
    includeFlag = ifelse(round(lower,digit)<= 0 & round(upper,digit) >= 0,1,0)}	
  power = 1-includeFlag
  
  power = length(power[power==1])/nrep
  if(sum(is.na(lower))>0){power=NA}
  return(power)
}	


coverage = function(lower, upper, truex, isRound=0,digit){
  
  nrep = length(lower)
  
  if (isRound==0){
    includeFlag2 = ifelse(lower <= truex & upper >= truex,1,0) #Does the interval include truex? If so, return a 1, else return a 0
  } else{
    includeFlag2 = ifelse(round(lower,digit)<= truex & round(upper,digit) >= truex,1,0)}  
  
  coverage = sum(includeFlag2,na.rm=T)/nrep*100
  return(coverage)
}

rope = function(lower, upper, isRound=0, digit){
  nrep = length(lower)
  if (isRound==0){
    includeFlag3 = ifelse(lower<= 0.05 & lower>=0,1,0) #Does the interval include (0,0.05)? If so, return a 1, else return a 0
  } else{
    includeFlag3 = ifelse(round(lower,digit)<= 0.05 & round(lower,digit)>= 0,1,0)}	
  rope = 1-includeFlag3
  
  rope = length(rope[rope==1])/nrep
  return(rope)
}



Deviance=resulttable[grep("deviance",rownames(resulttable)), "mean"]
# Get entropy
pr2 = resulttable[grep("^midx", rownames(resulttable)), "mean"] - 1
pr1 = 1-pr2
Entropy = 1 + (sum(pr2*log(pr2),na.rm=T)+sum(pr1*log(pr1),na.rm=T))/N/O/log(2)
sBIC = npar*log(N*O*(N*O+2)/24) + Deviance
BIC = npar*log(N*O) + Deviance
AIC = 2*npar + Deviance

# parameterlist <-c("beta00","sd_beta0","beta10",#"sd_beta1",
#                   "AR0","sd_AR",
#                   "gamma0","sd_gamma",#"shift",
#                   "alpha","IIV_0",
#                   #"betaIIV00",#"sd_betaIIV0",#"betaIIV10","sd_betaIIV1",
#                   "midx",
#                   "deviance")


# Get parameter estimates
temp = resulttable[c(grep("beta00",rownames(resulttable)),
                     grep("sd_beta0",rownames(resulttable)),
                     grep("beta10",rownames(resulttable)),
                     #grep("sd_beta1",rownames(resulttable)),
                     #grep("betaIIV10",rownames(resulttable)),
                     #grep("sd_betaIIV1",rownames(resulttable)),
                     grep("AR0",rownames(resulttable)),
                     grep("sd_AR",rownames(resulttable)),
                     grep("gamma0_2",rownames(resulttable)),
                     grep("sd_gamma_2",rownames(resulttable)),
                     grep("IIV_0",rownames(resulttable)),
                     #grep("^gamma0",rownames(resulttable)),
                     #grep("shift",rownames(resulttable)),
                     #grep("^sd_gamma",rownames(resulttable)),
                     #grep("betaIIV00",rownames(resulttable)),
                     #grep("sd_betaIIV0",rownames(resulttable)),
                     grep("^alpha\\[.*2,1\\]",rownames(resulttable)), #alpha[1,2,1],alpha[2,2,1]
                     grep("^alpha\\[.*1,2\\]",rownames(resulttable))), #alpha[1,1,2],alpha[2,1,2]
                   #c(1,4,5,6)
                   ]
temp
Entropy
sBIC
BIC
AIC