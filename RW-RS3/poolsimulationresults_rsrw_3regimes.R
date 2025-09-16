#rm(list=ls())
#args=(commandArgs(TRUE))
#if(length(args)==0){
#  print("No arguments supplied.")
#  ##supply default values
#}else{
#  for(i1 in 1:length(args)){
#    eval(parse(text=args[[i1]]))
#  }
#}


N=100
O=200
shift=0.6

library(plotROC)
library(MCMCvis)
#setwd("/gpfs/group/quc16/default/yxl823/Dissertation")


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

load(paste0("RSTVP_gcm_simulateddata_shift",shift,"_N",N,"O",O,"_1.Rdata"))

truevalue = c(beta00,sd_beta0,gamma0[2],NA,sd_gamma[2],NA,
              AR,sd_noise,NA,
              alpha[1,2,1],alpha[1,1,2],NA,NA,NA,NA)

parameter = c("beta00", "sd_beta0","gamma0[2]", "gamma0[3]", "sd_gamma[2]", "sd_gamma[3]", 
              "AR", "sd_noise","sd_int",
              "alpha21", "alpha12","alpha31","alpha13","alpha32","alpha23")

###################################################
#load(paste0("Result_RSTVP_gcm_shift",shift,"_N",N,"O",O,"_",r1,".Rdata"))
K = length(truevalue)
npar = K
myfiles = list.files(pattern=paste0("^Result_RSTVP_rw_3regimes_shift",shift,"_N",N,"O",O,"_.*\\.Rdata"),full.names = TRUE)
Nrep = length(myfiles)
myfiles_new = c()
for (i in 1:Nrep){
  load(myfiles[i])
  if(max(resulttable$RHAT,na.rm=T) < 1.2){
    myfiles_new = c(myfiles_new, myfiles[i])
  }
}

nrep = length(myfiles_new)
poolresult = array(data = rep(NA,K*nrep*4),dim = c(K,nrep,4))

AIC=rep(NA,nrep)
BIC=rep(NA,nrep)
sBIC=rep(NA,nrep)
Entropy=rep(NA,nrep)
ACC = rep(NA,nrep)
Recall = rep(NA,nrep)
Precision = rep(NA,nrep)
AUC = rep(NA,nrep)
Bias = rep(NA,nrep)
RMSE = rep(NA,nrep)

for (i in 1:nrep){
    load(myfiles_new[i])

    # Get IC measures
    Deviance=resulttable[grep("deviance",rownames(resulttable)), "mean"]
    #npar=xxx
    sBIC[i] = npar*log(N*O*(N*O+2)/24) + Deviance
    BIC[i] = npar*log(N*O) + Deviance
    AIC[i] = 2*npar + Deviance
    
    # Get entropy
    #midx_pos = data.frame(MCMCchains(codaSamples,params="midx"))
niter=1000
pr = array(NA, dim=c(3,N,O))
Midx = array(NA, dim=c(niter,N,O))
for(pp in 1:N){
  for(oo in 1:O){
    #Midx[,pp,oo] = as.numeric(midx_pos[1:niter,paste0("midx.",pp,".",oo,".")])
    Midx[,pp,oo] = MCMCchains(codaSamples,params=paste0("midx[",pp,",",oo,"]"),ISB = FALSE, exact = TRUE)[1:niter]
    prtable = table(Midx[,pp,oo])
    prtable = data.frame(prtable)
    prtable$Var1 = as.numeric(prtable$Var1)
    prtable2 = data.frame(matrix(NA,ncol=2,nrow=3))
    colnames(prtable2)=c("Var1","Freq")
    prtable2$Var1 = 1:3
    if(1 %in% prtable$Var1){
      prtable2$Freq[1] = prtable$Freq[prtable$Var1==1]
    }
    if(2 %in% prtable$Var1){
      prtable2$Freq[2] = prtable$Freq[prtable$Var1==2]
    }
    if(3 %in% prtable$Var1){
      prtable2$Freq[3] = prtable$Freq[prtable$Var1==3]
    }
    pr[1,pp,oo] = prtable2$Freq[1]/niter
    pr[2,pp,oo] = prtable2$Freq[2]/niter
    pr[3,pp,oo] = prtable2$Freq[3]/niter
  }}

    Entropy[i] = 1+sum(pr*log(pr),na.rm=T)/N/O/log(3)

    
    load(gsub("Result_RSTVP_rw_3regimes", "RSTVP_gcm_simulateddata", myfiles_new[i]))
    
    # check MAE and RMSE
    #Bias[i] = abs(int_est - int)/N/O
    #RMSE[i] = sqrt((int_est - int)^2/N/O)
    
    # check AUC
    Regime = as.vector(midx) - 1
    pr2 = 1 - as.vector(pr[1,,])    
    df = data.frame(cbind(Regime,pr2))
    rocplot = ggplot(df, aes(m = pr2, d = Regime))+ geom_roc(n.cuts=10,labels=TRUE)
    AUC[i] = calc_auc(rocplot)$AUC

    # Get parameter estimates
      temp = resulttable[c(grep("beta00",rownames(resulttable)),
                         grep("sd_beta0",rownames(resulttable)),
                         grep("^gamma0",rownames(resulttable)),
                         grep("^sd_gamma",rownames(resulttable)),
                         grep("AR",rownames(resulttable)),
                         grep("sd_noise",rownames(resulttable)),
			 grep("sd_int",rownames(resulttable)),
                        grep("^alpha\\[.*2,1\\]",rownames(resulttable)),
                        grep("^alpha\\[.*1,2\\]",rownames(resulttable)),
                        grep("^alpha\\[.*3,1\\]",rownames(resulttable)),
                        grep("^alpha\\[.*1,3\\]",rownames(resulttable)),
                        grep("^alpha\\[.*3,2\\]",rownames(resulttable)),
                        grep("^alpha\\[.*2,3\\]",rownames(resulttable))), c(1,4,5,6)]
    for(k in 1:K){
    for(j in 1:4){
      poolresult[k,i,j] = temp[k,j]
    }
  }
}


print(paste0("AIC: ", "Mean ", round(mean(AIC,na.rm=T),2), " SD ", round(sd(AIC,na.rm=T),2), " low ", round(as.numeric(quantile(AIC, 0.025, na.rm=T)),2)," high ", round(as.numeric(quantile(AIC, 0.975, na.rm=T)),2)))
print(paste0("BIC: ", "Mean ", round(mean(BIC,na.rm=T),2), " SD ", round(sd(BIC,na.rm=T),2), " low ", round(as.numeric(quantile(BIC, 0.025, na.rm=T)),2)," high ", round(as.numeric(quantile(BIC, 0.975, na.rm=T)),2)))
print(paste0("sBIC: ", "Mean ", round(mean(sBIC,na.rm=T),2), " SD ", round(sd(sBIC,na.rm=T),2), " low ", round(as.numeric(quantile(sBIC, 0.025, na.rm=T)),2)," high ", round(as.numeric(quantile(sBIC, 0.975, na.rm=T)),2)))
print(paste0("Entropy: ", "Mean ", round(mean(Entropy,na.rm=T),2), " SD ", round(sd(Entropy,na.rm=T),2), " low ", round(as.numeric(quantile(Entropy, 0.025, na.rm=T)),2)," high ", round(as.numeric(quantile(Entropy, 0.975, na.rm=T)),2)))
#print(paste0("Accuracy: ", "Mean ", round(mean(ACC,na.rm=T),2), " SD ", round(sd(ACC,na.rm=T),2), " low ", round(as.numeric(quantile(ACC, 0.025, na.rm=T)),2)," high ", round(as.numeric(quantile(ACC, 0.975, na.rm=T)),2)))
#print(paste0("Recall: ", "Mean ", round(mean(Recall,na.rm=T),2), " SD ", round(sd(Recall,na.rm=T),2), " low ", round(as.numeric(quantile(Recall, 0.025, na.rm=T)),2)," high ", round(as.numeric(quantile(Recall, 0.975, na.rm=T)),2)))
#print(paste0("Precision: ", "Mean ", round(mean(Precision,na.rm=T),2), " SD ", round(sd(Precision,na.rm=T),2), " low ", round(as.numeric(quantile(Precision, 0.025, na.rm=T)),2)," high ", round(as.numeric(quantile(Precision, 0.975, na.rm=T)),2)))
print(paste0("AUC: ", "Mean ", round(mean(AUC,na.rm=T),2), " SD ", round(sd(AUC,na.rm=T),2), " low ", round(as.numeric(quantile(AUC, 0.025, na.rm=T)),2)," high ", round(as.numeric(quantile(AUC, 0.975, na.rm=T)),2)))
#print(paste0("Bias: ", "Mean ", round(mean(Bias,na.rm=T),2), " SD ", round(sd(Bias,na.rm=T),2), " low ", round(as.numeric(quantile(Bias, 0.025, na.rm=T)),2)," high ", round(as.numeric(quantile(Bias, 0.975, na.rm=T)),2)))
#print(paste0("RMSE: ", "Mean ", round(mean(RMSE,na.rm=T),2), " SD ", round(sd(RMSE,na.rm=T),2), " low ", round(as.numeric(quantile(RMSE, 0.025, na.rm=T)),2)," high ", round(as.numeric(quantile(RMSE, 0.975, na.rm=T)),2)))



#save(ACC, Recall, Precision, AUC, Bias, RMSE, file=paste0("SimulationResult_REZIPX_",condition,".Rdata"))


#pdf("checkdistributions.pdf")
#for(k in 1:K){
#hist(poolresult[k,,1],main=parameter[k])
#}
#dev.off()

Est = c()
Bias = c()
Relbias = c()
SE = c()
SE_MC = c()
RMSE = c()
Power = c()
Coverage = c()
ROPE = c()

for(k in 1:K){
  Est[k] = mean(poolresult[k,,1])
     Bias[k] = bias(poolresult[k,,1], truevalue[k])
    Relbias[k] = relbias(poolresult[k,,1], truevalue[k])
    SE_MC[k] = sd(poolresult[k,,1])
    SE[k] = mean(poolresult[k,,2])
  RMSE[k] = rmse(poolresult[k,,1], truevalue[k])
  Power[k] = power(poolresult[k,,3], poolresult[k,,4], isRound=0,digit)
  Coverage[k] = coverage(poolresult[k,,3], poolresult[k,,4], truevalue[k],isRound=0,digit)
   ROPE[k] = rope(poolresult[k,,3], poolresult[k,,4], isRound=0,digit)

 }

sumdat=data.frame(cbind(truevalue, Est, Bias, Relbias, SE, SE_MC, RMSE, Power, ROPE, Coverage))
rownames(sumdat) = parameter
colnames(sumdat) = c("True","Est","Bias","RBias","SE","MCSE","RMSE","Power","ROPE","CR(%)")
sumdat[,2:(ncol(sumdat)-1)] = round(sumdat[,2:(ncol(sumdat)-1)], 2)
sumdat[,ncol(sumdat)] = round(sumdat[,ncol(sumdat)], 0)

write.csv(sumdat, file=paste0("SummaryMeasures_RSTVP_rw_3regimes_shift",shift,"_N",N,"O",O,".csv"))
sumdat
nrep








  
    



  



