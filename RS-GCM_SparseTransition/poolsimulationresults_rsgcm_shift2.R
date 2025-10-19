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
shift=2

library(plotROC)

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


data_path = "/storage/work/xjx5093/GOHIAR_RS/DataGeneration_SparseTransition/"
load(paste0(data_path, "RSTVP_gcm_simulateddata_shift",shift,"_N",N,"O",O,"_1.Rdata"))

truevalue = c(beta00,sd_beta0,beta10,sd_beta1,
              AR0,sd_AR, # gamma0[1], sd_gamma[1],
              gamma0[2], shift, sd_gamma[2],
              alpha[1,2,1],alpha[2,2,1],alpha[1,1,2],alpha[2,1,2])

parameter = c("beta00", "sd_beta0", "beta10", "sd_beta1", 
              "AR0", "sd_AR", #"gamma0[1]","sd_gamma[1]",
              "gamma0[2]", "shift", "sd_gamma[2]", 
              "alpha121", "alpha221", "alpha112", "alpha212")

###################################################
#load(paste0("Result_RSTVP_gcm_shift",shift,"_N",N,"O",O,"_",r1,".Rdata"))
K = length(truevalue)
npar = K
myfiles = list.files(pattern=paste0("^Result_RSTVP_gcm_shift",shift,"_N",N,"O",O,"_.*\\.Rdata"),full.names = TRUE)
Nrep = length(myfiles)
myfiles_new = c()
for (i in 1:Nrep){
  load(myfiles[i])
  if(max(resulttable$RHAT,na.rm=T) < 1.05){
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
    pr2 = resulttable[grep("^midx", rownames(resulttable)), "mean"] - 1
    pr1 = 1-pr2
    Entropy[i] = 1 + (sum(pr2*log(pr2),na.rm=T)+sum(pr1*log(pr1),na.rm=T))/N/O/log(2)
    
    
    load(paste0(data_path, gsub("Result_RSTVP_gcm", "RSTVP_gcm_simulateddata", myfiles_new[i])))
    
    # check MAE and RMSE
    #Bias[i] = abs(int_est - int)/N/O
    #RMSE[i] = sqrt((int_est - int)^2/N/O)
    
    # check ACC, Recall, Precision, and AUC
    Midx = resulttable[grep("^midx", rownames(resulttable)), "median"] - 1
    Regime = as.vector(midx) - 1    
    diff = data.frame(cbind(Regime,Midx))
    diff$TP = ifelse(diff$Regime==1&diff$Midx==1,1,0)
    diff$FP = ifelse(diff$Regime==0&diff$Midx==1,1,0)
    diff$TN = ifelse(diff$Regime==0&diff$Midx==0,1,0)
    diff$FN = ifelse(diff$Regime==1&diff$Midx==0,1,0)
    
    ACC[i] = (sum(diff$TP)+sum(diff$TN))/N/O
    Recall[i] = sum(diff$TP)/(sum(diff$TP)+sum(diff$FN))
    Precision[i] = sum(diff$TP)/(sum(diff$TP)+sum(diff$FP))
    df = data.frame(cbind(Regime,pr2))
    rocplot = ggplot(df, aes(m = pr2, d = Regime))+ geom_roc(n.cuts=10,labels=TRUE)
    AUC[i] = calc_auc(rocplot)$AUC
    
    # Get parameter estimates
    temp = resulttable[c(grep("beta00",rownames(resulttable)),
                         grep("sd_beta0",rownames(resulttable)),
                         grep("beta10",rownames(resulttable)),
                         grep("sd_beta1",rownames(resulttable)),
                         grep("AR0",rownames(resulttable)),
                         grep("sd_AR",rownames(resulttable)),
                         grep("^gamma0",rownames(resulttable))[2],
                         grep("shift",rownames(resulttable)),
                         grep("^sd_gamma",rownames(resulttable))[2],
                        grep("^alpha\\[.*2,1\\]",rownames(resulttable)), #alpha[1,2,1],alpha[2,2,1]
                        grep("^alpha\\[.*1,2\\]",rownames(resulttable))), #alpha[1,1,2],alpha[2,1,2]
                       c(1,4,5,6)]
    for(k in 1:K){
    for(j in 1:4){
      poolresult[k,i,j] = temp[k,j]
    }
  }
}


# 1) Put your metrics into a named list
metrics <- list(
  AIC       = AIC,
  BIC       = BIC,
  sBIC      = sBIC,
  Entropy   = Entropy,
  Accuracy  = ACC,
  Recall    = Recall,
  Precision = Precision,
  AUC       = AUC
)

# 2) Build the summary table
summary_table <- data.frame(
  Metric = names(metrics),
  Mean   = sapply(metrics, function(x) round(mean(x,    na.rm=TRUE), 2)),
  SD     = sapply(metrics, function(x) round(sd(x,      na.rm=TRUE), 2)),
  Low    = sapply(metrics, function(x) round(quantile(x, .025, na.rm=TRUE), 2)),
  High   = sapply(metrics, function(x) round(quantile(x, .975, na.rm=TRUE), 2)),
  row.names = NULL,
  stringsAsFactors = FALSE
)

# 3) Print to console
print(summary_table)

# 4) Write to a csv file
write.csv(
  summary_table,
  file = paste0("ModelSelectionMetrics_RSTVP_gcm_shift",shift,"_N",N,"O",O,".csv")
)


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

write.csv(sumdat, file=paste0("SummaryMeasures_RSTVP_gcm_shift",shift,"_N",N,"O",O,".csv"))
sumdat
print(paste("Number of replications:", Nrep))
print(paste("Number of converged replications:", nrep))
#print(myfiles_new)







  
    



  



