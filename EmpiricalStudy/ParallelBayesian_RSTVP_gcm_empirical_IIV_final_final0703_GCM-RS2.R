rm(list=ls())
#library(rjags)
library(R2jags)

load.module("dic")
library(coda)
#library(rstudioapi)
#current_path <- getSourceEditorContext()$path
#setwd(dirname(current_path))
#setwd("D:/XXYDATAanalysis/GP_Yanling")
source("postcalc.R")

set.seed(123)
load("./data/datblock_Meaning.Rdata")
load("./data/datblock_relationship.Rdata")
load("./data/datblock_accomplishment.Rdata")
load("./data/datCov.Rdata")

datblock_relationship$Relationship_scale=scale(datblock_relationship$Relationship)
colnames(datblock_relationship)[ncol(datblock_relationship)]="Relationship_scale"
datblock_accomplishment$Accomplishment_scale=scale(datblock_accomplishment$Accomplishment)
colnames(datblock_accomplishment)[ncol(datblock_accomplishment)]="Accomplishment_scale"

O = 250
timepoint = 0:3
Y = matrix(NA, N, O)
Cov = array(NA, dim=c(2, N, O))
G3 = rep(NA,N)
Gender = rep(NA,N)
Age = rep(NA,N)
maxT = rep(NA,N)
Period2 = matrix(NA,N,O)
Period3 = matrix(NA,N,O)
Period4 = matrix(NA,N,O)
time = matrix(NA,N,O)
nrObs = matrix(NA,N,5)
for(pp in 1:N){
  tmp = datblock[datblock$PID==pInd[pp],]
  tmp_relationship = datblock_relationship[datblock_relationship$PID==pInd[pp],]
  tmp_accomplishment = datblock_accomplishment[datblock_accomplishment$PID==pInd[pp],]
  maxT[pp] = nrow(tmp)
  time[pp,1:O] = c(rep(0, 56),rep(1, 60),rep(2, 52),rep(3, O-56-60-52))
  nrObs[pp,1] = 1
  nrObs[pp,2] = 56
  nrObs[pp,3] = 56+60
  nrObs[pp,4] = 56+60+52  
  nrObs[pp,5] = maxT[pp]
  group = unique(dat$Group[dat$PID == pInd[pp]])
  G3[pp] = ifelse(group==3,1,0)
  Gender[pp] = datCov[datCov$PID==pInd[pp],"Gender"]
  Age[pp] = datCov[datCov$PID==pInd[pp],"Age"] - 3
  for(oo in 1:maxT[pp]){
    Y[pp,oo] = tmp$Meaning[oo]
    Cov[1,pp,oo] = tmp_relationship$Relationship_scale[oo]
    Cov[2,pp,oo] = tmp_accomplishment$Accomplishment_scale[oo]# pre-intervention
    Period2[pp,oo] = ifelse(tmp$DayOfStudy[oo] > 14 & tmp$DayOfStudy[oo] <= 29, 1, 0)#intervention period
    Period3[pp,oo] = ifelse(tmp$DayOfStudy[oo] > 29 & tmp$DayOfStudy[oo] <= 42, 1, 0)
    Period4[pp,oo] = ifelse(tmp$DayOfStudy[oo] > 42, 1, 0)
  }
}


######################## Model Fitting


jags_data <- list(Y = Y, 
                  P = dim(Y)[1],
                  maxT = maxT,
                  Period2 = Period2,
                  #Period3 = Period3,
                  #Period4 = Period4,
                  time = time,
                  timepoint = timepoint, 
                  G3 = G3,
                  #Gender = Gender,
                  #Age = Age,
                  #Cov = Cov,
                  nrPeriod = 4,
                  nrObs = nrObs
                  )
parameterlist <-c("beta00","sd_beta0","beta10",#"sd_beta1",
                  "AR0","sd_AR",
                  "gamma0_2","sd_gamma_2",#"shift",
                  "alpha","IIV_0",
                  #"betaIIV00",#"sd_betaIIV0",#"betaIIV10","sd_betaIIV1",
                  "midx",
                  "deviance")
# ---- Using rjags ----
#jagsModel <- jags.model(file = "Cleanedrsonly_empirical_IIV_final0703_GCM-RS2.txt", 
#                        data = jags_data, 
#                        n.chains = 2, n.adapt = 4000) 
#update(jagsModel, n.iter = 2000)

#codaSamples <- coda.samples(jagsModel, 
#                            variable.names = parameterlist,
#                            n.iter = 30000)
#resulttable <- zcalc(codaSamples)


# ---- using R2jags to run in parallel ----
fit <- jags.parallel(
  data = jags_data,
  parameters.to.save = parameterlist,
  model.file = "Cleanedrsonly_empirical_IIV_final0703_GCM-RS2.txt",
  n.chains = 2,
  n.iter = 30000,
  n.burnin = 10000, #Previously 2000
  DIC = TRUE
)

mcmc_samples <- fit$BUGSoutput$sims.list
str(mcmc_samples)

# Or: get an `mcmc.list` object for direct use with coda
mcmc_list <- as.mcmc(fit)
summary(mcmc_list)
plot(mcmc_list)
resulttable <- zcalc(mcmc_list)


save(resulttable, file="Updated_resulttable_rstvar_gcm_empirical_IIV_final_0703_GCM-RS2.Rdata")




