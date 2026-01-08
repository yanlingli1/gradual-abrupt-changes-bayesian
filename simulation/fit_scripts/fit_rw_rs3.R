rm(list=ls())
args=(commandArgs(TRUE))
if(length(args)==0){
  print("No arguments supplied.")
  ##supply default values
}else{
  for(i1 in 1:length(args)){
    eval(parse(text=args[[i1]]))
  }
}


library(rjags)
load.module("dic")
library(coda)
library(MCMCvis)
#library(rstudioapi)
#current_path <- getSourceEditorContext()$path
#setwd(dirname(current_path))
source("postcalc.R")

load(paste0("RSTVP_gcm_simulateddata_shift",shift,"_N",N,"O",O,"_",r1,".Rdata"))

jags_data <- list(Y = Y, 
                  P = dim(Y)[1],
                  maxT = dim(Y)[2],
                  Period2 = Period2,
                  time = time, 
                  nrObs = nrObs,
                  nrPeriod = length(time),
                  G3 = G3
                  )

fixedinits<- list(list(.RNG.seed=1,.RNG.name="base::Mersenne-Twister"),
                  list(.RNG.seed=2,.RNG.name="base::Mersenne-Twister"))

jagsModel <- jags.model(file = "model_rw_rs3.txt", data = jags_data, inits = fixedinits, 
                        n.chains = 2, n.adapt = 4000) 
update(jagsModel, n.iter = 1000)

parameterlist <- c("sd_int","beta00","sd_beta0",
                   "gamma0","sd_gamma","alpha","AR","sd_noise",
                   #"gamma","midx","mu_int","delta","int",
                   "midx","deviance")

codaSamples <- coda.samples(jagsModel, 
                            variable.names = parameterlist,
                            n.iter = 20000)
resulttable <- zcalc(codaSamples)


save(codaSamples, resulttable, file=paste0("Result_RSTVP_rw_3regimes_shift",shift,"_N",N,"O",O,"_",r1,".Rdata"))




