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

# Check if R_LIBRARY_LOAD_PATH exists in the environment
library_path <- Sys.getenv("R_LIBRARY_LOAD_PATH")

# If the environment variable is not set, use the default path
if (library_path == "") {
  library_path <- "/storage/work/xjx5093/.conda/envs/jagstest/lib/R/library"
}

# Set the library path
.libPaths(library_path)
# Check the library path
.libPaths()

library(rjags)
load.module("dic")
library(coda)
#library(rstudioapi)
#current_path <- getSourceEditorContext()$path
#setwd(dirname(current_path))
source("postcalc.R")

data_path = "/storage/work/xjx5093/GOHIAR_RS/DataGeneration_SparseTransition/"
load(paste0(data_path, "RSTVP_gcm_simulateddata_shift",shift,"_N",N,"O",O,"_",r1,".Rdata"))

jags_data <- list(Y = Y, 
                  P = dim(Y)[1],
                  maxT = dim(Y)[2],
                  #Period2 = Period2,
                  time = time, 
                  nrObs = nrObs,
                  nrPeriod = length(time),
                  X = X
                  )

fixedinits<- list(list(.RNG.seed=1,.RNG.name="base::Mersenne-Twister"),
                  list(.RNG.seed=2,.RNG.name="base::Mersenne-Twister"))

jagsModel <- jags.model(file = "model_gcm_rs2_iiv_sparse.txt", data = jags_data, inits = fixedinits, n.chains = 2, n.adapt = 4000) 
update(jagsModel, n.iter = 2000)

parameterlist <- c("beta00","sd_beta0","beta10","sd_beta1",
                   "AR0","sd_AR",
                   "gamma0","shift","sd_gamma",
                   "alpha",
                   "midx","deviance")

codaSamples <- coda.samples(jagsModel, 
                            variable.names = parameterlist,
                            n.iter = 25000)
resulttable <- zcalc(codaSamples)


save(resulttable, file=paste0("Result_RSTVP_gcm_shift",shift,"_N",N,"O",O,"_",r1,".Rdata"))

