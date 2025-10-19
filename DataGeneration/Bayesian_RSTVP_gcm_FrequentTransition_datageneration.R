# "Frequent transition" condition

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

set.seed(r1)
#N = 100 # number of persons 
#O = 200 # number of time points
shift = as.numeric(shift)  # shift should be changed in the job submission script!
#G3 = c(rep(0, N/2),rep(1, N/2)) # two groups
time = 0:3
nrPeriod = length(time) # 4 periods
#Period2 = c(rep(0, O/4), rep(1, O/4), rep(0, O/4), rep(0, O/4))
nrObs=c(1, seq(O/4, O, O/4))
X = matrix(runif(N*O, -3, 3), nrow=N, ncol=O)

# parameters related to intercepts
beta00 = 0 # intercept in GCM
sd_beta0 = 0.2
beta10 = 0.5
sd_beta1 = 0.2


# parameters in the AR model
AR0 = 0.3  # AR parameter
sd_AR = 0.1

# parameters in the IIV model
#gamma0 = c(0, shift) # levels in two regimes #@@@ change it back, fix the first gamma01=0
gamma01 = 0 # fix gamma0[1] at 0
gamma0 = c(gamma01, gamma01 + shift)
sd_gamma = c(0, 0.2) # remove random effects for the shift level 
sd_gamma[1] = 0 #0.2
if(shift == 0){
  sd_gamma = c(0, 0)
}

# parameters in the RS model
# e.g., LO(S_{i,t}=r|S_{i,t-1}=s) = alpha[1,r,s]+alpha[2,r,s]*G3[pp]
alpha = array(rep(NA, 2*2*2), dim=c(2,2,2))
# fix alpha11 and alpha22 to 0
for(i in 1:2){
  alpha[i,1,1] = 0 
  alpha[i,2,2] = 0
}

#alpha[1,2,1] = -3 #empirical
#alpha[1,1,2] = -2 #empirical

#alpha[2,2,1] = .5 #empirical
#alpha[2,1,2] = -.5 #empirical

alpha[1,2,1] = -3  #Yanling
alpha[1,1,2] = -2  #Yanling

alpha[2,2,1] = 1  #Yanling
alpha[2,1,2] = -1  #Yanling

# simulate data

odds = array(rep(NA,N*O*2),dim = c(N,O,2))
midx = matrix(NA,N,O) # regime indicator
p_2s = matrix(NA,N,O) # probability of being in regime 2

Y=matrix(NA,N,O) # observed variable
gamma=matrix(NA,2,N) 
beta0=rep(NA,N)
beta1=rep(NA,N)
AR=rep(NA,N)
#intercept=matrix(NA,N,nrPeriod) 
intercept=matrix(NA,N,O) 
logIIV=matrix(NA,N,O) 
for(pp in 1:N){
  odds[pp,1,1] = 1
  odds[pp,1,2] = 0
  p_2s[pp,1] = odds[pp,1,2]/(odds[pp,1,1]+odds[pp,1,2])
  midx[pp,1] = sample.int(n=2,size=1,prob=odds[pp,1,])
  gamma[1,pp] = rnorm(1, gamma0[1], sd_gamma[1])
  gamma[2,pp] = rnorm(1, gamma0[2], sd_gamma[2])
  beta0[pp] = rnorm(1, beta00, sd_beta0)
  beta1[pp] = rnorm(1, beta10, sd_beta1)
  AR[pp] = rnorm(1, AR0, sd_AR)
  intercept[pp,1] = beta0[pp]
  logIIV[pp,1] = gamma[1,pp]
  Y[pp,1] = rnorm(1, intercept[pp,1], sqrt(exp(logIIV[pp,1])))
  
  for (bb in 1:nrPeriod){

   for(oo in (nrObs[bb]+1):nrObs[bb+1]){
    # from regime midx[pp,oo-1] to regime 1
    odds[pp,oo,1] = exp(alpha[1,1,midx[pp,oo-1]]+alpha[2,1,midx[pp,oo-1]]*X[pp,oo]) # remove Period2
    # from regime midx[pp,oo-1] to regime 2 
    odds[pp,oo,2] = exp(alpha[1,2,midx[pp,oo-1]]+alpha[2,2,midx[pp,oo-1]]*X[pp,oo]) # remove Period2
    midx[pp,oo]=sample.int(n=2,size=1,prob=odds[pp,oo,])
    p_2s[pp,oo] = odds[pp,oo,2]/(odds[pp,oo,1]+odds[pp,oo,2])

    intercept[pp,oo] = beta0[pp] + beta1[pp]*time[bb]    
    logIIV[pp,oo] = gamma[midx[pp,oo], pp] 
    Y[pp,oo] = rnorm(1, intercept[pp,oo] + AR[pp]*(Y[pp,oo-1] - intercept[pp,oo-1]), sqrt(exp(logIIV[pp,oo])))
   }
  }
}


#plot(Y[1,], type = "l",ylim=c(min(Y),max(Y)))
#for(i in 2:N){
#  lines(Y[i,])
#}


save(Y,N,O,X,#G3,
     time,nrPeriod,#Period2,
     nrObs,
     beta00, beta10, sd_beta0, sd_beta1,AR0,sd_AR,
     gamma0,sd_gamma,
     alpha,midx,
     file=paste0("RSTVP_gcm_simulateddata_shift",shift,"_N",N,"O",O,"_",r1,".Rdata"))
