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
shift = as.numeric(shift)
G3 = c(rep(0, N/2),rep(1, N/2)) # two groups
time = 0:3
nrPeriod = length(time) # 4 periods
Period2 = c(rep(0, O/4), rep(1, O/4), rep(0, O/4), rep(0, O/4))
nrObs=c(1, seq(O/4, O, O/4))


# parameters related to intercepts
beta00 = 0 # intercept in GCM
sd_beta0 = 0.2
#beta1 = 0.5 # slope in GCM
beta10 = 0.5
sd_beta1 = 0.2
gamma0 = c(0, shift) # levels in two regimes 
sd_gamma = c(0, 0.2)

# parameters in the AR model
AR0 = 0.3  # AR parameter
sd_AR = 0.1
sd_noise0 = 0.5 # process noise standard deviation
sd_sd_noise = 0.1

# parameters in the RS model
# e.g., p(S_{i,t}=r|S_{i,t}=s) = alpha[1,r,s]+alpha[2,r,s]*Period2[tt]*G3[pp]
alpha = array(rep(NA, 2*2*2), dim=c(2,2,2))
# fix alpha11 and alpha22 to 0
for(i in 1:2){
  alpha[i,1,1] = 0 
  alpha[i,2,2] = 0
}

alpha[1,2,1] = -0.5   
alpha[1,1,2] = -0.5 

alpha[2,2,1] = 1 # G3 is more likely to switch to the high-level regime in P2
alpha[2,1,2] = -1 # G3 is less likely to switch to the low-level regime in P2


# simulate data

odds = array(rep(NA,N*O*2),dim = c(N,O,2))
midx = matrix(NA,N,O) # regime indicator
p_2s = matrix(NA,N,O) # probability of being in regime 2

Y=matrix(NA,N,O) # observed variable
gamma=matrix(NA,2,N) 
beta0=rep(NA,N)
beta1=rep(NA,N)
AR=rep(NA,N)
sd_noise=rep(NA,N)
intercept=matrix(NA,N,nrPeriod) 
for(pp in 1:N){
  odds[pp,1,1] = 0.9 # more likely to stay in regime 1 at the first time point
  odds[pp,1,2] = 0.1
  p_2s[pp,1] = odds[pp,1,2]/(odds[pp,1,1]+odds[pp,1,2])
  midx[pp,1] = sample.int(n=2,size=1,prob=odds[pp,1,])
  gamma[1,pp] = rnorm(1, gamma0[1], sd_gamma[1])
  gamma[2,pp] = rnorm(1, gamma0[2], sd_gamma[2])
  beta0[pp] = rnorm(1, beta00, sd_beta0)
  beta1[pp] = rnorm(1, beta10, sd_beta1)
  AR[pp] = rnorm(1, AR0, sd_AR)
  sd_noise[pp] = rnorm(1, sd_noise0, sd_sd_noise)
  Y[pp,1] = rnorm(1, 0, sd_noise[pp])
  
  for (bb in 1:nrPeriod){
   #intercept[pp,bb] = beta0[pp] + beta1*time[bb]
   intercept[pp,bb] = beta0[pp] + beta1[pp]*time[bb]
   for(oo in (nrObs[bb]+1):nrObs[bb+1]){
    # from regime midx[pp,oo-1] to regime 1
    odds[pp,oo,1] = exp(alpha[1,1,midx[pp,oo-1]]+alpha[2,1,midx[pp,oo-1]]*Period2[oo]*G3[pp])
    # from regime midx[pp,oo-1] to regime 2
    odds[pp,oo,2] = exp(alpha[1,2,midx[pp,oo-1]]+alpha[2,2,midx[pp,oo-1]]*Period2[oo]*G3[pp])
    midx[pp,oo]=sample.int(n=2,size=1,prob=odds[pp,oo,])
    p_2s[pp,oo] = odds[pp,oo,2]/(odds[pp,oo,1]+odds[pp,oo,2])
    
    Y[pp,oo] = rnorm(1, gamma[midx[pp,oo],pp] + intercept[pp,bb] + AR[pp]*(Y[pp,oo-1] - intercept[pp,bb] - gamma[midx[pp,oo-1],pp]), sd_noise[pp])
   }
  }
}


#plot(Y[1,], type = "l",ylim=c(min(Y),max(Y)))
#for(i in 2:N){
#  lines(Y[i,])
#}


save(Y,N,O,G3,time,nrPeriod,Period2,nrObs,
     beta00,sd_beta0,beta10,sd_beta1,AR0,sd_AR,sd_noise0,sd_sd_noise,gamma0,sd_gamma,AR,sd_noise,alpha,
     odds,midx,gamma,beta0,intercept, 
     file=paste0("RSTVP_gcm_simulateddata_shift",shift,"_N",N,"O",O,"_",r1,".Rdata"))
