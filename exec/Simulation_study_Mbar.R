## Script for simulating data

## First extract parameters from real data runs

library(AgeHerring)
library(rstan)
library(ggplot2)
require(gridExtra)
library(dirmult)

## Analysis using stan
data(HerringSurvey)
data(Ring.readErr)
df = HerringSurvey

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
Amin=3;Amax=12
A = Amax-Amin+1
Astar = as.integer(mean(Amin,Amax))
Astar = 7
alpha0=1
R1 = dim(Ring.readErr)[1]
R2=length(unique(df$unitreaderindex))
R=R1+R2
Y=length(unique(df$stratayearindex))
S=dim(df$d)[1]


set.seed(456)

load("fitcomb0.year100000.RData")              #reader dependent age errors
xi1 = apply(rstan::extract(fitcomb0.year,pars="xi1")$xi1,c(2,3),mean)
xi2 = apply(rstan::extract(fitcomb0.year,pars="xi2")$xi2,c(2,3),mean)
delta = apply(rstan::extract(fitcomb0.year,pars="delta")$delta,c(2),mean)

## Reader-dependent error matrices

M1 = apply(rstan::extract(fitcomb0.year,pars="M1")$M1,c(2,3,4),mean)

#Number of simulations
Nsim=10
Nmcmc=100000
for(m in 1:Nsim)
{
  ## First simulate xi2
  S = dim(xi2)[1]
  for(s in 1:S){
    xi2[s,] = rdirichlet(1,delta[2]*xi1[df$stationstrataindex[s],]);
  }

  ## simulate ring data
  Sim.Ring.readErr = Ring.readErr
  for(r in 1:R1)
  {
    for(a in 1:A)
    {
      N = sum(Ring.readErr[r,a,])
      Sim.Ring.readErr[r,a,] = table(c(sample(1:A,N,prob=M1[r,a,],replace=TRUE),1:A))-1
    }
  }
  usethis::use_data(S2Sim.Ring.readErr,overwrite=TRUE)

  ## S2Simulate Survey data
  Sim.HerringSurvey = HerringSurvey
  for(s in 1:nrow(HerringSurvey$d))
  {
    N =  sum(HerringSurvey$d[s,])
    u = HerringSurvey$unitstationindex[s]
    r = R1 + HerringSurvey$unitreaderindex[s]
    Sim.HerringSurvey$d[s,] = table(c(sample(1:A,N,prob=M1[r,,]%*%xi2[u,],replace=TRUE),1:A))-1
  }
  usethis::use_data(S2Sim.HerringSurvey,overwrite=TRUE)

  #Perform inference
  df = Sim.HerringSurvey

  ## Analysis using stan
  rstan_options(auto_write = TRUE)
  options(mc.cores = parallel::detectCores())
  Amin=3;Amax=12
  A = Amax-Amin+1
  Astar = as.integer(mean(Amin,Amax))
  Astar = 7
  alpha0=1
  R1 = dim(dataRing)[1]
  R2=length(unique(df$unitreaderindex))
  R=R1+R2
  Y=length(unique(df$stratayearindex))
  S=dim(df$d)[1]

  ## Initial values by running only on ring data
  d.comb0 = list(R=R1,Amin=Amin,Amax=Amax,A=Amax-Amin+1,Astar=Astar,
                 readErr=Sim.Ring.readErr,
                 alpharep=rep(alpha0,A),
                 taupar=c(20,2,5,1,20),
                 deltapar=c(1,1),
                 eps=0.001)

  simfitring0.year = stan(file="exec/agereader_year_strat0.stan",data=d.comb0,
                              iter=Nmcmc/2,chains=4,thin=25)
  save(simfitring0.year,file=paste0("Simfitring0Mbar.",m,".",Nmcmc/2,".RData"))

  simext = rstan::extract(simfitring0.year)
  M = length(simext$Palpha0)
  init1 = function(id)
  {
    ind = sample(1:M,1)
    list(
      Palpha0=simext$Palpha0[ind],Palpha1=simext$Palpha1[ind],Pbeta0=simext$Pbeta0[ind],Pbeta1=simext$Pbeta1[ind],
      Ptrphi=simext$Ptrphi[ind],delta=simext$delta[ind,],
      tau=simext$tau[ind,],
      alpha0=c(simext$alpha0[ind,],rnorm(R2,simext$Palpha0[ind],sqrt(simext$tau[ind,1]))),
      beta0=c(simext$beta0[ind,],rnorm(R2,simext$Pbeta0[ind],sqrt(simext$tau[ind,2]))),
      alpha1=c(simext$alpha1[ind,],rnorm(R2,simext$Palpha1[ind],sqrt(simext$tau[ind,3]))),
      beta1=c(simext$beta1[ind,],rnorm(R2,simext$Pbeta1[ind],sqrt(simext$tau[ind,4]))),
      #beta1=c(simext$beta1[ind,],rtnorm(R2,simext$Pbeta1[ind],sqrt(simext$tau[ind,4]),ub=0.0)),
      trphi=c(simext$trphi[ind,],rnorm(R2,simext$Ptrphi[ind],sqrt(simext$tau[ind,5]))))
  }

  init2 = lapply(1:4,function(id) init1(id=id))

  d.comb = list(Y=Y,S=df$S,K=df$K,U=df$U,
                R1=R1,R2=R2,Amin=Amin,Amax=Amax,A=Amax-Amin+1,Astar=Astar,
                readErr=Sim.Ring.readErr,
                D=df$d,N=df$N,
                stratayearindex=df$stratayearindex,
                stationstrataindex=df$stationstrataindex,
                unitstationindex=df$unitstationindex,
                unitreaderindex=df$unitreaderindex,
                alpharep=rep(alpha0,A),
                taupar=c(1,1,1,1,1),
                deltapar=c(1,1),
                eps=0.001)

  simfitcomb0.year = stan(file="exec/agereader_comb_year_strat0.stan",data=d.comb,
                            iter=Nmcmc,chains=4,thin=100,init=init2)
  save(simfitcomb0.year,file=paste0("simfitcomb0Mbar",m,".",Nmcmc,".RData"))
}

#Summarize results
library(scoringRules)
Mxi1 = array(xi1,c(1,15,10))
Truef = apply(Make_Numfish_sample(Mxi1,Herring_D_boot,HerringSurvey$stratayearindex),2:3,mean)
fMbar = array(NA,c(Nsim,3,3,10))
crpsMbar = array(NA,c(Y,A,Nsim))
for(m in 1:Nsim)
{
  print(m)
  load(paste0("simfitcomb0Mbar",m,".",Nmcmc,".RData"))              #reader dependent age errors
  foo = rstan::extract(simfitcomb0.year)$xi1
  res0 = Make_Numfish_sample(foo,Herring_D_boot,HerringSurvey$stratayearindex)
  fMbar[m,,,] = apply(res0,c(2,3),quantile,probs=c(0.025,0.5,0.975))
  for(y in 1:3)
    crpsMbar[y,,m] = crps_sample(Truef[y,],t(res0[,y,]))
}
save(crpsMbar,file="crpsMbar.RData")
ages = matrix(rep(2+1:A,Nsim),nrow=A,byrow=T)
for(y in 1:3)
{
  pdf(paste0("Sim_Mbar_N_year",2020+y,".pdf"),height=5,width=10)
  matplot(t(ages),t(fMbar[,2,y,1:A]),type="l",lty=1,col=2,log="y",xlab="",
          ylab="",
          #main=paste0("Year",2020+y),
          cex.axis=2)
  lines(2+1:A,Truef[y,],lwd=3)
  dev.off()
}

withinfMbar = matrix(nrow=3,ncol=A)
for(y in 1:3)
  for(a in 1:A)
    withinfMbar[y,a] = mean(Truef[y,a] %[]% fMbar[,c(1,3),y,a])
