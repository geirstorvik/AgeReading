#Script for analyzing age-readings from annual Herring spawning survey in February/March
#Only use data from 2022
library(AgeHerring)
library(rstan)
library(TruncatedNormal)
logit = function(p){log(p/(1-p))}
invlogit = function(x){1/(1+exp(-x))}

## Choose dataset  (HerringSurvey/HerringSurvey202?)
data(SimR5.Ring.readErr)
dataRing = SimR5.Ring.readErr
data(SimR5.HerringSurvey)
df = SimR5.HerringSurvey

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
#Number of strata
set.seed(3453)

## Initial values by running only on ring data
d.comb0 = list(R=R1,Amin=Amin,Amax=Amax,A=Amax-Amin+1,Astar=Astar,
               readErr=dataRing,
               alpharep=rep(alpha0,A),
               taupar=c(20,2,5,1,20),
               deltapar=c(1,1),
               eps=0.001)

runinit = 1
if(runinit)
{
SimR5fitring0.year = stan(file="exec/agereader_year_strat0.stan",data=d.comb0,
                        iter=50000,chains=4,thin=25)
save(SimR5fitring0.year,file="SimR5fitring0.year50000.RData")
}
if(!runinit)
 load("SimR5fitring0.year50000.RData")

simext = rstan::extract(SimR5fitring0.year)
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
       beta1=c(simext$beta1[ind,],rtnorm(R2,simext$Pbeta1[ind],sqrt(simext$tau[ind,4]),ub=0.0)),
       trphi=c(simext$trphi[ind,],rnorm(R2,simext$Ptrphi[ind],sqrt(simext$tau[ind,5]))))
}

init2 = lapply(1:4,function(id) init1(id=id))

d.comb = list(Y=Y,S=df$S,K=df$K,U=df$U,
              R1=R1,R2=R2,Amin=Amin,Amax=Amax,A=Amax-Amin+1,Astar=Astar,
              readErr=dataRing,
              D=df$d,N=df$N,
              stratayearindex=df$stratayearindex,
              stationstrataindex=df$stationstrataindex,
              unitstationindex=df$unitstationindex,
              unitreaderindex=df$unitreaderindex,
              alpharep=rep(alpha0,A),
              taupar=c(20,2,5,1,20),
              deltapar=c(1,1),
              eps=0.001)

SimR5fitcomb0.year = stan(file="exec/agereader_comb_year_strat0.stan",data=d.comb,
                    iter=100000,chains=4,thin=100,init=init2)
save(SimR5fitcomb0.year,file="SimR5fitcomb0.year100000.RData")


