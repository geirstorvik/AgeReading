#Script for analyzing age-readings from annual Herring spawning survey in February/March
#Only use data from 2022
library(AgeHerring)
library(rstan)
logit = function(p){log(p/(1-p))}
invlogit = function(x){1/(1+exp(-x))}

## Choose dataset  (HerringSurvey/HerringSurvey202?)
data(SimR1.Ring.readErr)
dataRing = SimR1.Ring.readErr
data(SimR1.HerringSurvey)
df = SimR1.HerringSurvey

## Analysis using stan
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
Amin=3;Amax=12
A = Amax-Amin+1
Astar = as.integer(mean(Amin,Amax))
Astar = 5
alpha0=1
R1 = dim(dataRing)[1]
R2=length(unique(df$unitreaderindex))
R=R1+R2
Y=length(unique(df$stratayearindex))
S=dim(df$d)[1]
#Number of strata
d.comb = list(Y=Y,S=df$S,K=df$K,U=df$U,
              R1=R1,R2=R2,Amin=Amin,Amax=Amax,A=Amax-Amin+1,Astar=Astar,
              readErr=dataRing,
              D=df$d,N=df$N,
              stratayearindex=df$stratayearindex,
              stationstrataindex=df$stationstrataindex,
              unitstationindex=df$unitstationindex,
              unitreaderindex=df$unitreaderindex,
              alpharep=rep(alpha0,A),
              taupar=c(0,10),
              deltapar=c(1,1))

SimR1fitcomb1.year = stan(file="exec/agereader_comb_year_strat1.stan",data=d.comb,
                    iter=100000,chains=4,thin=100)
save(SimR1fitcomb1.year,file="SimR1fitcomb1.year100000.RData")
