#Script for analyzing age-readings from annual Herring spawning survey in February/March
#Only use data from 2022
library(AgeHerring)
library(rstan)
logit = function(p){log(p/(1-p))}
invlogit = function(x){1/(1+exp(-x))}
data(Ring.readErr)
data(Herring.Stratum)

## Choose dataset  (HerringSurvey/HerringSurvey202?)
data(HerringSurvey)
df = HerringSurvey

## Analysis using stan
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
Amin=3;Amax=12
A = Amax-Amin+1
Astar = as.integer(mean(Amin,Amax))
Astar = 5
alpha0=1
R1 = dim(Ring.readErr)[1]
R2=length(unique(df$unitreaderindex))
R=R1+R2
Y=length(unique(df$stratayearindex))
S=dim(df$d)[1]
#Number of strata

d.comb2 = list(Y=Y,S=df$S,K=df$K,U=df$U,
              Amin=Amin,Amax=Amax,A=Amax-Amin+1,Astar=Astar,
              D=df$d,N=df$N,
              stratayearindex=df$stratayearindex,
              stationstrataindex=df$stationstrataindex,
              unitstationindex=df$unitstationindex,
              alpharep=rep(alpha0,A),
              deltapar=c(1,1))
fitcomb2.year = stan(file="exec/agereader_comb_year_strat2.stan",data=d.comb2,
                     iter=100000,chains=4,thin=100)
save(fitcomb2.year,file="fitcomb2.year100000.RData")
