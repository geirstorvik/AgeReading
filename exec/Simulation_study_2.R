## Script for simulating data

## First extract parameters from real data runs

library(AgeHerring)
library(rstan)
library(ggplot2)
require(gridExtra)
library(dirmult)

## Analysis using stan
data(HerringSurvey)
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
#Number of strata

d.comb = list(Y=Y,S=df$S,K=df$K,U=df$U,
              R1=R1,R2=R2,Amin=Amin,Amax=Amax,A=Amax-Amin+1,Astar=Astar,
              readErr=Ring.readErr,
              D=df$d,N=df$N,
              stratayearindex=df$stratayearindex,
              stationstrataindex=df$stationstrataindex,
              unitstationindex=df$unitstationindex,
              unitreaderindex=df$unitreaderindex,
              alpharep=rep(alpha0,A),
              taupar=c(3,3),
              deltapar=c(1,1),
              eps=0.001)


set.seed(456)

load("fitcomb0.year100000.RData")              #reader dependent age errors
xi1 = apply(rstan::extract(fitcomb0.year,pars="xi1")$xi1,c(2,3),mean)
xi2 = apply(rstan::extract(fitcomb0.year,pars="xi2")$xi2,c(2,3),mean)
delta = apply(rstan::extract(fitcomb0.year,pars="delta")$delta,c(2),mean)

## First simulate xi2
S = dim(xi2)[1]
for(s in 1:S){
  xi2[s,] = rdirichlet(1,delta[2]*xi1[df$stationstrataindex[s],]);
}

## First using reader-dependent error matrices

M1 = apply(rstan::extract(fitcomb0.year,pars="M1")$M1,c(2,3,4),mean)

## simulate ring data
S2Sim.Ring.readErr = Ring.readErr
for(r in 1:R1)
{
  for(a in 1:A)
  {
    N = sum(Ring.readErr[r,a,])
    S2Sim.Ring.readErr[r,a,] = table(c(sample(1:A,N,prob=M1[r,a,],replace=TRUE),1:A))-1
  }
}
usethis::use_data(S2Sim.Ring.readErr,overwrite=TRUE)

## S2Simulate Survey data
S2Sim.HerringSurvey = HerringSurvey
for(s in 1:nrow(HerringSurvey$d))
{
  N =  sum(HerringSurvey$d[s,])
  u = HerringSurvey$unitstationindex[s]
  r = R1 + HerringSurvey$unitreaderindex[s]
  S2Sim.HerringSurvey$d[s,] = table(c(sample(1:A,N,prob=M1[r,,]%*%xi2[u,],replace=TRUE),1:A))-1
}
usethis::use_data(S2Sim.HerringSurvey,overwrite=TRUE)


## Then using reader-independent error matrices

load("fitcomb1.year100000.RData")              #reader dependent age errors
M2 = apply(rstan::extract(fitcomb1.year,pars="M1")$M1,c(2,3),mean)

## simulate ring data
S2Sim2.Ring.readErr = Ring.readErr
for(r in 1:R1)
{
  for(a in 1:A)
  {
    N = sum(Ring.readErr[r,a,])
    S2Sim2.Ring.readErr[r,a,] = table(c(sample(1:A,N,prob=M2[a,],replace=TRUE),1:A))-1
  }
}
usethis::use_data(S2Sim2.Ring.readErr,overwrite=TRUE)

## S2Simulate Survey data
S2Sim2.HerringSurvey = HerringSurvey
for(s in 1:nrow(HerringSurvey$d))
{
  N =  sum(HerringSurvey$d[s,])
  u = HerringSurvey$unitstationindex[s]
  S2Sim2.HerringSurvey$d[s,] = table(c(sample(1:A,N,prob=M2%*%xi2[u,],replace=TRUE),1:A))-1
}
usethis::use_data(S2Sim2.HerringSurvey,overwrite=TRUE)

## Using reader 1

M1 = apply(rstan::extract(fitcomb0.year,pars="M1")$M1,c(2,3,4),mean)
M2 = M1[1,,]

## simulate ring data
S2SimR1.Ring.readErr = Ring.readErr
for(r in 1:R1)
{
  for(a in 1:A)
  {
    N = sum(Ring.readErr[r,a,])
    S2SimR1.Ring.readErr[r,a,] = table(c(sample(1:A,N,prob=M2[a,],replace=TRUE),1:A))-1
  }
}
usethis::use_data(S2SimR1.Ring.readErr,overwrite=TRUE)

## S2Simulate Survey data
S2SimR1.HerringSurvey = HerringSurvey
for(s in 1:nrow(HerringSurvey$d))
{
  N =  sum(HerringSurvey$d[s,])
  u = HerringSurvey$unitstationindex[s]
  S2SimR1.HerringSurvey$d[s,] = table(c(sample(1:A,N,prob=M2%*%xi2[u,],replace=TRUE),1:A))-1
}
usethis::use_data(S2SimR1.HerringSurvey,overwrite=TRUE)

## Using reader 3

M1 = apply(rstan::extract(fitcomb0.year,pars="M1")$M1,c(2,3,4),mean)
M2 = M1[3,,]

## simulate ring data
S2SimR3.Ring.readErr = Ring.readErr
for(r in 1:R1)
{
  for(a in 1:A)
  {
    N = sum(Ring.readErr[r,a,])
    S2SimR3.Ring.readErr[r,a,] = table(c(sample(1:A,N,prob=M2[a,],replace=TRUE),1:A))-1
  }
}
usethis::use_data(S2SimR3.Ring.readErr,overwrite=TRUE)

## S2Simulate Survey data
S2SimR3.HerringSurvey = HerringSurvey
for(s in 1:nrow(HerringSurvey$d))
{
  N =  sum(HerringSurvey$d[s,])
  u = HerringSurvey$unitstationindex[s]
  S2SimR3.HerringSurvey$d[s,] = table(c(sample(1:A,N,prob=M2%*%xi2[u,],replace=TRUE),1:A))-1
}
usethis::use_data(S2SimR3.HerringSurvey,overwrite=TRUE)

## Using reader 5

M1 = apply(rstan::extract(fitcomb0.year,pars="M1")$M1,c(2,3,4),mean)
M2 = M1[5,,]

## simulate ring data
S2SimR5.Ring.readErr = Ring.readErr
for(r in 1:R1)
{
  for(a in 1:A)
  {
    N = sum(Ring.readErr[r,a,])
    S2SimR5.Ring.readErr[r,a,] = table(c(sample(1:A,N,prob=M2[a,],replace=TRUE),1:A))-1
  }
}
usethis::use_data(S2SimR5.Ring.readErr,overwrite=TRUE)

## S2Simulate Survey data
S2SimR5.HerringSurvey = HerringSurvey
for(s in 1:nrow(HerringSurvey$d))
{
  N =  sum(HerringSurvey$d[s,])
  u = HerringSurvey$unitstationindex[s]
  S2SimR5.HerringSurvey$d[s,] = table(c(sample(1:A,N,prob=M2%*%xi2[u,],replace=TRUE),1:A))-1
}
usethis::use_data(S2SimR5.HerringSurvey,overwrite=TRUE)



## Using reader 6

M1 = apply(rstan::extract(fitcomb0.year,pars="M1")$M1,c(2,3,4),mean)
M2 = M1[6,,]

## simulate ring data
S2SimR6.Ring.readErr = Ring.readErr
for(r in 1:R1)
{
  for(a in 1:A)
  {
    N = sum(Ring.readErr[r,a,])
    S2SimR6.Ring.readErr[r,a,] = table(c(sample(1:A,N,prob=M2[a,],replace=TRUE),1:A))-1
  }
}
usethis::use_data(S2SimR6.Ring.readErr,overwrite=TRUE)

## S2Simulate Survey data
S2SimR6.HerringSurvey = HerringSurvey
for(s in 1:nrow(HerringSurvey$d))
{
  N =  sum(HerringSurvey$d[s,])
  u = HerringSurvey$unitstationindex[s]
  S2SimR6.HerringSurvey$d[s,] = table(c(sample(1:A,N,prob=M2%*%xi2[u,],replace=TRUE),1:A))-1
}
usethis::use_data(S2SimR6.HerringSurvey,overwrite=TRUE)

