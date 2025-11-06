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


set.seed(25423)

load("fitcomb0.year100000.RData")              #reader dependent age errors
xi1 = apply(extract(fitcomb0.year,pars="xi1")$xi1,c(2,3),mean)
xi2 = apply(extract(fitcomb0.year,pars="xi2")$xi2,c(2,3),mean)
delta = apply(extract(fitcomb0.year,pars="delta")$delta,c(2),mean)

## First simulate xi2
S = dim(xi2)[1]
for(s in 1:S){
  xi2[s,] = rdirichlet(1,delta[2]*xi1[df$stationstrataindex[s],]);
}

## First using reader-dependent error matrices

M1 = apply(extract(fitcomb0.year,pars="M1")$M1,c(2,3,4),mean)

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
usethis::use_data(Sim.Ring.readErr,overwrite=TRUE)

## Simulate Survey data
Sim.HerringSurvey = HerringSurvey
for(s in 1:nrow(HerringSurvey$d))
{
  N =  sum(HerringSurvey$d[s,])
  u = HerringSurvey$unitstationindex[s]
  r = R1 + HerringSurvey$unitreaderindex[s]
  Sim.HerringSurvey$d[s,] = table(c(sample(1:A,N,prob=M1[r,,]%*%xi2[u,],replace=TRUE),1:A))-1
}
usethis::use_data(Sim.HerringSurvey,overwrite=TRUE)


## Then using reader-independent error matrices

load("fitcomb1.year100000.RData")              #reader dependent age errors
M2 = apply(extract(fitcomb1.year,pars="M1")$M1,c(2,3),mean)

## simulate ring data
Sim2.Ring.readErr = Ring.readErr
for(r in 1:R1)
{
  for(a in 1:A)
  {
    N = sum(Ring.readErr[r,a,])
    Sim2.Ring.readErr[r,a,] = table(c(sample(1:A,N,prob=M2[a,],replace=TRUE),1:A))-1
  }
}
usethis::use_data(Sim2.Ring.readErr,overwrite=TRUE)

## Simulate Survey data
Sim2.HerringSurvey = HerringSurvey
for(s in 1:nrow(HerringSurvey$d))
{
  N =  sum(HerringSurvey$d[s,])
  u = HerringSurvey$unitstationindex[s]
  Sim2.HerringSurvey$d[s,] = table(c(sample(1:A,N,prob=M2%*%xi2[u,],replace=TRUE),1:A))-1
}
usethis::use_data(Sim2.HerringSurvey,overwrite=TRUE)

## Using reader 1

M1 = apply(extract(fitcomb0.year,pars="M1")$M1,c(2,3,4),mean)
M2 = M1[1,,]

## simulate ring data
SimR1.Ring.readErr = Ring.readErr
for(r in 1:R1)
{
  for(a in 1:A)
  {
    N = sum(Ring.readErr[r,a,])
    SimR1.Ring.readErr[r,a,] = table(c(sample(1:A,N,prob=M2[a,],replace=TRUE),1:A))-1
  }
}
usethis::use_data(SimR1.Ring.readErr,overwrite=TRUE)

## Simulate Survey data
SimR1.HerringSurvey = HerringSurvey
for(s in 1:nrow(HerringSurvey$d))
{
  N =  sum(HerringSurvey$d[s,])
  u = HerringSurvey$unitstationindex[s]
  SimR1.HerringSurvey$d[s,] = table(c(sample(1:A,N,prob=M2%*%xi2[u,],replace=TRUE),1:A))-1
}
usethis::use_data(SimR1.HerringSurvey,overwrite=TRUE)

## Using reader 3

M1 = apply(extract(fitcomb0.year,pars="M1")$M1,c(2,3,4),mean)
M2 = M1[3,,]

## simulate ring data
SimR3.Ring.readErr = Ring.readErr
for(r in 1:R1)
{
  for(a in 1:A)
  {
    N = sum(Ring.readErr[r,a,])
    SimR3.Ring.readErr[r,a,] = table(c(sample(1:A,N,prob=M2[a,],replace=TRUE),1:A))-1
  }
}
usethis::use_data(SimR3.Ring.readErr,overwrite=TRUE)

## Simulate Survey data
SimR3.HerringSurvey = HerringSurvey
for(s in 1:nrow(HerringSurvey$d))
{
  N =  sum(HerringSurvey$d[s,])
  u = HerringSurvey$unitstationindex[s]
  SimR3.HerringSurvey$d[s,] = table(c(sample(1:A,N,prob=M2%*%xi2[u,],replace=TRUE),1:A))-1
}
usethis::use_data(SimR3.HerringSurvey,overwrite=TRUE)

## Using reader 5

M1 = apply(extract(fitcomb0.year,pars="M1")$M1,c(2,3,4),mean)
M2 = M1[5,,]

## simulate ring data
SimR5.Ring.readErr = Ring.readErr
for(r in 1:R1)
{
  for(a in 1:A)
  {
    N = sum(Ring.readErr[r,a,])
    SimR5.Ring.readErr[r,a,] = table(c(sample(1:A,N,prob=M2[a,],replace=TRUE),1:A))-1
  }
}
usethis::use_data(SimR5.Ring.readErr,overwrite=TRUE)

## Simulate Survey data
SimR5.HerringSurvey = HerringSurvey
for(s in 1:nrow(HerringSurvey$d))
{
  N =  sum(HerringSurvey$d[s,])
  u = HerringSurvey$unitstationindex[s]
  SimR5.HerringSurvey$d[s,] = table(c(sample(1:A,N,prob=M2%*%xi2[u,],replace=TRUE),1:A))-1
}
usethis::use_data(SimR5.HerringSurvey,overwrite=TRUE)



## Using reader 6

M1 = apply(extract(fitcomb0.year,pars="M1")$M1,c(2,3,4),mean)
M2 = M1[6,,]

## simulate ring data
SimR6.Ring.readErr = Ring.readErr
for(r in 1:R1)
{
  for(a in 1:A)
  {
    N = sum(Ring.readErr[r,a,])
    SimR6.Ring.readErr[r,a,] = table(c(sample(1:A,N,prob=M2[a,],replace=TRUE),1:A))-1
  }
}
usethis::use_data(SimR6.Ring.readErr,overwrite=TRUE)

## Simulate Survey data
SimR6.HerringSurvey = HerringSurvey
for(s in 1:nrow(HerringSurvey$d))
{
  N =  sum(HerringSurvey$d[s,])
  u = HerringSurvey$unitstationindex[s]
  SimR6.HerringSurvey$d[s,] = table(c(sample(1:A,N,prob=M2%*%xi2[u,],replace=TRUE),1:A))-1
}
usethis::use_data(SimR6.HerringSurvey,overwrite=TRUE)

