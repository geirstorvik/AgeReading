#Script for analyzing age-readings from annual Herring spawning survey in February/March
#Only use data from 2022
library(AgeHerring)
library(rstan)
library(bridgesampling)
logit = function(p){log(p/(1-p))}
invlogit = function(x){1/(1+exp(-x))}
data(Herring.survey.year)
data(Ring.readErr)
data(Herring.Stratum)

## Analysis using stan
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
Amin=3;Amax=15
A = Amax-Amin+1
Astar = as.integer(mean(Amin,Amax))
alpha0=1
R1 = dim(Ring.readErr)[1];R2=dim(Herring.survey.year$d)[3];R=R1+R2
Y=dim(Herring.survey.year$d)[1]
S=dim(Herring.survey.year$d)[2]
#Number of strata
K = dim(Herring.Stratum$INDSUM)[2]
w = rep(NA,K)
d.comb = list(Y=Y,S=S,K=dim(Herring.Stratum$w)[2],R1=R1,R2=R2,Amin=Amin,Amax=Amax,A=Amax-Amin+1,Astar=Astar,
              readErr=Ring.readErr,
              D=Herring.survey.year$d,N=Herring.survey.year$N,alpharep=rep(alpha0,A),
              taupar=c(10,10),
              deltapar=c(2,15))
#              w=Herring.Stratum$w,INDSUM=Herring.Stratum$INDSUM)
##Make initial values

fitcomb.year = stan(file="exec/agereader_comb_year.stan",data=d.comb,iter=10000,chains=2,thin=25)

save(fitcomb.year,file="fitcomb.year10000.RData")
load("fitcomb.year10000.RData")

## Calculate stratum estimates
Nest = calc_stratum_estimate(fitcomb.year)
M1 = bridge_sampler(fitcomb.year,silent=TRUE)
estpar = summary(fitcomb.year)$summary[,1]
pars = c(paste("alpha0[",1:R,"]",sep=""),
         paste("beta0[",1:R,"]",sep=""))
plot(fitcomb.year,pars=pars)

pars = c(paste("alpha0[",1:R,"]",sep=""))
plot(fitcomb.year,pars=pars)

pars = c(paste("beta0[",1:R,"]",sep=""))
plot(fitcomb.year,pars=pars)


pars = c("Palpha0","Pbeta0","Pbeta1","Pphi")
plot(fitcomb.year,pars=pars)

pars = c("delta")
stan_dens(fitcomb.year,pars=pars)
traceplot(fitcomb.year,pars=pars)

ind = expand.grid(1:Y,1:A)
pars=paste("xi0[",ind[,1],",",ind[,2],"]",sep="")
plot(fitcomb.year,pars=pars)
xi0 = matrix(unlist(extract(fitcomb.year,pars=pars)),ncol=39)
traceplot(fitcomb.year,pars=pars)
pairs(fitcomb.year,pars=pars[1:5])

pars=c("xi0")
plot(fitcomb.year,pars=pars)

ind = expand.grid(1:2,2,1:A)
pars2=paste("xis[",ind[,1],",",ind[,2],",",ind[,3],"]",sep="")
plot(fitcomb.year,pars=pars2)
xis3 = matrix(unlist(extract(fitcomb.year,pars=pars)),ncol=78)

foo = summary(fitcomb.year)$summary[,1:2]
foo1 = cbind(foo[paste("alpha0[",1:5,"]",sep=""),1:2],
#             foo[paste("alpha1[",1:5,"]",sep=""),1:2],
             foo[paste("beta0[",1:5,"]",sep=""),1:2],
             foo[paste("Pbeta1",sep=""),1:2],
             foo[paste("Pphi",sep=""),1:2])
rownames(foo1) = 1:5
xtable(foo1,digits=3)
foo2 = cbind(foo[4+1:6,],foo[10+1:6,])
rownames(foo2) = 1:6
xtable(foo2,digits=3)

