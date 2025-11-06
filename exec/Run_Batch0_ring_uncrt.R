#Script for analyzing age-readings from annual Herring spawning survey in February/March
#Only use data from 2022
library(AgeHerring)
library(rstan)
logit = function(p){log(p/(1-p))}
invlogit = function(x){1/(1+exp(-x))}

## Analysis using stan
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
Amin=3;Amax=12
A = Amax-Amin+1
Astar = as.integer(mean(Amin,Amax))
Astar = 7
alpha0=1
R1 = 6
U=nrow(Ring.full)
#Number of strata
d.comb = list(U=U,R=R1,Amin=Amin,Amax=Amax,A=Amax-Amin+1,Astar=Astar,
              D=Ring.full,
              alpharep=rep(alpha0,A),
              taupar=c(3,3),
              deltapar=c(1,1),
              true_age=rep(NA,U),
              eps=0.001)

fitringuncrt0.year = stan(file="exec/agereader_year_strat0_uncrt.stan",data=d.comb,
                    iter=5000,chains=4,thin=100)
save(fitringuncrt0.year,file="fitring0.year50000.RData")

pars = "trphi"
foo = extract(fitring0.year,pars)[[1]]
for(i in 1:5){
  for(j in (i+1):6){
    show(c(i,j,ks.test(foo[,i],foo[,j])$p.value))
    }}
