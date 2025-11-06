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
Astar = 7
alpha0=1
R1 = dim(Ring.readErr)[1]
#Number of strata
d.comb = list(R=R1,Amin=Amin,Amax=Amax,A=Amax-Amin+1,Astar=Astar,
              readErr=Ring.readErr,
              alpharep=rep(alpha0,A),
              taupar=c(0,10),
              deltapar=c(1,1),
              eps=0.001)

fitring0.year = stan(file="exec/agereader_year_strat0.stan",data=d.comb,
                    iter=5000,chains=4,thin=10)
save(fitring0.year,file="fitring0.year50000.RData")

if(0)
{
pars = "trphi"
foo = extract(fitring0.year,pars)[[1]]
for(i in 1:5){
  for(j in (i+1):6){
    show(c(i,j,ks.test(foo[,i],foo[,j])$p.value))
  }}


M1 = apply(extract(fitring0.year,pars="M1")$M1,c(2,3,4),mean)
M1diag = diag(M1[1,,])
for(r in 2:R1)
  M1diag = cbind(M1diag,diag(M1[r,,]))
matplot(M1diag,type="l",lty=1)
legend("bottomleft",1:R1,1:R1,lty=1,col=1:R1)

M1samp = extract(fitring0.year,pars="M1")$M1
ks.test(M1samp[,2,9,9],M1samp[,5,9,9])

d = list()
a = 4
for(i in 1:R1)
  d[[i]] = density(M1samp[,i,a,a])

matplot(cbind(d[[1]]$x,d[[2]]$x,d[[3]]$x,d[[4]]$x,d[[5]]$x,d[[6]]$x),
        cbind(d[[1]]$y,d[[2]]$y,d[[3]]$y,d[[4]]$y,d[[5]]$y,d[[6]]$y),
        lty=1,col=1:6,type="l",lwd=2)
}
