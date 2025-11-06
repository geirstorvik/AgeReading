#Script for analyzing age-readings from annual Herring spawning survey in February/March
#Only use data from 2022
library(AgeHerring)
library(rstan)
#library(bridgesampling)
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
Astar = as.integer(mean(c(Amin,Amax)))
Astar = 5
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
              taupar=c(10,10),
              deltapar=c(1,1))
#              w=Herring.Stratum$w,INDSUM=Herring.Stratum$INDSUM)
##Make initial values

fitcomb.year = stan(file="exec/agereader_comb_year_strat.stan",data=d.comb,
                    iter=10000,chains=4,thin=250)
save(fitcomb.year,file="fitcomb.year10000.RData")
load("fitcomb0.year100000.RData")

data("Herring_D_boot")
foo = rstan::extract(fitcomb.year)$xi1
res = Make_Numfish_sample(foo,Herring_D_boot,df$stratayearindex)

d.comb2 = list(Y=Y,S=df$S,K=df$K,U=df$U,
              Amin=Amin,Amax=Amax,A=Amax-Amin+1,Astar=Astar,
              D=df$d,N=df$N,
              stratayearindex=df$stratayearindex,
              stationstrataindex=df$stationstrataindex,
              unitstationindex=df$unitstationindex,
              alpharep=rep(alpha0,A),
              deltapar=c(1,1))
fitcombnoErr.year = stan(file="exec/agereader_comb_year_strat_noAgeerr.stan",data=d.comb2,
                     iter=50000,chains=4,thin=250)
save(fitcombnoErr.year,file="fitcombnoErr.year50000.RData")
load("fitcombnoErr.year50000.RData")


foo = rstan::extract(fitcomb1.year)$xi1
resM = Make_Numfish_sample(foo,Herring_D_boot,df$stratayearindex)


foo = rstan::extract(fitcombnoErr.year)$xi1
resnoErr = Make_Numfish_sample(foo,Herring_D_boot,df$stratayearindex)
Dya = matrix(nrow=Y,ncol=A)
par(mfrow=c(2,5))
for(y in 1:3)
  for(a in 1:A)
  {
    Dya[y,a] = ks.test(res[,y,a],resnoErr[,y,a])$statistic
    Dya[y,a] = ks.test(res[,y,a],resM[,y,a])$statistic
    if(Dya[y,a]>0.0)
    {
      xlim = range(c(res[,y,a],resM[,y,a]))
      hist(resnoErr[,y,a],xlim=xlim,,breaks=30,
           xlab="Numbers at age",main=paste("Year=",2020+y,"Age=",a+2),
           probability=TRUE,col=rgb(1,0,0,1/4))
      hist(resM[,y,a],xlim=xlim,breaks=30,
           probability=TRUE,col=rgb(0,1,1,1/4),add=T)
      hist(res[,y,a],xlim=xlim,breaks=30,
           probability=TRUE,col=rgb(0,0,1,1/4),add=T)
    }
  }

require(gridExtra)
library(ggplot2)
l = list()
for(a in 1:A)
{
 df = data.frame(y=c(res[,y,a],resM[,y,a],resnoErr[,y,a]),
                cl=c(rep(1,dim(res)[1]),rep(2,dim(resM)[1]),
                     rep(3,dim(resnoErr)[1])))
 df$cl = as.factor(df$cl)
 l[[a]] = ggplot(df,aes(x=y,color=cl,fill=cl)) +
   geom_density(alpha=0.5,show.legend=FALSE)
}
grid.arrange(l[[1]],l[[2]],l[[3]],l[[4]],
             l[[5]],l[[6]],l[[7]],l[[8]],
             l[[9]],l[[10]],ncol=5,nrow=2)
par(mfrow=c(2,))
for(a in c(2,3,4,5,6,7))
{
  hist(res[,1,a],xlim=c(0,2e+10),xlab="Numbers at age",main=paste("Age=",a+2),col=rgb(0,0,1,1/4))
  hist(resnoErr[,1,a],xlim=c(0,2e+10),main=paste("Age=",a+3),col=rgb(1,0,0,1/4),add=T)
}
par(mfrow=c(1,1))
a=3
hist(res[,1,a],xlim=c(0,2e+10),xlab="Numbers at age",main=paste("Age=",a+2),col=rgb(0,0,1,1/4))
hist(resnoErr[,1,a],xlim=c(0,2e+10),main=paste("Age=",a+3),col=rgb(1,0,0,1/4),add=T)



## Calculate stratum estimates
Nest = calc_stratum_estimate(fitcomb.year,Herring.Stratum$w,Herring.Stratum$INDSUM,A)
M1 = bridge_sampler(fitcomb.year,silent=TRUE)
estpar = summary(fitcomb.year)$summary[,1]
pars = c(paste("alpha0[",1:R,"]",sep=""),
         paste("beta0[",1:R,"]",sep=""))
plot(fitcomb.year,pars=pars)

pars = c(paste("alpha0[",1:R,"]",sep=""))
plot(fitcomb.year,pars=pars)

pars = c(paste("beta0[",1:R,"]",sep=""))
plot(fitcomb.year,pars=pars)


pars = c("Palpha0","Palpha1","Pbeta0","Pbeta1","Pphi")
plot(fitcomb0.year,pars=pars)

pars = c("delta")
plot(fitcomb.year,pars=pars)
stan_dens(fitcomb.year,pars=pars)
traceplot(fitcomb.year,pars=pars)

pdf("xi0_2021.pdf",height=10,width=4)
ind = expand.grid(1,1:A)
pars=paste("xi0[",ind[,1],",",ind[,2],"]",sep="")
plot(fitcomb.year,pars=pars)
dev.off()
pdf("xi0_2022.pdf",height=10,width=4)
ind = expand.grid(2,1:A)
pars=paste("xi0[",ind[,1],",",ind[,2],"]",sep="")
plot(fitcomb.year,pars=pars)
dev.off()
pdf("xi0_2023.pdf",height=10,width=4)
ind = expand.grid(3,1:A)
pars=paste("xi0[",ind[,1],",",ind[,2],"]",sep="")
plot(fitcomb.year,pars=pars)
dev.off()

xi0 = matrix(unlist(extract(fitcomb.year,pars=pars)),ncol=39)
traceplot(fitcomb.year,pars=pars)
pairs(fitcomb.year,pars=pars[1:5])

library(fields)
ind = expand.grid(1:A,1:A)
pars=paste("M1[6,",ind[,1],",",ind[,2],"]",sep="")
foo = matrix(summary(fitcomb.year,pars)$summary[,1],ncol=A)
image.plot(2+1:A,2+1:A,foo,xlab="True age",ylab="Read age")

ind = expand.grid(1:14,1:A)
pars=paste("M1[",ind[,1],",",ind[,2],",",ind[,2],"]",sep="")
foo = matrix(summary(fitcomb0.year,pars)$summary[,1],ncol=A)
matplot(t(foo),type="l")
##Extract probabilities for correct reading
Mdiag = NULL
for(r in 1:R)



ind = expand.grid(1:6,1:A)
pars=paste("xi1[",ind[,1],",",ind[,2],"]",sep="")
xi1year1 = matrix(summary(fitcombnoErr.year,pars)$summary[,1],ncol=A)
matplot(t(xi1year1),type="l")

ind = expand.grid(7:11,1:A)
pars=paste("xi1[",ind[,1],",",ind[,2],"]",sep="")
xi1year2 = matrix(summary(fitcombnoErr.year,pars)$summary[,1],ncol=A)
matplot(t(xi1year2),type="l")

ind = expand.grid(12:15,1:A)
pars=paste("xi1[",ind[,1],",",ind[,2],"]",sep="")
xi1year3 = matrix(summary(fitcombnoErr.year,pars)$summary[,1],ncol=A)
matplot(t(xi1year3),type="l")


ind = expand.grid(1,1:A)
pars=c("xi0")
pars=paste("xi0[",ind[,1],",",ind[,2],"]",sep="")
plot(fitcomb.year,pars=pars)

pdf("alpha0.pdf",height=10,width=8)
pars=c("alpha0")
plot(fitcomb.year,pars=pars)
dev.off()

pdf("beta0.pdf",height=10,width=8)
pars=c("beta0")
plot(fitcomb.year,pars=pars)
dev.off()

pars=c("xi1[61,4]")
traceplot(fitcomb.year,pars=pars)

pars=c("lp__")
traceplot(fitcomb.year,pars=pars)

ind = expand.grid(1,1:A)
pars2=paste("xi1[",ind[,1],",",ind[,2],"]",sep="")
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

