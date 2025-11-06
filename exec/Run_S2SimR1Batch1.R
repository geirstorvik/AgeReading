#Script for analyzing age-readings from annual Herring spawning survey in February/March
#Only use data from 2022
library(AgeHerring)
library(rstan)
logit = function(p){log(p/(1-p))}
invlogit = function(x){1/(1+exp(-x))}

## Choose dataset  (HerringSurvey/HerringSurvey202?)
data(S2SimR1.Ring.readErr)
dataRing = S2SimR1.Ring.readErr
data(S2SimR1.HerringSurvey)
df = S2SimR1.HerringSurvey

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

S2SimR1fitcomb1.year = stan(file="exec/agereader_comb_year_strat1.stan",data=d.comb,
                    iter=100000,chains=4,thin=100)
save(S2SimR1fitcomb1.year,file="S2SimR1fitcomb1.year100000.RData")

l = list()
ind=0
sled = c(FALSE,FALSE,TRUE,FALSE,FALSE,TRUE,FALSE,FALSE,TRUE)
for(y in 1:3)
{
  for(a in (1+y):(3+y))
  {
    #df = data.frame(y=c(res0[,y,a],resR1[,y,a],resR3[,y,a],resR6[,y,a])/10^9,
    #                reader=rep(c(0,1,3,6),each=dim(resR1)[1]))
    #                  cl=c(rep(1,dim(res0)[1]),rep(2,dim(resR1)[1]),rep(2,dim(resR3)[1]),rep(2,dim(resR5)[1]),rep(3,dim(resR6)[1])))
    df = data.frame(y=c(resR1[,y,a],resR3[,y,a],resR6[,y,a])/10^9,
                    reader=rep(c(1,3,6),each=dim(resR1)[1]))
    df$reader = as.factor(df$reader)
    ind = ind+1
    l[[ind]] = ggplot(df,aes(x=y,color=reader,fill=reader)) +
      geom_density(alpha=0.5,show.legend=sled[ind]) +
      ggtitle(paste0("Year=",2020+y,",Age=",a+Amin-1)) + xlab('') + ylab('') +
      geom_vline(xintercept=Ntrue[y,a]/10^9)
  }
}
#pdf("NestlargeA_extrm_0_0_sim_true.pdf",height=5,width=7)
grid.arrange(l[[1]],l[[2]],l[[3]],l[[4]],l[[5]],l[[6]],l[[7]],l[[8]],l[[9]],ncol=3,nrow=3)
dev.off()
