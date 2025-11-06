#Script for analyzing age-readings from annual Herring spawning survey in February/March
#Only use data from 2022
library(AgeHerring)
library(rstan)
logit = function(p){log(p/(1-p))}
invlogit = function(x){1/(1+exp(-x))}

## Choose dataset  (HerringSurvey/HerringSurvey202?)
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

#Extract true M
load("fitcomb0.year100000.RData")              #reader dependent age errors
M1 = apply(rstan::extract(fitcomb0.year,pars="M1")$M1,c(2,3,4),mean)
M1sim = M1[3,,]


#Number of strata
d.comb = list(Y=Y,S=df$S,K=df$K,U=df$U,
              R2=R2,A=Amax-Amin+1,
              D=df$d,N=df$N,
              stratayearindex=df$stratayearindex,
              stationstrataindex=df$stationstrataindex,
              unitstationindex=df$unitstationindex,
              unitreaderindex=df$unitreaderindex,
              alpharep=rep(alpha0,A),
              M1=M1sim)

SimR1fitcomb1.trueM.year = stan(file="exec/agereader_comb_year_strat1_trueM.stan",data=d.comb,
                    iter=100000,chains=4,thin=100)
save(SimR1fitcomb1.trueM.year,file="SimR1fitcomb1.trueM.year100000.RData")

foo = rstan::extract(SimR1fitcomb1.year)$xi1
resR1 = Make_Numfish_sample(foo,Herring_D_boot,HerringSurvey$stratayearindex)

xi1 = apply(rstan::extract(fitcomb0.year,pars="xi1")$xi1,c(2,3),mean)
Dmean = aggregate(Herring_D_boot$Abundance,list(Herring_D_boot$unitstrat),FUN=mean)$x
Ntrue = array(NA,c(3,A))
for(y in 1:3)
{
  Dy = Herring_D_boot[Herring_D_boot$Year==(2020+y),]
  Dmean = aggregate(Dy$Abundance,list(Dy$unitstrat),FUN=mean)$x
  for(a in 1:A)
    Ntrue[y,a] = sum(Dmean*xi1[yearstrat[[y]],a])
}

foo = rstan::extract(SimR1fitcomb1.trueM.year)$xi1
resR1 = Make_Numfish_sample(foo,Herring_D_boot,HerringSurvey$stratayearindex)

##Similar with densities for largest agegroup
l = list()
ind=0
sled = c(FALSE,FALSE,TRUE,FALSE,FALSE,TRUE,FALSE,FALSE,TRUE)
for(y in 1:3)
{
  for(a in (1+y):(3+y))
  {
    df = data.frame(y=c(resR1[,y,a])/10^9,
                    reader=rep(c(1),each=dim(resR1)[1]))
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
#dev.off()
