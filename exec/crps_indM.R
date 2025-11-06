## Analysis using stan
data(HerringSurvey)
data(Ring.readErr)
df = HerringSurvey


#Summarize results
library(scoringRules)
load("../save_results/fitcomb0.year100000.RData")              #reader dependent age errors
Amin=3;Amax=12
A = Amax-Amin+1
Y=length(unique(df$stratayearindex))
xi1 = apply(rstan::extract(fitcomb0.year,pars="xi1")$xi1,c(2,3),mean)
Mxi1 = array(xi1,c(1,15,10))
source("R/Make_Numfish_sample.R")
Truef = apply(Make_Numfish_sample(Mxi1,Herring_D_boot,HerringSurvey$stratayearindex),2:3,mean)
Nsim = 10
Nmcmc=100000

if(NULL)
{

  fM1 = array(NA,c(Nsim,3,3,10))
  crpsM1 = array(NA,c(Y,A,Nsim))
  for(m in 1:Nsim)
  {
    fil = paste0("../save_results/simfitcomb1M1_indM",m,".",Nmcmc,".RData")
    if(file.exists(fil))
    {
      print(m)
      load(fil)              #reader dependent age errors
      foo = rstan::extract(simfitcomb1.year)$xi1
      res0 = Make_Numfish_sample(foo,Herring_D_boot,HerringSurvey$stratayearindex)
      fM1[m,,,] = apply(res0,c(2,3),quantile,probs=c(0.025,0.5,0.975))
      for(y in 1:Y)
        crpsM1[y,,m] = crps_sample(Truef[y,],t(res0[,y,]))
    }
  }
  save(crpsM1,file="../save_results/crpsM1_indM.RData")


  fM3 = array(NA,c(Nsim,3,3,10))
  crpsM3 = array(NA,c(Y,A,Nsim))
  for(m in 1:Nsim)
  {
    fil = paste0("../save_results/simfitcomb1M3_indM",m,".",Nmcmc,".RData")
    if(file.exists(fil))
    {
      print(m)
      load(fil)              #reader dependent age errors
      foo = rstan::extract(simfitcomb1.year)$xi1
      res0 = Make_Numfish_sample(foo,Herring_D_boot,HerringSurvey$stratayearindex)
      fM3[m,,,] = apply(res0,c(2,3),quantile,probs=c(0.025,0.5,0.975))
      for(y in 1:Y)
        crpsM3[y,,m] = crps_sample(Truef[y,],t(res0[,y,]))
    }
  }
  save(crpsM3,file="../save_results/crpsM3_indM.RData")

  fM6 = array(NA,c(Nsim,3,3,10))
  crpsM6 = array(NA,c(Y,A,Nsim))
  for(m in 1:Nsim)
  {
    fil = paste0("../save_results/simfitcomb1M6_indM",m,".",Nmcmc,".RData")
    if(file.exists(fil))
    {
      print(m)
      load(fil)              #reader dependent age errors
      foo = rstan::extract(simfitcomb1.year)$xi1
      res0 = Make_Numfish_sample(foo,Herring_D_boot,HerringSurvey$stratayearindex)
      fM6[m,,,] = apply(res0,c(2,3),quantile,probs=c(0.025,0.5,0.975))
      for(y in 1:Y)
        crpsM6[y,,m] = crps_sample(Truef[y,],t(res0[,y,]))
    }
  }
  save(crpsM6,file="../save_results/crpsM6_indM.RData")

  fMbar = array(NA,c(Nsim,3,3,10))
  crpsMbar = array(NA,c(Y,A,Nsim))
  for(m in 1:Nsim)
  {
    fil = paste0("../save_results/simfitcomb1Mbar_indM",m,".",format(Nmcmc,scientific=F),".RData")
    if(file.exists(fil))
    {
      print(m)
      load(fil)              #reader dependent age errors
      foo = rstan::extract(simfitcomb1.year)$xi1
      res0 = Make_Numfish_sample(foo,Herring_D_boot,HerringSurvey$stratayearindex)
      fMbar[m,,,] = apply(res0,c(2,3),quantile,probs=c(0.025,0.5,0.975))
      for(y in 1:3)
        crpsMbar[y,,m] = crps_sample(Truef[y,],t(res0[,y,]))
    }
  }
  save(crpsMbar,file="../save_results/crpsMbar_indM.RData")
}


load("../save_results/crpsMbar_indM.RData")
load("../save_results/crpsM1_indM.RData")
load("../save_results/crpsM3_indM.RData")
load("../save_results/crpsM6_indM.RData")
for(y in 1:3)
{
  crpsMbar[y,,] = t(t(crpsMbar[y,,])/Truef[y,])
  crpsM1[y,,] = t(t(crpsM1[y,,])/Truef[y,])
  crpsM3[y,,] = t(t(crpsM3[y,,])/Truef[y,])
  crpsM6[y,,] = t(t(crpsM6[y,,])/Truef[y,])
}
crpscomb = rbind(crpsMbar[1,,],crpsM1[1,,],crpsM3[1,,],crpsM6[1,,],
                 crpsMbar[2,,],crpsM1[2,,],crpsM3[2,,],crpsM6[2,,],
                 crpsMbar[3,,],crpsM1[3,,],crpsM3[3,,],crpsM6[3,,])
df = data.frame(crps=as.vector(crpscomb),age=2+rep(rep(1:A,4),3),M=rep(rep(c("Mbar","M1","M3","M6"),each=A),3))
df$year = as.factor(rep(c(rep(2021,100),rep(2022,100),rep(2023,100)),4))
df$age = as.factor(df$age)
df$M = as.factor(df$M)
year = c(2021,2022,2023)
M = c("M1","M3","M6","Mbar")
age = 3:12
dfnew = NULL
for(y in 1:3)
  for(a in 1:A)
  {
    dfnew = rbind(dfnew,data.frame(crps=crpsMbar[y,a,],year=rep(year[y],Nsim),age=rep(age[a],Nsim),M=rep("Mbar",Nsim)))
    dfnew = rbind(dfnew,data.frame(crps=crpsM1[y,a,],year=rep(year[y],Nsim),age=rep(age[a],Nsim),M=rep("M1",Nsim)))
    dfnew = rbind(dfnew,data.frame(crps=crpsM3[y,a,],year=rep(year[y],Nsim),age=rep(age[a],Nsim),M=rep("M3",Nsim)))
    dfnew = rbind(dfnew,data.frame(crps=crpsM6[y,a,],year=rep(year[y],Nsim),age=rep(age[a],Nsim),M=rep("M6",Nsim)))
  }
df = dfnew
df$year = as.factor(df$year)
df$age = as.factor(df$age)
df$M = as.factor(df$M)
library(ggplot2)
df1 = df[df$year==2021,]
pdf("Sim_Mfixed_crps_2021_indM.pdf",height=5,width=8)
ggplot(df1,aes(age,crps))+geom_boxplot(aes(fill=M))+
  theme(axis.text=element_text(size=20),axis.title=element_text(size=20))
dev.off()
df2 = df[df$year==2022,]
pdf("Sim_Mfixed_crps_2022_indM.pdf",height=5,width=8)
ggplot(df2,aes(age,crps))+geom_boxplot(aes(fill=M))+
  theme(axis.text=element_text(size=20),axis.title=element_text(size=20))
dev.off()
df3 = df[df$year==2023,]
pdf("Sim_Mfixed_crps_2023_indM.pdf",height=5,width=8)
ggplot(df3,aes(age,crps))+geom_boxplot(aes(fill=M))+
  theme(axis.text=element_text(size=20),axis.title=element_text(size=20))
dev.off()

df = within(df,M<-relevel(M,ref="Mbar"))
contrasts(df) = contr.treatment()
fit1 = lm(log(crps)~age+M+year,data=df)
fit0 = lm(log(crps)~age+year,data=df)
fit2 = lm(log(crps)~age*M*year,data=df)
fit3 = lm(log(crps)~age*year,data=df)
anova(fit0,fit1,fit2,fit3)

df2 = within(df,age=relevel(age,ref=4),year=relevel(year,ref=2))
fit02 = lm(log(crps)~relevel(age,ref="4")+relevel(year,ref="2022"),data=df)
