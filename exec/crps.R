load("../save_results/crpsMbar.RData")
load("../save_results/crpsM1.RData")
load("../save_results/crpsM3.RData")
load("../save_results/crpsM6.RData")
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
df$year = as.factor(c(rep(2021,400),rep(2022,400),rep(2023,400)))
df$age = as.factor(df$age)
df$M = as.factor(df$M)
library(ggplot2)
df1 = df[df$year==2021,]
pdf("Sim_Mfixed_crps_2021.pdf",height=5,width=8)
ggplot(df1,aes(age,crps))+geom_boxplot(aes(fill=M))+
  theme(axis.text=element_text(size=20),axis.title=element_text(size=20))
dev.off()
df2 = df[df$year==2022,]
pdf("Sim_Mfixed_crps_2022.pdf",height=5,width=8)
ggplot(df2,aes(age,crps))+geom_boxplot(aes(fill=M))+
  theme(axis.text=element_text(size=20),axis.title=element_text(size=20))
dev.off()
df3 = df[df$year==2023,]
pdf("Sim_Mfixed_crps_2023.pdf",height=5,width=8)
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
