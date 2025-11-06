library(AgeHerring)
library(rstan)
library(ggplot2)
require(gridExtra)

Amin=3;Amax=12
A = Amax-Amin+1

##Comparing estimated xi^1 with truth
load("fitcomb0.year100000.RData")              #reader dependent age errors
xi1 = apply(extract(fitcomb0.year,pars="xi1")$xi1,c(2,3),mean)

yearstrat = list(y1=c(1:6),y2=c(7:11),y3=c(12:15))
##CRSP:
data(Herring_D_boot)
library(verification)
K=15
crpsxi = array(NA,c(K,A,4))
Dmean = aggregate(Herring_D_boot$Abundance,list(Herring_D_boot$unitstrat),FUN=mean)$x
Ntrue = array(NA,c(3,A))
for(y in 1:3)
{
  Dy = Herring_D_boot[Herring_D_boot$Year==(2020+y),]
  Dmean = aggregate(Dy$Abundance,list(Dy$unitstrat),FUN=mean)$x
  for(a in 1:A)
   Ntrue[y,a] = sum(Dmean*xi1[yearstrat[[y]],a])
}

#Data: M_r: Model:M_r
data(Herring_D_boot)
data(HerringSurvey)
load("S2simfitcomb0.year100000.RData")              #reader dependent age errors
foo = rstan::extract(S2simfitcomb0.year)$xi1
res0 = Make_Numfish_sample(foo,Herring_D_boot,HerringSurvey$stratayearindex)

load("S2simfitcomb1.year100000.RData")              #reader dependent age errors
foo = rstan::extract(S2simfitcomb1.year)$xi1
res1 = Make_Numfish_sample(foo,Herring_D_boot,HerringSurvey$stratayearindex)

#load("S2sim2fitcomb1.year100000.RData")              #reader dependent age errors
#foo = rstan::extract(S2sim2fitcomb1.year)$xi1
#res2 = Make_Numfish_sample(foo,Herring_D_boot,HerringSurvey$stratayearindex)

load("S2SimR1fitcomb1.year100000.RData")              #reader dependent age errors
foo = rstan::extract(S2SimR1fitcomb1.year)$xi1
resR1 = Make_Numfish_sample(foo,Herring_D_boot,HerringSurvey$stratayearindex)

load("S2SimR3fitcomb1.year100000.RData")              #reader dependent age errors
foo = rstan::extract(S2SimR3fitcomb1.year)$xi1
resR3 = Make_Numfish_sample(foo,Herring_D_boot,HerringSurvey$stratayearindex)

load("S2SimR6fitcomb1.year100000.RData")              #reader dependent age errors
foo = rstan::extract(S2SimR6fitcomb1.year)$xi1
resR6 = Make_Numfish_sample(foo,Herring_D_boot,HerringSurvey$stratayearindex)


pdf("S2Sim_fit0_xi1est.pdf",height=5,width=10)
par(mfrow=c(5,3),mar=c(1,1,1,1),oma=c(4,5,0,0))
xi1.est = extract(S2simfitcomb0.year,pars="xi1")$xi1
for(k in 1:15)
{
  q = apply(xi1.est[,k,],2,quantile,probs=c(0.25,0.5,0.975))
  print(max(q))
  matplot(2+cbind(1:A,1:A,1:A),t(q),type="l",lty=c(2,1,2),col=2,xlab="",ylab="",lwd=1.5,ylim=c(0,0.85))
  lines(2+1:A,xi1[k,])
  crpsxi[k,,1] = crps(xi1[k,],t(xi1.est[,1,]))$crps
}
mtext("Age (years)",side=1,line=1.5,outer=TRUE,cex=1.3)
mtext(expression(xi[k]^1),side=2,line=2,outer=TRUE,cex=1.3,las=0)
dev.off()

pdf("S2Sim_fit0_xi1est_2.pdf",height=5,width=10)
l = list()
for(k in 1:15)
{
  q = apply(xi1.est[,k,],2,quantile,probs=c(0.25,0.5,0.975))
  df = data.frame(age=2+1:A,q1=q[1,],q2=q[2,],q3=q[3,],xi1true=xi1[k,])
  l[[k]] = ggplot(df,aes(age,q2)) + geom_line(colour="red")+
    geom_ribbon(aes(ymin=q1,ymax=q3),fill="red",alpha=0.2) +
    theme(axis.title.x=element_blank(),axis.title.y=element_blank()) + ylim(0,0.85) +
    geom_line(aes(y=xi1true),show.legend=FALSE)
}
grid.arrange(l[[1]],l[[2]],l[[3]],l[[4]],l[[5]],l[[6]],
             l[[7]],l[[8]],l[[9]],l[[10]],l[[11]],l[[12]],
             l[[13]],l[[14]],l[[15]],ncol=5,nrow=3)
dev.off()



#Data: M_r: Model:M
pdf("S2Sim_fit1_xi1est.pdf",height=5,width=10)
par(mfrow=c(5,3),mar=c(1,1,1,1),oma=c(4,5,0,0))
xi1.est = extract(S2simfitcomb1.year,pars="xi1")$xi1
for(k in 1:15)
{
  q = apply(xi1.est[,k,],2,quantile,probs=c(0.25,0.5,0.975))
  matplot(2+cbind(1:A,1:A,1:A),t(q),type="l",lty=c(2,1,2),col=2,xlab="",ylab="",lwd=1.5,ylim=c(0,0.85))
  lines(2+1:A,xi1[k,])
  crpsxi[k,,2] = crps(xi1[k,],t(xi1.est[,1,]))$crps
}
mtext("Age (years)",side=1,line=1.5,outer=TRUE,cex=1.3)
mtext(expression(xi[k]^1),side=2,line=2,outer=TRUE,cex=1.3,las=0)
dev.off()

#Data: M, Model: M_r
pdf("S2Sim2_fit0_xi1est.pdf",height=5,width=10)
par(mfrow=c(5,3),mar=c(1,1,1,1),oma=c(4,5,0,0))
xi1.est = extract(sim2fitcomb0.year,pars="xi1")$xi1
for(k in 1:15)
{
  q = apply(xi1.est[,k,],2,quantile,probs=c(0.25,0.5,0.975))
  matplot(2+cbind(1:A,1:A,1:A),t(q),type="l",lty=c(2,1,2),col=2,xlab="",ylab="",lwd=1.5,ylim=c(0,0.85))
  lines(2+1:A,xi1[k,])
  crpsxi[k,,3] = crps(xi1[k,],t(xi1.est[,1,]))$crps
}
mtext("Age (years)",side=1,line=1.5,outer=TRUE,cex=1.3)
mtext(expression(xi[k]^1),side=2,line=2,outer=TRUE,cex=1.3,las=0)
dev.off()

#Data: M, Model: M
pdf("S2Sim2_fit1_xi1est.pdf",height=5,width=10)
par(mfrow=c(5,3),mar=c(1,1,1,1),oma=c(4,5,0,0))
xi1.est = extract(sim2fitcomb1.year,pars="xi1")$xi1
for(k in 1:15)
{
  q = apply(xi1.est[,k,],2,quantile,probs=c(0.25,0.5,0.975))
  matplot(2+cbind(1:A,1:A,1:A),t(q),type="l",lty=c(2,1,2),col=2,xlab="",ylab="",lwd=1.5,ylim=c(0,0.85))
  lines(2+1:A,xi1[k,])
  crpsxi[k,,4] = crps(xi1[k,],t(xi1.est[,1,]))$crps
}
mtext("Age (years)",side=1,line=1.5,outer=TRUE,cex=1.3)
mtext(expression(xi[k]^1),side=2,line=2,outer=TRUE,cex=1.3,las=0)
dev.off()

##Compare estimated against "true" N(y,a)
f = apply(res0,c(2,3),quantile,probs=c(0.025,0.5,0.975))
f2 = apply(res2,c(2,3),quantile,probs=c(0.025,0.5,0.975))
l = list()
crpsN0 = matrix(nrow=3,ncol=A)
for(y in 1:3)
{
  #matplot(cbind(1:A,1:A,1:A,1:A),cbind(t(f[,y,]),Ntrue[y,]),type="l",lty=c(2,1,2,1),lwd=2,col=c(1,1,1,2))
  df = data.frame(age=c(1:A,1:A,1:A,1:A,1:A,1:A,1:A),
                  N=c(f[1,y,],f[2,y,],f[3,y,],f2[1,y,],f2[2,y,],f2[3,y,],Ntrue[y,]),
                  cl=c(rep(1,A),rep(2,A),rep(3,A),rep(4,A),rep(5,A),rep(6,A),rep(7,A)),
                  col=c(rep(1,3*A),rep(2,3*A),rep(3,A)),
                  lty=c(rep(2,A),rep(1,A),rep(2,A),rep(2,A),rep(1,A),rep(2,A),rep(1,A)))
  df$N = df$N/10^9
  df$cl = as.factor(df$cl)
  df$col = as.factor(df$col)
  df$lty = as.factor(df$lty)
  l[[y]] = ggplot(df,aes(x=age,y=N,group=cl,color=col)) + geom_line(aes(lty=lty),show.legend=TRUE)
    crpsN0[y,] = crps(Ntrue[y,],t(res0[,y,]))$crps

}
pdf("S2Nest_comb_0_0_sim_true.pdf",height=5,width=7)
grid.arrange(l[[1]],l[[2]],l[[3]],ncol=1,nrow=3)
dev.off()

##Compare against extreme readers
f3 = apply(resR3,c(2,3),quantile,probs=c(0.025,0.5,0.975))
f6 = apply(resR6,c(2,3),quantile,probs=c(0.025,0.5,0.975))
l = list()
crpsN0 = matrix(nrow=3,ncol=A)
for(y in 1:3)
{
  #matplot(cbind(1:A,1:A,1:A,1:A),cbind(t(f[,y,]),Ntrue[y,]),type="l",lty=c(2,1,2,1),lwd=2,col=c(1,1,1,2))
  df = data.frame(age=c(1:A,1:A,1:A,1:A,1:A,1:A,1:A,1:A,1:A,1:A),
                  N=c(f[1,y,],f[2,y,],f[3,y,],f3[1,y,],f3[2,y,],f3[3,y,],f6[1,y,],f6[2,y,],f6[3,y,],Ntrue[y,]),
                  cl=c(rep(1,A),rep(2,A),rep(3,A),rep(4,A),rep(5,A),rep(6,A),rep(7,A),rep(8,A),rep(9,A),rep(10,A)),
                  col=c(rep(1,3*A),rep(2,3*A),rep(3,3*A),rep(4,A)),
                  lty=c(rep(2,A),rep(1,A),rep(2,A),rep(2,A),rep(1,A),rep(2,A),rep(2,A),rep(1,A),rep(2,A),rep(1,A)))
  df$N = df$N/10^9
  df$cl = as.factor(df$cl)
  df$col = as.factor(df$col)
  df$lty = as.factor(df$lty)
  l[[y]] = ggplot(df,aes(x=age,y=N,group=cl,color=col)) + geom_line(aes(lty=lty),show.legend=TRUE)
  crpsN0[y,] = crps(Ntrue[y,],t(res0[,y,]))$crps

}
pdf("S2Nest_extrm_0_0_sim_true.pdf",height=5,width=7)
grid.arrange(l[[1]],l[[2]],l[[3]],ncol=1,nrow=3)
dev.off()

##Compare against extreme readers
f = apply(res0,c(2,3),quantile,probs=c(0.5))
f1 = apply(resR1,c(2,3),quantile,probs=c(0.5))
f3 = apply(resR3,c(2,3),quantile,probs=c(0.5))
f5 = apply(resR5,c(2,3),quantile,probs=c(0.5))
f6 = apply(resR6,c(2,3),quantile,probs=c(0.5))
l = list()
crpsN0 = array(NA,c(3,5,A))
for(y in 1:3)
{
  #matplot(cbind(1:A,1:A,1:A,1:A),cbind(t(f[,y,]),Ntrue[y,]),type="l",lty=c(2,1,2,1),lwd=2,col=c(1,1,1,2))
  df = data.frame(age=c(1:A,1:A,1:A,1:A,1:A,1:A,1:A,1:A,1:A,1:A),
                  N=c(f[y,],f1[y,],f3[y,],f6[y,],Ntrue[y,]),
                  reader=c(rep(0,A),rep(1,A),rep(3,A),rep(5,A),rep("TRUE",A)))
  df$N = df$N/10^9
  df$cl = as.factor(df$reader)
  l[[y]] = ggplot(df,aes(x=age,y=N,group=reader,color=reader)) + geom_line(show.legend=TRUE)
  crpsN0[y,1,] = crps(Ntrue[y,],t(res0[,y,]))$crps
  crpsN0[y,2,] = crps(Ntrue[y,],t(resR1[,y,]))$crps
  crpsN0[y,3,] = crps(Ntrue[y,],t(resR3[,y,]))$crps
  crpsN0[y,4,] = crps(Ntrue[y,],t(resR5[,y,]))$crps
  crpsN0[y,5,] = crps(Ntrue[y,],t(resR6[,y,]))$crps

}
pdf("S2Nest_extrm_0_0_sim_true.pdf",height=5,width=7)
grid.arrange(l[[1]],l[[2]],l[[3]],ncol=1,nrow=3)
dev.off()

##Similar with densities for largest agegroup
l = list()
ind=0
sled = c(FALSE,FALSE,TRUE,FALSE,FALSE,TRUE,FALSE,FALSE,TRUE)
for(y in 1:3)
{
  for(a in (1+y):(3+y))
  {
    #df = data.frame(y=c(res0[,y,a],resR1[,y,a],resR3[,y,a],resR3[,y,a],resR6[,y,a])/10^9,
    #                cl=c(rep(1,dim(res0)[1]),rep(2,dim(resR1)[1]),,rep(2,dim(resR3)[1]),,rep(2,dim(resR5)[1]),rep(3,dim(resR6)[1])))
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
pdf("S2NestlargeA_extrm_0_0_sim_true.pdf",height=5,width=7)
grid.arrange(l[[1]],l[[2]],l[[3]],l[[4]],l[[5]],l[[6]],l[[7]],l[[8]],l[[9]],ncol=3,nrow=3)
dev.off()

l = list()
ind=0
sled = c(FALSE,FALSE,TRUE,FALSE,FALSE,TRUE,FALSE,FALSE,TRUE)
for(y in 1:3)
{
  for(a in (4+y):(6+y))
  {
    #df = data.frame(y=c(res0[,y,a],resR1[,y,a],resR3[,y,a],resR3[,y,a],resR6[,y,a])/10^9,
    #                cl=c(rep(1,dim(res0)[1]),rep(2,dim(resR1)[1]),,rep(2,dim(resR3)[1]),,rep(2,dim(resR5)[1]),rep(3,dim(resR6)[1])))
    df = data.frame(y=c(res0[,y,a],resR1[,y,a],resR3[,y,a],resR6[,y,a])/10^9,
                    reader=rep(c(0,1,3,6),each=dim(resR1)[1]))
    df$reader = as.factor(df$reader)
    ind = ind+1
    l[[ind]] = ggplot(df,aes(x=y,color=reader,fill=reader)) +
      geom_density(alpha=0.5,show.legend=sled[ind]) +
      ggtitle(paste0("Year=",2020+y,",Age=",a+Amin-1)) + xlab('') + ylab('') +
      geom_vline(xintercept=Ntrue[y,a]/10^9)
  }
}
pdf("S2NestlargeA_extrm_0_0_sim_true_2.pdf",height=5,width=7)
grid.arrange(l[[1]],l[[2]],l[[3]],l[[4]],l[[5]],l[[6]],l[[7]],l[[8]],l[[9]],ncol=3,nrow=3)
dev.off()


z_values <- round(t(M1[6, , ]), 3)
df <- expand.grid(x = 2 + 1:A, y = 2 + 1:A)
df$Prob <- c(t(z_values))
pdf("Heatmap_M6.pdf",height=5,width=5)
ggplot(df, aes(x = x, y = y, fill = Prob)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "black") +
  theme_minimal() +
  labs(title = "", x = "True age", y = "Read age", fill = "Prob") +
  scale_x_continuous(breaks = seq(min(df$x), max(df$x), by = 2)) +
  scale_y_continuous(breaks = seq(min(df$y), max(df$y), by = 2))
dev.off()


crpsN1 = matrix(nrow=3,ncol=A)
f = apply(res2,c(2,3),quantile,probs=c(0.025,0.5,0.975))
l = list()
for(y in 1:3)
{
  #matplot(cbind(1:A,1:A,1:A,1:A),cbind(t(f[,y,]),Ntrue[y,]),type="l",lty=c(2,1,2,1),lwd=2,col=c(1,1,1,2))
  df = data.frame(age=c(1:A,1:A,1:A,1:A),N=c(f[1,y,],f[2,y,],f[3,y,],Ntrue[y,]),
                  cl=c(rep(1,A),rep(2,A),rep(3,A),rep(4,A)),
                  col=c(rep(1,3*A),rep(2,A)),lty=c(rep(2,A),rep(1,A),rep(2,A),rep(1,A)))
  df$N = df$N/10^9
  df$cl = as.factor(df$cl)
  df$col = as.factor(df$col)
  df$lty = as.factor(df$lty)
  l[[y]] = ggplot(df,aes(x=age,y=N,group=cl,color=col)) + geom_line(aes(lty=lty),show.legend=FALSE)
  crpsN1[y,] = crps(Ntrue[y,],t(res2[,y,]))$crps
}
pdf("S2Nest_1_1_sim_true.pdf",height=5,width=5)
grid.arrange(l[[1]],l[[2]],l[[3]],ncol=1,nrow=3)
dev.off()


## Compare two extrems
l = list()
ind = 0
for(y in 1:3)
{
  for(a in (1+y):(3+y))
  {
    ind = ind + 1
    df = data.frame(y=c(res1[,y,a],res3[,y,a],res5[,y,a],res6[,y,a]),
                    cl=c(rep(1,dim(res1)[1]),
                         rep(2,dim(res1)[1])))
    df$y = df$y/10^9
    df$cl = as.factor(df$cl)
    l[[ind]] = ggplot(df,aes(x=y,color=cl,fill=cl)) +
      geom_density(alpha=0.5,show.legend=FALSE) +
      ggtitle(paste0("Age=",a+Amin-1)) + xlab('') + ylab('')
  }
}
pdf(paste0("NestSim_compare_extrems.pdf"),height=5,width=10)
grid.arrange(l[[1]],l[[2]],l[[3]],
             l[[4]],l[[5]],l[[6]],
             l[[7]],l[[8]],l[[9]],ncol=3,nrow=3)
dev.off()

##Compare reader dependent data/model with reader independent data/model

l = list()
y=1
for(a in 1:A)
{
  df = data.frame(y=c(res0[,y,a],res1[,y,a],res2[,y,a])/10^9,
                  cl=c(rep(1,dim(res0)[1]),rep(2,dim(res1)[1]),rep(3,dim(res2)[1])))
  df$cl = as.factor(df$cl)
  l[[a]] = ggplot(df,aes(x=y,color=cl,fill=cl)) +
    geom_density(alpha=0.5,show.legend=TRUE) +
    ggtitle(paste0("Age=",a+Amin-1)) + xlab('') + ylab('')

}
grid.arrange(l[[1]],l[[2]],l[[3]],l[[4]],
             l[[5]],l[[6]],l[[7]],l[[8]],
             l[[9]],l[[10]],ncol=5,nrow=2)


pars="alpha0"
plot(S2simfitcomb0.year,pars=pars)




## Make density plots
for(y in 1:3)
{
  pdf(paste0("Density_year",2020+y,".pdf"),height=5,width=10)
  y=1
  l = list()
  for(a in 1:A)
  {
    df = data.frame(y=c(res0[,y,a],res1[,y,a],res2[,y,a]),
                    cl=c(rep(1,dim(res0)[1]),rep(2,dim(res1)[1]),
                         rep(3,dim(res2)[1])))
    df$cl = as.factor(df$cl)
    l[[a]] = ggplot(df,aes(x=y,color=cl,fill=cl)) +
      geom_density(alpha=0.5,show.legend=FALSE)
  }
  grid.arrange(l[[1]],l[[2]],l[[3]],l[[4]],
               l[[5]],l[[6]],l[[7]],l[[8]],
               l[[9]],l[[10]],ncol=5,nrow=2)
  dev.off()
}

