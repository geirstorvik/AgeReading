library(AgeHerring)
library(rstan)
library(ggplot2)
require(gridExtra)
library(loo)
library(bayesplot)

Amin=3;Amax=12
A = Amax-Amin+1
R = 14

dir = "../save_results/"
load(paste0(dir,"fitcomb0.year100000.RData"))              #reader dependent age errors
load(paste0(dir,"fitcomb1.year100000.RData"))             #reader independent age errors
load(paste0(dir,"fitcomb2.year100000.RData"))             #assuming modal age is true age.

#Harmonic mean for marginal likelihood
foo0 = as.vector(extract(fitcomb0.year)$loglik)
ML0 = max(foo0)-log(mean(exp(-foo0+max(foo0))))
foo1 = as.vector(extract(fitcomb1.year)$loglik)
ML1 = max(foo1)-log(mean(exp(-foo1+max(foo1))))
foo2 = extract(fitcomb2.year)$lp__
ML2 = max(foo2)-log(mean(exp(-foo2+max(foo2))))

##Calculating number at age
data(Herring_D_boot)
data(HerringSurvey)

Bootmean = Herring_D_boot[Herring_D_boot$BootstrapID<2,]
nunit = length(unique(Bootmean$unitstrat))
#Only using mean of bootstrap samples
for(i in 1:nunit)
{
  Bootmean$Abundance[Bootmean$unitstrat==i] = mean(Herring_D_boot$Abundance[Herring_D_boot$unitstrat==i])
  Bootmean$Biomass[Bootmean$unitstrat==i] = mean(Herring_D_boot$Biomass[Herring_D_boot$unitstrat==i])
}


foo = rstan::extract(fitcomb0.year)$xi1
res0 = Make_Numfish_sample(foo,Herring_D_boot,HerringSurvey$stratayearindex)
res0mean = Make_Numfish_sample(foo,Bootmean,HerringSurvey$stratayearindex)
foo2 = array(apply(foo,2:3,mean),c(1,15,10))
res0mean2 = Make_Numfish_sample(foo2,Herring_D_boot,HerringSurvey$stratayearindex)

foo = rstan::extract(fitcomb1.year)$xi1
res1 = Make_Numfish_sample(foo,Herring_D_boot,HerringSurvey$stratayearindex)
res1mean = Make_Numfish_sample(foo,Bootmean,HerringSurvey$stratayearindex)
foo2 = array(apply(foo,2:3,mean),c(1,15,10))
res1mean2 = Make_Numfish_sample(foo2,Herring_D_boot,HerringSurvey$stratayearindex)

foo = rstan::extract(fitcomb2.year)$xi1
res2 = Make_Numfish_sample(foo,Herring_D_boot,HerringSurvey$stratayearindex)
res2mean = Make_Numfish_sample(foo,Bootmean,HerringSurvey$stratayearindex)
foo2 = array(apply(foo,2:3,mean),c(1,15,10))
res2mean2 = Make_Numfish_sample(foo2,Herring_D_boot,HerringSurvey$stratayearindex)



##Rename xi0 so ages become right
library(stringr)
ind = grep("xi0",names(fitcomb0.year))
foo = names(fitcomb0.year)[ind]
for(i in 10:1)
{
  foo = str_replace_all(foo,paste0(",",i,"]"),paste0(",",i+2,"]"))
}
names(fitcomb0.year)[ind] = foo
foo2 = list()
foo3 = substring(foo,first=5,last=5)
foo4 = substring(foo,first=7,last=8)
foo4[1:15] = substring(foo4[1:15],first=1,last=1)
for(i in 1:length(foo))
  foo2[[i]] = substitute(expression(x[i,j]^0),list(i=foo3[i],j=foo4[i]))

##Rename M1 so ages become right
ind = grep("M1",names(fitcomb0.year))
ind2 = matrix(nrow=length(ind),ncol=3)
i=0
for(a in 1:A)
  for(a2 in 1:A)
    for(r in 1:R)
    {
      i = i+1
      ind2[i,] = c(r,a2+Amin-1,a+Amin-1)
    }
names(fitcomb0.year)[ind] = paste0("M1[",ind2[,1],",",ind2[,2],",",ind2[,3],"]")

M1 = extract(fitcomb0.year,pars="M1")$M1
M1mean = apply(M1,c(2,3,4),mean)
##Heatmat av M matriser
for(r in c(1,3,5,6))
{
  print(r)
  z_values <- round(t(M1mean[r, , ]), 3)
  df <- expand.grid(x = 2 + 1:A, y = 2 + 1:A)
  df$Prob <- c(t(z_values))
  pdf(paste0("Heatmap_M",r,".pdf"),height=5,width=5)
  ggplot(df, aes(x = x, y = y, fill = Prob)) +
    geom_tile() +
    scale_fill_gradient(low = "white", high = "black") +
    theme_minimal() +
    labs(title = "", x = "True age", y = "Read age", fill = "Prob") +
    scale_x_continuous(breaks = seq(min(df$x), max(df$x), by = 2)) +
    scale_y_continuous(breaks = seq(min(df$y), max(df$y), by = 2))
  dev.off()
}


pdf("fit0_alpha0.pdf",height=6,width=5)
pars="alpha0"
xx =as.matrix(fitcomb0.year,pars=pars)
mcmc_intervals(xx,prob=0.8,prob_outer=0.95)+scale_y_discrete(labels=3:12)+labs(y="Reader",x="Proportion")+
  theme(axis.text.x = element_text(size = 14),axis.text.y = element_text(size = 14),
        axis.title.x=element_text(size=22),axis.title.y=element_text(size=22))
#plot(fitcomb0.year,pars=pars)
dev.off()

pdf("fit0_beta0.pdf",height=6,width=5)
pars="beta0"
xx =as.matrix(fitcomb0.year,pars=pars)
mcmc_intervals(xx,prob=0.8,prob_outer=0.95)+scale_y_discrete(labels=3:12)+labs(y="Reader",x="Proportion")+
  theme(axis.text.x = element_text(size = 14),axis.text.y = element_text(size = 14),
        axis.title.x=element_text(size=22),axis.title.y=element_text(size=22))
#plot(fitcomb0.year,pars=pars)
dev.off()

pdf("fit0_alpha1.pdf",height=6,width=5)
pars="alpha1"
xx =as.matrix(fitcomb0.year,pars=pars)
mcmc_intervals(xx,prob=0.8,prob_outer=0.95)+scale_y_discrete(labels=3:12)+labs(y="Reader",x="Proportion")+
  theme(axis.text.x = element_text(size = 14),axis.text.y = element_text(size = 14),
        axis.title.x=element_text(size=22),axis.title.y=element_text(size=22))
#plot(fitcomb0.year,pars=pars)
dev.off()

pdf("fit0_beta1.pdf",height=6,width=5)
pars="beta1"
xx =as.matrix(fitcomb0.year,pars=pars)
mcmc_intervals(xx,prob=0.8,prob_outer=0.95)+scale_y_discrete(labels=3:12)+labs(y="Reader",x="Proportion")+
  theme(axis.text.x = element_text(size = 14),axis.text.y = element_text(size = 14),
        axis.title.x=element_text(size=22),axis.title.y=element_text(size=22))
#plot(fitcomb0.year,pars=pars)
dev.off()

pdf("fit0_trphi.pdf",height=6,width=5)
pars="trphi"
xx =as.matrix(fitcomb0.year,pars=pars)
mcmc_intervals(xx,prob=0.8,prob_outer=0.95)+scale_y_discrete(labels=3:12)+labs(y="Reader",x="Proportion")+
  theme(axis.text.x = element_text(size = 14),axis.text.y = element_text(size = 14),
        axis.title.x=element_text(size=22),axis.title.y=element_text(size=22))
#plot(fitcomb0.year,pars=pars)
dev.off()

pdf("fit0_global.pdf",height=6,width=5)
pars = c("Palpha0","Palpha1","Pbeta0","Pbeta1","Ptrphi")
xx =as.matrix(fitcomb0.year,pars=pars)
mcmc_intervals(xx,prob=0.8,prob_outer=0.95)+scale_y_discrete(labels=c(expression(alpha[0]),expression(alpha[1]),expression(beta[0]),expression(beta[1]),"Ptrphi"))+
  theme(axis.text.x = element_text(size = 14),axis.text.y = element_text(size = 14),
        axis.title.x=element_text(size=22),axis.title.y=element_text(size=22))
#plot(fitcomb0.year,pars=pars)
dev.off()

pdf("fit0_xi0_2021.pdf",height=7,width=4)
pars = paste0("xi0[1,",2+1:A,"]")
xx =as.matrix(fitcomb0.year,pars=pars)
mcmc_intervals(xx,prob=0.8,prob_outer=0.95)+scale_y_discrete(labels=3:12)+labs(y="Age",x="Proportion")+
  theme(axis.text.x = element_text(size = 14),axis.text.y = element_text(size = 14),
        axis.title.x=element_text(size=22),axis.title.y=element_text(size=22))
#plot(fitcomb0.year,pars=pars)+labs(y="Age",x="Proportion")
dev.off()

pdf("fit0_xi0_2022.pdf",height=7,width=4)
xx =as.matrix(fitcomb0.year,pars=pars)
mcmc_intervals(xx,prob=0.8,prob_outer=0.95)+scale_y_discrete(labels=3:12)+labs(y="Age",x="Proportion")+
  theme(axis.text.x = element_text(size = 14),axis.text.y = element_text(size = 14),
        axis.title.x=element_text(size=22),axis.title.y=element_text(size=22))
#pars = paste0("xi0[2,",2+1:A,"]")
plot(fitcomb0.year,pars=pars)
dev.off()

pdf("fit0_xi0_2023.pdf",height=7,width=4)
pars = paste0("xi0[3,",2+1:A,"]")
xx =as.matrix(fitcomb0.year,pars=pars)
mcmc_intervals(xx,prob=0.8,prob_outer=0.95)+scale_y_discrete(labels=3:12)+labs(y="Age",x="Proportion")+
  theme(axis.text.x = element_text(size = 14),axis.text.y = element_text(size = 14),
        axis.title.x=element_text(size=22),axis.title.y=element_text(size=22))
#plot(fitcomb0.year,pars=pars)
dev.off()

pdf("fit0_M1_11_11.pdf",height=7,width=4)
ind = 1:R
pars = paste0("M1[",ind,",11,11]")
xx =as.matrix(fitcomb0.year,pars=pars)
mcmc_intervals(xx,prob=0.8,prob_outer=0.95)+scale_y_discrete(labels=1:14)+labs(y="Reader",x="Proportion")+
  theme(axis.text.x = element_text(size = 14),axis.text.y = element_text(size = 14),
        axis.title.x=element_text(size=22),axis.title.y=element_text(size=22))
#plot(fitcomb0.year,pars=pars)
dev.off()

pdf("fit0_M1_45.pdf",height=5,width=4)
ind = 1:R
pars = paste0("M1[",ind,",4,5]")
xx =as.matrix(fitcomb0.year,pars=pars)
mcmc_intervals(xx,prob=0.8,prob_outer=0.95)+scale_y_discrete(labels=1:14)+labs(y="Reader",x="Proportion")+
  theme(axis.text.x = element_text(size = 14),axis.text.y = element_text(size = 14),
        axis.title.x=element_text(size=22),axis.title.y=element_text(size=22))
#plot(fitcomb0.year,pars=pars)
dev.off()

pdf("fit0_M1_76.pdf",height=5,width=4)
ind = 1:R
pars = paste0("M1[",ind,",7,6]")
xx =as.matrix(fitcomb0.year,pars=pars)
mcmc_intervals(xx,prob=0.8,prob_outer=0.95)+scale_y_discrete(labels=1:14)+labs(y="Reader",x="Proportion")+
  theme(axis.text.x = element_text(size = 14),axis.text.y = element_text(size = 14),
        axis.title.x=element_text(size=22),axis.title.y=element_text(size=22))
#plot(fitcomb0.year,pars=pars)
dev.off()

Rp = 14
M1mean = apply(rstan::extract(fitcomb0.year,pars="M1")$M1,c(2,3,4),mean)
M1L = apply(rstan::extract(fitcomb0.year,pars="M1")$M1,c(2,3,4),quantile,probs=0.025)
M1U = apply(rstan::extract(fitcomb0.year,pars="M1")$M1,c(2,3,4),quantile,probs=0.975)
DM = array(NA,c(Rp,3,A))
for(i in 1:Rp)
{
  DM[i,1,] = diag(M1mean[i,,])
  DM[i,2,] = diag(M1L[i,,])
  DM[i,3,] = diag(M1U[i,,])
}
pdf("fit0_M1diag.pdf",height=7,width=6)
matplot(matrix(2+rep(1:A,Rp),nrow=A),t(DM[,1,]),
        type="l",lwd=1.75,lty=c(rep(2,6),rep(1,8)),ylim=c(0.5,1),
        xlab="Age",ylab="Probability of correct age reading",cex.lab=1.5)
legend("bottomleft",paste0("Reader ",1:Rp),lwd=1.75,lty=c(rep(2,6),rep(1,8)),col=1:14)

dev.off()

foo = extract(fitcomb0.year,pars="M1")$M1
x = NULL
cl = NULL
M = dim(foo)[1]
for(r in 1:14)
{
  x = c(x,foo[,r,9,9])
  cl = c(cl,rep(r,M))
}
df = data.frame(x=x,reader=cl)
df$reader = as.factor(df$reader)
pdf("fit0_M1_11_11.pdf",height=5,width=6)
ggplot(df,aes(x=x,color=reader,fill=reader)) +
  geom_density(alpha=0.5,show.legend=TRUE) + scale_fill_viridis_d() +
  ggtitle(paste0("Age=",11)) + xlab('') + ylab('')  + theme_minimal()
dev.off()

matplot(matrix(2+rep(1:A,3*Rp),nrow=A),cbind(t(DM[,1,]),t(DM[,2,]),t(DM[,3,])),
        type="l",lty=c(rep(1,Rp),rep(2,2*Rp)),ylim=c(0,1),
        xlab="Age",ylab="Prob of correct age reading")

## Test if diagonal elements are different
foo = extract(fitcomb0.year,pars="M1")$M1
RT = array(1,c(A,R,R))
for(a in 1:A)
{
  for(r in 1:(R-1))
    for(r2 in (r+1):R)
    {
       p = mean(foo[,r,a,a]>foo[,r2,a,a])
       RT[a,r,r2] = RT[a,r2,r] = min(p,1-p)
    }
}


for(y in 1:3)
{
  l = list()
  for(a in 1:A)
  {
    df = data.frame(y=c(res0[,y,a],res1[,y,a],res2[,y,a]),
                    cl=c(rep(1,dim(res0)[1]),rep(2,dim(res1)[1]),
                         rep(3,dim(res2)[1])))
    df$y = df$y/10^9
    df$cl = as.factor(df$cl)
    l[[a]] = ggplot(df,aes(x=y,color=cl,fill=cl)) +
      geom_density(alpha=0.5,show.legend=FALSE) +
      ggtitle(paste0("Age=",a+Amin-1)) + xlab('') + ylab('')
  }
  pdf(paste0("Nest_202",y,".pdf"),height=5,width=10)
  grid.arrange(l[[1]],l[[2]],l[[3]],l[[4]],
               l[[5]],l[[6]],l[[7]],l[[8]],
               l[[9]],l[[10]],ncol=5,nrow=2)
  dev.off()
}

#Only using mean of bootstrap samples
for(y in 1:3)
{
  l = list()
  for(a in 1:A)
  {
    df = data.frame(y=c(res0mean[,y,a],res1mean[,y,a],res2mean[,y,a]),
                    cl=c(rep(1,dim(res0mean)[1]),rep(2,dim(res1mean)[1]),
                         rep(3,dim(res2mean)[1])))
    df$y = df$y/10^9
    df$cl = as.factor(df$cl)
    l[[a]] = ggplot(df,aes(x=y,color=cl,fill=cl)) +
      geom_density(alpha=0.5,show.legend=FALSE) +
      ggtitle(paste0("Age=",a+Amin-1)) + xlab('') + ylab('')
  }
  pdf(paste0("Nest_mean_202",y,".pdf"),height=5,width=10)
  grid.arrange(l[[1]],l[[2]],l[[3]],l[[4]],
               l[[5]],l[[6]],l[[7]],l[[8]],
               l[[9]],l[[10]],ncol=5,nrow=2)
  dev.off()
}

#Only using mean of xi1
for(y in 1:3)
{
  l = list()
  for(a in 1:A)
  {
    df = data.frame(y=c(res0mean2[,y,a],res1mean2[,y,a],res2mean2[,y,a]),
                    cl=c(rep(1,dim(res0mean2)[1]),rep(2,dim(res1mean2)[1]),
                         rep(3,dim(res2mean2)[1])))
    df$y = df$y/10^9
    df$cl = as.factor(df$cl)
    l[[a]] = ggplot(df,aes(x=y,color=cl,fill=cl)) +
      geom_density(alpha=0.5,show.legend=FALSE) +
      ggtitle(paste0("Age=",a+Amin-1)) + xlab('') + ylab('')
  }
  pdf(paste0("Nest_mean2_202",y,".pdf"),height=5,width=10)
  grid.arrange(l[[1]],l[[2]],l[[3]],l[[4]],
               l[[5]],l[[6]],l[[7]],l[[8]],
               l[[9]],l[[10]],ncol=5,nrow=2)
  dev.off()
}

## Comparing different variances
for(y in 1:3)
{
  l = list()
  for(a in 1:A)
  {
    df = data.frame(y=c(res0[,y,a],res1[,y,a],res2[,y,a]),
                    cl=c(rep(1,dim(res0)[1]),
                         rep(2,dim(res1)[1]),rep(3,dim(res2)[1])))
    df$y = df$y/10^9
    df$cl = as.factor(df$cl)
    l[[a]] = ggplot(df,aes(x=y,color=cl,fill=cl)) +
      geom_density(alpha=0.5,show.legend=FALSE) +
      ggtitle(paste0("Age=",a+Amin-1)) + xlab('') + ylab('')
  }
  pdf(paste0("Nest_meanvar_202",y,".pdf"),height=5,width=10)
  grid.arrange(l[[1]],l[[2]],l[[3]],l[[4]],
               l[[5]],l[[6]],l[[7]],l[[8]],
               l[[9]],l[[10]],ncol=5,nrow=2)
  dev.off()
}

l = list()
ind = 0
for(y in 1:3)
{
  for(a in (1+y):(3+y))
  {
    ind = ind + 1
    df = data.frame(y=c(res0[,y,a],res1[,y,a],res2[,y,a]),
                    cl=c(rep(1,dim(res0)[1]),
                         rep(2,dim(res1)[1]),rep(3,dim(res2)[1])))
    df$y = df$y/10^9
    df$cl = as.factor(df$cl)
    l[[ind]] = ggplot(df,aes(x=y,color=cl,fill=cl)) +
      geom_density(alpha=0.5,show.legend=FALSE) +
      ggtitle(paste0("Age=",a+Amin-1)) + xlab('') + ylab('')
  }
}
pdf(paste0("Nest_meanvar_allyear.pdf"),height=5,width=10)
grid.arrange(l[[1]],l[[2]],l[[3]],
             l[[4]],l[[5]],l[[6]],
             l[[7]],l[[8]],l[[9]],ncol=3,nrow=3)
dev.off()


##Comparing with no uncertainty in D
for(y in 1:3)
{
  l = list()
  for(a in 1:A)
  {
    df = data.frame(y=c(res0[,y,a],res0mean[,y,a]),
                    cl=c(rep(1,dim(res0)[1]),
                         rep(2,dim(res0mean)[1])))
   df$y = df$y/10^9
    df$cl = as.factor(df$cl)
    l[[a]] = ggplot(df,aes(x=y,color=cl,fill=cl)) +
      geom_density(alpha=0.5,show.legend=FALSE) +
      ggtitle(paste0("Age=",a+Amin-1)) + xlab('') + ylab('')
  }
  pdf(paste0("Nest_novarD_202",y,".pdf"),height=5,width=10)
  grid.arrange(l[[1]],l[[2]],l[[3]],l[[4]],
               l[[5]],l[[6]],l[[7]],l[[8]],
               l[[9]],l[[10]],ncol=5,nrow=2)
  dev.off()
}


l = list()
ind=0
for(y in 1:3)
{
  for(a in (1+y):(3+y))
  {
    ind = ind+1
  df = data.frame(y=c(res0[,y,a],res0mean[,y,a]),
                    cl=c(rep(1,dim(res0)[1]),
                         rep(2,dim(res0mean)[1])))
    df$y = df$y/10^9
    df$cl = as.factor(df$cl)
    l[[ind]] = ggplot(df,aes(x=y,color=cl,fill=cl)) +
      geom_density(alpha=0.5,show.legend=FALSE) +
      ggtitle(paste0("Age=",a+Amin-1)) + xlab('') + ylab('')
  }}
pdf(paste0("Nest_nvarD_allyear.pdf"),height=5,width=10)
grid.arrange(l[[1]],l[[2]],l[[3]],l[[4]],l[[5]],l[[6]],l[[7]],l[[8]],l[[9]],ncol=3,nrow=3)
dev.off()


for(y in 1:3)
{
  l = list()
  for(a in 1:A)
  {
    df = data.frame(y=c(res0[,y,a],res1[,y,a]),
                    cl=c(rep(1,dim(res0)[1]),
                         rep(3,dim(res1)[1])))
    df$y = df$y/10^9
    df$cl = as.factor(df$cl)
    l[[a]] = ggplot(df,aes(x=y,color=cl,fill=cl)) +
      geom_density(alpha=0.5,show.legend=FALSE) +
      ggtitle(paste0("Age=",a+Amin-1)) + xlab('') + ylab('')
  }
  grid.arrange(l[[1]],l[[2]],l[[3]],l[[4]],
               l[[5]],l[[6]],l[[7]],l[[8]],
               l[[9]],l[[10]],ncol=5,nrow=2)
}

for(y in 1:3)
{
  dx = NULL
  dy = NULL
  for(a in 1:A)
  {
    print(c(y,a))
    dx = rbind(dx,rep(a,9))
    dy = rbind(dy,c(quantile(res0[,y,a]/10^9,probs=c(0.25,0.5,0.75)),
              quantile(res0mean[,y,a]/10^9,probs=c(0.25,0.5,0.75)),
              quantile(res1[,y,a]/10^9,probs=c(0.25,0.5,0.75))))
  }
  pdf(paste0("Nest_meanvar_202",y,".pdf"),height=5,width=10)
  grid.arrange(l[[1]],l[[2]],l[[3]],l[[4]],
               l[[5]],l[[6]],l[[7]],l[[8]],
               l[[9]],l[[10]],ncol=5,nrow=2)
  dev.off()
}



## Comparing diagonal elements
foo = rstan::extract(fitcomb0.year)$M1
l = list()
for(a in 2:A)
{
  y=NULL;cl=NULL
  for(r in (R1+1):R)
  {
    y = c(y,foo[,r,a,a-1])
    cl = c(cl,rep(r,dim(foo)[1]))
  }
 df  = data.frame(y=y,cl=as.factor(cl))
  l[[a]] = ggplot(df,aes(x=y,color=cl,fill=cl)) +
    geom_density(alpha=0.5,show.legend=FALSE) +
    ggtitle(paste0("Age=",a+Amin-1)) + xlab('') + ylab('')
}
pdf("Mdiag.pdf",height=5,width=10)
grid.arrange(l[[1]],l[[2]],l[[3]],l[[4]],
             l[[5]],l[[6]],l[[7]],l[[8]],
             l[[9]],l[[10]],ncol=5,nrow=2)
dev.off()
grid.arrange(l[[2]],l[[3]],l[[4]],
             l[[5]],l[[6]],l[[7]],l[[8]],
             l[[9]],l[[10]],ncol=5,nrow=2)


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

Dya = matrix(nrow=Y,ncol=A)
eps = rep(NA,A)
#par(mfrow=c(2,5))
for(a in 1:A){
  eps[a] = mean(res0[,,a])*0.1
  for(y in 1:3)
  {
    Dya[y,a] = mean(abs(res1[,y,a]-res0[,y,a])>eps[a])
    if(Dya[y,a]>10.0)
    {
      xlim = range(c(res0[,y,a],res1[,y,a]))
      hist(res0[,y,a],xlim=xlim,,breaks=30,
           xlab="Numbers at age",main=paste("Year=",2020+y,"Age=",a+2),
           probability=TRUE,col=rgb(1,0,0,1/4))
      hist(res1[,y,a],xlim=xlim,breaks=30,
           probability=TRUE,col=rgb(0,1,1,1/4),add=T)
      hist(res2[,y,a],xlim=xlim,breaks=30,
           probability=TRUE,col=rgb(0,0,1,1/4),add=T)
    }
  }
}

matplot(cbind(2+1:A,2+1:A,2+1:A),t(Dya),type="l",lwd=2,lty=1,col=1:3,ylim=c(0,1),
        xlab="Ages",ylab="P-value")
abline(h=0.05,lty=3,col=4)

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


pars = c("Palpha0","Pbeta0","Pbeta1","Pphi")
plot(fitcomb.year,pars=pars)

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
foo = matrix(summary(fitcomb.year,pars)$summary[,1],ncol=A)


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


##Compare models
library(bridgesampling)
M2 =  bridge_sampler(fitcomb2.year,silent=TRUE)

##Comparing models with WAIC
loglik0 = extract_log_lik(fitcomb0.year,parameter_name="lp__")
waic0 = loo::waic(loglik0)
loglik1 = extract_log_lik(fitcomb1.year,parameter_name="lp__")
waic1 = loo::waic(loglik1)
loglik2 = extract_log_lik(fitcomb2.year,parameter_name="lp__")
waic2 = loo::waic(loglik2)


