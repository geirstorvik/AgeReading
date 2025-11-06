library(xtable)
foo1 = HerringSet$age[HerringSet$startyear==2021]
foo1 = pmax(3,pmin(foo1,12))
foo2 = HerringSet$age[HerringSet$startyear==2022]
foo2 = pmax(3,pmin(foo2,12))
foo3 = HerringSet$age[HerringSet$startyear==2023]
foo3 = pmax(3,pmin(foo3,12))

tab = cbind(table(c(3:12,foo1))-1,
            table(c(3:12,foo2))-1,
            table(c(3:12,foo3))-1)

for(i in 1:3)
  tab[,i] = tab[,i]/sum(tab[,i])

age = cbind(3:12,3:12,3:12)
matplot(age,tab,type="l",lwd=2,lty=1,col=1:3,xlab="Age",ylab="Proportion",cex.axis=1.5,cex.lab=1.5)
legend("topright",legend=2021:2023,lwd=2,lty=1,col=1:3,cex=1.5)

foo = table(HerringSet$startyear,substring(HerringSet$SPATIAL_STRATUM,first=5))
xtable(foo[,c(3:6,1:2)],digits=0)

foo2 = table(HerringSet$startyear,HerringSet$AGEREADER)
colnames(foo2) = 1:8
xtable(foo2)

foo3 = cbind(foo[,c(3:6,1:2)],foo2)

foo = table(HerringSet$SPATIAL_STRATUM,HerringSet$AGEREADER)
foo = foo[rowSums(foo)>0,]

foo = foo[c(3:6,1:2,7:15),]

foo = cbind(foo,rowSums(foo))
foo = rbind(foo,colSums(foo))
colnames(foo) = c(1:8,"Sum")
xtable(foo)
