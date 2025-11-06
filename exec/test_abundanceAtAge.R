library(AgeHerring)

#herringSet <- data.table::fread(paste(dir,"Data/agereadingset_herring_spawning_survey.csv",sep=""))
data(herringSet)
agerange <- 0:max(herringSet$age, na.rm = T)
ageAbundance <- sapply(agerange, function(x){abundanceAtAge(herringSet, 2021, x)})
names(ageAbundance) <- agerange

meanAge <- sum(agerange %*% ageAbundance) / sum(ageAbundance)

print(ageAbundance)
plot(agerange,ageAbundance)

foo = abundanceAtReader(herringSet,2022)
reader = unique(herringSet$AGEREADER)
reader = sort(reader[!is.na(reader)])
par(mfrow=c(2,3))
x = NULL
y = NULL
for(r in reader)
{
  dataH = herringSet[AGEREADER==r]
  tab = table(dataH$serialnumber)
  ##Change agesamplecount to the numbers related to the current reader
  for(i in 1:length(tab))
  {
    dataH$agesamplecount[dataH$serialnumber==names(tab)[i]]=tab[i]
  }
  show(dim(dataH))
  ageAbundance <- sapply(agerange, function(x){abundanceAtAge(dataH, 2022, x)})
  names(ageAbundance) <- agerange
  x = cbind(x,agerange)
  y = cbind(y,ageAbundance)
  plot(agerange,ageAbundance)

}
par(mfrow=c(1,1))
matplot(x,y,type="l",col=1:5,lty=1,lwd=2,xlab="age",ylab="Mean abundance")
legend("topright",legend=reader,col=1:5,lty=1,lwd=2)
