## Script for making Herring.surveydata
## Counted within unit, where unit  falpha0is defined by the combinationof year, serialnumber and agereader
## d is an S times A matrix giving counts to different agegroups for each unit
## N is a vector of length S giving the total numbers for each unit
## strataindex is a vector of length S giving the strata
## yearindex is a vector of length S giving the year
## readerindex is a vector of length S giving the reader  (note: NA is given reader number 9 now!!!)

library(AgeHerring)
library(tidyr)

dir = "/nr/project/stat/HI-alderstolkning/"
HerringSet <- data.table::fread(paste(dir,"Data/agereadingset_herring_spawning_survey_2023.csv",sep=""))
#x <- data.table::fread(dir,paste("Data/SmartDots_Event_448_2022 NSS herring exchange (scales).csv",sep=""))
#x <- data.table::fread(paste("data/SmartDots_Event_448_2022_NSS_herring_exchange_scales.csv",sep=""))

#Removing NA´s for age
HerringSet = HerringSet[!is.na(age)]

# Ad hoc solution - setting NA to reader 9
#HerringSet$AGEREADER[is.na(HerringSet$AGEREADER)] = "reader 9"

#Remove observations with missing AGEREADER
HerringSet = HerringSet[!is.na(AGEREADER)]

year = sort(unique(HerringSet$startyear))
Y = length(year)
Amin = min(HerringSet$age)
Amax = max(HerringSet$age)
Amin = 3
Amax = 12
HerringSet$age = pmin(pmax(HerringSet$age,Amin),Amax)-Amin+1
Amin = 1
Amax = 10
A = Amax-Amin+1

#Make indices for strata
strat.u = unique(HerringSet$SPATIAL_STRATUM)
K = length(strat.u)
year = c(2020,2022,2023)
stratayearindex = rep(NA,K)
i = 0;yy=0
for(y in year)
{
  yy = yy+1
  strat.y = sort(strat.u[substring(strat.u,first=1,last=4)==y])
  for(j in 1:length(strat.y))
  {
    i = i + 1
    stratayearindex[i] = yy
  }
}

#Make indices for stations
stations = paste0(HerringSet$serialnumber,".",HerringSet$startyear)
stations.unique = unique(stations)
S = length(stations.unique)
stationstrataindex = rep(NA,S)
for(s in 1:S)
{
  data.s = HerringSet[stations==stations.unique[s],]
  stationstrataindex[s] = pmatch(data.s$SPATIAL_STRATUM[1],strat.u)
}


#Make table of readings within unit
#Unit based on combination of year, serialnumber and agereader
serialn = paste0(HerringSet$startyear,".",HerringSet$serialnumber,".",HerringSet$AGEREADER)
Unique.serialn = unique(serialn)
U = length(Unique.serialn)
reader = sort(unique(HerringSet$AGEREADER))
R = length(reader)
d = array(0,c(U,A))
unitstationindex = rep(NA,U)
unitreaderindex = rep(NA,U)
for(u in 1:U)
  {
    dataS = HerringSet[serialn==Unique.serialn[u]]
    #Check if unit has several spatial locations
    if(length(unique(dataS$longitudestart))>1)
      cat("Non-unique longitude within unit",u,"\n")
    #Check if units has different years
    if(length(unique(dataS$startyear))>1)
      cat("Non-unique year within unit",u,"\n")
    #Check if unit has different readers
    if(length(unique(dataS$AGEREADER))>1)
      cat("Non-unique reader within unit",u,length(unique(dataS$AGEREADER)),"\n")
    #Check if unit has different strata
    if(length(unique(dataS$SPATIAL_STRATUM))>1)
      cat("Non-unique stratum within unit",u,length(unique(dataS$SPATIAL_STRATUM)),"\n")
    station = paste0(dataS$serialnumber[1],".",dataS$startyear[1])
    unitstationindex[u] = pmatch(station,stations.unique)
    unitreaderindex[u] = pmatch(dataS$AGEREADER[1],reader)
    for(r in 1:R)
    {
      dataR = dataS[AGEREADER==reader[r]]
      if(nrow(dataR)>0)
      {
        d[u,] = table(c(dataR$age,Amin:Amax))-1
      }
    }
  }
N = rowSums(d)
HerringSurvey =
 list(d=d,N=N,K=K,S=S,U=U,
      stratayearindex=stratayearindex,stationstrataindex=stationstrataindex,
      unitstationindex=unitstationindex,
      unitreaderindex=unitreaderindex)
usethis::use_data(HerringSurvey,overwrite=TRUE)


if(0)
{
## Making strata info
HerringSet$longitudestart = round(HerringSet$longitudestart,2)
loc = unique(HerringSet[,c("longitudestart")])
station = rep(NA,nrow(HerringSet))
n = length(station)
for(i in 1:nrow(loc))
#for(i in 1:3)
  {
  #print(i)
  ind = which(HerringSet$longitudestart %in% loc[i])
  #print(ind)
  station[ind] = i
}
HerringSet$station = station
Herring0 = HerringSet[HerringSet$startyear==2021]
strat = substring(HerringSet$SPATIAL_STRATUM,first=6)
strat0 = strat[HerringSet$startyear==2021]
loc = matrix(nrow=length(unique(strat0)),ncol=2)
data0 = HerringSet[HeringSet$startyear==2021,]

Y = length(HeringSet$startyear)
K = length(unique(strat))
INDSUM = table(HeringSet$startyear,strat,HeringSet$serialnumber)
lab = labels(INDSUM)
w = array(0,dim(INDSUM))
for(i in 1:dim(INDSUM)[1])
  for(j in 1:dim(INDSUM)[2])
    for(k in 1:dim(INDSUM)[3])
    {
      f = HeringSet[(HeringSet$startyear==lab[[1]][i]) &
                       (strat==lab[[2]][j]) &
                       (HeringSet$serialnumber==lab[[3]][k]),]
      if(nrow(f)>0)
      {
       show(c(i,j,k,length(unique(f$SPATIAL_STRATUM))))
       w[i,j,k] = f$REL_ABUNDANCE[1]*f$SPATIAL_STRATUM_AREA[1]/nrow(f)
      }
    }
wstrat = apply(w,c(1,3),sum)
Herring.Stratum = list(INDSUM=INDSUM,w=apply(w,c(1,3),sum))
usethis::use_data(Herring.Stratum,overwrite=TRUE)
}
