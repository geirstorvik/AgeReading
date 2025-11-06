# Making ring data for stan processing
# Truncates readings to the interval {3,15} which is shifted to {1,13}
# Only readers 1-6 (the experts) are considered
library(readxl)
library(tidyr)
library(stringr)

#Reading data, skipping last row
dir = "/nr/project/stat/HI-alderstolkning/"
data_s = read_xlsx(paste(dir,"Data/NSS Exchange All readers.xlsx",sep=""), sheet = 1,n_max=255,na="NA")
dir = "/nr/samba/user/storvik/prj/agereading_herring/data/"
modaleAge = read.table(paste0(dir,"modale_aldre.txt"),header=TRUE)
n = nrow(data_s)
data_s$Modal2 = NA
data_s$Modal3 = NA
set.seed(123)
for(i in 1:n)
{
  ind = grep(data_s$`Fish ID`[i],modaleAge$FishID)
  if(is.null(ind))
    print(i)
  data_s$Modal2[i] = sample(modaleAge$modal_age_closest[c(ind,ind)],1)
  data_s$Modal3[i] = mean(as.numeric(data_s[i,7:12]),na.rm=TRUE)
}

data_slong = gather(data_s, reader, age, `R02 NO (Adam)`:`R32 FO`, factor_key=TRUE)
data_slong$age = as.numeric(data_slong$age)
data_slong$exp = 1
data_slong$exp[data_slong$reader %in% c("R16 NO (Timo)", "R30 FO", "R32 FO")] = 0
data_slong$country = sapply(as.character(data_slong$reader), function(x) strsplit(x, " ")[[1]][2], USE.NAMES=FALSE)
data_slong$reader = as.factor(data_slong$reader)
data_slong = data_slong[!is.na(data_slong$age),]
data_slong = data_slong[data_slong$exp==1,]
ID = unique(data_slong$`Fish ID`)
data_slong$fishID = match(data_slong$`Fish ID`,ID)
Reader = unique(data_slong$reader)
data_slong$ReaderN = match(data_slong$reader,Reader)
Ring.full = data_slong[,c("fishID","age","ReaderN")]
Ring.full = Ring.full[!is.na(Ring.full$age),]
Ring.full$age = pmin(pmax(Ring.full$age,Amin),Amax)-Amin+1
usethis::use_data(Ring.full,overwrite=TRUE)

Amin = 3
Amax = 12
age = pmin(pmax(data_slong$age,Amin),Amax)-Amin+1
agetrue = pmin(pmax(data_slong$`Modal age (exp)`,Amin),Amax)-Amin+1
agetrue = pmin(pmax(data_slong$Modal2,Amin),Amax)-Amin+1
Amin = 1
Amax = 10
A = Amax-Amin+1

readErr = table(data_slong$reader,agetrue,age)
readErr = readErr[1:6,,]
Ring.readErr = readErr
usethis::use_data(Ring.readErr,overwrite=TRUE)
