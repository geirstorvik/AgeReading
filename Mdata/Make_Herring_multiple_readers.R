## Script for making data set 1 (several agereaders on same fish)
library(readxl)
library(AgeHerring)
library(data.table)
#Read data
dir = "/nr/project/stat/HI-alderstolkning/"
data_s = read_xlsx(paste(dir,"Data/NSS Exchange All readers.xlsx",sep=""), sheet = 1)

#Extract relevant info
data_slong = gather(data_s, reader, age, `R02 NO (Adam)`:`R32 FO`, factor_key=TRUE)
data_slong$age = as.numeric(data_slong$age)
#Only consider first 6 readers
data_slong$exp = 1
data_slong$exp[data_slong$reader %in% c("R16 NO (Timo)", "R30 FO", "R32 FO")] = 0
data_slong = data_slong[data_slong$exp==1,]

data_slong$country = sapply(as.character(data_slong$reader), function(x) strsplit(x, " ")[[1]][2], USE.NAMES=FALSE)
data_slong$reader = as.factor(data_slong$reader)
data_exp = data_slong[data_slong$reader%in%names(data_s)[7:12] & !is.na(data_slong$age), ]
true_age = pmin(data_exp$`Modal age (exp)`-2,13)
age = pmin(data_exp$age-2,13)
read = as.numeric(data_exp$reader)
Herring_multiple_readers = data.frame(ID=data_exp$`Fish ID`, trueage=true_age,readage=age,reader=data_exp$reader)
usethis::use_data(Herring_multiple_readers,overwrite=TRUE)

