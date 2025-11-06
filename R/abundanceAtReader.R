#' Estimate of abundance for all ages for each reader
#' Essentially extracting data corresponding to each reader, but agesamplecount needs to
#' be corrected with respect to how many that are read for the specific reader.
#' @param data data.table formatted as explained above
#' @param year survey year to calculate abundance for
#' @return estimated abundance for age for each reader.
#' @import data.table
#' @export
abundanceAtReader <- function(data, year,Plot=TRUE){
  agerange <- 0:max(data$age, na.rm = T)
  reader = unique(data$AGEREADER)
  reader = sort(reader[!is.na(reader)])
  x = NULL
  y = NULL
  for(r in reader)
  {
    dataH = data[AGEREADER==r]
    tab = table(dataH$serialnumber)
    ##Change agesamplecount to the numbers related to the current reader
    for(i in 1:length(tab))
    {
      dataH$agesamplecount[dataH$serialnumber==names(tab)[i]]=tab[i]
    }
    ageAbundance <- sapply(agerange, function(x){abundanceAtAge(dataH, year, x)})
    names(ageAbundance) <- agerange
    x = cbind(x,agerange)
    y = cbind(y,ageAbundance)

  }
  if(Plot)
  {
    matplot(x,y,type="l",col=1:5,lty=1,lwd=2,xlab="age",ylab="Mean abundance")
    legend("topright",legend=reader,col=1:5,lty=1,lwd=2)
  }
  list(agerange=x,ageAbundance=y)
}
