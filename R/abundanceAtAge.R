#' Estimate of abundance at age
#' @param data data.table formatted as explained above
#' @param year survey year to calculate abundance for
#' @param age to calculate abundance for
#' @return estimated abundance for age.
#' @import data.table
#' @export
abundanceAtAge <- function(data, year, age){

  require(data.table)
  #filter by year
  data <- data[data$startyear == year,]

  #get totals for each station / haul
  data$atAge <- !is.na(data$age) & data$age==age
  data$atAgeCs <- data$atAge * data$catchcount / data$agesamplecount
  stationTotal <- data[,list(totalFishH=sum(catchcount), atAgeH=sum(atAgeCs), abundance=REL_ABUNDANCE), by=c("startyear", "serialnumber", "SPATIAL_STRATUM", "SPATIAL_STRATUM_AREA")]

  #get proportion at station
  stationTotal$propAtAgeH <- stationTotal$atAgeH/stationTotal$totalFishH

  #get mean abundance by spatial strata
  strataMean <- stationTotal[,list(strataMean=mean(propAtAgeH*abundance)), by=c("startyear", "SPATIAL_STRATUM", "SPATIAL_STRATUM_AREA")]

  #get mean abundance for total area
  return(sum(strataMean$strataMean/strataMean$SPATIAL_STRATUM_AREA)*sum(strataMean$SPATIAL_STRATUM_AREA))
}
