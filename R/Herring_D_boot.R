#' Herring D boot dataset
#'
#' The data set  is obtained from data recorded at the annual spawning survey
#' in February / March (cruise series ‘Norwegian Sea NOR Norwegian spring-spawning
#' herring spawning cruise in Feb_Mar’).
#'
#' @format ## `HerringSurvey`
#' A list containing:
#' \describe{
#' \item{d}{A matrix of dimension 128x13 containing the number of fish read in the A=13 different agegroups within $S=128$ units
#' where units are defined as combinations of year, serialnumber and reader.}
#' \item{N}{A vector of length 128 giving the total number of fish read within each unit}
#' \item{K}{The number of strata (total over all years)}
#' \item{S}{The number of stations (total over all years)}
#' \item{U}{The number of observation units (combination of station/reader,total over all years)}
#' \item{stratayearindex}{a vector of length K giving the year for each strata}
#' \item{stationstrataindex}{a vector of length  a vector of length S giving the strata for each obs unit}
#' \item{unitstationindex}{a vector of length U giving the station for each obs unit}
#' \item{unitreaderindex}{a vector of length U giving the reader for each obs unit}
#' }
#' @source HI
"HerringSurvey"
