#' Herring data set with multiple readers
#' Three non-expert readers are removed from the dataset
#' Ages are first truncated to the interval {3,15} and thereafter shifted to {1,13}
#'
#' @docType data
#'
#' @usage data(Herring_multiple_readers)
#'
#' @format A data frame with 1473 rows and 4 columns:
#' \describe{
#' \item{trueage}{Age of fish, estimated by taking the mode of all readers}
#' \item{readage}{Read age}
#' \item{ID}{Fish ID}
#' \item{reader}{reader, converted to numeric values 1-6 where 1=R02 NO (Adam), 2=R04 NO (Justine), 3=R08 IS R14 NO (Stine), 4=R18 IS R24 NO (Bjørn Vidar), 5=R16 NO (Timo) R30 FO and 6=R32 FO}
#' }
#' @source HI
"Herring_multiple_readers"
