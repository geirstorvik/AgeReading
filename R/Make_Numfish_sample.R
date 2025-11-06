#' Calculate samples of number of fish based on posterior samples of xi1
#' and boostrap samples of abundance
#'
#' @param xi1 Array of dimension MxKxA giving posterior samples of xi1, the age distribution per strata
#' @param bootD bootstrap samples of abundance within each strata
#' @return estimated abundance for age for each reader.
#' @import data.table
#' @export
Make_Numfish_sample = function(xi1,bootD,stratayearindex)
{
  M = dim(xi1)[1]
  K = dim(xi1)[2]
  A = dim(xi1)[3]
  Y = length(unique(stratayearindex))
  B = max(bootD[,1])
  Numfish = array(NA,c(M*B,Y,A))
 for(b in 1:B)
  {
    for(y in 1:Y)
    {
    Db = bootD[(bootD[,1]==b) & (stratayearindex[bootD[,"unitstrat"]]==y),"Abundance"]
    for(m in 1:M)
      Numfish[(b-1)*M+m,y,] = colSums(diag(Db)%*%xi1[m,stratayearindex==y,])
    }
  }
 Numfish
}
