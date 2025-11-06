#' Calculate within each stratum estimate of number of fish within each age group
#'
#' @param fitobj object fitted by stan
#' @param weights Weight within each stratum (D*A/n)
#' @param INDSUM indsum
#' @return Matrix containing probabilities for error readings
#' @export
calc_stratum_estimate = function(fitobj,weights,INDSUM,A)
{
  Y = dim(INDSUM)[1]
  K = dim(INDSUM)[2]
  S = dim(INDSUM)[3]
  ind = expand.grid(1:Y,1:S,1:A)
  pars=paste("xis[",ind[,1],",",ind[,2],",",ind[,3],"]",sep="")
  xis = matrix(unlist(extract(fitobj,pars=pars)),ncol=Y*S*A)
  M = dim(xis)[1]
  Mxis = array(xis,c(M,Y,S,A))

  ind = expand.grid(1:Y,1:A)
  pars=paste("xi0[",ind[,1],",",ind[,2],",",ind[,3],"]",sep="")
  xi0 = matrix(unlist(extract(fitobj,pars=pars)),ncol=Y*S*A)
  M0 = dim(xi0)[1]
  Mxi0 = array(xi0,c(M,Y,A))


  Nhat = array(NA,c(M,Y,A))
  Nhat0 = array(NA,c(M,Y,A))
  for(y in 1:Y){
    for(a in 1:A){
      for(m in 1:M){
       Nhat[m,y,a]=0
       Nhat0[m,y,a]=0
       for(k in 1:K){
         Nhat0[m,y,a] = Nhat0[m,y,a] + weights[y,k]*Mxi0[m,y,a]
         for(s in 1:S){
           Nhat[m,y,a] = Nhat[m,y,a] + weights[y,s]*Mxis[m,y,s,a]
         }}}}}
  Nhat
}
