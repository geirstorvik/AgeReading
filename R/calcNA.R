#' Calculation of transition probabilities for age error readings.
#'
#' Default values on parameters are from fitted model to ...
#' @param w An y times K matrix giving
#' @param Amin mininum age
#' @param Amax maximum age
#' @return Transition matrix
#' @export
calNA = function(w,xiks,stratum,Y,K,S,A){
  for(y in 1:Y){
    for(a in 1:A){
      for(k in 1:K){
        Nhat[y,k,a] = 0
      }
      for(s in 1:S)
      {
        k = stratum[s];
        Nhat[y,k,a] = Nhat[y,k,a] + w[y,s]*xiks[y,k,s,a]

      }
    }
  }
  Nhat
}
