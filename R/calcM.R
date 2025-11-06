#' Calculation of transition probabilities for age error readings.
#'
#' Default values on parameters are from fitted model to ...
#' @param theta (alpha0,alpha1,beta0,beta1,log(phi))
#' @param Amin mininum age
#' @param Amax maximum age
#' @return Transition matrix
#' @export
calcM = function(theta,Amin,Amax){
  inv_logit = function(x){exp(x)/(1+exp(x))}
  A = Amax-Amin+1
  Amed = mean(c(Amin,Amax))
  #Amed=0
  M = matrix(NA,A,A)
  q = rep(NA,A)
  phi = inv_logit(theta[5]);
  pA = inv_logit(theta[3]+theta[4]*(Amin-Amed));
  q[1] = inv_logit(theta[1]+theta[2]*(Amin-Amed));
  M[1,1] = 1-q[1]+pA*q[1];
  l = 2:(A-1)
  M[1,l] = (1-pA)*q[1]*(1-phi)*phi^(l-2);
  M[1,A] = (1-pA)*q[1]*phi^(A-2);
  for(k in 2:(A-1))
  {
    a = Amin+k-1;
    q[k] = inv_logit(theta[1]+theta[2]*(a-Amed));
    M[k,k] = inv_logit(theta[3]+theta[4]*(a-Amed));
    konst = (1-M[k,k])*(1-q[k])*(1-phi)/(1-phi^(k-1));
    l = 1:(k-1)
    M[k,l] = konst*phi^(k-l-1);
    konst = (1-M[k,k])*q[k]*(1-phi);
    l =  (k+1):A
    M[k,l] = konst*phi^(l-k-1);
    M[k,A] = (1-M[k,k])*q[k]*phi^(A-k-1);
  }
  q[A] = inv_logit(theta[1]+theta[2]*(Amax-Amed));
  pA = inv_logit(theta[3]+theta[4]*(Amax-Amed));
  M[A,A] = 1-(1-pA)*(1-q[A]);
  konst = (1-pA)*(1-q[A])*(1-phi)/(1-phi^(Amax-Amin));
  l = 1:(A-1)
  M[A,l] = konst*phi^(A-l-1);
  return(M);
}
