//
// (alpha0,alpha1)=theta[1:2]
// (beta0,beta1)=theta[3:4]
// phi=exp(trphi)
// Changed M to its transpose
functions {
  matrix calcM(int Amin,int A,int Astar,real alpha0,real alpha1,real beta0,real beta1,real trphi){
    int Amax = Amin+A-1;
    int a;
    real konst;
    real pA;
    matrix[A,A] M;
    vector[A] Msum;
    vector[A] q;
    real phi=inv_logit(trphi);

    phi = inv_logit(trphi);
    pA = inv_logit(beta0+beta1*(Amin-Astar));
    q[1] = inv_logit(alpha0+alpha1*(Amin-Astar));
    M[1,1] = 1-q[1]+pA*q[1];
    for(l in 2:(A-1))
     M[l,1] = (1-pA)*q[1]*(1-phi)*phi^(l-2);
    M[A,1] = (1-pA)*q[1]*phi^(A-2);
    for(k in 2:(A-1))
    {
      a = Amin+k-1-Astar;
      q[k] = inv_logit(alpha0+alpha1*a);
      M[k,k] = inv_logit(beta0+beta1*a);
      konst = (1-M[k,k])*(1-q[k])*(1-phi)/(1-phi^(k-1));
      //M[k,1:(k-1)] = konst*phi^((k-2):0);
      for(l in 1:(k-1))
        M[l,k] = konst*phi^(k-l-1);
      konst = (1-M[k,k])*q[k]*(1-phi);
      for(l in (k+1):A)
        M[l,k] = konst*phi^(l-k-1);
      M[A,k] = (1-M[k,k])*q[k]*phi^(A-k-1);
    }
    q[A] = inv_logit(alpha0+alpha1*(Amax-Astar));
    pA = inv_logit(beta0+beta1*(Amax-Astar));
    M[A,A] = 1-(1-pA)*(1-q[A]);
    konst = (1-pA)*(1-q[A])*(1-phi)/(1-phi^(A-1));
    for(l in 1:(A-1))
       M[l,A] = konst*phi^(A-l-1);
     return(M);
  }
  vector calcr(int K,int S,int[,] INDSUM,real[] xis){
    vector[K] r;
    for(k in 1:K)
    {
      r[k] = 0.0;
        for(s in 1:S){
          r[k] = r[k] + INDSUM[k,s]*xis[s];
        }
    }
    return(r);
  }
  real[,,] calcNA(real[,] w,real[,,,] xiks,int[] stratum,int Y,int K,int S,int A){
  real Nhat[Y,K,A];
  int k;
  for(y in 1:Y){
    for(a in 1:A){
      for(l in 1:K){
        Nhat[y,l,a] = 0;
      }
      for(s in 1:S)
      {
        k = stratum[s];
        Nhat[y,k,a] = Nhat[y,k,a] + w[y,s]*xiks[y,k,s,a];
      }
    }
  }
  return(Nhat);
 }
}
