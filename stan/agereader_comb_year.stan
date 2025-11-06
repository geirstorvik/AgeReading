#include calcM.stan

data {
  int<lower=1> Amin;    //Minimum age
  int<lower=1> A;       //Number of age groups
  int<lower=1> Astar;   //A* for standardization
  int<lower=1> R1;       //Number of readers for dataset 1
  int<lower=0> readErr[R1,A,A];  //Estimate ages for dataset 1
  int<lower=1> Y;       //Number of years for dataset 2
  int<lower=1> S;       //Number of stations for dataset 2
  int<lower=1> R2;       //Number of readers for dataset 2
  int<lower=0> D[Y,S,R2,A];//Number of fish within agegroup a within station s from reader r
  vector<lower=0>[A] alpharep;
  real deltapar[2];
  real taupar[2];
}


transformed data{
  int<lower=1> R=R1+R2;
}

parameters { //Overall relative abundance
real Palpha0;
real Pbeta0;
real Pbeta1;
real Ptrphi;
real<lower=0> tau[2];
real<lower=0> delta_std;           //Variance parameter for variability between stations
simplex[A] xi0[Y];
real alpha0_std[R];
real beta0_std[R];
simplex[A] xis[Y,S];
//real r[K,A];                      //Sum of proportions over all stations within stratum
}

transformed parameters {
  //Using standardized variables in order to get reasonable initial values
  real alpha0[R];
  real beta0[R];
  real delta;
  real Pphi;
  matrix[A,A] M1[R1+R2];
  Pphi = inv_logit(Ptrphi);
  for(j in 1:R){
    alpha0[j] = Palpha0+tau[1]*alpha0_std[j];
    beta0[j] = Pbeta0+tau[2]*beta0_std[j];
  }
  delta = ((deltapar[2]-deltapar[1])*delta_std+2.0*(deltapar[1]+deltapar[2]))/4.0;
  for(j in 1:R1){
    M1[j] = calcM(Amin,A,Astar,alpha0[j],0.0,beta0[j],Pbeta1,Ptrphi);
  }
  for(j in 1:R2)
    M1[R1+j] = calcM(Amin,A,Astar,alpha0[R1+j],0.0,beta0[R1+j],Pbeta1,Ptrphi);
}

model {
  //Global parameters
  target += normal_lpdf(Palpha0|0,10);
  target += normal_lpdf(Pbeta0|0,10);
  target += normal_lpdf(Pbeta1|0,10);
  target += normal_lpdf(Ptrphi|0,10);
  target += uniform_lpdf(delta_std|-2.0,2.0);
  for(j in 1:R){
    target += std_normal_lpdf(alpha0_std[j]);
    target += std_normal_lpdf(beta0_std[j]);
  }
  target += inv_gamma_lpdf(tau|taupar[1],taupar[2]);
  for(y in 1:Y){
    target += dirichlet_lpdf(xi0[y]|alpharep);
  }
  for(y in 1:Y) {
    for(s in 1:S){
      target += dirichlet_lpdf(xis[y,s]|delta*xi0[y]);
    }
  }
  //Dataset 1
  for(j in 1:R1){
    vector[A] sum1;
    for(a in 1:A) {
      sum1[a] = multinomial_lpmf(readErr[j,a,1:A]|M1[j,1:A,a]);
    }
  target += sum(sum1);
  }
  //Dataset 2
  for(j in 1:R2){
    for(y in 1:Y) {
      vector[S] sum2;
      for(s in 1:S) {
        vector[A] pi_ = M1[R1+j]*xis[y,s];
        sum2[s] = multinomial_lpmf(D[y,s,j]|pi_);
        //target += multinomial_lpmf(D[y,s,j]|pi_);
      }
      target += sum(sum2);
    }
  }
}
