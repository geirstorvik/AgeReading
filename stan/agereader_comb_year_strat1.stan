#include calcM.stan

data {
  int<lower=1> Amin;     //Minimum age
  int<lower=1> A;        //Number of age groups
  int<lower=1> Astar;    //A* for standardization
  int<lower=1> R1;       //Number of readers for dataset 1
  int<lower=0> readErr[R1,A,A];  //Estimate ages for dataset 1
  int<lower=1> Y;        //Number of years for dataset 2
  int<lower=1> S;        //Number of stations for dataset 2, within each year
  int<lower=1> U;        //Number of units for dataset 2, within each year
  int<lower=1> R2;       //Number of readers for dataset 2
  int<lower=0> D[U,A];  //Number of fish within agegroup a within station s from reader r
  int<lower=0> K;       // Number of strata
  int<lower=1,upper=Y>  stratayearindex[K];  //Index for strata for a specific station
  int<lower=1,upper=K>  stationstrataindex[S];    //Index for year for a specific station
  int<lower=1,upper=S>  unitstationindex[U];
  int<lower=1,upper=R2> unitreaderindex[U];
  vector<lower=0>[A] alpharep;
}


transformed data{
  int<lower=1> R=R1+R2;
}

parameters { //Overall relative abundance
real Palpha0;
real Palpha1;
real Pbeta0;
real<upper=0> Pbeta1;
real Ptrphi;
real<lower=0> delta[2];           //Variance parameter for variability between stations
simplex[A] xi0[Y];
simplex[A] xi1[K];
simplex[A] xi2[S];
}

transformed parameters {
  //Using standardized variables in order to get reasonable initial values
  real Pphi;
  real loglik;
  matrix[A,A] M1;
  Pphi = inv_logit(Ptrphi);

  M1 = calcM(Amin,A,Astar,Palpha0,Palpha1,Pbeta0,Pbeta1,Ptrphi);

  loglik = 0.0;

  for(j in 1:R1){
    for(a in 1:A) {
      loglik += multinomial_lpmf(readErr[j,a,1:A]|M1[1:A,a]);
     }
  }
  //Dataset 2
  for(u in 1:U) {
    vector[A] pi_ = M1*xi2[unitstationindex[u]];
    loglik += multinomial_lpmf(D[u]|pi_);
  }
}

model {
  //Global parameters
  target += normal_lpdf(Palpha0|0,3);
  target += normal_lpdf(Palpha1|0,3);
  target += normal_lpdf(Pbeta0|0,3);
  target += normal_lpdf(Pbeta1|0,3);
  target += normal_lpdf(Ptrphi|0,3);

  target += uniform_lpdf(delta |0.0,20.0);
  for(y in 1:Y){
    target += dirichlet_lpdf(xi0[y]|alpharep);
  }
  for(k in 1:K){
      target += dirichlet_lpdf(xi1[k]|delta[1]*xi0[stratayearindex[k]]);
    }
  for(s in 1:S){
    target += dirichlet_lpdf(xi2[s]|delta[2]*xi1[stationstrataindex[s]]);
  }
  //Dataset 1
  for(j in 1:R1){
    vector[A] sum1;
    for(a in 1:A) {
      sum1[a] = multinomial_lpmf(readErr[j,a,1:A]|M1[1:A,a]);
     }
    target += sum(sum1);
  }
  //Dataset 2
  for(u in 1:U) {
    vector[A] pi_ = M1*xi2[unitstationindex[u]];
    target += multinomial_lpmf(D[u]|pi_);
  }
}
