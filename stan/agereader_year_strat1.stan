#include calcM.stan

data {
  int<lower=1> Amin;     //Minimum age
  int<lower=1> A;        //Number of age groups
  int<lower=1> Astar;    //A* for standardization
  int<lower=1> R;       //Number of readers for dataset 1
  int<lower=0> readErr[R,A,A];  //Estimate ages for dataset 1
  vector<lower=0>[A] alpharep;
  real taupar[5];
}

parameters { //Overall relative abundance
real Palpha0;
real Palpha1;
real Pbeta0;
real Pbeta1;
real Ptrphi;
real<lower=0,upper=1> tau[5];
}

transformed parameters {
  //Using standardized variables in order to get reasonable initial values
  real Pphi;
  matrix[A,A] M1;

  Pphi = inv_logit(Ptrphi);
  M1 = calcM(Amin,A,Astar,Palpha0,Palpha1,Pbeta0,Pbeta1,Ptrphi);
}

model {
  //Global parameters
  target += normal_lpdf(Palpha0|0,3);
  target += normal_lpdf(Palpha1|0,3);
  target += normal_lpdf(Pbeta0|0,3);
  target += normal_lpdf(Pbeta1|0,3);
  target += normal_lpdf(Ptrphi|0,3);

  //Dataset 1
  for(j in 1:R){
    vector[A] sum1;
    for(a in 1:A) {
      sum1[a] = multinomial_lpmf(readErr[j,a,1:A]|M1[1:A,a]);
     }
    target += sum(sum1);
  }
}
