#include calcM.stan

data {
  int<lower=1> Amin;     //Minimum age
  int<lower=1> A;        //Number of age groups
  int<lower=1> Astar;    //A* for standardization
  int<lower=1> R;       //Number of readers for dataset 1
  int<lower=0> readErr[R,A,A];  //Estimate ages for dataset 1
  vector<lower=0>[A] alpharep;
  real deltapar[2];
  real taupar[2];
}

parameters { //Overall relative abundance
real Palpha0;
real Palpha1;
real Pbeta0;
real Pbeta1;
real Ptrphi;
real<lower=0> tau[5];
real<lower=0> delta[2];           //Variance parameter for variability between stations
real alpha0[R];
real beta0[R];
real alpha1[R];
real beta1[R];
real trphi[R];
}

transformed parameters {
  //Using standardized variables in order to get reasonable initial values
  matrix[A,A] M1[R];


  for(j in 1:R){
    M1[j] = calcM(Amin,A,Astar,alpha0[j],alpha1[j],beta0[j],beta1[j],trphi[j]);
  }
}

model {
  //Global parameters
  target += normal_lpdf(Palpha0|0,3);
  target += normal_lpdf(Palpha1|0,3);
  target += normal_lpdf(Pbeta0|0,3);
  target += normal_lpdf(Pbeta1|0,3);
  target += normal_lpdf(Ptrphi|0,3);

  target += uniform_lpdf(delta|0.0,20.0);
  target += normal_lpdf(alpha0|Palpha0,sqrt(tau[1]));
  target += normal_lpdf(beta0|Pbeta0,sqrt(tau[2]));
  target += normal_lpdf(alpha1|Palpha1,sqrt(tau[3]));
  target += normal_lpdf(beta1|Pbeta1,sqrt(tau[4]));
  target += normal_lpdf(trphi|Ptrphi,sqrt(tau[5]));
  target += uniform_lpdf(tau|taupar[1],taupar[2]);
  //Dataset 1
  for(j in 1:R){
    vector[A] sum1;
    for(a in 1:A) {
      sum1[a] = multinomial_lpmf(readErr[j,a,1:A]|M1[j,1:A,a]);
    }
    target += sum(sum1);
  }
}
