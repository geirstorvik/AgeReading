data {
  int<lower=1> Amin;     //Minimum age
  int<lower=1> A;        //Number of age groups
  int<lower=1> Astar;    //A* for standardization
  int<lower=1> Y;        //Number of years for dataset 2
  int<lower=1> S;        //Number of stations for dataset 2, within each year
  int<lower=1> U;        //Number of units for dataset 2, within each year
  int<lower=0> D[U,A];  //Number of fish within agegroup a within station s from reader r
  int<lower=0> K;       // Number of strata
  int<lower=1,upper=Y>  stratayearindex[K];  //Index for strata for a specific station
  int<lower=1,upper=K>  stationstrataindex[S];    //Index for year for a specific station
  int<lower=1,upper=S>  unitstationindex[U];
  vector<lower=0>[A] alpharep;
  real deltapar[2];
}

parameters { //Overall relative abundance
real<lower=0> delta[2];           //Variance parameter for variability between stations
simplex[A] xi0[Y];
simplex[A] xi1[K];
simplex[A] xi2[S];
}


model {
  //Global parameters
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
  //Dataset 2
  for(u in 1:U) {
    vector[A] pi_ = xi2[unitstationindex[u]];
    target += multinomial_lpmf(D[u]|pi_);
  }
}
