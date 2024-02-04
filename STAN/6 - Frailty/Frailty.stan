data{
  int n;
  int J;
  int Nbetas;
  matrix[n, J] time;
  matrix[n, J] delta;
  matrix[n, Nbetas] X;
}

parameters{  
  vector[Nbetas] beta;
  real<lower=0> alpha;
  real<lower=0> psi;
  real<lower=0> w[n];
}

transformed parameters{
  real lambda = exp( beta[1] );
}

model{
  matrix[n, J] logHaz;
  matrix[n, J] cumHaz;

  // FRAILTY SPECIFICATION
  for(i in 1:n){
     for(j in 1:J){
        // Log-hazard function
        logHaz[i,j] = log(alpha) + (alpha-1) * log(time[i,j]) + X[i,] * beta + log(w[i]) ;
        // Cumulative hazard function
        cumHaz[i,j] = (time[i,j]^alpha) * exp( X[i,] * beta + log(w[i]) );

        target += delta[i,j] * logHaz[i,j] - cumHaz[i,j];
     }
  }

  // LOG-PRIORS             
  // Coefficients
  target += normal_lpdf(beta | 0, sqrt(1000));  

  // Weibull shape parameter
  target += uniform_lpdf(alpha | 0, 10);

  // Frailty parameter
  target += gamma_lpdf(psi | 0.01, 0.01);

  // Frailty terms
  target += gamma_lpdf(w | psi, psi);
}
