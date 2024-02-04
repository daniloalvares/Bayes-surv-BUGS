data{
  int n;
  int Ntrans;
  int Nbetas;
  matrix[n, Ntrans] time;
  matrix[n, Ntrans] delta;
  matrix[n, Ntrans] X;
}

parameters{
  matrix[Nbetas, Ntrans] beta;
  real<lower=0> lambda[Ntrans];
  real<lower=0> alpha[Ntrans];
}

model{
  matrix[n, Ntrans] logHaz;
  matrix[n, Ntrans] cumHaz;

  // MULTI-STATE SPECIFICATION
  for(i in 1:n){
     for(k in 1:Ntrans){
        // Log-hazard function
        logHaz[i,k] = log(lambda[k]) + log(alpha[k]) + (alpha[k] - 1) * log(time[i,k]) + X[i,] * beta[,k];
        // Cumulative hazard function
        cumHaz[i,k] = lambda[k] * pow(time[i,k], alpha[k]) * exp( X[i,] * beta[,k] );
     }
     target += dot_product(delta[i,], logHaz[i,]) - sum(cumHaz[i,]);
  }

  // LOG-PRIORS             
  // Coefficients
  target += normal_lpdf(to_vector(beta) | 0, sqrt(1000));

  // Weibull scale and shape parameters
  target += gamma_lpdf(lambda | 0.01, 0.01);
  target += uniform_lpdf(alpha | 0, 10);
}
