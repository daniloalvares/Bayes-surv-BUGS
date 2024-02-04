data{
  int n;
  int Nrisks;
  int Nbetas;
  vector[n] time;
  matrix[n, Nrisks] delta;
  matrix[n, Nbetas] X;
}

parameters{
  matrix[Nbetas, Nrisks] beta;
  real<lower=0> lambda[Nrisks];
  real<lower=0> alpha[Nrisks];
}

model{
  matrix[n, Nrisks] logHaz;
  matrix[n, Nrisks] cumHaz;

  // COMPETING RISKS SPECIFICATION
  for(i in 1:n){
     for(k in 1:Nrisks){
        // Log-hazard function
        logHaz[i,k] = log(lambda[k]) + log(alpha[k]) + (alpha[k] - 1) * log(time[i]) + X[i,] * beta[,k];
        // Cumulative hazard function
        cumHaz[i,k] = lambda[k] * pow(time[i], alpha[k]) * exp( X[i,] * beta[,k] );
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
