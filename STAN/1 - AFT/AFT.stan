data{
  int n; 
  int Nbetas;
  vector[n] time;
  vector[n] delta;
  matrix[n, Nbetas] X;
}

parameters{  
  vector[Nbetas] beta;
  real<lower=0> alpha;
}

transformed parameters{
  real<lower=0> sigma = 1/alpha;
}

model{
  vector[n] logHaz;
  vector[n] cumHaz;

  // GUMBEL AFT SPECIFICATION
  for(i in 1:n){
     // Log-hazard function
     logHaz[i] = log(alpha) + (alpha-1) * log(time[i]) - X[i,] * beta * alpha;
     // Cumulative hazard function
     cumHaz[i] = (time[i]^alpha) * exp( - X[i,] * beta * alpha );

     target += delta[i] * logHaz[i] - cumHaz[i];
  }

  // LOG-PRIORS             
  // Coefficients
  target += normal_lpdf(beta | 0, sqrt(1000));  

  // Weibull shape parameter
  target += uniform_lpdf(alpha | 0, 10);
}

generated quantities{
  // Relative median
  real RM34 = exp( beta[3] - beta[4] );
}
