// Legend:
// C: cured/immune subpopulation
// U: uncured/susceptible subpopulation

data{
  int n;
  int Nbetas;
  vector[n] time;
  vector[n] delta;
  matrix[n, Nbetas] XC;
  matrix[n, Nbetas-1] XU;
}

parameters{  
  vector[Nbetas] betaC;
  vector[Nbetas-1] betaU;
  real<lower=0> lambda;
  real<lower=0> alpha;
}

model{
  vector[n] logHaz;
  vector[n] cumHaz;
  vector[n] eta; 

  // MIXTURE CURE SPECIFICATION
  for(i in 1:n){
     // Logistic regression model (cured subpopulation)
     eta[i] = inv_logit( XC[i,] * betaC );

     // Proportional hazard model (uncured subpopulation)
     // Log-hazard function
     logHaz[i] = log(lambda) + log(alpha) + (alpha - 1) * log(time[i]) + XU[i,] * betaU;
     // Cumulative hazard function
     cumHaz[i] = lambda * pow(time[i], alpha) * exp( XU[i,] * betaU );

     target += delta[i] * (log(1 - eta[i]) + logHaz[i] - cumHaz[i]) + 
               (1 - delta[i]) * log(eta[i] + (1 - eta[i]) * exp( - cumHaz[i] ));
  }

  // LOG-PRIORS             
  // Coefficients
  target += normal_lpdf(betaC | 0, sqrt(1000));
  target += normal_lpdf(betaU | 0, sqrt(1000));

  // Weibull scale and shape parameters
  target += gamma_lpdf(lambda | 0.01, 0.01);
  target += uniform_lpdf(alpha | 0, 10);
}
