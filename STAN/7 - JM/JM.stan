functions{
  vector linear_predictor(matrix X, vector tt, int[] ID, vector betaL, matrix bi){
    int N = num_elements(tt);
    vector[N] out = betaL[1] + bi[ID,1] + betaL[2] * tt + rows_dot_product(bi[ID,2],tt) + betaL[3] * X[ID,2];

    return out;
  } 
}

data{
  int N; 
  int n;
  int NbetasL;
  int NbetasS;
  int Nb;
  vector[N] y;
  vector[N] tt;
  int<lower=1, upper=n> ID[N];
  vector[n] time;
  vector[n] delta;
  matrix[n, NbetasS] X;
  int K;
  vector[K] xk;
  vector[K] wk;
}

parameters{
  vector[NbetasL] betaL;
  vector[NbetasS] betaS;
  real gamma;
  real<lower=0> alpha;
  real<lower=0> tau;
  cov_matrix[Nb] Sigma;
  matrix[n, Nb] bi;
}

transformed parameters{
  real lambda = exp( betaS[1] );
  real sigma = sqrt(1/tau);
}

model{
  vector[N] mu;
  matrix[n, K] HazK;
  vector[n] cumHaz;

  // LINEAR MIXED MODEL SPECIFICATION
  // Linear predictor
  mu = linear_predictor(X, tt, ID, betaL, bi);
  target += normal_lpdf(y | mu, sigma);

  // WEIBULL PROPORTIONAL HAZARD SPECIFICATION
  for(i in 1:n){
     // Hazard function evaluated at Gauss-Legendre quadrature integration points
     for(j in 1:K){
        HazK[i,j] = alpha * pow(time[i] / 2 * (xk[j]+1), alpha-1) * exp( X[i,] * betaS + gamma * ( bi[i,1] + bi[i,2] * time[i] / 2 * (xk[j]+1) ) );
     }
     // Cumulative hazard function with Gauss-Legendre quadrature
     cumHaz[i] = time[i] / 2 * dot_product(wk, HazK[i,]);
       
     target += delta[i] * log(HazK[i,K]) - cumHaz[i];
  }

  // LOG-PRIORS             
  // Coefficients
  target += normal_lpdf(betaL | 0, sqrt(1000));  
  target += normal_lpdf(betaS | 0, sqrt(1000));

  // Association parameter
  target += normal_lpdf(gamma | 0, sqrt(1000));

  // Weibull shape parameter
  target += gamma_lpdf(alpha | 0.01, 0.01);

  // Error precision
  target += gamma_lpdf(tau | 0.01, 0.01);

  // Random-effects variance-covariance matrix
  Sigma ~ inv_wishart(Nb, diag_matrix(rep_vector(1,Nb)));

  // Random-effects
  for(i in 1:n){ target += multi_normal_lpdf(bi[i,1:Nb] | rep_vector(0,Nb), Sigma); }
}   
