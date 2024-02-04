functions{
  vector cumulative_baseline_hazard(vector time, vector a, int[] intobs, real[] lambda){
    int n = num_elements(time);
    vector[n] out;

    for(i in 1:n){
       vector[intobs[i]] HH;
       for(k in 1:intobs[i]){
          real cond = step(time[i] - a[k + 1]);
          HH[k] = cond * (a[k + 1] - a[k]) * lambda[k] + (1 - cond) * (time[i] - a[k]) * lambda[k];
       }
       // Cumulative hazard function
       out[i] = sum(HH);
    }
    return out;
  } 
}

data{
  int n;
  int m;
  int Nbetas;
  vector[n] time;
  vector[n] delta;
  matrix[n, Nbetas] X;
  vector[m+1] a;
  int intobs[n];
}

parameters{  
  vector[Nbetas] beta;
  real<lower=0> lambda[m];
}

model{
  vector[n] logHaz;
  vector[n] cumHaz;
  vector[n] H = cumulative_baseline_hazard(time, a, intobs, lambda); 

  // PROPORTIONAL HAZARD WITH PIECEWISE CONSTANT SPECIFICATION
  for(i in 1:n){
     // Log-hazard function
     logHaz[i] = log(lambda[intobs[i]]) + X[i,] * beta;
     // Cumulative hazard function
     cumHaz[i] = H[i] * exp( X[i,] * beta );

     target += delta[i] * logHaz[i] - cumHaz[i];
  }

  // LOG-PRIORS             
  // Coefficients
  target += normal_lpdf(beta | 0, sqrt(1000));  

  // Piecewise constants
  target += gamma_lpdf(lambda | 0.01, 0.01);
}
