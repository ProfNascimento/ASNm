//
// This Stan program defines a Alpha-Skew Normal (ASN) model, with a
// vector of values 'y' modeled as modified normally distributed
// with mean 'mu', standard deviation 'sigma' and asymmetry 'alpha'.
//

functions{
  real ASN_lpdf(vector x, real mu0, real sigma0, real alpha0){
    vector[num_elements(x)] prob;
    real lprob;
    for (i in 1:num_elements(x)){
      prob[i] = log( (1-alpha0*((x[i]-mu0)/sigma0))^2+1 ) - log((2+alpha0^2)*sigma0) + normal_lpdf(x[i] | mu0, sigma0);
    }
    lprob = sum(prob);
    return lprob;
  }
}

// The input data is a vector 'y' of length 'N'.
data {
  int <lower=0> N;
  vector[N] Y;
}

// The parameters accepted by the model. Our model
// accepts three parameters 'mu', 'sigma' and 'alpha'.
parameters {
  real mu;
  real<lower=0> sigma;
  real alpha;
}

// The model to be estimated. We model the output
// 'y' to be ASN distributed with mean 'mu',
// standard deviation 'sigma' and asymmetry 'alpha'.
model {
  //PRIORI
  mu ~ normal(0,15);
  sigma ~ gamma(0.001,0.001);
  alpha ~ normal(0,15);
  //LIKELIHOOD
  Y ~ ASN(mu, sigma, alpha);
}

// GENERATY RANDOM NUMBERS
generated quantities {
  vector[N] y_sim;
  real x;
  for(i in 1:N){
   x = uniform_rng(min(Y), max(Y));  
   y_sim[i] = (( (1-alpha*((x-mu)/sigma))^2+1 )/((2+alpha^2)*sigma)) * exp(normal_lpdf(x | mu, sigma)); 
  }
}