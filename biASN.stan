// Bivariate ASN model
functions{
  real biASN_lpdf(matrix Y, vector mu0, matrix sigma0, vector alpha0){
    vector[rows(Y)] prob;
    real lprob;
    for (i in 1:rows(Y)){
      prob[i] = log( (1-alpha0[1]*((Y[i,1]-mu0[1])/sigma0[1,1])-alpha0[2]*((Y[i,2]-mu0[2])/sigma0[2,2]))^2 +1 ) 
      - log(2+alpha0[1]^2+alpha0[2]^2+2*alpha0[1]*alpha0[2]*sigma0[1,2]) + multi_normal_lpdf(Y[i,:] | mu0, sigma0);
    }
    lprob = sum(prob);
    return lprob;
  }
}

data { 
  int<lower=0> N;
  matrix[N,2] Y;
}
parameters {
  vector[2] mu;
  vector[2] alpha;
  vector<lower=0>[2] lambda;
  real<lower=-1,upper=1> r;
} 
transformed parameters {
  vector<lower=0>[2] sigma;
  cov_matrix[2] T;

  // Reparameterization
  sigma[1] = inv_sqrt(lambda[1]);
  sigma[2] = inv_sqrt(lambda[2]);
  T[1,1] = square(sigma[1]);
  T[1,2] = r * sigma[1] * sigma[2];
  T[2,1] = r * sigma[1] * sigma[2];
  T[2,2] = square(sigma[2]);
}
model {
  // Data
  Y ~ biASN(mu, T, alpha);
}
