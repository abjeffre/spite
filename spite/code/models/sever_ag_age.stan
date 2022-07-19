
functions{
  matrix cov_GPL2(matrix x, real sq_alpha, real sq_rho, real delta) {
    int N = dims(x)[1];
    matrix[N, N] K;
    for (i in 1:(N-1)) {
      K[i, i] = sq_alpha + delta;
      for (j in (i + 1):N) {
        K[i, j] = sq_alpha * exp(-sq_rho * square(x[i,j]) );
        K[j, i] = K[i, j];
      }
    }
    K[N, N] = sq_alpha + delta;
    return K;
  }
}
data{
  int N;
  int region[N];
  int y[N];
  int outgroup[N];
  int village[N];
  int C[N];
  int S[N];
  int LA;
  matrix[LA,LA] DmatAGE;
  int AGE[N];
  int G[N];
  vector[N] IN;
  int ET[N];
  int NET;
}

parameters{
  real a;
  real aV[6];
  real aX[6];
  vector[LA] zA[6];
  real bET[NET];    
  real<lower=0> Sigma_ET;
  real bIN[6];    
  real<lower=0> Sigma_IN;
  real bIN_mu;
  //non-centered cross classifed
  matrix[2,6] zG;
  vector<lower=0>[2] Sigma_GX;
  cholesky_factor_corr[2] L_RhoGX;
  real<lower =0> sigma_X;
  real<lower =0> sigma_V;
  matrix[6,2] zz;
  vector<lower=0>[6] Sigma_zX;
  cholesky_factor_corr[6] L_RhozX;
   //Multilevel Gaussian Parameters
  real<lower=0> etasq_tildeA[6];
  real<lower=0> rhosq_tildeA[6];
  real<lower=0> sigmagp_tildeA[6];
  //group level means
  real<lower=0> etasq_muA[6];
  real<lower=0> rhosq_muA[6];
  real<lower=0> sigmagp_muA[6];
  //linkage sd's
  real<lower=0> etasq_sdA[6];
  real<lower=0> rhosq_sdA[6];
  real<lower=0> sigmagp_sdA[6];
}

transformed parameters{
  matrix[6,2] bG;
  matrix[2,6] bz;
  vector[LA] kA[6];
  matrix[LA, LA] LA_sigma[6];
  matrix[LA,LA] SIGMAA[6];
  real<lower=0> etasqA[6];
  real<lower=0> rhosqA[6];
  real<lower=0> sigmagpA[6];
  // Building non-centered multi-level parameterization for GP
  for (i in 1:6){
    rhosqA[i] = exp(log(rhosq_muA[i]) + rhosq_sdA[i] * rhosq_tildeA[i]);
    etasqA[i] = exp(log(etasq_muA[i]) + etasq_sdA[i] * etasq_tildeA[i]);
    sigmagpA[i] = exp(log(sigmagp_muA[i]) + sigmagp_sdA[i] * sigmagp_tildeA[i]);

  }
  for (i in 1:6){
   SIGMAA[i] = cov_GPL2(DmatAGE, etasqA[i], rhosqA[i], sigmagpA[i]);
   LA_sigma[i] = cholesky_decompose(SIGMAA[i]);
   kA[i] = LA_sigma[i] * zA[i];
    }
 bz = (diag_pre_multiply(Sigma_zX, L_RhozX) * zz)';
 bG = (diag_pre_multiply(Sigma_GX, L_RhoGX) * zG)';
}


model{
  vector[N] p;
  sigma_X ~ exponential(2);
  sigma_V ~ exponential(2);
  a ~ normal(-2,.3);
  aV ~ normal(0,.5);
  aX ~ normal(0,.5);
  L_RhoGX ~ lkj_corr_cholesky( 2 );
  Sigma_GX ~ exponential( 2 );
  to_vector( zG ) ~ normal( 0 , .5 );
  bET ~ normal(0,.5);
  Sigma_ET ~ exponential(2);
  bIN ~ normal(0, .5);
  bIN_mu ~ normal(0, .5);
  Sigma_IN ~ exponential(2);
  //AGE
  for (i in 1:6){
  etasq_sdA[i] ~ normal(0,.25);
  rhosq_sdA[i] ~ normal(0,.5);
  sigmagp_sdA[i] ~ normal(0,.1);
  etasq_tildeA[i] ~ normal(0,.25);
  rhosq_tildeA[i] ~ normal(0,.5);
  sigmagp_tildeA[i] ~ normal(0,.25);
  etasq_muA[i] ~ normal(0, .25);
  rhosq_muA[i] ~ gamma(3,5);
  sigmagp_muA[i] ~ normal(0, .25);
  
  zA[i] ~ normal(0, 1);
  }
	//linear
  L_RhozX ~ lkj_corr_cholesky( 2 );
  Sigma_zX ~ exponential( 1 );
  to_vector( zz ) ~ normal( 0 , 1 );

 for(i in 1:N){
    p[i] = aV[village[i]]*sigma_V + aX[C[i]]*sigma_X +
    bz[S[i], C[i]] + kA[C[i]][AGE[i]] + 
    bG[C[i], G[i]]+ a +
    (bIN[C[i]]*Sigma_IN + bIN_mu)*IN[i] + bET[ET[i]]*Sigma_ET;
  }
  y ~ bernoulli_logit(p);
}
