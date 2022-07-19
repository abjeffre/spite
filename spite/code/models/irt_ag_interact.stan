
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
    int M1; //n questions drought
    int N;
    int Y[N,M1];//answers drought
    int region[N];
    int y[N];
    int outgroup[N];
    int village[N];
    int LA;
    matrix[LA,LA] DmatAGE;
    int AGE[N];
    int C[N];
    int G[N];
    vector[N] IN;
    int ET[N];
    int NET;
}

parameters{
    real av[6];
    vector[LA] zA[6];
    
    real mu;
    real bET[NET];    
    real<lower=0> Sigma_ET;
    
    real bIN[6];    
    real<lower=0> Sigma_IN;
    real bIN_mu;
    
    matrix[2,6] zG;
    vector<lower=0>[2] Sigma_GX;
    cholesky_factor_corr[2] L_RhoGX;
    
    
    real<lower = 0> sigma_V;
    vector<lower=0>[M1] a1;//discrimination drought
	  vector[M1] b1;//difficulty Drought
	  vector[N] theta;//latent trait
    matrix[2,3] zz;
    vector<lower=0>[2] Sigma_zX;
    cholesky_factor_corr[2] L_RhozX;
    matrix[2,3] zz2;
    vector<lower=0>[2] Sigma_zX2;
    cholesky_factor_corr[2] L_RhozX2;
    matrix[2,3] az ;
    vector<lower=0>[2] Sigma_aX;
    cholesky_factor_corr[2] L_RhoaX;
    
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
matrix[3,2] bz;
matrix[3,2] bz2;
matrix[3,2] aX;
vector[LA] kA[6];
matrix[LA, LA] LA_sigma[6];
matrix[LA,LA] SIGMAA[6];
matrix[6,2] bG;


real<lower=0> etasqA[6];
real<lower=0> rhosqA[6];
real<lower=0> sigmagpA[6];

 aX = (diag_pre_multiply(Sigma_aX, L_RhoaX) * az)';
 bz = (diag_pre_multiply(Sigma_zX, L_RhozX) * zz)';
 bz2 = (diag_pre_multiply(Sigma_zX2, L_RhozX2) * zz2)';
 bG = (diag_pre_multiply(Sigma_GX, L_RhoGX) * zG)';





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

}


model{
     vector[N] p;
     
     bET ~ normal(0,.5);
     Sigma_ET ~ exponential(2);
     bIN ~ normal(0, .25);
     bIN_mu ~ normal(0, .5);
     Sigma_IN ~ exponential(2);
     
      mu ~ normal(-2.5, .5);
      sigma_V ~ exponential(2);
      av ~ normal(0, .5);
       //irt
      a1 ~ lognormal(0, 0.5); //value constrained above zero
	    b1 ~ normal(0,1);
	    theta ~ normal(0,1);
	    //linear
      L_RhozX ~ lkj_corr_cholesky( 2 );
      Sigma_zX ~ exponential( 1 );
      to_vector( zz ) ~ normal( 0 , 1 );
      
      L_RhozX2 ~ lkj_corr_cholesky( 3 );
      Sigma_zX2 ~ exponential( 1.1 );
      to_vector( zz2 ) ~ normal( 0 , .5 );
    
     L_RhoaX ~ lkj_corr_cholesky( 2 );
     Sigma_aX ~ exponential( 1 );
     to_vector( az ) ~ normal( 0 , 1 );
     
             
    //Aroughts
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
    
     L_RhoGX ~ lkj_corr_cholesky( 2 );
     Sigma_GX ~ exponential( 2 );
     to_vector( zG ) ~ normal( 0 , .5 );
     
     //irt models
      for ( i in 1:N ) {
    		for (j in 1:M1 ) {
    			real p1 = inv_logit(a1[j]*(theta[i]-b1[j]));
    			Y[i,j] ~ bernoulli( p1 );
    		}
    	}

     for(i in 1:N){
        p[i] =  av[village[i]]*sigma_V + mu +
        aX[region[i], outgroup[i]] +
        bz[region[i], outgroup[i]]*theta[i] +
        bz2[region[i], outgroup[i]]*theta[i]*IN[i] +
        kA[C[i]][AGE[i]] +
        bG[C[i], G[i]]+
        (bIN[C[i]]*Sigma_IN + bIN_mu)*IN[i] +
        bET[ET[i]]*Sigma_ET;
      }

    y ~ bernoulli_logit(p);

}
