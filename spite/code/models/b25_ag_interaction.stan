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
    int y[N];
    int C[N];
    int village[N];
    int LD;
    matrix[LD,LD] DmatDROUGHT;
    int DROUGHT[N];
    int LA;
    matrix[LA,LA] DmatAGE;
    int AGE[N];
    int TA[N];
    int TB[N];
    int TC[N];
    int TD[N];
    vector[N]IN;
    int G[N];
    int NET;
    int ET[N];
  
}

parameters{
    real aV[6];    
    real village_mu;
    real<lower=0> Sigma_V;
    
    real aX[6];    
    real<lower=0> Sigma_CC;
    
    real bET[NET];    
    real<lower=0> Sigma_ET;
    
    real bIN[6];    
    real<lower=0> Sigma_IN;
    real bIN_mu;
    
    real bA[2];    
    real<lower=0> Sigma_A;
    
    real bB[2];    
    real<lower=0> Sigma_B;
    
    real bC[2];    
    real<lower=0> Sigma_C;
    

    real bD[2];    
    real<lower=0> Sigma_D;
    
    matrix[2,6] zG;
    vector<lower=0>[2] Sigma_GX;
    cholesky_factor_corr[2] L_RhoGX;
    
    
  //GAUSSIAN PROCESS 
  vector[LD] zD[6];
  vector[LA] zA[6];
  vector[LD] zD2[6];
  
  //Multilevel Gaussian Parameters
  
  real<lower=0> etasq_tildeD[6];
  real<lower=0> rhosq_tildeD[6];
  real<lower=0> sigmagp_tildeD[6];
    
  
  
  //group level means
  real<lower=0> etasq_muD[6];
  real<lower=0> rhosq_muD[6];
  real<lower=0> sigmagp_muD[6];
  

  
  //linkage sd's
  real<lower=0> etasq_sdD[6];
  real<lower=0> rhosq_sdD[6];
  real<lower=0> sigmagp_sdD[6];
     
  //Multilevel Gaussian Parameters
  
  real<lower=0> etasq_tildeD2[6];
  real<lower=0> rhosq_tildeD2[6];
  real<lower=0> sigmagp_tildeD2[6];
    
  
  
  //group level means
  real<lower=0> etasq_muD2[6];
  real<lower=0> rhosq_muD2[6];
  real<lower=0> sigmagp_muD2[6];
  

  
  //linkage sd's
  real<lower=0> etasq_sdD2[6];
  real<lower=0> rhosq_sdD2[6];
  real<lower=0> sigmagp_sdD2[6];
          


  
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

  // GAUSSIAN PROCESS
   matrix[6,2] bG;
   vector[LD] kD[6];
   vector[LD] kD2[6];
   vector[LA] kA[6];
   matrix[LD, LD] LD_sigma[6];
   matrix[LD, LD] LD_sigma2[6];
   matrix[LA, LA] LA_sigma[6];
   matrix[LD,LD] SIGMAD[6];
   matrix[LD,LD] SIGMAD2[6];
   matrix[LA,LA] SIGMAA[6];


   real<lower=0> etasqD[6];
   real<lower=0> rhosqD[6];
   real<lower=0> sigmagpD[6];
   
   real<lower=0> etasqD2[6];
   real<lower=0> rhosqD2[6];
   real<lower=0> sigmagpD2[6];

   real<lower=0> etasqA[6];
   real<lower=0> rhosqA[6];
   real<lower=0> sigmagpA[6];


  // Building non-centered multi-level parameterization for GP

    for (i in 1:6){
                rhosqD[i] = exp(log(rhosq_muD[i]) + rhosq_sdD[i] * rhosq_tildeD[i]);
                etasqD[i] = exp(log(etasq_muD[i]) + etasq_sdD[i] * etasq_tildeD[i]);
                sigmagpD[i] = exp(log(sigmagp_muD[i]) + sigmagp_sdD[i] * sigmagp_tildeD[i]);
                rhosqD2[i] = exp(log(rhosq_muD2[i]) + rhosq_sdD2[i] * rhosq_tildeD2[i]);
                etasqD2[i] = exp(log(etasq_muD2[i]) + etasq_sdD2[i] * etasq_tildeD2[i]);
                sigmagpD2[i] = exp(log(sigmagp_muD2[i]) + sigmagp_sdD2[i] * sigmagp_tildeD2[i]);
                rhosqA[i] = exp(log(rhosq_muA[i]) + rhosq_sdA[i] * rhosq_tildeA[i]);
                etasqA[i] = exp(log(etasq_muA[i]) + etasq_sdA[i] * etasq_tildeA[i]);
                sigmagpA[i] = exp(log(sigmagp_muA[i]) + sigmagp_sdA[i] * sigmagp_tildeA[i]);  
          }

    for (i in 1:6){
                   SIGMAD[i] = cov_GPL2(DmatDROUGHT, etasqD[i], rhosqD[i], sigmagpD[i]);
                   LD_sigma[i] = cholesky_decompose(SIGMAD[i]);
                   kD[i] = LD_sigma[i] * zD[i];
                   SIGMAD2[i] = cov_GPL2(DmatDROUGHT, etasqD2[i], rhosqD2[i], sigmagpD2[i]);
                   LD_sigma2[i] = cholesky_decompose(SIGMAD2[i]);
                   kD2[i] = LD_sigma2[i] * zD2[i];

                   SIGMAA[i] = cov_GPL2(DmatAGE, etasqA[i], rhosqA[i], sigmagpA[i]);
                   LA_sigma[i] = cholesky_decompose(SIGMAA[i]);
                   kA[i] = LA_sigma[i] * zA[i];
  }
  bG = (diag_pre_multiply(Sigma_GX, L_RhoGX) * zG)';

}


model{
     vector[N] p;
     
     
     bET ~ normal(0,.5);
     Sigma_ET ~ exponential(2);
     bIN ~ normal(0, .5);
     bIN_mu ~ normal(0, 1);
     Sigma_IN ~ exponential(2);
     aV ~ normal(0, .5);
     village_mu ~ normal(-2, .5);
     Sigma_V ~ exponential(2);
     aX ~ normal(0, 1);
     Sigma_CC ~ exponential(1);
     bA ~ normal(0, 1);
     Sigma_A ~ exponential(1);
     bB ~ normal(0, 1);
     Sigma_B ~ exponential(1);
     bC ~ normal(0, 1);
     Sigma_C ~ exponential(1);
     bD ~ normal(0, 1);
     Sigma_D ~ exponential(1);
     L_RhoGX ~ lkj_corr_cholesky( 2 );
     Sigma_GX ~ exponential( 1 );
     to_vector( zG ) ~ normal( 0 , 1 );
     
    //Gaussian Process Priors 
    for(i in 1:6){
        //Droughts
        etasq_sdD[i] ~ normal(0,.75);
        rhosq_sdD[i] ~ normal(0,.1);
        sigmagp_sdD[i] ~ normal(0,.1);
        
        etasq_tildeD[i] ~ normal(0,.75);
        rhosq_tildeD[i] ~ normal(0,.1);
        sigmagp_tildeD[i] ~ normal(0,.1);
        
        etasq_muD[i] ~ normal(0, .75);
        rhosq_muD[i] ~ gamma(4,5);
        sigmagp_muD[i] ~ normal(0, .1);
        
        // Interaction
        etasq_sdD2[i] ~ normal(0,.75);
        rhosq_sdD2[i] ~ normal(0,.1);
        sigmagp_sdD2[i] ~ normal(0,.1);
        
        etasq_tildeD2[i] ~ normal(0,.75);
        rhosq_tildeD2[i] ~ normal(0,.1);
        sigmagp_tildeD2[i] ~ normal(0,.1);
        
        etasq_muD2[i] ~ normal(0, .75);
        rhosq_muD2[i] ~ gamma(4,5);
        sigmagp_muD2[i] ~ normal(0, .1);
        
        //AGE
        etasq_sdA[i] ~ normal(0,.5);
        rhosq_sdA[i] ~ normal(0,.75);
        sigmagp_sdA[i] ~ normal(0,.2);
        
        etasq_tildeA[i] ~ normal(0,.25);
        rhosq_tildeA[i] ~ normal(0,.5);
        sigmagp_tildeA[i] ~ normal(0,.25);
        
        etasq_muA[i] ~ normal(0, .25);
        rhosq_muA[i] ~ gamma(3,5);
        sigmagp_muA[i] ~ normal(0, .25);

    zD[i] ~ normal(0, 1);
    zD2[i] ~ normal(0, 1);
    zA[i] ~ normal(0, 1);

    }
    

     for(i in 1:N){
        p[i] =  village_mu + 
        aV[village[i]]*Sigma_V +
        aX[C[i]]*Sigma_CC + 
        bG[C[i], G[i]]+
        (bIN[C[i]]*Sigma_IN + bIN_mu)*IN[i] +
        kD[C[i]][DROUGHT[i]] +
        kD2[C[i]][DROUGHT[i]]*IN[i] +
        kA[C[i]][AGE[i]] +
        bET[ET[i]]*Sigma_ET;// +
        //bA[TA[i]]*Sigma_A +
        //bB[TB[i]]*Sigma_B +
        //bC[TC[i]]*Sigma_C +
        //bD[TD[i]]*Sigma_D ;
      }

    y ~ bernoulli_logit(p);
}
