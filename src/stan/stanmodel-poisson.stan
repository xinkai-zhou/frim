functions {
  real partial_sum_lpmf(array[] int y, int start, int end, vector eta) {
    return poisson_lupmf(y | eta[start:end]);
  }
}
data {
  int<lower=0>            Nsubs;       // number of individuals (first level)
  int<lower=0>            Nfuns;       // number of observed functions
  int<lower=0>            S;           // number of observations per function
  int<lower=0>            L;           // number of level 1 eigenfunctions
  int<lower=0>            M;           // number of level 2 eigenfunctions
  int<lower=0>            ylen;        // Length of the outcome vector in long format
  int<lower=0>            p;           // number of columns in the fixed effect design matrix
  array[ylen]  int        visit;       // a vector of the visit variable
  array[Nsubs] int        nvisits;     // a vector of the number of repeated functions within individuals
  array[ylen]  int        non_na;      // indices of non-NA values in the original data.
  array[ylen]  int<lower=0> y;         // outcome (binary functional data)
  matrix[ylen, p]         X;           // design matrix for fixed effects
  matrix[p,S]             betaHat;     // smoothed fixed effect estimates.
  row_vector[S]           mu_0;        // population mean, from MFPCA
  matrix[S, max(nvisits)] mu_j;        // visit specific shift from overall mean, from MFPCA
  matrix[L,S]             efuncs_l1;   // level 1 eigen functions
  matrix[M,S]             efuncs_l2;   // level 2 eigen functions
  int<lower=1>            grainsize;
}
parameters {
  matrix[Nsubs, L] xil1;   // level 1 random effects
  matrix[Nfuns, M] xil2;   // level 2 random effects
  vector<lower=0>[L] sigmal1sq; // sd of the random intercept term
  vector<lower=0>[M] sigmal2sq; // sd of the random intercept term
}
transformed parameters {
  matrix[Nsubs, L] xil1_sc;   // level 1 random effects
  matrix[Nfuns, M] xil2_sc;   // level 2 random effects
  for(l in 1:L){
      xil1_sc[,l] = xil1[,l] * sqrt(sigmal1sq[l]);
  }
  for(m in 1:M){
      xil2_sc[,m] = xil2[,m] * sqrt(sigmal2sq[m]);
  }
}
model {
  int inx;
  vector[ylen] eta_vec;
  matrix[Nsubs, S] eta_mat_i;
  matrix[Nfuns, S] eta_mat_ij;

  // priors
  to_vector(xil1) ~ normal(0,1);
  to_vector(xil2) ~ normal(0,1);
  sigmal1sq       ~ inv_gamma(0.01, 0.01);
  sigmal2sq       ~ inv_gamma(0.01, 0.01);

  // New implementation
  eta_mat_i = xil1_sc * efuncs_l1;
  for(i in 1:Nsubs){
     eta_mat_i[i, ] = eta_mat_i[i, ] + mu_0;
  }
  eta_mat_ij = xil2_sc * efuncs_l2;
  inx = 1;
  for(i in 1:Nsubs){
          for(j in 1:nvisits[i]){
                fixef = X[inx, ] * betaHat; // fixed effect at all time points
                eta_mat_ij[inx,] =  fixef + (mu_j[, visit[inx]])' + eta_mat_ij[inx,] + eta_mat_i[i,];
                inx = inx + 1;
          }
  }
  eta_vec = to_vector((eta_mat_ij)')[non_na];

  // likelihood
  // target += reduce_sum(partial_sum_lpmf, y, grainsize, eta_vec + mu_0);
  target += reduce_sum(partial_sum_lpmf, y, grainsize, eta_vec);
}
