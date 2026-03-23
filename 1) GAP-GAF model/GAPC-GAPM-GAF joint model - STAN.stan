/*
  NOTATION:

  p = pressure
  f = frequency

  g = gulf
  y = year
  m = mother-calf

  X = design matrices (observation level)
  Z = prediction matrices (year level)
  V = prediction matrices (average level)

*/

data {

  int<lower=0> N_f;      // f for frequency
  int<lower=0> N_pred_f;
  int<lower=0> N_p;      // p for pressure (attack number)
  int<lower=0> N_pred_p;

  int with[N_f];    // intervals with attacks
  int total[N_f];  // intervals observed
  int gan[N_p];    // gull attack number (for pressure)

  vector[N_p] hours; // hours observed (to compute pressure)

  // Number of
  int N_g; // gulfs
  int N_gm; // gulf * mc combinations
  int N_gy; // gulf * year combinations
  int N_gmy; // gulf * year * mc combinations

  int gulf_id_gy[N_pred_f]; // integer gulf identifier [1:2], long as year * gulf

  // Design matrices
      // gulf matrices
  matrix[N_p, N_g] X_g_p;
  matrix[N_f, N_g] X_g_f;
    // gulf - year matrix
  matrix[N_p, N_gy] X_gy_p;
  matrix[N_f, N_gy] X_gy_f;
    // gulf-mc matrix (for gap only)
  matrix[N_p, N_gm] X_gm_p;
    // gulf-mc-year matrix (for gap only)
  matrix[N_p, N_gmy] X_gmy_p;

  // Prediction matrices (year level)
    // gap
  matrix[N_pred_p, N_gm] Z_gm_p; // gulf * mc
  matrix[N_pred_p, N_gy] Z_gy_p; // gulf * year
  matrix[N_pred_p, N_gmy] Z_gmy_p; // gulf * mc * year
    //gaf
  matrix[N_pred_f, N_g] Z_g_f;
  matrix[N_pred_f, N_gy] Z_gy_f;

  // Prediction matrices (average across years)
    // gap
  matrix[N_gm, N_gm] V_gm_p;
  matrix[N_gm, N_g] V_g_p;

    //gaf
  matrix[N_g, N_g] V_g_f;

  // priors
  real prior_mu_f_sd;
  real prior_mu_p_sd;
  real prior_mu_p_mean;
  real prior_sd_gy_f_sd;
  real prior_sd_gy_p_sd;
  real prior_sd_gmy_sd;
  real prior_sd_obs_sd;
  real prior_scale_sd;

}

parameters {

  /*
  NOTATION:
  mu for fixed effects (means at linear predictor scale)
  e for error terms (random effects)
  */

  // means
  vector[N_g] mu_f_raw;
  vector[N_gm] mu_p_raw;

  // year * gulf random effects
  matrix[N_gy, 2] e_gy_raw;
  /* these are correlated random effects. column 1 is for frequency (logit),
  column 2 is for pressure (log) */

  // year * mc * gulf random effect (only for gap)
  vector[N_gmy] e_gmy_raw;

  // observation-level random effect for gap
  //vector[N_f-2] e_obs_f_raw;
  vector[N_f] e_obs_f_raw;

  // standard deviations
  matrix<lower = 0>[N_g, 2] sd_gy_raw; // 2 for gulfs and 2 for gap and gaf
  vector<lower = 0>[N_gm] sd_gmy_raw;
  real<lower = 0> sd_obs_raw; // obs level random effect for gaf (logit)
  real<lower = 0> scale_raw;  // scale for negative binomial (gap)

  // correlation among gap and gap
  real<lower = -1, upper = 1> rho;

}

transformed parameters {

  // mus
  vector[N_g] mu_f = mu_f_raw * prior_mu_f_sd;
  vector[N_gm] mu_p = mu_p_raw * prior_mu_p_sd + prior_mu_p_mean;

  // sds
  matrix<lower = 0>[N_g, 2] sd_gy;//  N_g for gulfs and 2 for gap and gaf
  vector<lower = 0>[N_gm] sd_gmy = sd_gmy_raw * prior_sd_gmy_sd;
  real<lower = 0> sd_obs = sd_obs_raw * prior_sd_obs_sd; // obs level random effect for gaf (logit)
  real<lower = 0> scale = scale_raw * prior_scale_sd;
  real<lower = 0> phi = 1 / scale;

  // random effects

  // year * gulf random effects
  matrix[N_gy, 2] e_gy;
  // year * mc * gulf random effect (only for gap)
  vector[N_gmy] sd_gmy_long = Z_gm_p * sd_gmy; /* this lengthens the sd vectors
  from N_gm to N_gmy, to match e_gmy element-wise */
  vector[N_gmy] e_gmy = e_gmy_raw .* sd_gmy_long;

  // (Define variables to compute later)
  // observation-level random effect for gaf
  vector[N_f] e_obs_f;

  // correlation matrix for gaf and gap
  matrix[2, 2] rho_mat;
  matrix[2, 2] rho_chol;

  // expected values and sds
  vector[N_f] gaf_hat_logit;
  vector[N_p] gap_hat;
  vector[N_g] gaf_sds_logit;
  vector[N_g] gaf_means_logit;
  vector[N_gm] gap_sds;
  vector[N_gm] gap_means;

  // Predictions at year level
  vector[N_pred_f] gaf_pred; // at logit scale
  vector[N_pred_p] gap_pred; // at response scale

  // observation-level random effect for gaf
  /*
  e_obs_f[1:2] = rep_vector(0, 2); // first two years are taken as having only 1 day
  e_obs_f[3:N_f] = e_obs_f_raw * sd_obs;
  */
  e_obs_f = e_obs_f_raw * sd_obs;

  // sd among years by gulf
  sd_gy[, 1] = sd_gy_raw[, 1] * prior_sd_gy_f_sd; // gaf
  sd_gy[, 2] = sd_gy_raw[, 2] * prior_sd_gy_p_sd; // gap

  // correlated random effects
    // corr matrix
  for (i in 1:2)
    rho_mat[i, i] = 1;
  rho_mat[2, 1] = rho;
  rho_mat[1, 2] = rho;
  rho_chol = cholesky_decompose(rho_mat);

  // get correlated random effect for each gulf-year.
  for(i in 1:N_gy) {
    int g = gulf_id_gy[i];
    e_gy[i, ] = e_gy_raw[i, ] * diag_pre_multiply(sd_gy[g, ], rho_chol);
  }
  /*Pre-multiplying the -corr matrix- Cholesky factor by the scale produces the
  Cholesky factor of the final covariance matrix*/

  // Expected values

  // gaf for posterior predictive check (unconditional to observation)
  gaf_hat_logit = X_g_f * mu_f + X_gy_f * e_gy[, 1];

  // gap for pp check and likelihood
  gap_hat = exp(X_gm_p * mu_p +
                X_gy_p * e_gy[, 2] +
                X_gmy_p * e_gmy +
                log(hours)); // offset term

  // Standard deviations and means unconditional to year

  // gaf variability among years and among observations (to compute the gulf means)
  // at logit scale
  gaf_sds_logit = sqrt(sd_gy[, 1] .* sd_gy[, 1] + sd_obs ^ 2);
  gaf_means_logit = mu_f;

  // gap variability among years|gulf and among years|(gulf * mc) (to compute the gulf-mc means)
  {
    vector[N_gm] sd_gm = V_gm_p * sd_gmy;
    vector[N_gm] sd_g = V_g_p * sd_gy[, 2];
    real sigma2;

    for(i in 1:N_gm) {
      // mean at exp scale
      gap_means[i] = exp(mu_p[i] + 0.5 * (sd_gm[i] ^ 2 + sd_g[i] ^ 2));

      // variance in log scale is sum of variances
      sigma2 = sd_gm[i] ^ 2 + sd_g[i] ^ 2;
      // variance in exp scale is:
      gap_sds[i] = sqrt((exp(sigma2) - 1) * exp(2 * mu_p[i] + sigma2));
    }
  }

  // Predictions at year level
  gaf_pred = Z_g_f * mu_f + Z_gy_f * e_gy[, 1]; // logit scale
  gap_pred = exp(Z_gm_p * mu_p + Z_gy_p * e_gy[, 2] + Z_gmy_p * e_gmy); // response scale

}

model {

  // Priors
  mu_f_raw ~ std_normal(); // fixed effects for gaf
  mu_p_raw ~ std_normal(); // fixed effects for gap
  to_vector(e_gy_raw) ~ std_normal(); // gulf-year random effect for gaf and gap
  e_gmy_raw ~ std_normal(); // gulf-mc-year random effect for gap
  e_obs_f_raw ~ std_normal(); // obs level random effect for gaf (logit)
  to_vector(sd_gy_raw) ~ std_normal(); // gulf-year random effect for gaf and gap
  sd_gmy_raw ~ std_normal(); // gulf-mc-year random effect for gap
  sd_obs_raw ~ std_normal(); // obs level random effect for gaf (logit)
  scale_raw ~ std_normal();  // scale for negative binomial (gap)
  // rho is given an implicit flat prior on (-1, 1)

  // Likelihood for gaf
  {
    vector[N_f] gaf_hat = inv_logit(X_g_f * mu_f + X_gy_f * e_gy[, 1] + e_obs_f);
    with ~ binomial(total, gaf_hat);
  }

  // Likelihood for gap
  gan ~ neg_binomial_2(gap_hat, phi);

}