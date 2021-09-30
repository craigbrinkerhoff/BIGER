//Craig BIKER stan model for remotely sensing gas transfer velocity
//Fall 2021

functions {
  // Conversion from array to vector takes place by row.
  // Nested elements are appended in sequence.

  // Convert an array to a vector based on a binary matrix
  // indicating non-missing data
  vector ragged_vec(vector[] x, int[,] bin) {
    vector[num_elements(x)] out;
    int ind;

    ind = 1;
    for (i in 1:size(x)) {
      for (t in 1:num_elements(x[1])) {
        if (bin[i, t] == 1) {
          out[ind] = x[i, t];
          ind += 1;
        }
      }
    }
    // print(out);
    return(out[1:(ind - 1)]);
  }

  // Repeat elements of a "row" vector to match with 2-D array vectorization
  vector ragged_row(vector x, int[,] bin) {
    vector[num_elements(bin)] out;
    int ind;

    ind = 0;
    for (i in 1:size(bin)) {
      for (t in 1:num_elements(bin[1])) {
        if (bin[i, t] == 1) {
          ind += 1;
          out[ind] = x[t];
        }
      }
    }
    return(out[1:ind]);
  }

  // Repeat elements of a "column" vector to match with 2-D array vectorization
  vector ragged_col(vector x, int[,] bin) {
    vector[num_elements(bin)] out;
    int ind;

    ind = 0;
    for (i in 1:size(bin)) {
      for (t in 1:num_elements(bin[1])) {
        if (bin[i, t] == 1) {
          ind += 1;
          out[ind] = x[i];
        }
      }
    }
    return(out[1:ind]);
  }

  // indices of vectorized bin1 that are also in vectorized bin2
  int[] commoninds(int[,] bin1, int[,] bin2) {
    int out[num_elements(bin1)];
    int vecinds[size(bin1), num_elements(bin1[1])];
    int ctr;
    int ind;

    ctr = 0;
    for (i in 1:size(bin1)) {
      for (t in 1:num_elements(bin1[1])) {
        if (bin1[i, t] == 1) {
          ctr += 1;
          vecinds[i, t] = ctr;
        }
      }
    }

    ind = 0;
    for (i in 1:size(vecinds)) {
      for (t in 1:size(vecinds[1])) {
        if (bin2[i, t] == 1) {
          ind += 1;
          out[ind] = vecinds[i, t];
        }
      }
    }
    return(out[1:ind]);
  }
}

//Setup data object
data {

  // Options
  int<lower=0, upper=1> inc; // Include Manning? 0=no; 1=yes;
  int<lower=0, upper=1> meas_err; //0=no, 1=yes;

  // Dimensions
  int<lower=0> nx; // number of reaches
  int<lower=0> nt; // number of times
  int<lower=0> ntot_man; // total number of non-missing Manning observations

  // // Missing data
  int<lower=0,upper=1> hasdat_man[nx, nt]; // matrix of 0 (missing), 1 (not missing)

  // *Actual* data
  vector[nt] Wobs[nx]; // measured widths, including placeholders for missing
  vector[nt] Sobs[nx]; // measured slopes
  vector[nt] dAobs[nx]; // measured partial area
  vector[nx] dA_shift; // adjustment from min to median

  real<lower=0> Serr_sd; //measurement error slopes
  real<lower=0> Werr_sd;
  real<lower=0> dAerr_sd; //measurement error slopes

  // Hard bounds on parameters
  real lowerbound_A0; // These must be scalars, unfortunately.
  real upperbound_A0;
  real upperbound_logk;
  real lowerbound_logk;

  // *Known* likelihood parameters
  vector<lower=0>[nt] sigma_post[nx]; // Manning error standard deviation

  // Hyperparameters
  vector[nt] logk_hat; // prior mean on logF g/m2*dy
  real logA0_hat[nx]; //space-varying A0 prior m2

  vector<lower=0>[nt] logk_sd;
  real<lower=0> logA0_sd[nx];
}

//set up transformed parameters
transformed data {
  // Transformed data are *vectors*, not arrays. This to allow ragged structure

  vector[nt] dApos_array[nx];

  vector[ntot_man] Wobsvec_man;
  vector[ntot_man] Sobsvec_man;

  vector[ntot_man] logWobs_man;
  vector[ntot_man] logSobs_man;
  vector[ntot_man] dApos_obs;
  vector[ntot_man] sigmavec_man;

  int ntot_w; // how many widths in likelihood: equal to ntot_man unless inc_a
  ntot_w = ntot_man;

  for (i in 1:nx) {
    dApos_array[i] = dAobs[i] - min(dAobs[i]); // make all dA positive
  }

  // convert pseudo-ragged arrays to vectors
  Wobsvec_man = ragged_vec(Wobs, hasdat_man);
  Sobsvec_man = ragged_vec(Sobs, hasdat_man);
  dApos_obs = ragged_vec(dApos_array, hasdat_man);

  logWobs_man = log(Wobsvec_man);
  logSobs_man = log(Sobsvec_man);

  sigmavec_man = ragged_vec(sigma_post, hasdat_man);
}

//Set up surface-level parameters
parameters {
  vector<lower=lowerbound_logk,upper=upperbound_logk>[nt] logk;
  vector<lower=lowerbound_A0,upper=upperbound_A0>[nx] A0[inc];

  vector<lower=0>[ntot_man] Sact[meas_err * inc];
  vector<lower=0>[ntot_w] Wact[meas_err * inc];
  vector[ntot_man] dApos_act[meas_err * inc];
}

//Calculated transformed parameters
transformed parameters {

  vector[ntot_man] eq_lhs[inc]; // LHS for Manning likelihood
  vector[ntot_man] logA_man[inc]; // log area for Manning's equation
  vector[ntot_man] logk_man[inc]; // location-repeated logk
  vector[ntot_man] eq_rhs[inc]; // RHS for Manning likelihood

  // Params
  if (inc) {
    if (meas_err) { //Measurement error in slopes and heights
      logA_man[1] = log(ragged_col(A0[1], hasdat_man) + dApos_act[1]);
      logk_man[1] = ragged_row(logk, hasdat_man);

      //Brinkerhoff k600~Ustar model
      // man_lhs[1] = log(56.0294) + 0.5*log(9.8) + 0.5*log(Sact[1]) - 0.5*log(Wact[1]);
      //  man_rhs[1] = logk_man[1] - 0.5*logA_man[1];

      //man_lhs[1] = 0.3997133*logWobs_man - 0.899355*log(Sact[1]) - log(85.10025) - 0.59957*log(9.8);
      //man_rhs[1] = 0.3997133*(logA_man[1]) - 0.59957*logN_man[1] - logk_man[1];

      //Brinkehoff implementation of Moog & Jirka's 'chainsaw model'
      eq_lhs[1] = log(76.4) + (0.5625)*log(9.8) + (0.5625)*log(Sact[1]) - (0.6875)*log(Wact[1]);
      eq_rhs[1] = logk_man[1] - (0.6875)*logA_man[1];
    }

    else { //No measurement error in slopes and heights
      logk_man[1] = ragged_row(logk, hasdat_man);
      logA_man[1] = log(ragged_col(A0[1], hasdat_man) + dApos_obs);

      //Brinkerhoff k600~Ustar model
      // man_lhs[1] = log(56.0294) + 0.5*log(9.8) + 0.5*logSobs_man - 0.5*logWobs_man;
      // man_rhs[1] = logk_man[1] - 0.5*logA_man[1];

      //   man_lhs[1] = 0.3997133*logWobs_man - 0.899355*logSobs_man - log(85.10025) - 0.59957*log(9.8);
      //   man_rhs[1] = 0.3997133*(logA_man[1]) - 0.59957*logN_man[1] - logk_man[1];

      //Brinkehoff implementation of Moog & Jirka's 'chainsaw model'
      eq_lhs[1] = log(76.4) + (0.5625)*log(9.8) + (0.5625)*logSobs_man - (0.6875)*logWobs_man;
      eq_rhs[1] = logk_man[1] - (0.6875)*logA_man[1];
    }
  }
}

//Actual Bayesian inference model
model {
  // Priors
  if (inc) {
    A0[1] + dA_shift[1] ~ lognormal(logA0_hat, logA0_sd);
    logk[1] ~ normal(logk_hat, logk_sd);
  }

  //likelihood/sampling model
  if (inc){
    eq_lhs[1] ~ normal(eq_rhs[1], sigmavec_man);
  }

  //latent variables for measurement error
  if (meas_err){
    if (inc) {
      Wact[1] ~ normal(Wobsvec_man, Werr_sd); // W meas err
      Sact[1] ~ normal(Sobsvec_man, Serr_sd); // S meas err
      dApos_act[1] ~ normal(dApos_obs, dAerr_sd); // dA meas err

      target += -log(Wact[1]); // Jacobian adjustments
      target += -log(Sact[1]);
    }
  }
}
