//4 class stanmodel 'psma' coded as follows:
functions {
  real[] critical_value(real [] alpha_desc, int side) {
    int k = size(alpha_desc);
    real alpha_asc[k] = sort_asc(alpha_desc);
    real cutoff_pos[k];
    real cutoff_neg[k];
    real cutoffs[k * side]; 
    
    if (side == 1) {
      for (i in 1:k) {
        // calculate the one-sided cutoffs
        cutoffs[i] = inv_Phi(1 - alpha_desc[i]);
        // cutoffs[i] = inv_Phi(1 - alpha_desc[i]/2);
      }
    } else if (side == 2) {
      // for two-sided tests
      for (i in 1:k) {
        // calculate the two-sided cutoffs
        cutoff_pos[i] = inv_Phi(1 - alpha_desc[i]/2);
        cutoff_neg[i] = inv_Phi(alpha_asc[i]/2);
      }
      cutoffs = append_array(cutoff_neg, cutoff_pos);
    }
    
    return(cutoffs);
  } 
  
  real normal_lnorm(real theta, real tau, real sigma,
  real [] alpha, vector eta, int side, real [] cutoff) {
    int k = size(alpha);
    // real cutoff[k] = critical_value(alpha, 1);
    real cdf;
    real summands[k + 1];
    real se = sigma;
    real cutoff_pre;
    vector[k+1] omega_I = eta; //
    
    // summands[1] = eta[1];
    // 
    // for(i in 2:(k - 1)) {
    //   cutoff = inv_Phi(1 - alpha[i])*sigma;
    //   cdf = normal_cdf(cutoff, theta, sqrt(tau * tau + sigma * sigma));
    //   summands[i] = cdf*(eta[i] - eta[i - 1]);
    // }
    
    
    // calculate the denominator for one-sdied (positive) tests
    // similar to library(RoBMA)
    for (i in 1:k+1) {
          if (i == 1) {
            // cutoff = cutoffs[1] * se; 
            // cutoff = inv_Phi(1 - alpha[1])*se;
            summands[1] = normal_cdf(cutoff[1]*se, theta, sqrt(tau * tau + se * se)) * omega_I[1];
          } else if (i == k+1) {
            // cutoff_pre = inv_Phi(1 - alpha[k])*se; 
            summands[k+1] = (1-normal_cdf(cutoff[k]*se, theta, sqrt(tau * tau + se * se))) * omega_I[k+1];
          } else {
            // cutoff_pre = inv_Phi(1 - alpha[i-1])*se; 
            // cutoff = inv_Phi(1 - alpha[i])*se; 
            summands[i] = (normal_cdf(cutoff[i]*se, theta, sqrt(tau * tau + se * se)) -
            normal_cdf(cutoff[i-1]*se, theta, sqrt(tau * tau + se * se))) * omega_I[i];
          }
    }
    
    return(log(sum(summands)));
  }
  
  // Both the prior and likelihood make use of the same normalizing constant from
  // the likelihood, which is omitted in the 'mini' functions. 'Maxi' includes the
  // normalizing constant.
  
  // real psma_normal_prior_mini_lpdf(real theta, real theta0, real tau, real sigma,
  // real [] alpha, vector eta) {
  //   real y = normal_lpdf(theta | theta0, tau);
  //   // real normalizer = normal_lnorm(theta0, tau, sigma, alpha, eta);
  //   // return(y - normalizer);
  //   return(y);
  // }
  
  real psma_normal_mini_lpdf(real x, real theta, real theta0, real tau, real sigma,
  real [] alpha, vector eta, int side, real[] cutoff) {
    int k = size(alpha);
    real y = normal_lpdf(x | theta, sigma);
    real u = (1 - normal_cdf(x, 0, sigma));
    real alpha_desc[k] = alpha; // from large to small
    vector[k+1] omega = eta; //
    real normalizer = normal_lnorm(theta0, tau, sigma, alpha, eta, side, cutoff);
    // real normalizer = normal_lnorm(theta0, tau, sigma, alpha, eta);
    
    // for(i in 1:(k - 1)){
    //   if(alpha[i] < u && u <= alpha[i + 1]) {
    //     y += log(eta[i]);
    //     break;
    //   }
    // }
    
    // when alpha_asc is excluding 1 and 0
    if (u > alpha_desc[1]) {
      y += log(omega[1]);
    } else if(u <= alpha_desc[k]) {
      y += log(omega[k+1]);
    } else {
      for(i in 2:k){
        if(u > alpha_desc[i]  && u <= alpha_desc[i - 1]) {
          y += log(omega[i]);
          break;
        } 
      }
    }
    
    // return(y);
    return(y - normalizer);
  }
  
  // real psma_normal_maxi_lpdf(real x, real theta, real sigma,
  // real [] alpha, vector eta) {
  //   real y = psma_normal_mini_lpdf(x | theta, sigma, alpha, eta);
  //   real normalizer = normal_lnorm(theta, 0, sigma, alpha, eta);
  //   return(y - normalizer);
  // }
  
  // This is the marginal lpdf as in Hedges' paper.
  
  // real psma_normal_marginal_lpdf(real x, real theta0, real tau, real sigma,
  // real [] alpha, vector eta) {
  //   
  //   int k = size(alpha);
  //   real y = normal_lpdf(x | theta0, sqrt(tau * tau + sigma * sigma));
  //   real u = (1 - normal_cdf(x, 0, sigma));
  //   real normalizer = normal_lnorm(theta0, tau, sigma, alpha, eta);
  //   
  //   for(i in 1:(k - 1)){
  //     if(alpha[i] < u && u <= alpha[i + 1]) {
  //       y += log(eta[i]);
  //       break;
  //     }
  //   }
  //   
  //   return(y - normalizer);
  // }
}

data {
  
  // Input data.
  int<lower = 0> N;   // Number of observations.
  int<lower = 0> k;   // Length of alpha.
  real alpha[k];      // The vector of cuttoffs.
  real yi[N];         // The estimated effect sizes.
  real vi[N];         // The study-specific variances.
  
  // Prior parameters.
  vector[k + 1] eta0;
  // real tau_mean;
  // real <lower = 0> tau_sd;
  // real <lower = 0> u_min;
  // real <lower = 0> u_max;
  // real <lower = 0> shape;
  // real <lower = 0> scale;
  int tau_prior;
  
}

parameters {
  
  real theta0;
  real <lower = 0> tau;
  // positive_ordered[k - 1] weights;
  real theta[N];
  simplex[k+1] weights;
  
}

transformed parameters {
  vector[k + 1] eta;
  real cutoff[k] = critical_value(alpha, 1);
  // for(i in 1:(k - 1)) eta[i] = weights[k - i]/weights[k - 1];
  eta = cumulative_sum(weights);
  
}

model {
  
  theta0 ~ normal(0, 1);
  
  if(tau_prior == 1) {
    tau ~  normal(0, 3) T[0, ];
  } 
  // else if (tau_prior == 2) {
  //   tau ~ uniform(u_min, u_max);
  // } else if (tau_prior == 3) {
  //   tau ~ inv_gamma(shape, scale);
  // }
  
  // weights ~ gamma(eta0, 1);
  weights ~ dirichlet(eta0);
  
  for(n in 1:N) {
    theta[n] ~ normal(theta0, tau);
    // theta[n] ~ psma_normal_prior_mini(theta0, tau, sqrt(vi[n]), alpha, eta);
    yi[n] ~ psma_normal_mini(theta[n], theta0, tau, sqrt(vi[n]), alpha, eta, 1, cutoff);
  }
  
}

generated quantities {
  real cutoff_output[k] = cutoff;
  // vector[N] log_lik_marginal;
  // vector[N] log_lik;
  
  // for(n in 1:N)
  // log_lik[n] = psma_normal_maxi_lpdf(yi[n] | theta[n], sqrt(vi[n]), alpha, eta);
  // 
  // for(n in 1:N)
  // log_lik_marginal[n] = psma_normal_marginal_lpdf(yi[n] | theta0, tau, sqrt(vi[n]), alpha, eta);
  
} 
