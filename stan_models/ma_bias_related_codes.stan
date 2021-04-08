// generated with brms 2.15.0
functions {
  real[] critical_value(real [] alpha_d, int side) {
    int k = size(alpha_d);
    real alpha_a[k] = sort_asc(alpha_d);
    real cutoff_pos[k];
    real cutoff_neg[k];
    real cutoffs[k * side]; 
    
    if (side == 1) {
      for (i in 1:k) {
        // calculate the one-sided cutoffs
        cutoffs[i] = inv_Phi(1 - alpha_d[i]);
        // cutoffs[i] = inv_Phi(1 - alpha_d[i]/2);  // don't think should use this
      }
    } else if (side == 2) {
      // for two-sided tests
      for (i in 1:k) {
        // calculate the two-sided cutoffs
        cutoff_pos[i] = inv_Phi(1 - alpha_d[i]/2);
        cutoff_neg[i] = inv_Phi(alpha_a[i]/2);
      }
      cutoffs = append_array(cutoff_neg, cutoff_pos);
    }
    
    return(cutoffs);
  } 
  
  // Note: alpha is from largest (1) to smallest (0)
  real weighted_normal_lnorm(real mu, real tau, real se, real [] alpha_d, 
  vector omega, int side, real [] cutoffs) {
    // alpha is the critical alpha level (excluding 1 and 0) from large to small
    int k = size(alpha_d);  // length of alpha
    int I = k*side+1;  // I is the number of intervals
    // real cutoffs[I-1] = critical_value(alpha, side);
    real cutoff;
    real cutoff_pre;
    vector[k] omega_neg;
    vector[I] omega_I;
    real sumdenom[I];
    
    if (side == 2) {
      omega_neg = sort_desc(omega)[1:k];
      omega_I = append_row(omega_neg, omega);
    } else if (side == 1) {
      omega_I = omega;
    }
    
    // calculate the denominator for one-sdied (positive) tests
    // similar to library(RoBMA)
    for (i in 1:I) {
      if (i == 1) {
        cutoff = cutoffs[1] * se; 
        sumdenom[1] = normal_cdf(cutoff, mu, sqrt(tau * tau + se * se)) * omega_I[1];
      } else if (i == I) {
        cutoff_pre = cutoffs[I-1] * se; 
        sumdenom[I] = (1 - normal_cdf(cutoff_pre, mu, sqrt(tau * tau + se * se))) * omega_I[I];
      } else {
        cutoff_pre = cutoffs[i-1] * se; 
        cutoff = cutoffs[i] * se; 
        sumdenom[i] = (normal_cdf(cutoff, mu, sqrt(tau * tau + se * se)) -
        normal_cdf(cutoff_pre, mu, sqrt(tau * tau + se * se))) * omega_I[i];
      }
    }
    
    return(log(sum(sumdenom)));
  }
  
  real weighted_normal_lpdf(real x, real mu, real Intercept, real tau, 
  real se, real [] alpha_d, vector omega, int side, real [] cutoffs) {
    int k = size(alpha_d); // length of alpha
    real y = normal_lpdf(x | mu, se);
    real u;
    // there are two ways to calculate the normalizer; I prefer the second approach (at least now)
    // they are supposed to be equivalent?
    real normalizer = weighted_normal_lnorm(Intercept, tau, se, alpha_d, omega, side, cutoffs);
    // real normalizer = weighted_normal_lnorm(mu, 0, se, alpha_d, omega, side, cutoffs);
    
    if (side == 1) {
      u = (1 - normal_cdf(x, 0, se));
    } else if (side == 2) {
      u = (1 - normal_cdf(fabs(x), 0, se))*2;
    }
    
    // when alpha is excluding 1 and 0
    if (u > alpha_d[1]) {
      y += log(omega[1]);
    } else if(u <= alpha_d[k]) {
      y += log(omega[k+1]);
    } else {
      for(i in 2:k){
        if(u > alpha_d[i]  && u <= alpha_d[i - 1]) {
          y += log(omega[i]);
          break;
        } 
      }
    }
    
    return(y - normalizer);
  }
}
data {
  //...
  // (bias-related) 
  int<lower=1> N_alpha; // length of omega
  real<lower=0,upper=1> alpha[N_alpha]; // start from 1 and finish at 0
  int<lower=1, upper=2> side; // one-side or two-side tests
}
transformed data {
  //...
  real<lower=0,upper=1> alpha_d[N_alpha] = sort_desc(alpha); // order alpha from large to small
  real cutoffs[N_alpha*side] = critical_value(alpha_d, side); // critical cutoff values
}
parameters {
  // ...
  // (bias-related)
  simplex[N_alpha+1] eta; // 
}
transformed parameters {
  //...
  vector<lower=0>[N_alpha+1] omega;  // (bias-related) publication bias!!!
  //...
  // (bias-related) calculate omega based on eta from dirichlet
  omega = cumulative_sum(eta);
}
model {
  // likelihood including constants
  
  // initialize linear predictor term
  //...
  // (bias-related) For each study separately
  target += weighted_normal_lpdf(Y[n] | mu[n], Intercept, sd_1[1], se[n], 
      alpha_d, omega, side, cutoffs);
  // target += normal_lpdf(Y | mu, se);
  
  // ...
  // (bias-related) 
  target += dirichlet_lpdf(eta | rep_vector(1, N_alpha+1));
}
generated quantities {
  // ...
}
