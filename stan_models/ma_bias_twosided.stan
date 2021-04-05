// generated with brms 2.15.0
functions {
  
  real[] critical_value(real [] alpha, int side) {
    int k = size(alpha);
    real alpha_asc[k] = sort_asc(alpha);
    real alpha_desc[k] = sort_desc(alpha);
    real cutoff_pos[k];
    real cutoff_neg[k];
    real cutoffs[k * side]; 
    
    for (i in 1:k) {
      cutoff_pos[i] = inv_Phi(1 - alpha_desc[i]/2);
      cutoff_neg[k+1-i] = inv_Phi(alpha_asc[i]/2);
    }
    
    if (side == 1) {
      cutoffs = cutoff_pos;
    } else if (side == 2) {
      cutoffs = append_array(cutoff_neg, cutoff_pos);
    }
    
    return(cutoffs);
  } 
  
  // Note: alpha is from largest (1) to smallest (0)
  real bias_normal_lnorm(real mu, real tau, real se, real [] alpha, 
  vector omega, int side, real [] cutoffs) {
    // alpha is the critical alpha level (excluding 1 and 0) from large to small
    int k = size(alpha);  // length of alpha
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
  
  real bias_normal_mini_lpdf(real x, real mu, real Intercept, real tau, 
  real se, real [] alpha, vector omega, int side, real [] cutoffs) {
    int k = size(alpha)+1; // length of omega
    real y = normal_lpdf(x | mu, se);
    real u;
    // there are two ways to calculate the normalizer; I prefer the second approach (at least now)
    // they are supposed to be equivalent?
    // real normalizer = bias_normal_lnorm(Intercept, tau, se, alpha, omega);
    real normalizer = bias_normal_lnorm(mu, 0, se, alpha, omega, side, cutoffs);
    
    if (side == 1) {
      u = (1 - normal_cdf(x, 0, se));
    } else if (side == 2) {
      u = (1 - normal_cdf(fabs(x), 0, se))*2;
    }
    
    // when alpha is excluding 1 and 0
    if (u > alpha[1]) {
      y += log(omega[1]);
    } else if(u <= alpha[k-1]) {
      y += log(omega[k]);
    } else {
      for(i in 2:k-1){
        if(u > alpha[i]  && u <= alpha[i - 1]) {
          y += log(omega[i]);
          break;
        } 
      }
    }
    
    return(y - normalizer);
  }
}
data {
  int<lower=1> N;  // total number of observations
  vector[N] Y;  // response variable
  vector<lower=0>[N] se;  // known sampling error
  // data for group-level effects of ID 1
  int<lower=1> N_1;  // number of grouping levels
  int<lower=1> M_1;  // number of coefficients per level
  int<lower=1> J_1[N];  // grouping indicator per observation
  // group-level predictor values
  vector[N] Z_1_1;
  int prior_only;  // should the likelihood be ignored?
  // (bias-related) 
  int<lower=1> I; // length of omega
  real<lower=0,upper=1> alpha[I-1]; // start from 1 and finish at 0
  int<lower=1, upper=2> S; // one-side or two-side tests
}
transformed data {
  vector<lower=0>[N] se2 = square(se);
  real cutoffs[(I-1)*S] = critical_value(alpha, S); // critical cutoff values
}
parameters {
  real Intercept;  // temporary intercept for centered predictors
  vector<lower=0>[M_1] sd_1;  // group-level standard deviations
  vector[N_1] z_1[M_1];  // standardized group-level effects
  // (bias-related)
  simplex[I] theta; // 
}
transformed parameters {
  real<lower=0> sigma = 0;  // residual SD
  vector<lower=0,upper=1>[I] omega;  // (bias-related) publication bias
  vector[N_1] r_1_1;  // actual group-level effects
  r_1_1 = (sd_1[1] * (z_1[1]));
  // (bias-related) calculate omega based on theta from dirichlet
  omega = cumulative_sum(theta);
}
model {
  // likelihood including constants
  if (!prior_only) {
    // initialize linear predictor term
    vector[N] mu = Intercept + rep_vector(0.0, N);
    for (n in 1:N) {
      // add more terms to the linear predictor
      mu[n] += r_1_1[J_1[n]] * Z_1_1[n];
      // (bias-related)
      target += bias_normal_mini_lpdf(Y[n] | mu[n], Intercept, sd_1[1], se[n], 
      alpha, omega, S, cutoffs);
    }
    // target += normal_lpdf(Y | mu, se);
  }
  // priors including constants
  target += student_t_lpdf(Intercept | 3, 0.5, 2.5);
  target += student_t_lpdf(sd_1 | 3, 0, 2.5)
    - 1 * student_t_lccdf(0 | 3, 0, 2.5);
  target += std_normal_lpdf(z_1[1]);
  // (bias-related) 
  target += dirichlet_lpdf(theta | rep_vector(2, I));
}
generated quantities {
  // actual population-level intercept
  real b_Intercept = Intercept;
}

