// generated with brms 2.15.0
functions {
  // // function from library(publipha)
  // // Note: alpha is from smallest (0) to largest (1)
  // real normal_lnorm(real theta, real tau, real sigma,
  // real [] alpha, vector eta) {
  //   int k = size(alpha);
  //   real cutoff;
  //   real cdf;
  //   real summands[k - 1];
  // 
  //   summands[1] = eta[1];
  // 
  //   for(i in 2:(k - 1)) {
  //     cutoff = inv_Phi(1 - alpha[i])*sigma;
  //     cdf = normal_cdf(cutoff, theta, sqrt(tau * tau + sigma * sigma));
  //     summands[i] = cdf*(eta[i] - eta[i - 1]);
  //   }
  // 
  //   return(log(sum(summands)));
  // }

  // Note: alpha is from largest (1) to smallest (0)
  real bias_normal_lnorm(real mu, real tau, real se, real [] alpha, vector omega) {
    int k = size(alpha); // alpha is the critical alpha level (including 1 and 0) from large to small
    real cutoff;
    real cutoff_pre;
    real cdf;
    real sumdenom[k - 1];

    // // method from library(publipha) (publipha/src/stan_files/chunks/psma_likelihoods.stan)
    // for (i in 2:k-1) {
    //   cutoff = inv_Phi(1 - alpha[i]) * se;
    //   cdf = normal_cdf(cutoff, mu, sqrt(tau * tau + se * se));
    //   sumdenom[i-1] = cdf*(omega[i-1] - omega[i]);
    // }
    // sumdenom[k - 1] = omega[k - 1]; // always to be 1

    // another method (similar to library(RoBMA))
    for (i in 2:k) {
      if (i == 2) {
        cutoff = inv_Phi(1 - alpha[i]) * se;
        sumdenom[i-1] = normal_cdf(cutoff, mu, sqrt(tau * tau + se * se)) * omega[i-1];
      } else if (i == k) {
        cutoff_pre = inv_Phi(1 - alpha[i-1]) * se;
        sumdenom[i-1] = (1 - normal_cdf(cutoff_pre, mu, sqrt(tau * tau + se * se))) * omega[i-1];
      } else {
        cutoff_pre = inv_Phi(1 - alpha[i-1]) * se;
        cutoff = inv_Phi(1 - alpha[i]) * se;
        sumdenom[i-1] = (normal_cdf(cutoff, mu, sqrt(tau * tau + se * se)) -
        normal_cdf(cutoff_pre, mu, sqrt(tau * tau + se * se))) * omega[i-1];
      }
    }

    return(log(sum(sumdenom)));
  }
  
  // real bias_normal_prior_mini_lpdf(real mu, real Intercept, real tau, real se,
  // real [] alpha, vector omega) {
  //   real y = normal_lpdf(mu | Intercept, tau);
  //   real normalizer = normal_lnorm(Intercept, tau, se, alpha, omega);
  //   return(y - normalizer);
  // }
  
  real bias_normal_mini_lpdf(real x, real mu, real Intercept, real tau, 
  real se, real [] alpha, vector omega) {
    int k = size(alpha);
    real y = normal_lpdf(x | mu, se);
    real u = (1 - normal_cdf(x, 0, se));
    // real normalizer = normal_lnorm(Intercept, tau, se, alpha, omega);
    // real normalizer = normal_lnorm(mu, 0, se, alpha, omega);
    // real normalizer = bias_normal_lnorm(Intercept, tau, se, alpha, omega);
    real normalizer = bias_normal_lnorm(mu, 0, se, alpha, omega);

    // when alpha is from 1 to 0
    for(i in 1:(k - 1)){
      if(alpha[i] >= u && u > alpha[i + 1]) {
        y += log(omega[i]);
        break;
      }
    }
    
    // when alpha is from 0 to 1
    // for(i in 1:(k - 1)){
    // if(alpha[i] < u && u <= alpha[i + 1]) {
    //   y += log(omega[i]);
    //   break;
    // }
  // }
    
    return(y - normalizer);
  }
  
  // real bias_normal_maxi_lpdf(real x, real mu, real se,
  // real [] alpha, vector omega) {
  //   real y = bias_normal_mini_lpdf(x | mu, se, alpha, omega);
  //   real normalizer = bias_normal_lnorm(mu, 0, se, alpha, omega);
  //   return(y - normalizer);
  // }
  
  // // This is the marginal lpdf as in Hedges' paper.
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
  //   return(y - normalizer);
  // }
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
  // int<lower=1> I[N]; // index for intervals based on p-value
  real<lower=0,upper=1> alpha[3];
}
transformed data {
  vector<lower=0>[N] se2 = square(se);
}
parameters {
  real Intercept;  // temporary intercept for centered predictors
  vector<lower=0>[M_1] sd_1;  // group-level standard deviations
  vector[N_1] z_1[M_1];  // standardized group-level effects
  // (bias-related)
  simplex[2] theta; // 
}
transformed parameters {
  real<lower=0> sigma = 0;  // residual SD
  vector<lower=0,upper=1>[2] omega;  // (bias-related) publication bias
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
      target += bias_normal_mini_lpdf(Y[n] | mu[n], Intercept, sd_1[1], se[n], alpha, omega);

      // target += normal_lpdf(Y[n] | mu[n], se[n]);
      // target += log(omega[I[n]]);
      // // denominators
      // target += - log_sum_exp(
      //   normal_lcdf(1.96 | mu[n], sqrt(se[n] * se[n] + sd_1[1] * sd_1[1])) + log(omega[1]),
      //   normal_lccdf(1.96 | mu[n], sqrt(se[n] * se[n] + sd_1[1] * sd_1[1])) + log(omega[2])
      //   );
    }
    // target += normal_lpdf(Y | mu, se);
  }
  // priors including constants
  target += student_t_lpdf(Intercept | 3, 0.5, 2.5);
  target += student_t_lpdf(sd_1 | 3, 0, 2.5)
    - 1 * student_t_lccdf(0 | 3, 0, 2.5);
  target += std_normal_lpdf(z_1[1]);
  // (bias-related) 
  target += dirichlet_lpdf(theta | rep_vector(2, 2));
}
generated quantities {
  // actual population-level intercept
  real b_Intercept = Intercept;
}
