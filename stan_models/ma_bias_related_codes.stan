// generated with brms 2.15.0
functions {
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

  real bias_normal_mini_lpdf(real x, real mu, real Intercept, real tau, 
  real se, real [] alpha, vector omega) {
    int k = size(alpha);
    real y = normal_lpdf(x | mu, se);
    real u = (1 - normal_cdf(x, 0, se));
    // there are two ways to calculate the normalizer; I prefer the second approach (at least now)
    // they are supposed to be equivalent?
    real normalizer = bias_normal_lnorm(Intercept, tau, se, alpha, omega);
    // real normalizer = bias_normal_lnorm(mu, 0, se, alpha, omega);

    // when alpha is from 1 to 0
    for(i in 1:(k - 1)){
      if(alpha[i] >= u && u > alpha[i + 1]) {
        y += log(omega[i]);
        break;
      }
    }
    return(y - normalizer);
  }
}
data {
  //...
  // (bias-related) 
  int<lower=1> I; // number of critial alpha (not include 0 or 1)
  real<lower=0,upper=1> alpha[I+1]; // start from 1 and finish at 0
}
transformed data {
  //...
}
parameters {
  // ...
  // (bias-related)
  simplex[I] theta; // 
}
transformed parameters {
  //...
  vector<lower=0,upper=1>[I] omega;  // (bias-related) publication bias
  //...
  // (bias-related) calculate omega based on theta from dirichlet
  omega = cumulative_sum(theta);
}
model {
  // likelihood including constants
  
  // initialize linear predictor term
  //...
  // (bias-related) For each study separately
  target += bias_normal_mini_lpdf(Y[n] | mu[n], Intercept, sd_1[1], se[n], alpha, omega);
  // target += normal_lpdf(Y | mu, se);
  
  // ...
  // (bias-related) 
  target += dirichlet_lpdf(theta | rep_vector(2, I));
}
generated quantities {
  // ...
}