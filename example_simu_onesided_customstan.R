# Example X: one-sided selection models (mainly to compare all the custom and modified stan codes)
## Simulate data
library(tidyverse)
library(brms)
library(rstan)
library(MCMCvis)
library(publipha)
library(RoBMA)

rstan_options(auto_write = TRUE)

set.seed(2021)
N_exp <- 60            # number of studies
theta_exp <- 0.3       # the true effect

tau_e_exp <- 0.15      # the standard deviation of the true effects for studies
SE_exp <- 0.2          # the standard deviation of the sampling distribution

# the probablity to publish the study based the p-value
w_publish <- c(0.2, 0.5, 1)  # (0.10, 1] (0.05, 0.10] [0, 0.05]

alpha_e <- rnorm(N_exp, 0, tau_e_exp)      # studies

df_simu <- tibble(
    SE = extraDistr::rtnorm(N_exp, 0.1, SE_exp, a = 0)) %>%
    mutate(g = rnorm(n(), theta_exp+alpha_e, SE), # 
           experiment = paste0("E",1:n()),
           pvalue = round((1-pnorm(g, 0, SE))*2, 3), # still two-sided p-value
           # pvalue = round(1-pnorm(g, 0, SE), 3), # one-sdied p-value (should not use this)
           # p-value intervals
           K = case_when(pvalue <= .05 ~ 3,
                         (pvalue > .05 & pvalue <= .1) ~ 2,
                         pvalue > .1 ~ 1),
           published = extraDistr::rbern(n(), w_publish[K]))

dfx <- filter(df_simu, published==1)

df_simu %>%
    group_by(K) %>%
    summarize(mean(published), n())

dfx %>% 
    summarize(mean(g))

##################### Fit with brms #####################
### Fit the brm model (without publication bias)
# Note: brms still fit two-sided models
bfx <- bf(g | se(SE) ~ 1 + (1|experiment))
brmfit_one <- brm(bfx, data = dfx,
               chains = 6, cores = 6, seed = 12,
               control = list(adapt_delta = .9))


##################### Equivalent codes Set 1 ##################################
# 1.28 and 1.65
# Below functions essentially use one-sided p-value of 0.10 and 0.05 to set 
# intervals for selection models. 

## with custom Stan codes.
data_ls_one1 <- make_standata(bfx, data=dfx)
data_ls_one1$alpha <- c(0.10, 0.05)
data_ls_one1$N_alpha <- length(data_ls_one1$alpha)  # number of intervals
data_ls_one1$side <- 1
ex4_bias_one1 <- stan(file = 'stan_models/ma_bias_twosided.stan', #  'stan_bias.stan',
                      data = data_ls_one1,
                      iter = 6000, warmup = 2000, 
                      chains = 6, cores = 6, seed = 12,
                      control = list(adapt_delta = .9))
MCMCsummary(ex4_bias_one1, params=c("b_Intercept", "omega", "cutoff_output"))
# pairs(ex4_bias_one1, pars=c("omega", "b_Intercept"))

ex4_bias_one1_gamma_cumsum <- stan(file = 'stan_models/ma_bias_twosided_gamma_cumsum.stan', #  'stan_bias.stan',
                            data = data_ls_one1,
                            iter = 6000, warmup = 2000, 
                            chains = 6, cores = 6, seed = 12,
                            control = list(adapt_delta = .9))
MCMCsummary(ex4_bias_one1_gamma_cumsum, params=c("b_Intercept", "omega", "cutoff_output"))
# pairs(ex4_bias_one1_gamma_cumsum, pars=c("omega", "b_Intercept"))

ex4_bias_one1_gamma <- stan(file = 'stan_models/ma_bias_twosided_gamma.stan', #  'stan_bias.stan',
                            data = data_ls_one1,
                            iter = 6000, warmup = 2000, 
                            chains = 6, cores = 6, seed = 12,
                            control = list(adapt_delta = .9))
MCMCsummary(ex4_bias_one1_gamma, params=c("b_Intercept", "omega", "cutoff_output"))
# pairs(ex4_bias_one1_gamma, pars=c("omega", "b_Intercept"))


##### other packages

## With library(publipha)
psmafit1 <- psma(yi = dfx$g, vi = dfx$SE^2, alpha = c(0, 0.05, 0.1, 1),
                 iter = 4000, warmup = 2000, 
                 chains = 6, cores = 6, seed = 12,
                 control = list(adapt_delta = .9))
MCMCsummary(psmafit1, params=c("theta0", "eta"))


## with modified Stan codes from library(publipha)
psmafit_ls1 <- list(
    N = data_ls_one1$N,
    k = 2,
    alpha = c(0.1, 0.05),
    yi = data_ls_one1$Y,
    vi = data_ls_one1$se^2,
    eta0 = c(1,1,1),
    tau_prior = 1
)
# the following results should be the same as psmafit
psmafit1_stan_gamma <- stan(file = 'stan_models/psma_stan_gamma_prior.stan', #  'stan_bias.stan',
                      data = psmafit_ls1,
                      iter = 4000, warmup = 2000, 
                      chains = 6, cores = 6, seed = 12,
                      control = list(adapt_delta = .9))
MCMCsummary(psmafit1_stan_gamma, params=c("theta0", "eta", "cutoff_output"))

psmafit1_stan <- stan(file = 'stan_models/psma_stan.stan', #  'stan_bias.stan',
                      data = psmafit_ls1,
                      iter = 4000, warmup = 2000, 
                      chains = 6, cores = 6, seed = 12,
                      control = list(adapt_delta = .9))
MCMCsummary(psmafit1_stan, params=c("theta0", "eta", "cutoff_output"))
# Note: it seems that using gamma or dirichlet make differences. 

## With library(RoMBA)
robma_1 <- RoBMA(d = dfx$g, se = dfx$SE,
                   priors_omega = list(prior(distribution = "one.sided", 
                                             parameters = list(alpha = c(1, 1, 1), 
                                                               steps = c(.05, .10)), 
                                             prior_odds = 1)),
                   priors_mu_null = NULL,
                   priors_tau_null = NULL,
                   priors_omega_null = NULL)
robma_1$models


##################### Equivalent codes Set 2 ##################################
# 1.28, 1.65 and 1.96
# Below functions essentially use one-sided p-value of 0.1, 0.05 and 0.025 to set 
# intervals for selection models. This is equivalent to two-sided tests but only 
# the positive effects are easier to be published.

## Fit the model with the publication bias (with one-sided)
data_ls_one2 <- make_standata(bfx, data=dfx)
data_ls_one2$alpha <- c(0.1, 0.05, .025)
data_ls_one2$N_alpha <- length(data_ls_one2$alpha)  # number of intervals
data_ls_one2$side <- 1
ex4_bias_one2 <- stan(file = 'stan_models/ma_bias_twosided.stan', #  'stan_bias.stan',
                     data = data_ls_one2,
                     chains = 8, cores = 8, seed = 12,
                     iter = 4000, warmup = 2000, 
                     control = list(adapt_delta = .9))
MCMCsummary(ex4_bias_one2, params=c("b_Intercept", "omega", "cutoff_output"))
# pairs(ex4_bias_one2, pars=c("omega", "b_Intercept"))

ex4_bias_one2_gamma_cumsum <- stan(file = 'stan_models/ma_bias_twosided_gamma_cumsum.stan', #  'stan_bias.stan',
                                   data = data_ls_one2,
                                   iter = 4000, warmup = 2000, 
                                   chains = 6, cores = 6, seed = 12,
                                   control = list(adapt_delta = .9))
MCMCsummary(ex4_bias_one2_gamma_cumsum, params=c("b_Intercept", "omega", "cutoff_output"))
# pairs(ex4_bias_one2_gamma_cumsum, pars=c("omega", "b_Intercept"))

ex4_bias_one2_gamma <- stan(file = 'stan_models/ma_bias_twosided_gamma.stan', #  'stan_bias.stan',
                            data = data_ls_one2,
                            iter = 4000, warmup = 2000, 
                            chains = 6, cores = 6, seed = 12,
                            control = list(adapt_delta = .9))
MCMCsummary(ex4_bias_one2_gamma, params=c("b_Intercept", "omega", "cutoff_output"))
# pairs(ex4_bias_one2_gamma, pars=c("omega", "b_Intercept"))


## with library(publipha)
psmafit2 <- psma(yi = dfx$g, vi = dfx$SE^2, alpha = c(0, 0.025, 0.05, 1),
                 chains = 6, cores = 6, seed = 12,
                 iter = 4000, warmup = 2000, 
                 control = list(adapt_delta = .9))
MCMCsummary(psmafit2, params=c("theta0", "eta"))


## with Stan code modified from library(publipha)
psmafit_ls2 <- list(
    N = data_ls_one2$N,
    k = 3,
    alpha = c(0.1, 0.05, 0.025),
    yi = data_ls_one2$Y,
    vi = data_ls_one2$se^2,
    eta0 = c(1,1,1,1),
    tau_prior = 1
)
# the following results should be the same as psmafit
psmafit2_stan_gamma <- stan(file = 'stan_models/psma_stan_gamma_prior.stan', #  'stan_bias.stan',
                            data = psmafit_ls2,
                            iter = 4000, warmup = 2000, 
                            chains = 6, cores = 6, seed = 12,
                            control = list(adapt_delta = .9))
MCMCsummary(psmafit2_stan_gamma, params=c("theta0", "eta", "cutoff_output"))
# Note: it seems that using gamma or dirichlet make differences. 

psmafit2_stan <- stan(file = 'stan_models/psma_stan.stan', #  'stan_bias.stan',
                      data = psmafit_ls2,
                      iter = 4000, warmup = 2000, 
                      chains = 6, cores = 6, seed = 12,
                      control = list(adapt_delta = .9))
MCMCsummary(psmafit2_stan, params=c("theta0", "eta", "cutoff_output"))



## With library(RoMBA)
robma_2 <- RoBMA(d = dfx$g, se = dfx$SE,
                 priors_omega = list(prior(distribution = "one.sided", 
                                           parameters = list(alpha = c(1, 1, 1, 1), 
                                                             steps = c(.025, .05, .1)), 
                                           prior_odds = 1)),
                 priors_mu_null = NULL,
                 priors_tau_null = NULL,
                 priors_omega_null = NULL)
robma_2$models

