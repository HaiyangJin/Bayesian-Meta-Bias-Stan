# Example 1
library(brms)
library(rstan)
library(MCMCvis)
rstan_options(auto_write = TRUE)

# this dataset is from library(RoBMA)
df1 <- readr::read_csv("data/ExampleData.csv")

# meta-analysis with selection models

bf1 <- brmsformula(d | se(se) ~ 1 + (1|name))
# make_stancode(bf1, data=df1)
data_ls1 <- make_standata(bf1, data=df1) # obtain the stan data
data_ls1$alpha <- c(0.10, 0.05)
data_ls1$N_alpha <- length(data_ls1$alpha)  # number of intervals
data_ls1$side <- 2 # two-sided tests
ex1_bias <- stan(file = 'stan_models/ma_bias_twosided.stan',
                 data = data_ls1,
                 chains = 4, cores = 4, seed = 12)
MCMCsummary(ex1_bias, params=c("b_Intercept", "omega"))
# pairs(ex1_bias, pars=c("omega", "b_Intercept", "eta"))


###############################################################################
# Example 2
df2 <- metafor::dat.begg1989
bf2 <- bf(yi | se(sqrt(vi)) ~ 1 + (1|study))

# make_stancode(bf2, data = df2)
data_ls2 <- make_standata(bf2, data = df2)
data_ls2$alpha <- c(0.10, 0.05)
data_ls2$N_alpha <- length(data_ls2$alpha) # number of intervals
data_ls2$side <- 2 # two-sided tests
ex2_bias <- stan(file = 'stan_models/ma_bias_twosided.stan', #  'stan_bias.stan',
                 data = data_ls2,
                 chains = 4, cores = 4, seed = 1)
MCMCsummary(ex2_bias, params=c("b_Intercept", "omega"))
pairs(ex2_bias, pars=c("b_Intercept", "omega"))

