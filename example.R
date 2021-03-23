# Example 1
library(brms)
library(rstan)
rstan_options(auto_write = TRUE)

df1 <- readr::read_csv("data/ExampleData.csv")

# fit with brms
brmfit1 <- brm(d | se(se) ~ 1 + (1|name), data = df1,
              chains = 4, cores = 4, seed = 12)

# stancode(brmfit1)  # obtain the stancode
data_ls1 <- standata(brmfit1) # obtain the stan data

data_ls1$alpha <- c(1, 0.025, 0.05, 0)
data_ls1$I <- length(data_ls1$alpha)-1  # number of intervals
ex1_bias <- stan(file = 'stan_models/ma_bias.stan',
                 data = data_ls1,
                 chains = 4, cores = 4, seed = 12)
summary(ex1_bias, pars=c("omega", "b_Intercept"))
pairs(ex1_bias, pars=c("omega", "b_Intercept"))

###############################################################################
# Example 2
df2 <- metafor::dat.begg1989
brmfit2 <- brm(yi | se(sqrt(vi)) ~ 1 + (1|study), data = df2,
              chains = 4, cores = 4, seed = 12)

data_ls2 <- standata(brmfit2)
data_ls2$alpha <- c(1, 0.025, 0.05, 0)
data_ls2$I <- length(data_ls2$alpha)-1 # number of intervals
ex2_bias <- stan(file = 'stan_models/ma_bias.stan', #  'stan_bias.stan',
                 data = data_ls2,
                 chains = 4, cores = 4, seed = 1)
summary(ex2_bias, pars=c("omega", "b_Intercept"))
pairs(ex2_bias, pars=c("omega", "b_Intercept"))

