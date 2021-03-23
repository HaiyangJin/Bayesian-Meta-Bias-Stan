# Example 1
library(brms)
library(rstan)

df1 <- readr::read_csv("data/ExampleData.csv")

# fit with brms
brmfit1 <- brm(d | se(se) ~ 1 + (1|name), data = df1,
              chains = 4, cores = 4, seed = 12)

# stancode(brmfit1)  # obtain the stancode
data_ls1 <- standata(brmfit1) # obtain the stan data

data_ls1$alpha <- c(0, 0.025, 1)
ex1_bias <- stan(file = 'stan_models/ma_bias.stan',
                 data = data_ls1,
                 chains = 4, cores = 4, seed = 12)
summary(ex1_bias, pars=c("omega", "b_Intercept"))
pairs(ex1_bias, pars=c("omega", "b_Intercept"))

# stancode(brmfit)
###############################################################################
# Example 2
df2 <- metafor::dat.begg1989
brmfit2 <- brm(yi | se(sqrt(vi)) ~ 1 + (1|study), data = df2,
              chains = 4, cores = 4, seed = 12)

data_ls2 <- standata(brmfit2)
data_ls2$alpha <- c(0, 0.025, 1)
ex2_bias <- stan(file = 'stan_models/ma_bias.stan', #  'stan_bias.stan',
                 data = data_ls2,
                 chains = 4, cores = 4, seed = 1)
summary(ex2_bias, pars=c("omega", "b_Intercept"))
pairs(ex2_bias, pars=c("omega", "b_Intercept"))

stanfit_bias <- stan(file = 'stan_models/ma_bias.stan', #  'stan_bias.stan',
                     data = data_ls,
                     chains = 4, cores = 4)

summary(stanfit_bias, pars="omega")
