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
ex1_bias <- stan(file = '../stan_models/ma_bias_selections.stan',
                 data = data_ls1,
                 chains = 4, cores = 4, seed = 12)
summary(ex1_bias, pars=c("omega", "b_Intercept"))
pairs(ex1_bias, pars=c("omega", "b_Intercept"))

data_ls11 <- data_ls1
data_ls11$alpha <- c(1, 0.025, 0)
ex1_bias1 <- stan(file = '../stan_models/ma_bias.stan', #  'stan_bias.stan',
                     data = data_ls11,
                     chains = 4, cores = 4, seed = 12)
summary(ex1_bias1, pars=c("omega", "b_Intercept"))
pairs(ex1_bias1, pars=c("omega", "b_Intercept"))

###############################################################################
# Example 2
df2 <- metafor::dat.begg1989
brmfit2 <- brm(yi | se(sqrt(vi)) ~ 1 + (1|study), data = df2,
              chains = 4, cores = 4, seed = 12)

data_ls2 <- standata(brmfit2)
data_ls2$alpha <- c(0, 0.025, 1)
ex2_bias <- stan(file = '../stan_models/ma_bias_selections.stan', #  'stan_bias.stan',
                 data = data_ls2,
                 chains = 4, cores = 4, seed = 1)
summary(ex2_bias, pars=c("omega", "b_Intercept"))
pairs(ex2_bias, pars=c("omega", "b_Intercept"))


data_ls22 <- data_ls2
data_ls22$alpha <- c(1, 0.025, 0)
ex2_bias22 <- stan(file = '../stan_models/ma_bias.stan', #  'stan_bias.stan',
                      data = data_ls22,
                      chains = 4, cores = 4, seed = 1)
summary(ex2_bias22, pars=c("omega", "b_Intercept"))
pairs(ex2_bias22, pars=c("omega", "b_Intercept"))


