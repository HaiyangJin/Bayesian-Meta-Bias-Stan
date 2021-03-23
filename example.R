library(tidyverse)
library(brms)
library(rstan)

df <- read_csv("data/ExampleData.csv")

brmfit <- brm(d | se(se) ~ 1 + (1|name), data = df,
              chains = 4, cores = 4)

# stancode(brmfit)

data_ls <- standata(brmfit)
data_ls$I <- df$I

stanfit_bias <- stan(file = 'stan_models/ma_bias.stan', #  'stan_bias.stan',
                     data = data_ls,
                     chains = 4, cores = 4)

summary(stanfit_bias, pars="omega")
