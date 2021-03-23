# Example 3
## Simulate data
library(tidyverse)

set.seed(12)
N_pw <- 50           # number of studies for PW
N_pw_a <- 38         # number of articles
theta_pw <- 0.3      # the true effect for PW (small effects) 

tau_e_pw <- 0.15     # the standard deviation of the true effects for studies
tau_a_pw <- 0.1      # the standard deviation of the true effects for articles
SE_pw <- 0.2         # the standard deviation of the sampling distribution

# the probablity to publish the study based the p-value
w_publish <- c(0.2, 0.5, 1)  # (0.10, 1] (0.05, 0.10] [0, 0.05]

alpha_e <- rnorm(N_pw, 0, tau_e_pw)      # studies
alpha_a <- rnorm(N_pw_a, 0, tau_a_pw)    # articles

df_simu <- tibble(
    article = c(sample(1:N_pw_a, N_pw_a), 
                sample(1:N_pw_a, N_pw-N_pw_a, replace = T)),
    SE = extraDistr::rtnorm(N_pw, 0.1, SE_pw, a = 0)) %>%
    arrange(article) %>% 
    mutate(g = rnorm(n(), theta_pw+alpha_e+alpha_a[article], SE), # 
           experiment = paste0("E",1:n()),
           article = paste0("A",article),
           pvalue = round((1-pnorm(abs(g), 0, SE))*2, 3),
           # p-value intervals
           K = case_when(pvalue <= .05 ~ 3,
                         (pvalue > 0.05 & pvalue <= .1) ~ 2,
                         pvalue > .1 ~ 1),
           published = extraDistr::rbern(n(), w_publish[K]))

df3 <- filter(df_simu, published==1)

## Fit the model
library(brms)
library(rstan)
rstan_options(auto_write = TRUE)

### Fit the brm model (without publication bias)
brmfit3 <- brm(g | se(SE) ~ 1 + (1|experiment), data = df3,
               chains = 4, cores = 4, seed = 12)
data_ls3 <- standata(brmfit3)
data_ls3$alpha <- c(1, 0.05, 0.025, 0)
data_ls3$K <- length(data_ls3$alpha)-2

### Fit the model with the publication bias
ex3_bias <- stan(file = 'stan_models/ma_bias.stan', #  'stan_bias.stan',
                 data = data_ls3,
                 chains = 4, cores = 4, seed = 12)
summary(ex3_bias, pars=c("b_Intercept", "omega"))
pairs(ex3_bias, pars=c("omega", "b_Intercept"))
