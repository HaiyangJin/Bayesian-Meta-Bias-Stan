There are (at least) three kinds of selection models.
1. Smaller p-value intervals correspond to higher probability of publishing, no matter if the effects are positive or negative (I called it "two-sided selection models").
2. Smaller p-value intervals with positive effects correspond to higher probability of publishing.
   1. The publishing probability of studies with small p-value ("significant") but negative effects is the same as those non-significant results. (I called it "one-sided" selection models).
   2. The publishing probability of studies with small p-value ("significant") but negative effects is even smaller than non-significant results. (I do not really think this model will be used often).


---


Some notes for the different Stan model files.

`ma_bias_twosided.stan` is supposed to be the main code.
- It could fit both "two-sided" and "one-sided" models.
- This model uses `Dirichlet` priors for the probabilities (with `cumsum()`).


`ma_bias_twosided_gamma.stan` is modified from `ma_bias_twosided.stan`.
- The main difference is that `gamma` and re-weighting (from `library(publipha)`) is used for the priors of probabilities.


`ma_bias_twosided_gamma_cumsum.stan` is modified from `ma_bias_twosided.stan`.
- The main difference is that `gamma` and `cumsum` (from `library(RoBMA)`) is used for the priors of probabilities.


`psma_stan_gamma_prior.stan` is modified from codes of `library(publipha)`.
- It could only fit "one-sided" models.
- It uses `gamma` priors for probabilities (with re-weighting).


`psma_stan.stan` is also modified from codes of `library(publipha)`. More precisely, it is modified from `psma_stan_gamma_prior.stan`.
- It could only fit "one-sided" models.
- This model uses `Dirichlet` priors for the probabilities (with `cumsum()`).


`ma_bias_twosided_gamma.stan` is modified from `ma_bias_twosided.stan`.
- It should be able to fit both "two-sided" and "one-sided" models.
- It uses `gamma` priors for probabilities (with re-weighting).


`ma_bias_related_codes.stan` stores the additional codes needed for adding selection models to the `brms()` codes (mainly based on `ma_bias_twosided.stan`).
One caveat is that it seems different priors (i.e., with Dirichlet or gamma) lead to slightly different results.
