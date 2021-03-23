# Bayesian meta-analysis with selection models (publication bias)


Here are the Stan and R codes for performing Bayesian meta-analysis with selection models to take into account publication bias (selections). 

The **assumption** applied in the selection models (implemented here) is that the probability one experiment was published depends on its p-value. More specifically:
  
  1. The p-value range (i.e., from 0 to 1) could be divided into several intervals. And within the same interval, the probablities that experiments were published were the same;
  2. The probability was higher for the interval with smaller p-values.
  3. Experiments in which the p-value was smaller than 0.05 were all published (i.e., with the probability of 100%);  
   
Parts of the weighted probability density distribution is from [library(publipha)](https://github.com/JonasMoss/publipha). 
