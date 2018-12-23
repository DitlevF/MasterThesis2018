# MasterThesis2018

This page holds the code for the M.Sc. thesis in Economics at the University of Copenhagen "High Dimensional Sparse 
Econometrics - A Bayesian Perspective" by Ditlev Kiersgaard Frisch.

The files starting with "ss" are spike-and-slab R-functions to estimate high dimensional sparse models, and they cover:

1) Spike-and-slab models with Dirac delta spike and independence slab:
  1.1) ss_individual with an individual Beta prior on the inclusion indicators
  1.2) ss_mixture_MH with a mixture Beta prior on the inclusion indicators with Metropolis-Hastings updating of the parameters
  1.3) ss_mixture with a mixture Beta prior on the inclusion indicators with parameters fixed
  1.3) ss_common with a common Beta prior on the inclusion indicators
  
  1.4) ss_treatment with a two-step treatment prior on the inclusion indicators with a first-step approximation via a Lasso          with a 10-fold cross-validated penalty term.
  
2) Spike-and-slab model with a normal spike and NMIG prior
  2.1) ss_NMIG

       

