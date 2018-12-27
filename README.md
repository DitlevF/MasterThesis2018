# MasterThesis2018

This page holds the code for the M.Sc. thesis in Economics at the University of Copenhagen "High Dimensional 
Econometrics - A Bayesian Perspective" by Ditlev Kiersgaard Frisch.

The file "ss_example" is an example code for conducting spike-and-slab estimation of
the three types of HDS models discussed in the paper, which requires placing the files ss_common, ss_treatment and
ss_finite_mixture_reg_individual in your own working directory folder.

The files starting with "ss" are spike-and-slab R-functions to estimate high dimensional sparse models, and they cover:

Spike-and-slab models with Dirac delta spike and independence slab
  1) ss_individual with an individual Beta prior on the inclusion indicators
  2) ss_mixture_MH with a mixture Beta prior on the inclusion indicators with a Metropolis-Hastings step for sampling the parameters
  3) ss_mixture_fixed with a mixture Beta prior on the inclusion indicators with parameters fixed
  4) ss_common with a common Beta prior on the inclusion indicators
  
  5) ss_treatment with a two-step treatment prior on the inclusion indicators with a first-step approximation via a Lasso          with a 10-fold cross-validated penalty term.
  
  6) ss_finite_mixture_reg_individual for estimating finite mixture models of linear regressions with an individual beta
  prior on the inclusion indicators.
  
Spike-and-slab model with a normal spike and NMIG prior
  1) ss_NMIG

       
General note: To avoid confusion with the in-built R-function gamma(), the inclusion indicator variables labelled as gamma throughout the paper are labelled as delta in the code.
