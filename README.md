# Files and code for "Accounting for uncertainty in dormant life stages in stochastic demographic models"
[DOI: 10.1111/oik.03696](http://onlinelibrary.wiley.com/doi/10.1111/oik.03696/abstract)

**Data files**: 

`dataDroso.csv`: Demographic transitions of Drosophyllum lusitanicum populations recorded in annual censuses (from 2011to 2015) in five populations. These data are used to quantify vital rates of above-ground individuals. 

`dataDrosoSB.csv`: Seed fates (in a binary format) inferred from two experiments. These data are used to quantify the transitions related to the seed bank and associated parameter uncertainties.

In case the reader wishes to forego the step of fitting the Bayesian models, we provided a *mcmcOUT.csv* file with 1000 posterior parameter values for each of the parameters estimated with Bayesian models using uninformative priors.

**R code**:

`BayModel.R`:  Executes and saves the results of a Bayesian model quantifying all vital rates; illustrates basic diagnostics that can be run on the results of an MCMC run (i.e., the posterior parameter distribution) to check for model convergence and autocorrelation of the posterior samples.

`makeIPM.R`: Demonstrates how to construct IPMs including continuous and discrete (seed bank) transitions for (A) mean parameter values and (B) from the parameter distributions of the Bayesian models; saves IPMs for all parameters related to seed-bank ingression, stasis, and ingression. The code is based on the supporting material in Ellner and Rees (2006), Am. Nat., 167, 410-428. 

`perturbVR.R`: Demonstrates how to construct IPMs from perturbed vital rates. Each IPM is obtained by (a) perturbing a vital rate by its mean or standard deviation (see `makeVRmu.R` on constructing mean vital-rate kernels) and (b) constructing a new IPM kernel incorporating the perturbed vital rate. 

`makeIPMmu.R` and `makeVRmu.R`: functions to constructs IPMs and vital-rate kernels, respectively, for average environments. 

`sLambdaSimul.R`: Runs simulations, based on different fire return intervals, of the stochastic population growth rate using IPMs constructed (A) from mean parameter values, (B) from perturbed vital rates, and (C) for each posterior sample of the parameters describing seed-bank ingression (goSB), stasis (staySB) and egression (outSB); calculates the stochastic population growth rate, its elasticities, and the probability of quasi-extinction at time t. The structure of the code is based on Tuljapurkar et al. (2003), Am. Nat., 162, 489-502 and Trotter et al. (2013), Methods Ecol. Evol., 4, 290-298.

`sLambdaRmpi.R`: Implements the simulations of the stochastic population growth rate using parallel processing, where simulations are split into different processors of a supercomputer to greatly speed up computational time. 




