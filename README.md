This repository provides a framework to jointly model the bulk and tail data in a Bayesian approach.

#########################################################################
Simulation (Folder Simulation)
#########################################################################

Main Code 
--------------------
Result_Analysis.R 
- Code for diagnosing the goodness of fit in Scenario 1.1, 1.2, and 1.3 and comparing results from Scenario 2 with the model in André et al. (2024).

Auxiliaries
--------------------
Run_Simulation_1.1.R  
Run_Simulation_1.2.R  
Run_Simulation_1.3.R
- Code for obtaining MCMC results of scenarios 1.1, 1.2, and 1.3 in the simulation section.

Run_Simulation_2.R
- Code for fitting the Bivariate Extreme Mixture Model on the Gaussian distributed data, referred to as Scenario 2
  in the following. The result is compared with that in André et al. (2024).
  
BEMM_Functions.R  
BEMM_Functions_AFSlice.R
- Nimble code for running MCMC on the Bivariate Extreme Mixture Model. The first uses a blocked 
  random walk sampler, while the second uses the Automated Factor Sliced Sampler.
  
Posterior_Dependence.R 
- Code for calculating chi, chi bar, and Kendall's tau of the posterior predictive replicate in Scenario 2.
  
plots.R
- Code for illustrating the Bivariate Extreme Mixture Model.

Functions.R  
Lidia_model.R
- Code for reproducing results from André et al. (2024)

RevExp_U_Functions.r  
CommonFunctions.r  
- Code for fitting multivariate Generalized Pareto Distribution of a U representation and reverse exponential distribution as the generator. These two files are from the supporting material of Kiriliouk et al. (2019).

#########################################################################
Data Application (Folder UK_Temp)
#########################################################################
Main Code 
--------------------
Temp_Result_Analysis.R
- Code for reproducing the analysis results.

Auxiliaries
--------------------
east-sussex_oxfordshire.RData
- Data for analysis. The RData contains two data frames:
  1. temp.daily.max.all: raw daily maxixa air temperature data used in the analysis.
      column City 1 is the record from a station located in bishops lane, Ringmer,
      while column city 2 represents the record from a station located in Model Farm, Shirburn
  2. Y.all: residuals after removing the trend and periodic components.
  
East_sussex_Oxfordshire.R
- Code for fitting the Bivariate Extreme Mixture Model on the daily air temperature data.

BEMM_Function_Temp.R
- Nimble code for running the MCMC using the Automated Factor Sliced Sampler.

-------------------------------------------------------------------------
Reference 
-------------------------------------------------------------------------
André, L., Wadsworth, J. and O'Hagan, A. (2024) Joint modelling of the body and tail of
bivariate data. Computational Statistics & Data Analysis 189, 107841.

Kiriliouk, A., Rootzén, H., Segers, J. and Wadsworth, J. L. (2019) Peaks over thresholds
modeling with multivariate generalized Pareto distributions. Technometrics 61(1), 123–135.
