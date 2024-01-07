This project provides a framework to jointly model the bulk and tail data in a Bayesian approach.

-------------------------------------------------------------------------
Simulation
-------------------------------------------------------------------------
Run_Simulation_1.1.R  
Run_Simulation_1.2.R  
Run_Simulation_1.3.R
- Code for running simulations of scenarios 1.1, 1.2, and 1.3 in the simulation section.

Run_Simulation_2.R
- Code for fitting the Bivariate Extreme Mixture Model on the Gaussian distributed data, referred to as Scenario 2
  in the following. The result is compared with Lidia's model
  
BEMM_Functions.R  
BEMM_Functions_AFSlice.R
- Nimble code for running MCMC on the Bivariate Extreme Mixture Model. The first uses a blocked 
  random walk sampler, while the second uses the Automated Factor Sliced Sampler.
  
  
Posterior_Dependence.R 
- Code for calculating chi, chi bar, and Kendall's tau of the posterior predictive replicate in Scenario 2

Result_Analysis.R 
- Code for diagnosing the goodness of fit in Scenario 1.1, 1.2, and 1.3 and comparing results from Scenario 2 with Lidia's model
  
plots.R
- Code for illustrating the Bivariate Extreme Mixture Model

Functions.R  
Lidia_model.R
- Code for generating results from the paper Joint modeling of the body and tail of bivariate data.

-------------------------------------------------------------------------
UK_Temp
-------------------------------------------------------------------------
east-sussex_oxfordshire.RData
- The RData contains two data frames:
  1. temp.daily.max.all: daily maxixa air temperature data used in the analysis.
      column City 1 is the record from a station located in bishops lane, Ringmer,
      while column city 2 represents the record from a station located in Model Farm, Shirburn
  2. Y.all: residuals after removing the trend and periodic components.
  
East_sussex_Oxfordshire.R
- Code for fitting the Bivariate Extreme Mixture Model on the daily air temperature data

BEMM_Function_Temp.R
- Nimble code for running the MCMC using the Automated Factor Sliced Sampler

Temp_Result_Analysis.R
- Code for reproducing the analysis results.
  

-------------------------------------------------------------------------
KRSW
-------------------------------------------------------------------------
This folder contains the code for fitting multivariate GPD in the paper
Peaks over thresholds modeling with multivariate generalized Pareto distributions,
check the ReadMe file in the folder for details.


