This folder contains all code and data needed to reproduce the analyses in Sections 5 (UK Banks) and 6 (Landslides), and simulation study in the Supplementary Material.


-------------------------------------------------------------------------
Section 5
-------------------------------------------------------------------------

Top level
----------

BanksAnalysisFull.r

- Code for reproducing the UK Banks Analysis


Folder: Data
------------

BanknameW_Adj_Oct16.txt - adjusted weekly closing prices for Bankname
BanknameWDates_Adj_Oct16.txt - corresponding dates 

Folder: Functions
-----------------

Gumbel_T_Functions.r
Gumbel_U_Functions.r
RevExp_T_Functions.r
RevExp_U_Functions.r
MVGauss_T_Functions.r

- functions for simulating and fitting h_U, h_T for the Gumbel and reverse exponential generators, and h_T for the multivariate Gaussian generator. Front end functions (intended to be called directly by the user) have names that fit with notation in the paper; back end functions may not.

CommonFunctions.r

- Some functions used for multiple models

ModelDiagnosticsNewNames.r

- Functions for model diagnostics

-------------------------------------------------------------------------
Section 6
-------------------------------------------------------------------------

Top level
----------

LandslideAnalysisFull.r

- Code for reproducing the Landslide Analysis


Folder: Data
------------

precipitation.csv - daily precipitation amounts in Abisko between 01/01/1913 and 01/01/2015.

Folder: Functions
-----------------

Landslide-analysis-function.R

- Censored likelihood for the ordered components model.
- Functions for model diagnostics.
- Goodness of fit test based on chi.

-------------------------------------------------------------------------
Extra
-------------------------------------------------------------------------

Top level
---------

ExampleUse.r

- Ancillary code to help the interested reader understand the functionality of key simulation and fitting functions

SimulationForSM.r

- Code for the simulation results in the SM
