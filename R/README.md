# Description of files in this folder  
  
This folder contains all code to reproduce results and figures.  
  
## Main analysis  
  
-   functions  
    Helper functions used throughout the analysis. Placed in one file to make the other files more legible/ clean  

-   Data_preparation  
    Filtering and cleaning of fossil and temperature data. Merging of both into `final_data.RData` used in the subsequent analysis  

-   pal_int_comparison  
    Fitting of GLMMs. Calculation of difference in medians and effect sizes via bootstrapping and Bayesian estimation  

-   log_odds_phyla  
    Analysis of log-odds ratios per phyla and through time  

-   cont_fragmentation  
    New analysis with continental fragmentation as predictor variable instead of temperature  

-   environmental_differences  
    Comparison of origination responses between onshore and offshore taxa  
  
  
## Robustness tests  
  
-   bayesian_model_check  
    Convergence metrics and posterior predictive checks for the Bayesian estimation  

-   cross_correlation_test  
    Test for cross-correlation between paleo-temperature estimates and the continental fragmentation index  

-   extinction_dynamics    
    Tests for cross-correlation between a variety of metrics for extinction rates and the origination signal used in our analysis  

-   oxygen_isotope_test  
    Test of dependency of our results of choice of compilation and methodology of oxygen isotope data  

-   random_effect_test  
    Test for for the inclusion of a second random effect term in our analysis by comparing a model with both random effects to a model with genus-level variation and a second constant random effect by means of AIC  
  
  
## Additional analysis  
  
-   log_odds_period  
    log-odds ratios at finer time-scales  

-   model_comparison  
    Comparison between tradition models and model with paleoclimate interaction added via information criterion  


## Supplemental information  
  
-   supplemental_plots  
    Additional visualisation  

-   supplemental_tables  
    All tables produced via the `flextable` package and exported to a word file  
  
  
  
