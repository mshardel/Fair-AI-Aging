# Fair-AI-Aging
Predict post-discharge days at home after hospitalization for hip fracture among Medicare fee-for-service beneficiaries. Programs written by Rong (Jeffery) Zhao, MS.

1. internal_crossvalidation.R fits OLS and average constrained regression and validates the models with multiple metrics using data from study years 2010-2017.
2. external_validation.R fits OLS and average constrained regression using data from study years 2010-2017 and validates/tests the models with multiple metrics using data from study year 2018.
3. alldata_table_figures.R fits OLS and average constrained regression using data from all study years (2010-2018) and produces a table and figures of model coefficients.
