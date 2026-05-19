# Fair-AI-Aging
Predict post-discharge days at home after hospitalization for hip fracture among Medicare fee-for-service beneficiaries. Programs written by Rong (Jeffery) Zhao, MS.

1. internal_crossvalidation.R fits OLS and average constrained regression and validates the models with multiple metrics using data from study years 2010-2017.
2. external_validation.R fits OLS and average constrained regression using data from study years 2010-2017 and validates/tests the models with multiple metrics using data from study year 2018.
3. alldata_table_figures.R fits OLS and average constrained regression using data from all study years (2010-2018) and produces a table and figures of model coefficients.
4. external_validation_evaluation.R fits OLS and average constrained regression using data from study years 2010-2017 and validates/tests the models with additional metrics and fairness scenarios using data from study year 2018.
5. external_validation_evaluation_sensanalysis.R performs the same analyses and assessments as external_validation_evaluation.R, but includes participants with missing hospital-acquired conditions (HAC).
6. internal_crossvalidation_evaluation.R fits OLS and average constrained regression and validates the models with multiple metrics with additional metrics and fairness scenarios using data from study years 2010-2017.
7. internal_crossvalidation_evaluation.R performs the same analyses and assessments as internal_crossvalidation_evaluation.R, but includes participants with missing hospital-acquired conditions (HAC).