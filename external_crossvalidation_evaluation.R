###################################################################################################
# 4/3/26
# Author: Jeffery Zhao, MS
# Department of Biostatistics, Pennsylvania State University
#
# This program performs external validation: training using study years 2010-2017, 
# validation using 2018 data, 50 replications
# Predictions of 6-month days at home using OLS and Average Constrained Regression
# Evaluated using R2, MAE, RMSE predictive ratios by group, mean residuals by group, mean residual difference, 
# and fair covariance.
# Fairness metric quantiles also computed
# Fairness scenarios are Age/Sex, Age, and Sex
# Performance evaluated by demographics: Age, Sex, and Race
####################################################################################################

# Required Libraries
library(haven)
library(dplyr)
library(CVXR)
library(writexl)

# --- 1. HELPER FUNCTIONS ---

rsquared <- function(y, predy) {
  y     <- as.vector(y)
  predy <- as.vector(predy)
  SSR   <- sum((y - predy)^2)
  SST   <- sum((y - mean(y))^2)
  if (SST == 0) return(ifelse(SSR == 0, 1, 0))
  1 - SSR / SST
}

# Six fairness metrics for a group contrast (used for the 3 primary scenarios)
calculate_fairness_metrics <- function(y_true, y_pred, test_grp) {
  y_true <- as.vector(y_true)
  y_pred <- as.vector(y_pred)
  list(
    predictive_ratio_grp = sum(y_pred[test_grp == 1]) / sum(y_true[test_grp == 1]),
    predictive_ratio_ref = sum(y_pred[test_grp == 0]) / sum(y_true[test_grp == 0]),
    mean_residual_grp    = mean(y_pred[test_grp == 1] - y_true[test_grp == 1]),
    mean_residual_ref    = mean(y_pred[test_grp == 0] - y_true[test_grp == 0]),
    mean_residual_diff   = mean(y_pred[test_grp == 1] - y_true[test_grp == 1]) -
                           mean(y_pred[test_grp == 0] - y_true[test_grp == 0]),
    fair_covariance      = as.vector(cov(test_grp, y_true - y_pred))
  )
}

# Extended metrics: adds subgroup-specific MAE/RMSE (used for primary scenarios and RACE/SEX/AGECAT analysis)
calculate_subgroup_metrics <- function(y_true, y_pred, test_grp) {
  base          <- calculate_fairness_metrics(y_true, y_pred, test_grp)
  base$mae_grp  <- mean(abs(y_pred[test_grp == 1] - y_true[test_grp == 1]))
  base$mae_ref  <- mean(abs(y_pred[test_grp == 0] - y_true[test_grp == 0]))
  base$rmse_grp <- sqrt(mean((y_pred[test_grp == 1] - y_true[test_grp == 1])^2))
  base$rmse_ref <- sqrt(mean((y_pred[test_grp == 0] - y_true[test_grp == 0])^2))
  base
}

# Quantile-based metrics for YoungMen scenario (5 metrics x 3 quantiles = 15 values)
calculate_quantile_metrics <- function(y_true, y_pred, test_grp,
                                       probs = c(0.25, 0.5, 0.75)) {
  resid   <- y_pred - y_true
  metrics <- list()
  for (p in probs) {
    tag <- paste0("q", as.integer(p * 100))
    metrics[[paste0("pred_ratio_grp_", tag)]] <-
      quantile(y_pred[test_grp == 1], p) / quantile(y_true[test_grp == 1], p)
    metrics[[paste0("pred_ratio_ref_", tag)]] <-
      quantile(y_pred[test_grp == 0], p) / quantile(y_true[test_grp == 0], p)
    metrics[[paste0("residual_grp_", tag)]]  <- quantile(resid[test_grp == 1], p)
    metrics[[paste0("residual_ref_", tag)]]  <- quantile(resid[test_grp == 0], p)
    metrics[[paste0("residual_diff_", tag)]] <-
      quantile(resid[test_grp == 1], p) - quantile(resid[test_grp == 0], p)
  }
  metrics
}


# --- 2. CORE EVALUATION FUNCTION ---

run_one_replication <- function(data, train_ratio = 0.7) {

  ## Train-Test Split: Train on 70% of 2010-2017, test on all 2018 data
  pool_data  <- data[data$yearAdmission != 8, ]
  test_data  <- data[data$yearAdmission == 8, ]
  train_idx  <- sample(nrow(pool_data), floor(nrow(pool_data) * train_ratio))
  train_data <- pool_data[train_idx, ]

  ## Covariate matrices
  covariates <- c("ALZDMT_PRIOR", "HAC", "LOS_DAY_CNT", "PRE_DAYS",
                  "SURG_TYPE2", "SURG_TYPE3", "TEACH_HOSPITAL",
                  "TERT_BED_SIZE1", "TERT_BED_SIZE2", "yearAdmission",
                  grep("^ELX_GRP_", names(data), value = TRUE))

  X_train <- as.matrix(train_data[, covariates])
  y_train <- train_data$y
  X_test  <- as.matrix(test_data[, covariates])
  y_test  <- test_data$y

  all_results <- list()

  ## Fit OLS, clip to [0, 180], compute global metrics
  beta_ols      <- as.matrix(coef(lm(y_train ~ X_train + 1)))
  pred_ols_test <- pmax(0, pmin(180, cbind(1, X_test) %*% beta_ols))
  r2_ols   <- rsquared(y_test, pred_ols_test)
  mae_ols  <- mean(abs(y_test - pred_ols_test))
  rmse_ols <- sqrt(mean((y_test - pred_ols_test)^2))

  ## Pre-compute demographic subgroup indicators on test set
  # RACE: 1=White (ref), 2=Black (grp)
  sub_race <- ifelse(test_data$RACE == 1, 0L,
               ifelse(test_data$RACE == 2, 1L, NA_integer_))
  # SEX: 2=Female (ref), 1=Male (grp)
  sub_sex  <- ifelse(test_data$SEX == 2, 0L,
               ifelse(test_data$SEX == 1, 1L, NA_integer_))
  # AGECAT: 1=65-74 (ref), 2|3=75+ (grp)
  sub_age  <- ifelse(test_data$AGECAT == 1, 0L,
               ifelse(test_data$AGECAT %in% c(2, 3), 1L, NA_integer_))

  idx_race <- which(!is.na(sub_race))
  idx_sex  <- which(!is.na(sub_sex))
  idx_age  <- which(!is.na(sub_age))

  ## Scenario definitions
  scenarios <- list(
    list(name         = "YoungMen_vs_Others",
         train_def    = list(quote(AGECAT == 1 & SEX == 1)),
         eval_grp_def = quote(AGECAT == 1 & SEX == 1),
         eval_ref_def = "others"),
    list(name         = "Male_vs_Female",
         train_def    = list(quote(SEX == 1)),
         eval_grp_def = quote(SEX == 1),
         eval_ref_def = quote(SEX == 2)),
    list(name         = "Age65-74_vs_Age>=75",
         train_def    = list(quote(AGECAT == 1)),
         eval_grp_def = quote(AGECAT == 1),
         eval_ref_def = quote(AGECAT == 2 | AGECAT == 3))
  )

  ## CVXR variable setup (shared across scenarios)
  k    <- length(beta_ols)
  beta <- Variable(k)
  loss <- sum((y_train - cbind(1, X_train) %*% beta)^2)

  ## Loop over scenarios
  for (scenario in scenarios) {

    ## Fit ACR with scenario constraint, clip predictions
    constraints <- lapply(scenario$train_def, function(def) {
      grp_train <- with(train_data, eval(def))
      sum((cbind(1, X_train)[grp_train, ]) %*% beta) == sum(y_train[grp_train])
    })
    result        <- solve(Problem(Minimize(loss), constraints), solver = "ECOS")
    pred_acr_test <- pmax(0, pmin(180, cbind(1, X_test) %*% result$getValue(beta)))

    r2_acr   <- rsquared(y_test, pred_acr_test)
    mae_acr  <- mean(abs(y_test - pred_acr_test))
    rmse_acr <- sqrt(mean((y_test - pred_acr_test)^2))

    ## Scenario contrast group indicator
    grp_vec <- with(test_data, eval(scenario$eval_grp_def))
    if (is.character(scenario$eval_ref_def) && scenario$eval_ref_def == "others") {
      test_grp      <- as.integer(grp_vec)
      valid_indices <- seq_len(nrow(test_data))
    } else {
      ref_vec           <- with(test_data, eval(scenario$eval_ref_def))
      test_grp          <- rep(NA_integer_, nrow(test_data))
      test_grp[grp_vec] <- 1L
      test_grp[ref_vec] <- 0L
      valid_indices     <- which(!is.na(test_grp))
    }

    ## Primary scenario contrast metrics (r2, mae, rmse are global; mae/rmse _grp/_ref are contrast-specific)
    fair_ols <- calculate_subgroup_metrics(
      y_test[valid_indices], pred_ols_test[valid_indices], test_grp[valid_indices])
    fair_acr <- calculate_subgroup_metrics(
      y_test[valid_indices], pred_acr_test[valid_indices], test_grp[valid_indices])

    all_results[[scenario$name]] <- list(
      ols = c(r2 = r2_ols, mae = mae_ols, rmse = rmse_ols, fair_ols),
      acr = c(r2 = r2_acr, mae = mae_acr, rmse = rmse_acr, fair_acr)
    )

    ## Quantile metrics: YoungMen scenario only
    if (scenario$name == "YoungMen_vs_Others") {
      all_results[["YoungMen_quantiles"]] <- list(
        ols = calculate_quantile_metrics(
          y_test[valid_indices], pred_ols_test[valid_indices], test_grp[valid_indices]),
        acr = calculate_quantile_metrics(
          y_test[valid_indices], pred_acr_test[valid_indices], test_grp[valid_indices])
      )
    }

    ## Subgroup metrics by RACE, SEX, AGECAT using this scenario's model predictions
    sname <- scenario$name
    all_results[[paste0(sname, "_RACE")]] <- list(
      ols = calculate_subgroup_metrics(
        y_test[idx_race], pred_ols_test[idx_race], sub_race[idx_race]),
      acr = calculate_subgroup_metrics(
        y_test[idx_race], pred_acr_test[idx_race], sub_race[idx_race])
    )
    all_results[[paste0(sname, "_SEX")]] <- list(
      ols = calculate_subgroup_metrics(
        y_test[idx_sex], pred_ols_test[idx_sex], sub_sex[idx_sex]),
      acr = calculate_subgroup_metrics(
        y_test[idx_sex], pred_acr_test[idx_sex], sub_sex[idx_sex])
    )
    all_results[[paste0(sname, "_AGECAT")]] <- list(
      ols = calculate_subgroup_metrics(
        y_test[idx_age], pred_ols_test[idx_age], sub_age[idx_age]),
      acr = calculate_subgroup_metrics(
        y_test[idx_age], pred_acr_test[idx_age], sub_age[idx_age])
    )
  }

  all_results
}


# --- 3. AGGREGATION FUNCTION ---

aggregate_results <- function(results_list) {
  scenario_names <- names(results_list[[1]])
  final_sheets   <- list()

  for (scenario in scenario_names) {
    ols_runs     <- lapply(results_list, function(run) run[[scenario]][["ols"]])
    acr_runs     <- lapply(results_list, function(run) run[[scenario]][["acr"]])
    metric_names <- names(ols_runs[[1]])

    summary_data <- t(sapply(metric_names, function(metric) {
      ols_vals <- sapply(ols_runs, `[[`, metric)
      acr_vals <- sapply(acr_runs, `[[`, metric)
      c(Mean_OLS = mean(ols_vals, na.rm = TRUE),
        SD_OLS   = sd(ols_vals,   na.rm = TRUE),
        Mean_ACR = mean(acr_vals, na.rm = TRUE),
        SD_ACR   = sd(acr_vals,   na.rm = TRUE))
    }))

    final_df <- as.data.frame(summary_data)
    final_df <- cbind(Metric = rownames(final_df), final_df)
    rownames(final_df) <- NULL
    final_sheets[[scenario]] <- final_df
  }

  final_sheets
}


# --- 4. MAIN ---

main <- function() {
  data <- read_sas("adrd_analytic_wide_time_masked6.sas7bdat")
  data <- data[!(is.na(data$TERT_BED_SIZE) | is.na(data$HAC)), ]
  data <- data[data$AGECAT != 0, ]
  data$PRE_DAYS       <- rowSums(data[, paste0("PRE_DAYS_AT_HOME_UPD", 1:6)],        na.rm = TRUE)
  data$SURG_TYPE2     <- ifelse(data$SURG_TYPE == 2, 1, 0)
  data$SURG_TYPE3     <- ifelse(data$SURG_TYPE == 3, 1, 0)
  data$TERT_BED_SIZE1 <- ifelse(data$TERT_BED_SIZE == 1, 1, 0)
  data$TERT_BED_SIZE2 <- ifelse(data$TERT_BED_SIZE == 2, 1, 0)
  data$y              <- rowSums(data[, paste0("POST_DSCHRG_DAYS_AT_HOME_UPD", 1:6)], na.rm = TRUE)

  set.seed(133)
  results_list <- lapply(seq_len(50), function(i) run_one_replication(data))
  aggregate_results(results_list)
}

t0 <- proc.time()
final_sheets <- main()
elapsed <- proc.time() - t0
cat(sprintf("Elapsed time: %.1f seconds (%.2f minutes)\n",
            elapsed["elapsed"], elapsed["elapsed"] / 60))

# --- Save results (run manually when ready) ---
write_xlsx(final_sheets, "output/revision/Analysis3_v6_results.xlsx")
