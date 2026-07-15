# run all specified estimators a single time

# @param mc_id the monte-carlo id
# @param n the sample size
# @param data_args a list of arguments for data generation
# @param cached_datasets either NULL, if we're generating new datasets; or a list
#        of cached datasets
# @param estimators the estimator(s) to run
# @param data_only if TRUE, only generate and return a dataset. Otherwise, run estimators
# @param raking_args a list of arguments for raking
# @param tmle_args a list of arguments for TMLE (including Super Learner arguments)
# @param rescale rescale the IPC weights?
# @param rescale_ps rescale the treatment propensity score (stabilize)?
# @param weight_regressions use IPC weights for propensity and outcome regression?
# @param use_imputation, use_basis, projection_intercept, use_hal args controlling EIF projection estimate
# @param use_meat use meat for sandwich?
# @param parametric is this parametric?
# @param rake_on_outcome should we rake on the outcome?
# @param cgr_tx_type the type of treatment estimator for CGR, defaults to "ipw" ("adjust" is also accepted)
# @param stratified is this a stratified design or Poisson sampling?
# @param strata_form the formula for stratification
# @return the results; either a list of datasets or a list of estimation results
investigate_performance_once <- function(mc_id = 1, n = 100,
                                         data_args = list(scenario = 1,
                                                          gencor = 0),
                                         outcome_name = "Y", tx_name = "X",
                                         missing_indicator = "is.complete",
                                         fam = "binomial",
                                         cached_datasets = NULL,
                                         data_only = FALSE,
                                         estimators = c("raking_cond"),
                                         raking_args = list("n_imp" = 10,
                                                            "rake_on_outcome" = TRUE,
                                                            "rake_on_exposure" = FALSE,
                                                            "rake_on_covars" = FALSE,
                                                            "rake_on_eifs" = TRUE,
                                                            "calibrate_separate" = TRUE),
                                         tmle_args = list(
                                           "g_lib" = c("SL.glm"), "miss_lib" = c("SL.glm"),
                                           "q_lib" = c("SL.glm"), "K" = 5,
                                           "phase1_covars" = c("Y", "X", "Zs", "Zw"),
                                           "phase2_covars" = c("Ws", "Ww")
                                         ), use_imputation = FALSE, use_hal = TRUE,
                                         rescale = TRUE, rescale_ps = FALSE,
                                         weight_regressions = TRUE, use_basis = FALSE,
                                         projection_intercept = TRUE, use_meat = TRUE,
                                         parametric = TRUE, rake_on_outcome = TRUE,
                                         cgr_tx_type = "ipw", stratified = FALSE,
                                         strata_form = "") {
  # generate (or grab) a dataset
  if (isTRUE(is.null(cached_datasets))) {
    df <- gen_data(N = n, scenario = data_args$scenario, gencor = data_args$gencor)
  } else {
    df <- cached_datasets[[mc_id]]
  }
  # if data_only, return it
  if (isTRUE(data_only)) {
    return(df)
  } else {
    if (is.character(fam)) {
      if (fam == "gaussian") {
        fam <- gaussian()
      } else if (fam == "binomial") {
        fam <- binomial()
      }
    }
    # set up formulas
    formulas <- get_formulas(data_args$scenario)
    # set up observed variables
    obs_vars <- c("Y", "X", "Zs", "Zw", "As", "Aw", "Ws_obs", "Ww_obs", "Us", "Uw", "is.complete")
    # obs_vars <- c("Y", "X", "Zs", "Zw", "As", "Aw", "Ws_obs", "Ww_obs", "is.complete")
    # apply the estimators to the given dataset
    all_results <- lapply(as.list(estimators), function(est) {
      run_one_estimator(
        estimator = est, df = df, oracle_formula = formulas$oracle_outcome_formula,
        population_formula = formulas$census_outcome_formula,
        outcome_formula = formulas$working_outcome_formula,
        outcome_formula_factor = formulas$outcome_formula_factor,
        outcome_formula_confounded = formulas$confounded_outcome_formula,
        tx_formula = formulas$propensity_formula,
        tx_formula_factor = formulas$propensity_formula_factor,
        miss_formula = formulas$miss_formula,
        miss_formula_factor = formulas$miss_formula_factor,
        phase_1_formula = formulas$phase_1_formula,
        fam = fam, obs_vars = obs_vars, tmle_args = tmle_args, mi_args = mi_args,
        raking_args = raking_args, outcome_name = outcome_name,
        tx_name = tx_name, missing_indicator = missing_indicator,
        use_imputation = use_imputation, use_hal = use_hal,
        use_meat = use_meat, use_basis = use_basis,
        rescale = rescale, rescale_ps = rescale_ps, weight_regressions = weight_regressions,
        rake_on_outcome = rake_on_outcome, cgr_tx_type = cgr_tx_type,
        stratified = stratified, strata_form = strata_form
      )
    })
    # return results
    results_df <- do.call(rbind, lapply(all_results, function(x) x))
    results_df$mc_id <- mc_id
    return(results_df)
  }
}

# Run a single estimator -------------------------------------------------------

# @param estimator a character string specifying the estimator, e.g., "cc_oracle"
# @param df the dataset
# @param oracle_formula the true outcome regression formula
# @param population_formula the working regression formula (using unobserved variables)
# @param outcome_formula the working outcome regression formula
# @param outcome_formula_factor the working outcome regression formula with explicit factor variables, for TMLE-M
# @param outcome_formula_confounded the confounded outcome regression formula
# @param miss_formula the missing-data formula
# @param miss_formula_factor the missing-data formula with explicit factor variables, for TMLE
# @param tx_formula the treatment assignment formula
# @param tx_formula_factor the treatment assignment formula with explicit factor variables, for TMLE-M
# @param tmle_args a list of arguments for TMLE
# @param mi_args a list of arguments for MI methods
# @param raking_args a list of arguments for raking methods
# @param calOption the calibration option
# @param outcome_name the name of the outcome
# @param tx_name the name of the treatment variable
# @param missing_indicator the name of the missingness indicator
# @param plasmode logical, is this a plasmode simulation?
# @param fam the family, e.g., "binomial"
# @param rescale rescale the weights to have mean 1 before computing SEs?
# @param rescale_ps rescale the treatment propensity score (stabilize)?
# @param weight_regressions use IPC weights for propensity and outcome regression?
# @param use_hal logical, use HAL to estimate EIF projection
# @param use_basis logical, use basis expansion to estimate EIF projection
# @param use_projection_intercept logical, use an intercept in estimating EIF projection
# @param use_meat use the meat from the sandwich to estimate variances?
# @param parametric is this parametric?
# @param rake_on_outcome also rake on the outcome?
# @param cgr_tx_type the type of treatment estimator for CGR, defaults to "ipw" ("adjust" is also accepted)
# @param stratified is this a stratified design or Poisson sampling?
# @param strata_form the formula for stratification
# @return the results of running the estimator on the dataset
run_one_estimator <- function(estimator, df, oracle_formula,
                              population_formula, outcome_formula, outcome_formula_factor,
                              outcome_formula_confounded, miss_formula, miss_formula_factor,
                              tx_formula, tx_formula_factor, phase_1_formula,
                              obs_vars, tmle_args, mi_args,
                              raking_args, outcome_name = "Y", tx_name = "X",
                              missing_indicator = "is.complete", plasmode = FALSE,
                              fam = "binomial", use_hal = TRUE, rescale = TRUE,
                              rescale_ps = FALSE,
                              use_basis = FALSE, projection_intercept = TRUE,
                              use_imputation = FALSE,
                              weight_regressions = TRUE, use_meat = TRUE, parametric = TRUE,
                              rake_on_outcome = TRUE, cgr_tx_type = "ipw",
                              stratified = FALSE, strata_form = "") {
  if (isTRUE(parametric)) {
    this_model <- "parametric"
  } else {
    this_model <- "nonparametric"
  }
  obs_data <- df[, obs_vars]
  if (grepl("CC", estimator, ignore.case = TRUE) | grepl("population", estimator, ignore.case = TRUE)) {
    # complete-case estimators
    if (grepl("oracle", estimator)) {
      this_formula <- oracle_formula
      these_data <- df[, !grepl(missing_indicator, names(df))]
    } else if (grepl("population", estimator)) {
      this_formula <- population_formula
      these_data <- df[, !grepl(missing_indicator, names(df))]
    } else if (grepl("noW", estimator)) {
      this_formula <- outcome_formula_confounded
      these_data <- df[, !grepl(missing_indicator, names(df))]
    } else {
      this_formula <- outcome_formula
      these_data <- df[complete.cases(df), !grepl(missing_indicator, names(df))]
    }
    output <- run_cc(
      data = these_data, formula = this_formula, procedure = estimator,
      coefficient_of_interest = tx_name, fam = fam
    )
  } else if (grepl("ipw", estimator, ignore.case = TRUE) & grepl("marginal", estimator, ignore.case = TRUE) & !grepl("aipw", estimator, ignore.case = TRUE)) {
    # IPW
    output <- run_ipw_marginal(data = obs_data, formula = outcome_formula,
                               miss_formula = miss_formula,
                               tx_formula = tx_formula,
                               outcome_name = outcome_name,
                               family = fam,
                               coefficient_of_interest = tx_name,
                               missing_indicator = missing_indicator, model = this_model,
                               sl_args = list(q_lib = tmle_args$q_lib, g_lib = tmle_args$g_lib,
                                              miss_lib = tmle_args$miss_lib, K = tmle_args$K),
                               rescale = rescale,
                               outcome_tx_estimator = "aipw")
  } else if (grepl("CR", estimator, ignore.case = TRUE) |
             grepl("CGR", estimator, ignore.case = TRUE) |
             (grepl("rak", estimator, ignore.case = TRUE) & grepl("cond", estimator, ignore.case = TRUE))) {
    # CGR: conditional generalized raking
    phase_1_formula <- gsub("+ Aw", "",
                            gsub("As + ", "",
                                 gsub("Y + ", "",
                                      gsub("X + ", "", phase_1_formula, fixed = TRUE),
                                      fixed = TRUE),
                                 fixed = TRUE),
                            fixed = TRUE)
    # check to see what type of CGR it is
    if (grepl("-", estimator, fixed = TRUE)) {
      split_est <- unlist(strsplit(estimator, "-", fixed = TRUE))
      cgr_tx_type <- split_est[2]
    }
    if (grepl("ipw", cgr_tx_type, ignore.case = TRUE) | grepl("iptw", cgr_tx_type, ignore.case = TRUE)) {
      cgr_formula <- "Y ~ X"
    } else {
      cgr_formula <- outcome_formula
    }
    output <- run_raking_conditional(data = obs_data, outcome_formula = cgr_formula,
                                     tx_formula = tx_formula, miss_formula = miss_formula,
                                     phase_1_formula = phase_1_formula,
                                     coefficient_of_interest = tx_name,
                                     missing_indicator = missing_indicator,
                                     family = fam, tx_type = cgr_tx_type, procedure = "CGR",
                                     rake_on_outcome = rake_on_outcome,
                                     rescale_ps = rescale_ps,
                                     stratified = stratified,
                                     strata_form = strata_form)$result %>%
      rename(outcome_tx_estimator = outcome_tx_procedure,
             miss_estimator = miss_procedure)
  } else if (grepl("MR", estimator, ignore.case = TRUE) |
             grepl("MGR", estimator, ignore.case = TRUE) |
             (grepl("rak", estimator, ignore.case = TRUE) & grepl("marg", estimator, ignore.case = TRUE))) {
    # MGR: marginal generalized raking
    output <- run_marginal(data = obs_data, formula = outcome_formula,
                           miss_formula = miss_formula, tx_formula = tx_formula,
                           phase_1_formula = phase_1_formula,
                           outcome_name = outcome_name, family = fam,
                           coefficient_of_interest = tx_name,
                           missing_indicator = missing_indicator, model = this_model,
                           use_projection_intercept = projection_intercept,
                           rescale_miss = rescale, rescale_ps = rescale_ps, use_weights = weight_regressions,
                           use_meat = use_meat, use_imputation = use_imputation,
                           use_hal = use_hal, use_basis = use_basis,
                           rake_on_outcome = rake_on_outcome,
                           outcome_tx_estimator = "aipw", miss_estimator = "MR",
                           stratified = stratified, strata_form = strata_form)$result %>%
      mutate(miss_estimator = "MGR")
  } else if (grepl("AIPCW", estimator, ignore.case = TRUE) & grepl("AIPTW", estimator, ignore.case = TRUE)) {
    # AIPCW-AIPTW
    output <- run_marginal(data = obs_data, formula = outcome_formula,
                           miss_formula = miss_formula, tx_formula = tx_formula,
                           phase_1_formula = phase_1_formula,
                           outcome_name = outcome_name, family = fam,
                           coefficient_of_interest = tx_name,
                           missing_indicator = missing_indicator, model = this_model,
                           use_projection_intercept = projection_intercept,
                           rescale_miss = rescale, rescale_ps = rescale_ps, use_weights = weight_regressions,
                           use_meat = use_meat, use_imputation = use_imputation,
                           use_hal = use_hal, use_basis = use_basis,
                           rake_on_outcome = rake_on_outcome,
                           outcome_tx_estimator = "aipw", miss_estimator = "aipw",
                           stratified = stratified, strata_form = strata_form)$result %>%
      mutate(miss_estimator = "AIPCW", outcome_tx_estimator = "AIPTW")
  } else {
    output <- NULL
  }
  return(output %>% filter(estimand != "OR"))
}
