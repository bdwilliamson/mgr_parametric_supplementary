# analyze the VCCC data

# load required packages and functions -----------------------------------------
packages <- c("dplyr", "tidyr", "here", "rprojroot", "sandwich", "survey",
              "future", "future.apply", "marginaleffects")
lapply(packages, library, character.only = TRUE)

this_path <- normalizePath(".", mustWork = FALSE)
proj_root <- rprojroot::find_root_file(criterion = ".projectile", path = this_path)
# anonymized
path2 <- "path_to_VCCC_data"
data_dir <- paste0(path2, "/ade_death_mgr_analysis_dataset/")

if (!dir.exists(paste0(proj_root, "/results/data_analysis/"))) {
  dir.create(paste0(proj_root, "/results/data_analysis/"), recursive = TRUE)
}

source(paste0(proj_root, "/code/00_utils.R"))

# load original dataset --------------------------------------------------------
vccc_dat <- readRDS(file = paste0(data_dir, "vccc_data_clean.rds"))
# define new datasets with exposure set to 0 and 1
newdata_0 <- newdata_1 <- vccc_dat
newdata_0$contains_pi <- 0
newdata_1$contains_pi <- 1

# load sampled datasets --------------------------------------------------------
dataset_list <- readRDS(file = paste0(data_dir, "vccc_sampled_datasets.rds"))

# set up parallelization, for later --------------------------------------------
future::plan(multisession, workers = 20)
set.seed(20260130)
seeds <- future_lapply(as.list(seq_len(length(dataset_list))), FUN = function(x) .Random.seed,
                       future.chunk.size = Inf, future.seed = 20260130)

# define the outcome regression, propensity score, and missing-data model equations -------
outcome_regression_ideal <- "ade_death_bin ~ contains_pi + cd4 + year + age + race + risk"
outcome_regression_unvalidated <- "ade_death_bin_star ~ contains_pi_star + cd4_star + year_star + age + race + risk"

ps_ideal <- "contains_pi ~ cd4 + year + age + race + risk"
ps_unvalidated <- "contains_pi_star ~ cd4_star + year_star + age + race + risk"

miss_model <- "validated ~ contains_pi_star + ade_death_bin_star + cd4_star + year_star + sex + age + race + risk"

# phase 1 variables
phase_1_formula <- "~ ade_death_bin_star + cd4_star + contains_pi_star + year_star + sex + age + race + risk"
# outcome name, exposure name, missing-data indicator
outcome_name <- "ade_death_bin"
errorprone_outcome_name <- "ade_death_bin_star"
exposure_name <- "contains_pi"
missing_indicator <- "validated"
strata_form <- "~ ade_death_bin_star"
# validated and unvalidated estimates ------------------------------------------
# estimate the probability of being validated
# all_miss_fits <- lapply(dataset_list, function(dat) {
#   miss_fit <- glm(as.formula(miss_model), data = dataset_list[[1]], family = "binomial")
#   summary(miss_fit)
# })
# coefs <- unlist(lapply(all_miss_fits, function(x) x$coefficients[3,1]))
miss_fit <- glm("validated ~ contains_pi_star + ade_death_bin_star + cd4_star + year_star + age + race + risk", data = dataset_list[[1]], family = "binomial")
summary(predict(miss_fit, type = "response"))

# estimate the exposure propensity score
exposure_ps_validated <- glm(as.formula(ps_ideal), data = vccc_dat, family = "binomial")
exposure_ps_unvalidated <- glm(as.formula(ps_unvalidated), data = vccc_dat, family = "binomial")

# summary of true g
probs_validated <- predict(exposure_ps_validated, type = "response")
probs_unvalidated <- predict(exposure_ps_unvalidated, type = "response")
summary(probs_validated)
summary(probs_unvalidated)

# estimate the outcome regression
or_validated <- glm(as.formula(outcome_regression_ideal), data = vccc_dat, family = "binomial")
or_unvalidated <- glm(as.formula(outcome_regression_unvalidated), data = vccc_dat, family = "binomial")
# validated exposure-specific means; note these are the same as the g-computation versions
# and approximately the same as the IPTW version without rescaling the PS
exposure_mean_0_validated <- estimate_tx_specific_mean(
  y = vccc_dat$ade_death_bin, a = vccc_dat$contains_pi, a_level = 0, delta = rep(1, nrow(vccc_dat)),
  outcome_regression = predict(or_validated, type = "response"),
  outcome_regression_a = predict(or_validated, newdata = newdata_0, type = "response"),
  prob_ps = 1 - predict(exposure_ps_validated, type = "response"), type = "aipw",
  rescale_ps = FALSE
)
exposure_mean_1_validated <- estimate_tx_specific_mean(
  y = vccc_dat$ade_death_bin, a = vccc_dat$contains_pi, a_level = 1, delta = rep(1, nrow(vccc_dat)),
  outcome_regression = predict(or_validated, type = "response"),
  outcome_regression_a = predict(or_validated, newdata = newdata_1, type = "response"),
  prob_ps = predict(exposure_ps_validated, type = "response"), type = "aipw",
  rescale_ps = FALSE
)
exposure_mean_0_validated_g <- estimate_tx_specific_mean(
  y = vccc_dat$ade_death_bin, a = vccc_dat$contains_pi, a_level = 0, delta = rep(1, nrow(vccc_dat)),
  outcome_regression = predict(or_validated, type = "response"),
  outcome_regression_a = predict(or_validated, newdata = newdata_0, type = "response"),
  prob_ps = 1 - predict(exposure_ps_validated, type = "response"), type = "g_comp",
  rescale_ps = FALSE
)
exposure_mean_1_validated_g <- estimate_tx_specific_mean(
  y = vccc_dat$ade_death_bin, a = vccc_dat$contains_pi, a_level = 1, delta = rep(1, nrow(vccc_dat)),
  outcome_regression = predict(or_validated, type = "response"),
  outcome_regression_a = predict(or_validated, newdata = newdata_1, type = "response"),
  prob_ps = predict(exposure_ps_validated, type = "response"), type = "g_comp",
  rescale_ps = FALSE
)
exposure_mean_0_validated_iptw <- estimate_tx_specific_mean(
  y = vccc_dat$ade_death_bin, a = vccc_dat$contains_pi, a_level = 0, delta = rep(1, nrow(vccc_dat)),
  outcome_regression = predict(or_validated, type = "response"),
  outcome_regression_a = predict(or_validated, newdata = newdata_0, type = "response"),
  prob_ps = 1 - predict(exposure_ps_validated, type = "response"), type = "ipw",
  rescale_ps = FALSE
)
exposure_mean_1_validated_iptw <- estimate_tx_specific_mean(
  y = vccc_dat$ade_death_bin, a = vccc_dat$contains_pi, a_level = 1, delta = rep(1, nrow(vccc_dat)),
  outcome_regression = predict(or_validated, type = "response"),
  outcome_regression_a = predict(or_validated, newdata = newdata_1, type = "response"),
  prob_ps = predict(exposure_ps_validated, type = "response"), type = "ipw",
  rescale_ps = FALSE
)
data.frame(
  type = c("g_comp", "iptw", "aiptw"),
  mu0 = c(exposure_mean_0_validated_g, exposure_mean_0_validated_iptw, exposure_mean_0_validated),
  mu1 = c(exposure_mean_1_validated_g, exposure_mean_1_validated_iptw, exposure_mean_1_validated)
) %>%
  mutate(ate = mu1 - mu0,
         rr = mu1 / mu0)

# unvalidated exposure-specific means
exposure_mean_0_unvalidated <- estimate_tx_specific_mean(
  y = vccc_dat$ade_death_bin_star, a = vccc_dat$contains_pi_star, a_level = 0, delta = rep(1, nrow(vccc_dat)),
  outcome_regression = predict(or_unvalidated, type = "response"),
  outcome_regression_a = predict(or_unvalidated, newdata = newdata_0, type = "response"),
  prob_ps = 1 - predict(exposure_ps_unvalidated, type = "response"), type = "aipw",
  rescale_ps = FALSE
)
exposure_mean_1_unvalidated <- estimate_tx_specific_mean(
  y = vccc_dat$ade_death_bin_star, a = vccc_dat$contains_pi_star, a_level = 1, delta = rep(1, nrow(vccc_dat)),
  outcome_regression = predict(or_unvalidated, type = "response"),
  outcome_regression_a = predict(or_unvalidated, newdata = newdata_1, type = "response"),
  prob_ps = predict(exposure_ps_unvalidated, type = "response"), type = "aipw",
  rescale_ps = FALSE
)

# estimate the ATE and RR
ate_validated <- exposure_mean_1_validated - exposure_mean_0_validated
rr_validated <- exposure_mean_1_validated / exposure_mean_0_validated

ate_unvalidated <- exposure_mean_1_unvalidated - exposure_mean_0_unvalidated
rr_unvalidated <- exposure_mean_1_unvalidated / exposure_mean_0_unvalidated

full_results <- tibble::tibble(
  dataset = c(rep("unvalidated", 4), rep("validated", 4)),
  estimand = rep(c("mu0", "mu1", "ATE", "RR"), 2),
  value = c(exposure_mean_0_unvalidated, exposure_mean_1_unvalidated,
            ate_unvalidated, rr_unvalidated,
            exposure_mean_0_validated, exposure_mean_1_validated,
            ate_validated, rr_validated)
)
saveRDS(full_results, file = paste0(proj_root, "/results/data_analysis/full_data_results.rds"))

# summary of estimated treatment PS
all_ps_summaries <- lapply(as.list(seq_len(length(dataset_list))), FUN = function(i) {
  ps <- glm(as.formula(ps_ideal), data = dataset_list[[i]], family = "binomial",
            subset = (validated == 1))
  data.frame(t(unclass(summary(predict(ps, type = "response")))))
})
ps_summs <- data.table::rbindlist(all_ps_summaries)
summary(ps_summs[, "Min."])
summary(ps_summs[, "Max."])

# CGR --------------------------------------------------------------------------
cgr_results_iptw <- future_lapply(
  X = as.list(seq_len(length(dataset_list))), FUN = function(i) {
    tryCatch(
      run_raking_conditional(
        data = dataset_list[[i]],
        outcome_formula = "ade_death_bin ~ contains_pi",
        tx_formula = ps_ideal,
        miss_formula = miss_model,
        phase_1_formula = phase_1_formula,
        outcome_name = outcome_name,
        errorprone_outcome_name = errorprone_outcome_name,
        coefficient_of_interest = exposure_name,
        missing_indicator = missing_indicator,
        family = binomial(link = "logit"), tx_type = "ipw", procedure = "CGR",
        rake_on_outcome = TRUE, stratified = TRUE, strata_form = strata_form
      )$result %>%
        mutate(dataset = i, .before = "est") %>%
        rename(outcome_tx_estimator = outcome_tx_procedure,
               miss_estimator = miss_procedure),
      error = function(e) {
        message(conditionMessage(e))
        print(paste0("Error occurred when running Monte-Carlo iteration ", i))
        tibble::tibble(miss_estimator = "CGR", outcome_tx_estimator = "ipw",
                       estimand = c("mu0", "mu1", "ATE", "RR", "OR"),
                       dataset = i, est = rep(NA, 5), SE = rep(NA, 5))
      }
    )
  }, future.seed = seeds, future.stdout = structure(TRUE, drop = TRUE),
  future.conditions = "message"
)
saveRDS(cgr_results_iptw, file = paste0(proj_root, "/results/data_analysis/cgr_results_iptw.rds"))

cgr_results_ra <- future_lapply(
  X = as.list(seq_len(length(dataset_list))), FUN = function(i) {
    tryCatch(
      run_raking_conditional(
        data = dataset_list[[i]],
        outcome_formula = outcome_regression_ideal,
        tx_formula = ps_ideal,
        miss_formula = miss_model,
        phase_1_formula = phase_1_formula,
        outcome_name = outcome_name,
        errorprone_outcome_name = errorprone_outcome_name,
        coefficient_of_interest = exposure_name,
        missing_indicator = missing_indicator,
        family = binomial(link = "logit"), tx_type = "adjust", procedure = "CGR",
        rake_on_outcome = TRUE, stratified = TRUE, strata_form = strata_form
      )$result %>%
        mutate(dataset = i, .before = "est") %>%
        rename(outcome_tx_estimator = outcome_tx_procedure,
               miss_estimator = miss_procedure),
      error = function(e) {
        message(conditionMessage(e))
        print(paste0("Error occurred when running Monte-Carlo iteration ", i))
        tibble::tibble(miss_estimator = "CGR", outcome_tx_estimator = "adjust",
                       estimand = c("mu0", "mu1", "ATE", "RR", "OR"),
                       dataset = i, est = rep(NA, 5), SE = rep(NA, 5))
      }
    )
  }, future.seed = seeds, future.stdout = structure(TRUE, drop = TRUE),
  future.conditions = "message"
)
saveRDS(cgr_results_ra, file = paste0(proj_root, "/results/data_analysis/cgr_results_ra.rds"))

# MGR --------------------------------------------------------------------------
mgr_results <- future_lapply(
  X = as.list(seq_len(length(dataset_list))), FUN = function(i) {
    tryCatch(
      run_marginal(
        data = dataset_list[[i]],
        formula = outcome_regression_ideal,
        tx_formula = ps_ideal,
        phase_1_formula = phase_1_formula,
        miss_formula = miss_model,
        outcome_name = outcome_name,
        errorprone_outcome_name = errorprone_outcome_name,
        coefficient_of_interest = exposure_name,
        missing_indicator = missing_indicator,
        family = binomial(link = "logit"),
        use_imputation = TRUE, use_meat = FALSE, rake_on_outcome = TRUE,
        outcome_tx_estimator = "aipw", miss_estimator = "MGR",
        rescale_miss = TRUE, rescale_ps = FALSE, stratified = TRUE, strata_form = strata_form
      )$result %>%
        mutate(dataset = i, .before = "est"),
      error = function(e) {
        message(conditionMessage(e))
        print(paste0("Error occurred when running Monte-Carlo iteration ", i))
        tibble::tibble(miss_estimator = "MGR", outcome_tx_estimator = "aipw",
                       estimand = c("mu0", "mu1", "ATE", "RR", "OR"),
                       dataset = i, est = rep(NA, 5), SE = rep(NA, 5))
      }
    )
  }, future.seed = seeds, future.stdout = structure(TRUE, drop = TRUE),
  future.conditions = "message"
)
saveRDS(mgr_results, file = paste0(proj_root, "/results/data_analysis/mgr_results.rds"))

