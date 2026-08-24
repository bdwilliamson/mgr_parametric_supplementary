# generally useful functions for running GR, AIPW, TMLE

# SL function for two-stage GLM ------------------------------------------------
SL.glm.twostage <- function(Y, X, newX, family, obsWeights, model = TRUE, ...) {
  contains_Y <- (names(X) == "Y")
  if (any(contains_Y)) {
    names(X)[contains_Y] <- "Y.orig"
  }
  SL.glm(Y = Y, X = X, newX = newX, family = family, obsWeights = obsWeights, model = model, ...)
}


# efficient influence functions (EIFs) for mRD, mRR, mOR -----------------------

# Logistic regression influence function  - for Raking/survey calibration
# Tong's binary influence function for binary data
# @param fit a glm
inf.fun.logit <- function(fit) {
  dm <- model.matrix(fit)
  Ihat <- (t(dm) %*% (dm * fit$fitted.values * (1 - fit$fitted.values))) / nrow(dm)
  ## influence function
  infl <- (dm * resid(fit, type = "response")) %*% solve(Ihat)
  infl
}

# general EIFs for a regression model
# @param fit a glm
# @param use_meat should we use the meat to estimate the inverse information?
get_glm_eif <- function(fit, use_meat = FALSE, HC = "HC3") {
  score <- sandwich::estfun(fit)
  if (use_meat) {
    inv_information <- solve(sandwich::meatHC(fit, type = HC))
  } else {
    inv_information <- sandwich::bread(fit)
  }
  eif <- t(inv_information %*% t(score))
  return(eif)
}

# the ideal-data influence functions for g-computation, IPW, and AIPW estimators of the treatment-specific mean
# note that each has a difference influence function, not necessarily the efficient one
if_tx_mean_nonparametric <- function(a = NULL, a_level = 1, y = NULL, w = NULL,
                                     Q = NULL, Q_a = NULL, g = NULL,
                                     psi = NULL,
                                     type = "aipw") {
  if (type == "ipw") {
    eif <- (as.numeric(a == a_level) / g) * y - psi
  } else if (type == "g_comp") {
    eif <- Q_a - psi
  } else {
    eif <- (as.numeric(a == a_level) / g) * (y - Q) + Q_a - psi
  }
  return(eif)
}
# the full-data parametric EIF for the exposure-specific mean, under a logistic
# outcome regression model
# @param nonparametric_eif the nonparametric EIF for the treatment-specific mean,
#                          from a call to eif_tx_mean_nonparametric()
# @param glm_fit a fitted glm object
eif_tx_mean_parametric <- function(nonparametric_eif = NULL, glm_fit = NULL,
                                   use_meat = TRUE, HC = "HC3") {
  # extract the score, an n x k matrix
  score <- sandwich::estfun(glm_fit)
  # estimate the inverse Fisher information, a k by k matrix; bread is the inverse Fisher information
  if (use_meat) {
    inv_information <- solve(sandwich::meatHC(glm_fit, type = HC))
  } else {
    inv_information <- sandwich::bread(glm_fit)
  }
  # inv_information <- solve(sandwich::meatHC(glm_fit, type = "HC0"))
  # parametric EIF is n x k matrix
  parametric_eif <- t(inv_information %*% t(score))
  # get the EIF. note nonparametric_eif is an n x 1 vector, score is an n x k matrix, information
  eif <- as.numeric(t(colMeans(sweep(score, 1, nonparametric_eif, "*"))) %*% t(parametric_eif))
  # eif <- as.numeric(rowMeans(t(sweep(score, 1, nonparametric_eif, "*"))) %*% parametric_eif)
  return(eif)
}

# the full-data EIF for the marginal risk difference (RD), odds ratio (OR), or
# relative risk (RR)
# @param mu1 the estimated treatment-specific mean under A = 1
# @param mu0 the estimated treatment-specific mean under A = 0
# @param eif1 the estimated full-data EIF for the treatment-specific mean under A = 1
# @param eif0 the estimated full-data EIF for the treatment-specific mean under A = 0
# @param estimand either "mRD", "mOR", or "mRR"
# @return the EIF values. Note that for mOR and mRR, the EIF is for the log mRR and log mRR
full_data_eif <- function(mu1 = NULL, mu0 = NULL, eif1 = NULL, eif0 = NULL, estimand = "mRD") {
  if (grepl("RD", estimand, ignore.case = TRUE)) {
    eif <- eif1 - eif0
  } else if (grepl("RR", estimand, ignore.case = TRUE) & !grepl("nolog", estimand, ignore.case = TRUE)) {
    eif <- eif1 / mu1 - eif0 / mu0
  } else if (grepl("OR", estimand, ignore.case = TRUE)) {
    # eif <- (eif1 + mu1) / (mu1 * (1 - mu1)) - (eif0 + mu0) / (mu0 * (1 - mu0))
    eif <- (1 / (mu1 * (1 - mu1))) * eif1 - (1 / (mu0 * (1 - mu0))) * eif0
  } else if (grepl("RR", estimand, ignore.case = TRUE) & grepl("nolog", estimand, ignore.case = TRUE)) {
    eif <- eif1 / mu0 - eif0 * mu1 / (mu0 ^ 2)
  } else {
    stop("The requested estimand has not yet been implemented.")
  }
  return(eif)
}

# observed-data influence function
# @param delta indicator of being observed (1) or missing
# @param prob_obs probability of being observed
# @param full_eif the full-data EIF
# @param eif_projection the projection of the full-data EIF onto the observed data
observed_data_eif <- function(delta = NULL, prob_obs = NULL, full_eif = NULL,
                              eif_projection = NULL) {
  obs_eif <- vector("numeric", length = length(delta))
  obs_eif[delta == 1] <- 1 / prob_obs[delta == 1] * (full_eif - eif_projection[delta == 1]) + eif_projection[delta == 1]
  obs_eif[delta == 0] <- eif_projection[delta == 0]
  # obs_eif <- as.numeric(delta == 1) / prob_obs * full_eif -
  #   (as.numeric(delta == 1) - prob_obs) / prob_obs * eif_projection
  return(obs_eif)
}

# for SL: get the correct number of CV folds -----------------------------------
# get the effective sample size
# @param n the sample size
# @param n_rare the number of events in the minority class (if outcome_type == "binary")
# @param outcome_type the type of outcome ("continuous" or "binary")
# @return the effective sample size
get_n_eff <- function(n = 100, n_rare = 100, outcome_type = "continuous") {
  if (outcome_type == "continuous") {
    n_eff <- n
  } else {
    n_eff <- min(n, 5 * n_rare)
  }
  return(n_eff)
}

# Set the number of cross-validation
# folds as a function of n.effective
# See Phillips 2023 doi.org/10.1093/ije/dyad023
# @param n.effective the effective sample size
# @return the number of CV folds for SL
setV <- function(n.effective = 100) {
  if (n.effective <= 30) {
    V <- n.effective
  } else if (n.effective <= 500) {
    V <- 20
  } else if (n.effective <= 1000) {
    V <- 10
  } else if (n.effective <= 10000) {
    V <- 5
  } else {
    V <- 2
  }
  return(V)
}


# estimate propensity score, etc. ----------------------------------------------

# estimate a propensity score
# @param data the dataset
# @param formula the missing-data formula
# @param variable the variable name of interest
# @param model "parametric" to use a glm, otherwise use SL
# @param sl_args arguments to pass to superlearner
get_propensity_score <- function(data = NULL, formula = NULL,
                                 variable = "Delta", model = "parametric",
                                 sl_args = list(), use_weights = TRUE,
                                 weights = rep(1, nrow(data))) {
  if (use_weights) {
    ipw_weights <- weights
  } else {
    ipw_weights <- rep(1, nrow(data))
  }
  if (model == "parametric" | (length(sl_args$lib) == 1 & all(sl_args$lib == "SL.glm"))) {
    data$ipw_weights <- ipw_weights
    fit <- glm(formula = as.formula(formula), family = "quasibinomial", data = data,
               weights = ipw_weights)
    propensity_score_1 <- predict(fit, type = "response")
    propensity_score_0 <- 1 - propensity_score_1
  } else {
    these_terms <- terms.formula(as.formula(formula))
    these_data <- data[, gsub(" ", "", rownames(attr(these_terms, "factors")))]
    if (length(sl_args$lib) == 1 & all(sl_args$lib == "SL.glm")) {
      family <- "quasibinomial"
    } else {
      family <- "binomial"
    }
    fit <- SuperLearner::SuperLearner(
      Y = these_data[[variable]], X = these_data[, !grepl(variable, names(these_data))],
      SL.library = sl_args$lib, cvControl = list(V = sl_args$K,
                                                 stratifyCV = sl_args$stratifyCV),
      family = family, method = "method.CC_nloglik", obsWeights = ipw_weights
    )
    pred <- predict(fit, onlySL = TRUE)$pred
    propensity_score_1 <- pred
    propensity_score_0 <- (1 - pred)
    num_coefs <- NA
  }
  ip_weights <- 1 / propensity_score_1
  return(list("ps_1" = propensity_score_1, "ps_0" = propensity_score_0,
              "ip_weights" = ip_weights, "fit" = fit))
}

get_outcome_regression <- function(data = NULL, formula = NULL,
                                   coefficient_of_interest = "A",
                                   model = "parametric", family = "gaussian",
                                   use_weights = TRUE,
                                   weights = rep(1, nrow(data))) {
  data$wts <- rep(1, nrow(data))
  if (use_weights) {
    data$wts <- weights
  }
  if (model == "parametric") {
    if (family$family == "binomial") {
      if (family$link == "logit") {
        fam2 <- quasibinomial()
      } else {
        fam2 <- quasibinomial(link = "probit")
      }
    } else {
      fam2 <- gaussian()
    }
    fit <- glm(formula = as.formula(formula), data = data, family = fam2,
               weights = wts)
    preds <- predict(fit, type = "response")
    data0 <- data1 <- data
    data0[[coefficient_of_interest]] <- 0
    data1[[coefficient_of_interest]] <- 1
    preds_0 <- predict(fit, newdata = data0, type = "response")
    preds_1 <- predict(fit, newdata = data1, type = "response")
  }
  return(list("preds" = preds, "preds_0" = preds_0, "preds_1" = preds_1, "fit" = fit))
}

# get the full-data EIF for the treatment-specific mean and
# outcome regression predictions, using the specified model (parametric or SL)
# @param data the dataset
# @param formula the outcome regression formula
# @param coefficient_of_interest the name of the exposure variable
# @param outcome_name the name of the outcome variable
# @param ps_1 the propensity of being assigned exposure level 1
# @param ps_0 the propensity of being assigned exposure level 0
# @param miss_prob the estimated probability of being observed
# @param model "parametric" or "SL"
# @param sl_args list of arguments for SL
# @param use_weights logical, should we use IPW weights in estimating the outcome regression, defaults to true
get_eifs <- function(data = NULL, formula = NULL, coefficient_of_interest = "A",
                     outcome_name = "Y", ps_1 = NULL, ps_0 = NULL, miss_prob = rep(1, length(ps_1)),
                     model = "parametric", sl_args = list(), use_weights = TRUE) {
  # inverse probability weights
  if (use_weights) {
    ipw_weights <- 1 / miss_prob
  } else {
    ipw_weights <- rep(1, length(miss_prob))
  }
  scaled_ipw_weights <- ipw_weights / mean(ipw_weights)
  # subset to only the variables of interest
  these_terms <- terms.formula(as.formula(formula))
  these_data <- data[, gsub(" ", "", rownames(attr(these_terms, "factors")))]

  newdata_0 <- these_data
  newdata_0[[coefficient_of_interest]] <- 0
  newdata_1 <- these_data
  newdata_1[[coefficient_of_interest]] <- 1
  outcome <- these_data[[outcome_name]]
  if (model == "parametric") {
    if (length(unique(outcome)) == 2) {
      family <- "quasibinomial"
    } else {
      family <- "gaussian"
    }
    outcome_regression <- glm(formula = as.formula(formula),
                              family = family, data = these_data,
                              weights = scaled_ipw_weights)
    outcome_regression_unweighted <- glm(formula = as.formula(formula),
                                         family = family, data = these_data)
    preds <- predict(outcome_regression, type = "response")
    preds_0 <- predict(outcome_regression, newdata = newdata_0, type = "response")
    preds_1 <- predict(outcome_regression, newdata = newdata_1, type = "response")
  } else {
    covariates <- these_data[, -1] # always first
    if (length(unique(outcome)) == 2) {
      if (length(sl_args$q_lib) == 1 & all(sl_args$q_lib == "SL.glm")) {
        family <- "quasibinomial"
      } else {
        family <- "binomial"
      }
      method <- "method.CC_nloglik"
    } else {
      family <- "gaussian"
      method <- "method.CC_LS"
    }

    outcome_regression <- SuperLearner::SuperLearner(Y = outcome,
                                                     X = covariates,
                                                     SL.library = sl_args$q_lib,
                                                     cvControl = list(V = sl_args$K,
                                                                      stratifyCV = sl_args$stratifyCV),
                                                     family = family,
                                                     method = method,
                                                     obsWeights = scaled_ipw_weights)
    preds <- predict(outcome_regression, onlySL = TRUE)$pred
    preds_0 <- predict(outcome_regression, newdata = newdata_0[, -1], onlySL = TRUE)$pred
    preds_1 <- predict(outcome_regression, newdata = newdata_1[, -1], onlySL = TRUE)$pred
  }
  est_eif_nonparametric_1 <- eif_tx_mean_nonparametric(a = these_data[[coefficient_of_interest]],
                                                       a_level = 1, y = outcome,
                                                       Q = preds, Q_a = preds_1,
                                                       g = ps_1)
  est_eif_nonparametric_0 <- eif_tx_mean_nonparametric(a = these_data[[coefficient_of_interest]],
                                                       a_level = 0, y = outcome,
                                                       Q = preds, Q_a = preds_0,
                                                       g = ps_0)
  # if assuming a parametric model, get the EIF in that model
  if (model == "parametric") {
    eif_1 <- eif_tx_mean_parametric(nonparametric_eif = est_eif_nonparametric_1,
                                    glm_fit = outcome_regression_unweighted)
    eif_0 <- eif_tx_mean_parametric(nonparametric_eif = est_eif_nonparametric_0,
                                    glm_fit = outcome_regression_unweighted)
  } else {
    eif_1 <- est_eif_nonparametric_1
    eif_0 <- est_eif_nonparametric_0
  }
  return(list("eif_1" = eif_1, "eif_0" = eif_0,
              "or_preds_1" = preds_1, "or_preds_0" = preds_0))
}

# estimate the optimal raking variable using MI
# (expectation of full-data EIF on the entire dataset conditional on observed varaibles)
# @param formula the outcome regression formula
# @param tx_formula the treatment-assignment formula
# @param phase_1_formula the variables fully observed in phase 1
# @param data the dataset
# @param missing_indicator the missing-data id column
# @param n_imp number of imputations
# @param conditional are we estimating a conditional parameter first?
# @param use_sl logical, should we use SL (true) or MI to estimate projection
# @param model "parametric" or "SL", for outcome regression
# @param sl_args list of args for SL
# @param miss_prob the estimated probability of being observed
# @param tx_type the type of treatment-propensity adjustment (adjust, IPW, etc.)
# @return the optimal raking variable
estimate_eif_projection_mi <- function(outcome_formula = NULL, tx_formula = NULL,
                                       phase_1_formula = NULL,
                                       data = NULL, missing_indicator = "Delta",
                                       coefficient_of_interest = "A",
                                       outcome_name = "Y",
                                       n_imp = 10, conditional = TRUE, use_sl = FALSE,
                                       model = "parametric", sl_args = list(),
                                       family = "binomial",
                                       miss_prob = NULL, tx_type = "adjust",
                                       rescale_miss = TRUE, rescale_ps = FALSE, use_meat = TRUE) {
  if (is.character(family)) {
    family <- match.fun(family)()
  }
  init_fit <- glm(formula = as.formula(outcome_formula), family = family, data = data)
  num_coefs <- length(coefficients(init_fit))
  # get optimal raking variable: expected value of influence function by averaging over MI datasets
  # include any variables present in either "formula" or "tx_formula"
  if (tx_type == "adjust") {
    combined_formula <- outcome_formula
  } else {
    combined_formula <- paste0(outcome_formula, gsub(paste0(coefficient_of_interest, " ~ "), " + ", tx_formula))
  }
  if (!is.null(phase_1_formula)) {
    combined_formula <- paste0(combined_formula, gsub("~", " + ", phase_1_formula))
  }
  data_mi <- model.frame(as.formula(combined_formula), data = data, na.action = NULL)
  phase_2_indx <- (data[[missing_indicator]] == 1)
  # add on fake phase 2 data with missing values; this is for averaging
  miss_phase2 <- data_mi[data[[missing_indicator]] == 1, ]
  miss_cols <- which(colSums(is.na(data_mi)) > 0)
  miss_phase2[, miss_cols] <- NA
  data_mi2 <- rbind(data_mi, miss_phase2)
  fake_phase_2 <- ((1:nrow(data_mi2)) %in% (nrow(data_mi) + 1):nrow(data_mi2))
  data_mi <- data_mi2
  # start the imputation
  init_mi <- mice::mice(data_mi, maxit = 0)
  pred_matrix <- init_mi$predictorMatrix
  data_imputed <- mice::mice(data_mi, predictorMatrix = pred_matrix, m = n_imp,
                             maxit = 20, print = FALSE)
  if (conditional) {
    infMat_all <- array(data = 0, dim = c(nrow(data), num_coefs, n_imp))
    for (i in seq_len(n_imp)) {
      imp_init_i <- mice::complete(data_imputed, i)
      imputation_i <- imp_init_i[1:nrow(data), ]
      imputation_i[phase_2_indx, ] <- imp_init_i[fake_phase_2, ]
      # if not doing adjustment, fit the treatment propensity model
      if (!grepl("adjust", tx_type)) {
        tx_ps_list <- get_propensity_score(
          data = imputation_i, formula = tx_formula, variable = coefficient_of_interest,
          model = "parametric", use_weights = FALSE
        )
        tx_ps_weights <- rep(NA, length(tx_ps_list$ip_weights))
        p_1 <- mean(imputation_i[[coefficient_of_interest]] == 1)
        p_0 <- mean(imputation_i[[coefficient_of_interest]] == 0)
        g_1 <- switch(
          as.numeric(rescale_ps) + 1, tx_ps_list$ps_1, tx_ps_list$ps_1 / p_1
        )
        g_0 <- switch(
          as.numeric(rescale_ps) + 1, tx_ps_list$ps_0, tx_ps_list$ps_0 / p_0
        )
        tx_ps_weights[imputation_i[[coefficient_of_interest]] == 1] <- 1 / g_1[imputation_i[[coefficient_of_interest]] == 1]
        tx_ps_weights[imputation_i[[coefficient_of_interest]] == 0] <- 1 / g_0[imputation_i[[coefficient_of_interest]] == 0]
        this_fam <- switch(as.numeric(family$family != "gaussian") + 1,
                           gaussian(),
                           switch(as.numeric(family$link != "logit") + 1,
                                  quasibinomial(link = "logit"),
                                  quasibinomial(link = "probit")))
      } else {
        tx_ps_weights <- rep(1, nrow(imputation_i))
        this_fam <- family
      }
      imputation_i$weights <- tx_ps_weights
      mi_fit <- glm(formula = as.formula(outcome_formula), family = this_fam, data = imputation_i,
                    weights = weights)
      infMat_all[,,i] <- get_glm_eif(mi_fit, use_meat = use_meat)
    }
    obs_data_eif <- rowMeans(infMat_all, dims = 2)
    obs_data_eif_rd <- obs_data_eif_rr <- obs_data_eif_rr_nolog <- obs_data_eif_or <- obs_data_eif_1 <- obs_data_eif_0 <- obs_data_eif
  } else {
    inf_mat_1 <- matrix(0, nrow = nrow(data), ncol = n_imp)
    inf_mat_0 <- matrix(0, nrow = nrow(data), ncol = n_imp)
    inf_mat_rd <- matrix(0, nrow = nrow(data), ncol = n_imp)
    inf_mat_rr <- matrix(0, nrow = nrow(data), ncol = n_imp)
    inf_mat_rr_nolog <- matrix(0, nrow = nrow(data), ncol = n_imp)
    inf_mat_or <- matrix(0, nrow = nrow(data), ncol = n_imp)
    outcome_name <- gsub(" ", "", strsplit(outcome_formula, "~", fixed = TRUE)[[1]][1])
    delta <- rep(1, nrow(data))
    or_terms <- terms.formula(as.formula(outcome_formula))
    or_varnames <- gsub(" ", "", rownames(attr(or_terms, "factors")))
    for (i in seq_len(n_imp)) {
      # impute and use predicted values for missing-data columns
      imp_init_i <- mice::complete(data_imputed, i)
      imputation_i <- imp_init_i[1:nrow(data), ]
      imputation_i[phase_2_indx, ] <- imp_init_i[fake_phase_2, ]
      imputation_i_with_indicator <- cbind.data.frame(imputation_i, data[[missing_indicator]])
      names(imputation_i_with_indicator)[ncol(imputation_i_with_indicator)] <- missing_indicator
      or_data <- imputation_i_with_indicator[, or_varnames]
      fully_observed_cols <- (colSums(is.na(or_data)) == 0)
      fully_observed_data <- or_data[complete.cases(or_data), fully_observed_cols]
      # get probability of treatment assignment
      prop_score <- get_propensity_score(data = imputation_i_with_indicator,
                                         formula = tx_formula,
                                         model = model,
                                         variable = coefficient_of_interest,
                                         sl_args = list(lib = sl_args$g_lib, K = sl_args$K,
                                                        stratifyCV = sl_args$stratifyCV),
                                         use_weights = FALSE)
      # get outcome regression model
      or_list <- get_outcome_regression(data = imputation_i_with_indicator,
                                        formula = as.formula(outcome_formula),
                                        coefficient_of_interest = coefficient_of_interest,
                                        model = model, family = family, use_weights = FALSE)
      mu_1 <- estimate_tx_specific_mean(a = imputation_i_with_indicator[[coefficient_of_interest]],
                                        a_level = 1, y = imputation_i_with_indicator[[outcome_name]],
                                        delta = delta, outcome_regression = or_list$preds,
                                        outcome_regression_a = or_list$preds_1,
                                        prob_obs = delta, prob_ps = prop_score$ps_1,
                                        eif = rep(0, length(delta)), type = tx_type,
                                        type_miss = "ipw", rescale_miss = rescale_miss, rescale_ps = rescale_ps)
      mu_0 <- estimate_tx_specific_mean(a = imputation_i_with_indicator[[coefficient_of_interest]],
                                        a_level = 0, y = imputation_i_with_indicator[[outcome_name]],
                                        delta = delta, outcome_regression = or_list$preds,
                                        outcome_regression_a = or_list$preds_0,
                                        prob_obs = delta, prob_ps = prop_score$ps_0,
                                        eif = rep(0, length(delta)), type = tx_type,
                                        type_miss = "ipw", rescale_miss = rescale_miss, rescale_ps = rescale_ps)
      # get EIFs
      eifs <- get_ideal_data_eifs(data = imputation_i_with_indicator,
                                  Q = or_list$preds, Q_1 = or_list$preds_1,
                                  Q_0 = or_list$preds_0, g_1 = prop_score$ps_1,
                                  g_0 = prop_score$ps_0,
                                  coefficient_of_interest = coefficient_of_interest,
                                  type = "ipw", psi_1 = mu_1, psi_0 = mu_0,
                                  outcome_name = outcome_name, glm_fit = or_list$fit,
                                  use_meat = use_meat)
      eif_1 <- eifs$if_1
      inf_mat_1[, i] <- eif_1
      eif_0 <- eifs$if_0
      inf_mat_0[, i] <- eif_0
      inf_mat_rd[, i] <- full_data_eif(mu1 = mu_1, mu0 = mu_0,
                                       eif1 = eif_1, eif0 = eif_0,
                                       estimand = "RD")
      inf_mat_rr[, i] <- full_data_eif(mu1 = mu_1, mu0 = mu_0,
                                       eif1 = eif_1, eif0 = eif_0,
                                       estimand = "RR")
      inf_mat_rr_nolog[, i] <- full_data_eif(mu1 = mu_1, mu0 = mu_0,
                                             eif1 = eif_1, eif0 = eif_0,
                                             estimand = "nologRR")
      inf_mat_or[, i] <- full_data_eif(mu1 = mu_1, mu0 = mu_0,
                                       eif1 = eif_1, eif0 = eif_0,
                                       estimand = "OR")
    }
    obs_data_eif <- NULL
    obs_data_eif_1 <- rowMeans(inf_mat_1)
    obs_data_eif_0 <- rowMeans(inf_mat_0)
    obs_data_eif_rd <- rowMeans(inf_mat_rd)
    obs_data_eif_rr <- rowMeans(inf_mat_rr)
    obs_data_eif_rr_nolog <- rowMeans(inf_mat_rr_nolog)
    obs_data_eif_or <- rowMeans(inf_mat_or)
  }
  return(list("eif" = obs_data_eif, "eif_rr" = obs_data_eif_rr,
              "eif_rd" = obs_data_eif_rd, "eif_or" = obs_data_eif_or,
              "eif_1" = obs_data_eif_1, "eif_0" = obs_data_eif_0,
              "eif_rr_nolog" = obs_data_eif_rr_nolog))
}

# estimate the optimal raking variable using regression
# (expectation of full-data EIF on the entire dataset conditional on observed varaibles)
# @param eif_1, eif_0, eif_rd, eif_rr, eif_or EIFs for mu1, mu0, RD, RR, OR
# @param fully_observed_data rows and columns for people with fully-observed data
# @param fully_observed_cols columns that are always observed
# @param or_data dataset (complete cases) for the outcome regression
# @param use_hal logical, use HAL to estimate EIF projection
# @param use_basis logical, use basis expansion to estimate EIF projection
# @param use_projection_intercept logical, use an intercept in estimating EIF projection
# @return the optimal raking variable
estimate_eif_projection_regression <- function(eif_1 = NULL, eif_0 = NULL,
                                               eif_rd = NULL, eif_rr = NULL, eif_or = NULL,
                                               fully_observed_data = NULL,
                                               fully_observed_cols = NULL, or_data = NULL,
                                               use_hal = FALSE, use_basis = FALSE,
                                               use_projection_intercept = TRUE) {
  if (!use_hal) {
    if (use_basis) {
      projection_fit_1 <- glm(
        stats::as.formula(paste("eif_1 ~", paste("splines::bs(", names(fully_observed_data), ", df = 3)", collapse = " + "))),
        data = fully_observed_data
      )
      projection_fit_0 <- glm(
        stats::as.formula(paste("eif_0 ~", paste("splines::bs(", names(fully_observed_data), ", df = 3)", collapse = " + "))),
        data = fully_observed_data
      )
      projection_fit_rd <- glm(
        stats::as.formula(paste("eif_rd ~", paste("splines::bs(", names(fully_observed_data), ", df = 3)", collapse = " + "))),
        data = fully_observed_data
      )
      projection_fit_rr <- glm(
        stats::as.formula(paste("eif_rr ~", paste("splines::bs(", names(fully_observed_data), ", df = 3)", collapse = " + "))),
        data = fully_observed_data
      )
      projection_fit_or <- glm(
        stats::as.formula(paste("eif_or ~", paste("splines::bs(", names(fully_observed_data), ", df = 3)", collapse = " + "))),
        data = fully_observed_data
      )
    } else {
      if (use_projection_intercept) {
        projection_fit_1 <- glm(eif_1 ~ ., data = fully_observed_data)
        projection_fit_0 <- glm(eif_0 ~ ., data = fully_observed_data)
        projection_fit_rd <- glm(eif_rd ~ ., data = fully_observed_data)
        projection_fit_rr <- glm(eif_rr ~ ., data = fully_observed_data)
        projection_fit_or <- glm(eif_or ~ ., data = fully_observed_data)
      } else {
        projection_fit_1 <- glm(eif_1 ~ . - 1, data = fully_observed_data)
        projection_fit_0 <- glm(eif_0 ~ . - 1, data = fully_observed_data)
        projection_fit_rd <- glm(eif_rd ~ . - 1, data = fully_observed_data)
        projection_fit_rr <- glm(eif_rr ~ . - 1, data = fully_observed_data)
        projection_fit_or <- glm(eif_or ~ . - 1, data = fully_observed_data)
      }
    }
    if_projections_1 <- predict(projection_fit_1, newdata = or_data[, fully_observed_cols])
    if_projections_0 <- predict(projection_fit_0, newdata = or_data[, fully_observed_cols])
    if_projections_rd <- predict(projection_fit_rd, newdata = or_data[, fully_observed_cols])
    if_projections_rr <- predict(projection_fit_rr, newdata = or_data[, fully_observed_cols])
    if_projections_or <- predict(projection_fit_or, newdata = or_data[, fully_observed_cols])
  } else {
    projection_fit_1 <- hal9001::fit_hal(
      X = fully_observed_data, Y = eif_1,
      fit_control = list(cv_select = TRUE, use_min = TRUE, lambda.min.ratio = 1e-4),
      family = "gaussian"
    )
    projection_fit_0 <- hal9001::fit_hal(
      X = fully_observed_data, Y = eif_0,
      fit_control = list(cv_select = TRUE, use_min = TRUE, lambda.min.ratio = 1e-4),
      family = "gaussian"
    )
    projection_fit_rd <- hal9001::fit_hal(
      X = fully_observed_data, Y = eif_rd,
      fit_control = list(cv_select = TRUE, use_min = TRUE, lambda.min.ratio = 1e-4),
      family = "gaussian"
    )
    projection_fit_rr <- hal9001::fit_hal(
      X = fully_observed_data, Y = eif_rr,
      fit_control = list(cv_select = TRUE, use_min = TRUE, lambda.min.ratio = 1e-4),
      family = "gaussian"
    )
    projection_fit_or <- hal9001::fit_hal(
      X = fully_observed_data, Y = eif_or,
      fit_control = list(cv_select = TRUE, use_min = TRUE, lambda.min.ratio = 1e-4),
      family = "gaussian"
    )
    if_projections_1 <- predict(projection_fit_1, new_data = or_data[, fully_observed_cols])
    if_projections_0 <- predict(projection_fit_0, new_data = or_data[, fully_observed_cols])
    if_projections_rd <- predict(projection_fit_rd, new_data = or_data[, fully_observed_cols])
    if_projections_rr <- predict(projection_fit_rr, new_data = or_data[, fully_observed_cols])
    if_projections_or <- predict(projection_fit_or, new_data = or_data[, fully_observed_cols])
  }
  return(list("projection_1" = if_projections_1, "projection_0" = if_projections_0,
              "projection_rd" = if_projections_rd, "projection_rr" = if_projections_rr,
              "projection_or" = if_projections_or))
}

# get calibrated weights using survey
# @param data the dataset, with a column called "ip_weights"
# @param missing_indicator the name of the missing-data id column
# @param raking_formula what to rake on
# @param eif the estimated EIF
# @param ip_weights the initial inverse probability of missing-data weights
# @param iptw_weights IPTW weights (if necessary)
# @param tx_type how to adjust for treatment assignment ("adjust" or "iptw"/"aiptw")
# @param stratified is this a stratified design or Poisson sampling?
# @param strata_form the formula for stratification
get_calibrated_design <- function(data = NULL, missing_indicator = "Delta",
                                  raking_formula = NULL,
                                  eif = NULL, ip_weights = NULL,
                                  iptw_weights = rep(1, length(ip_weights)), tx_type = "adjust",
                                  stratified = FALSE, strata_form = "") {
  final_weights <- ip_weights
  if (length(iptw_weights) == sum(data[[missing_indicator]])) {
    the_iptw_weights <- iptw_weights
  } else {
    the_iptw_weights <- iptw_weights[data[[missing_indicator]] == 1]
  }
  final_weights[data[[missing_indicator]] == 1] <- ip_weights[data[[missing_indicator]] == 1] * the_iptw_weights
  data2 <- data
  data2$eif <- eif
  # data2$weights <- final_weights
  data2$weights <- ip_weights
  data2$ip_weights <- ip_weights
  if (stratified) {
    this_design <- survey::twophase(
      id = list(~ 1, ~ 1), subset = ~ I(data2[[missing_indicator]] == 1),
      probs = list(NULL, ~ I(1 / weights)),
      strata = list(NULL, as.formula(strata_form)),
      data = data2
    )
  } else {
    this_design <- survey::twophase(
      id = list(~ 1, ~ 1), subset = ~ I(data2[[missing_indicator]] == 1),
      probs = list(NULL, ~ I(1 / weights)),
      pps = list(NULL, survey::poisson_sampling(1 / data2$ip_weights[data2[[missing_indicator]] == 1])),
      data = data2
    )
  }

  calibrated_design <- survey::calibrate(this_design,
                                         formula = as.formula(raking_formula), phase = 2,
                                         calfun = "raking")
  # final weights, including IPTW
  calibrated_design2 <- calibrated_design
  calibrated_design2$prob <- calibrated_design$prob / the_iptw_weights
  return(calibrated_design2)
}

# get all calibrated weights using survey
# @param raking_variables a list of the raking variables
# @param rake_on_outcome, rake_on_exposure, rake_on_covars calibrate on outcome and/or covariates
# @param rescale should we rescale the weights to have mean 1?
# @param calibrate_separate should we calibrate the weights for each parameter separately?
# @param tx_type how we are adjusting for the propensity (ipw, aipw, adjusting)
# @return calibrated weights for treatment-specific means, RD, RR, and OR
# @param stratified is this a stratified design or Poisson sampling?
# @param strata_form the formula for stratification
get_all_calibrated_weights <- function(data = NULL, missing_indicator = "Delta",
                                       ip_weights = NULL,
                                       raking_variables = NULL,
                                       rake_on_outcome = TRUE,
                                       rake_on_exposure = TRUE,
                                       rake_on_covars = FALSE,
                                       rescale_miss = FALSE, family = "binomial",
                                       calibrate_separate = TRUE,
                                       conditional = FALSE,
                                       tx_type = "ipw",
                                       outcome_name = "Y",
                                       errorprone_outcome_name = "Y",
                                       exposure_name = "X",
                                       stratified = FALSE,
                                       strata_form = "") {
  if (is.character(family)) {
    family <- match.fun(family)()
  }
  raking_outcome_name <- outcome_name
  if (errorprone_outcome_name != outcome_name) {
    raking_outcome_name <- errorprone_outcome_name
  }
  if (isTRUE(calibrate_separate)) {
    # get the raking formula
    if (conditional) {
      raking_eifs <- paste0("~ ", paste0(paste0("inf", 1:ncol(raking_variables$eif)), collapse = " + "))
    } else {
      raking_eifs <- "~ eif"
    }
    raking_formula <- paste0(raking_eifs,
                             switch(as.numeric(rake_on_outcome) + 1, "", paste0(" + ", raking_outcome_name)),
                             switch(as.numeric(rake_on_exposure) + 1, "", paste0(" + ", exposure_name)),
                             switch(as.numeric(rake_on_covars) + 1, "", paste0(" +", unlist(strsplit(formula, "~")[[1]])[2])))
    # calibrate the weights for treatment-specific mean
    calibrated_design_1 <- get_calibrated_design(
      data = data, missing_indicator = missing_indicator, raking_formula = raking_formula,
      eif = raking_variables$eif_1, ip_weights = ip_weights, stratified = stratified,
      strata_form = strata_form
    )
    calibrated_weights_1 <- weights(calibrated_design_1)
    # calibrate the weights for control-specific mean
    calibrated_design_0 <- get_calibrated_design(
      data = data, missing_indicator = missing_indicator, raking_formula = raking_formula,
      eif = raking_variables$eif_0, ip_weights = ip_weights, stratified = stratified,
      strata_form = strata_form
    )
    calibrated_weights_0 <- weights(calibrated_design_0)
    # calibrate the weights for RD
    calibrated_design_rd <- get_calibrated_design(
      data = data, missing_indicator = missing_indicator, raking_formula = raking_formula,
      eif = raking_variables$eif_rd, ip_weights = ip_weights, stratified = stratified,
      strata_form = strata_form
    )
    calibrated_weights_rd <- weights(calibrated_design_rd)
    if (family$family != "gaussian") {
      # calibrate the weights for RR
      calibrated_design_rr <- get_calibrated_design(
        data = data, missing_indicator = missing_indicator, raking_formula = raking_formula,
        eif = raking_variables$eif_rr, ip_weights = ip_weights, stratified = stratified,
        strata_form = strata_form
      )
      calibrated_weights_rr <- weights(calibrated_design_rr)
      # calibrate the weights for OR
      calibrated_design_or <- get_calibrated_design(
        data = data, missing_indicator = missing_indicator, raking_formula = raking_formula,
        eif = raking_variables$eif_or, ip_weights = ip_weights, stratified = stratified,
        strata_form = strata_form
      )
      calibrated_weights_or <- weights(calibrated_design_or)
    } else {
      calibrated_weights_rr <- calibrated_weights_or <- calibrated_weights_rd
    }
  } else {
    # get the raking formula
    if (conditional) {
      raking_eifs <- paste0("~ ", paste0(paste0("inf", 1:ncol(raking_variables$eif)), collapse = " + "))
    } else {
      if (family$family == "gaussian") {
        raking_eifs <- c("0", "1", "rd")
      } else {
        # raking_eifs <- c("0", "1", "rd", "rr", "or")
        raking_eifs <- c("0", "1")
      }
      raking_eifs <- paste0("~ ", paste(paste0("eif_", raking_eifs), collapse = " + "))
    }
    raking_formula <- paste0(raking_eifs,
                             switch(as.numeric(rake_on_outcome) + 1, "", paste0(" + ", raking_outcome_name)),
                             switch(as.numeric(rake_on_exposure) + 1, "", paste0(" + ", exposure_name)),
                             switch(as.numeric(rake_on_covars) + 1, "", paste0(" +", unlist(strsplit(formula, "~")[[1]])[2])))
    if (conditional) {
      all_eifs <- raking_variables$eif
      data2 <- cbind.data.frame(data, all_eifs)
      names(data2)[(ncol(data) + 1):ncol(data2)] <- paste0("inf", 1:ncol(all_eifs))
    } else {
      # all_eifs <- do.call(cbind.data.frame, raking_variables[-1])
      all_eifs <- do.call(cbind.data.frame, raking_variables)
      data2 <- cbind.data.frame(data, all_eifs)
    }
    calibrated_design <- get_calibrated_design(
      data = data2, missing_indicator = missing_indicator, raking_formula = raking_formula,
      ip_weights = ip_weights, stratified = stratified,
      strata_form = strata_form
    )
    calibrated_weights <- weights(calibrated_design)
    calibrated_weights_1 <- calibrated_weights_0 <- calibrated_weights_rd <- calibrated_weights_rr <- calibrated_weights_or <- calibrated_weights
  }

  if (rescale_miss) {
    output_list <- list("tx_mean_1" = calibrated_weights_1 / mean(calibrated_weights_1),
                        "tx_mean_0" = calibrated_weights_0 / mean(calibrated_weights_0),
                        "rd" = calibrated_weights_rd / mean(calibrated_weights_rd),
                        "rr" = calibrated_weights_rr / mean(calibrated_weights_rr),
                        "or" = calibrated_weights_or / mean(calibrated_weights_or))
  } else {
    output_list <- list("tx_mean_1" = calibrated_weights_1, "tx_mean_0" = calibrated_weights_0,
                        "rd" = calibrated_weights_rd, "rr" = calibrated_weights_rr,
                        "or" = calibrated_weights_or)
  }
  return(output_list)
}

# marginalize a conditional estimator
# @param fit the fitted model object for estimating conditional quantities
# @param design the design object (for weights)
# @param coefficient_of_interest the variable corresponding to the coef of interest
# @param family the outcome family; if gaussian, only return the ATE. if not gaussian, return RD, RR, OR
marginalize_conditional_est <- function(fit = NULL, design = NULL, coefficient_of_interest = "A",
                                        family = "binomial") {
  these_weights <- weights(design)
  avg_preds <- avg_predictions(fit, variables = coefficient_of_interest,
                               wts = these_weights)
  rd_comparison <- avg_comparisons(fit, comparison = "differenceavg",
                                   variables = coefficient_of_interest,
                                   wts = these_weights)
  rr_comparison <- avg_comparisons(fit, comparison = "lnratioavg",
                                   variables = coefficient_of_interest,
                                   wts = these_weights)
  rr_comparison_nolog <- avg_comparisons(fit, comparison = "ratioavg",
                                   variables = coefficient_of_interest,
                                   wts = these_weights)
  or_comparison <- avg_comparisons(fit, comparison = "lnoravg",
                                   variables = coefficient_of_interest,
                                   wts = these_weights)
  results <- rbind(c(avg_preds[1, "estimate"], avg_preds[1, "std.error"]),
                   c(avg_preds[2, "estimate"], avg_preds[2, "std.error"]),
                   c(rd_comparison[1,"estimate"], rd_comparison[1,"std.error"]),
                   c(exp(rr_comparison[1,"estimate"]), rr_comparison[1,"std.error"]),
                   # c(rr_comparison_nolog[1, "estimate"], rr_comparison_nolog[1, "std.error"]),
                   c(exp(or_comparison[1,"estimate"]), or_comparison[1,"std.error"]))
  # return
  estimand <- c("mu0", "mu1", "ATE",
                # "logRR",
                "RR", "OR")

  results <- cbind.data.frame(estimand, results)
  colnames(results) <- c("estimand","est","SE")
  return(results)
}

# complete-case estimators, as a benchmark -------------------------------------
run_cc <- function(data = NULL,
                   outcome_formula = "Y ~ A + W",
                   tx_formula = "A ~ W",
                   coefficient_of_interest = "A",
                   outcome_name = "Y", tx_type = "adjust", use_meat = FALSE,
                   rescale_ps = FALSE, family = "binomial") {
  if (is.character(family)) {
    family <- switch(as.numeric(family == "gaussian"), gaussian(), binomial())
  }
  # estimate OR
  est_or <- get_outcome_regression(
    data = data, formula = outcome_formula, coefficient_of_interest = coefficient_of_interest,
    model = "parametric",  family = family, use_weights = FALSE
  )

  # estimate PS
  est_ps <- get_propensity_score(
    data = data, formula = tx_formula, variable = coefficient_of_interest,
    model = "parametric", use_weights = FALSE
  )

  # marginalize conditional estimator
  fit <- glm(as.formula(outcome_formula), data = data, family = family)
  avg_preds <- avg_predictions(fit, variables = coefficient_of_interest)
  rd_comparison <- avg_comparisons(fit, comparison = "differenceavg",
                                   variables = coefficient_of_interest)
  rr_comparison <- avg_comparisons(fit, comparison = "lnratioavg",
                                   variables = coefficient_of_interest)
  rr_comparison_nolog <- avg_comparisons(fit, comparison = "ratioavg",
                                         variables = coefficient_of_interest)
  results <- rbind(c(avg_preds[1, "estimate"], avg_preds[1, "std.error"]),
                   c(avg_preds[2, "estimate"], avg_preds[2, "std.error"]),
                   c(rd_comparison[1,"estimate"], rd_comparison[1,"std.error"]),
                   c(exp(rr_comparison[1,"estimate"]), rr_comparison[1,"std.error"]))
  # return
  estimand <- c("mu0", "mu1", "ATE", "RR")

  cond_results <- cbind.data.frame(estimand, results)
  colnames(cond_results) <- c("estimand","est","SE")
  cond_results$outcome_tx_estimator <- tx_type
  cond_results$miss_estimator <- "CGR"

  # use AIPTW
  mu_1 <- estimate_tx_specific_mean(
    a = data[[coefficient_of_interest]], a_level = 1, y = data[[outcome_name]],
    delta = rep(1, nrow(data)),
    outcome_regression = est_or$preds, outcome_regression_a = est_or$preds_1,
    prob_ps = est_ps$ps_1
  )
  mu_0 <- estimate_tx_specific_mean(
    a = data[[coefficient_of_interest]], a_level = 0, y = data[[outcome_name]],
    delta = rep(1, nrow(data)),
    outcome_regression = est_or$preds, outcome_regression_a = est_or$preds_0,
    prob_ps = est_ps$ps_0
  )
  all_eifs <- get_ideal_data_eifs(
    data = data[, 1:10], Q = est_or$preds, Q_1 = est_or$preds_1, Q_0 = est_or$preds_0,
    g_1 = est_ps$ps_1, g_0 = est_ps$ps_0, type = "ipw", psi_1 = mu_1, psi_0 = mu_0,
    coefficient_of_interest = coefficient_of_interest, outcome_name = outcome_name,
    glm_fit = est_or$fit, use_meat = use_meat
  )
  se_mu1 <- as.numeric(sqrt(crossprod(matrix(all_eifs$if_1))) / (nrow(data)))
  se_mu0 <- as.numeric(sqrt(crossprod(matrix(all_eifs$if_0))) / (nrow(data)))
  se_ate <- as.numeric(sqrt(crossprod(matrix(all_eifs$if_rd))) / (nrow(data)))
  se_rr <- as.numeric(sqrt(crossprod(matrix(all_eifs$if_rr))) / (nrow(data)))

  marg_result <- tibble::tibble(
    "estimand" = c("mu1", "mu0", "ATE", "RR"),
    "outcome_tx_estimator" = "aipw",
    "miss_estimator" = "MGR",
    "est" = c(mu_1, mu_0, mu_1 - mu_0, mu_1 / mu_0),
    "SE" = c(se_mu1, se_mu0, se_ate, se_rr)
  )
  return(dplyr::bind_rows(cond_results, marg_result))
}

# GR for conditional parameter, then marginalizing -----------------------------

# run raking for a conditional parameter, then marginalize it
# @param data the dataset
# @param outcome_formula the outcome regression formula
# @param tx_formula the propensity score formula (only used if tx_type != "adjust")
# @param miss_formula the missing-data model formula
# @param n_imp the number of imputations for the raking variable
# @param coefficient_of_interest the coefficient of interest, either a character string (e.g., "X") or its position in a formula (e.g., 2, behind the intercept)
# @param missing_indicator the name of the missingness indicator
# @param rake_on_outcome rake on the outcome in addition to EIFs?
# @param rake_on_exposure rake on the exposure variable in addition to EIFs?
# @param rake_on_covars rake on the covariates in addition to EIFs?
# @param rake_on_eifs rake on the EIFs?
# @param procedure the name of the procedure, e.g., "GR_conditional"
# @param tx_type use adjustment ("adjust"), IPW ("ipw"), or AIPW ("aipw") to handle the propensity score?
# @param use_meat logical, use meat for EIF calculation?
# @param rescale_ps logical, stabilize IPTW weights (for treatment propensity)?
# @param stratified is this a stratified design or Poisson sampling?
# @param strata_form the formula for stratification
# @return a list with the parameter estimate and standard error (in a data.frame) and the fitted model
run_raking_conditional <- function(data = NULL,
                                   outcome_formula = "Y ~ A + W",
                                   tx_formula = "A ~ W",
                                   miss_formula = "Delta ~ A + Y",
                                   phase_1_formula = NULL,
                                   coefficient_of_interest = "A",
                                   outcome_name = "Y",
                                   errorprone_outcome_name = "Y",
                                   n_imp = 10,
                                   missing_indicator = "Delta",
                                   rake_on_outcome = FALSE,
                                   rake_on_exposure = FALSE,
                                   rake_on_covars = FALSE,
                                   rake_on_eifs = TRUE,
                                   family = "binomial",
                                   procedure = NULL,
                                   tx_type = "adjust",
                                   use_meat = FALSE,
                                   rescale_ps = TRUE,
                                   stratified = FALSE,
                                   strata_form = "") {
  if (is.character(family)) {
    family <- match.fun(family)()
  }
  # get number of coefficients in outcome model
  init_fit <- glm(formula = as.formula(outcome_formula), family = family, data = data)
  num_coefs <- length(coefficients(init_fit))

  # get initial inverse probability of missing weights
  ps_list <- get_propensity_score(data = data, formula = miss_formula,
                                    variable = missing_indicator, model = "parametric",
                                    use_weights = FALSE)
  ip_weights <- ps_list$ip_weights
  # if using IPW for the tx propensity score, fit it
  if (!grepl("adjust", tx_type)) {
    tx_ps_list <- get_propensity_score(
      data = data, formula = tx_formula, variable = coefficient_of_interest,
      model = "parametric", use_weights = TRUE, weights = ip_weights
    )
    tx_ps_weights <- rep(NA, length(tx_ps_list$ip_weights))
    p_1 <- mean(data[[coefficient_of_interest]] == 1)
    p_0 <- mean(data[[coefficient_of_interest]] == 0)
    # g_1 <- switch(
    #   as.numeric(rescale_ps) + 1, tx_ps_list$ps_1, tx_ps_list$ps_1 / p_1
    # )
    # g_0 <- switch(
    #   as.numeric(rescale_ps) + 1, tx_ps_list$ps_0, tx_ps_list$ps_0 / p_0
    # )
    g_1 <- tx_ps_list$ps_1
    g_0 <- tx_ps_list$ps_0
    tx_ps_weights[data[data[[missing_indicator]] == 1, ][[coefficient_of_interest]] == 1] <- 1 / g_1[data[data[[missing_indicator]] == 1, ][[coefficient_of_interest]] == 1]
    tx_ps_weights[data[data[[missing_indicator]] == 1, ][[coefficient_of_interest]] == 0] <- 1 / g_0[data[data[[missing_indicator]] == 1, ][[coefficient_of_interest]] == 0]
  } else {
    tx_ps_weights <- rep(1, sum(complete.cases(data)))
  }
  # the full-data EIF, using MI
  obs_data_eif <- estimate_eif_projection_mi(outcome_formula = outcome_formula,
                                             tx_formula = tx_formula,
                                             phase_1_formula = phase_1_formula,
                                             coefficient_of_interest = coefficient_of_interest,
                                             data = data,
                                             missing_indicator = missing_indicator,
                                             n_imp = n_imp, conditional = TRUE,
                                             family = family,
                                             tx_type = tx_type, use_meat = use_meat,
                                             rescale_ps = rescale_ps)
  inf_mat <- obs_data_eif$eif
  # always rake at least on the EIF for treatment (add 1 for intercept)
  i <- grep(coefficient_of_interest, strsplit(outcome_formula, "+", fixed = TRUE)[[1]]) + 1
  raking_formula <- paste0("~ inf", i)
  data$inf <- inf_mat[, i]
  names(data)[names(data) == "inf"] <- paste0("inf", i)
  # run raking for the conditional parameter; rake on all influence functions if requested
  if (isTRUE(rake_on_eifs)) {
    raking_formula <- "~ inf1"
    for (i in seq_len(num_coefs)) {
      varname <- paste0("inf", i)
      data$inf <- inf_mat[, i]
      names(data)[names(data) == "inf"] <- varname
      if (i > 1) {
        raking_formula <- paste0(raking_formula, " + ", varname)
      }
    }
  }
  raking_outcome_name <- outcome_name
  if (errorprone_outcome_name != outcome_name) {
    raking_outcome_name <- errorprone_outcome_name
  }
  raking_formula <- paste0(raking_formula,
                           switch(as.numeric(rake_on_outcome) + 1, "", paste0(" + ", raking_outcome_name)),
                           switch(as.numeric(rake_on_exposure) + 1, "", paste0(" + ", coefficient_of_interest)),
                           switch(as.numeric(rake_on_covars) + 1, "", paste0(" +", unlist(strsplit(formula, "~")[[1]])[2])))
  # now handle IPTW, if we have done it
  if (!grepl("adjust", tx_type)) {
    raking_formula <- paste0(raking_formula, paste0(" +", unlist(strsplit(phase_1_formula, "~")[[1]])[2]))
  }
  # final weights include IPTW if that is what we have done. can include PS variables but only phase 1
  calibrated_design <- get_calibrated_design(
    data = data, missing_indicator = missing_indicator, raking_formula = raking_formula,
    ip_weights = ip_weights, iptw_weights = tx_ps_weights, tx_type = tx_type, stratified = stratified,
    strata_form = strata_form
  )
  # another way to get IPTW: do the projection for the propensity score explicitly/manually
  # tx_mat <- model.matrix(as.formula(tx_formula), data = data)
  # xz <- (data[[coefficient_of_interest]] - mean(data[[coefficient_of_interest]]))[data[[missing_indicator]] == 1] * tx_mat * ip_weights[data[[missing_indicator]] == 1]
  # qrxz <- qr(xz)
  # data2 <- data
  # data2$weights <- weights(calibrated_design)
  # or_fit <- glm(formula = as.formula(outcome_formula), family = family, data = data,
  #               weights = weights)
  # U <-
  # V <- qr.resid(qrxz, data[[outcome_name]][data[[missing_indicator]] == 1] - xz %*% solve(t(xz) %*% xz) %*% t(xz) %*% as.matrix(data[[outcome_name]][data[[missing_indicator]] == 1]))
  # svytotal(V, calibrated_design)
  # need to get the weights now by multiplying by IPW weights for tx propensity
  # if (!grepl("adjust", tx_type)) {
  #   these_weights <- weights(calibrated_design) * tx_ps_weights
  #   calibrated_design$prob <- 1 / these_weights
  # }
  prob_obs_long <- ps_list$ps_1
  prob_obs_long[data[[missing_indicator]] == 1] <- calibrated_design$prob
  if (family$family == "gaussian") {
    fit <- survey::svyglm(outcome_formula, design = calibrated_design, family = gaussian)
  } else if (family$family == "binomial") {
    if (family$link == "logit") {
      fit <- survey::svyglm(outcome_formula, design = calibrated_design, family = quasibinomial)
    } else {
      fit <- survey::svyglm(outcome_formula, design = calibrated_design, family = quasibinomial(link = "probit"))
    }
  } else if (family$family == "poisson") {
    fit <- survey::svyglm(outcome_formula, design = calibrated_design, family = family)
  } else {
    stop("Currently only implemented for gaussian and binomial families.")
  }

  # marginalize
  results <- marginalize_conditional_est(fit = fit, design = calibrated_design,
                                         coefficient_of_interest = coefficient_of_interest,
                                         family = family)
  if (is.null(procedure)) {
    procedure <- "GR_conditional"
  }
  return(list(
    "result" = dplyr::bind_cols(tibble("miss_procedure" = procedure, "outcome_tx_procedure" = tx_type), as_tibble(results)),
    "calibrated_weights" = prob_obs_long
  ))
}

# scale a set of weights by its mean
# @param weights the weights to scale
# @return the scaled weights
scale_weights <- function(weights) {
  return(weights / mean(weights))
}

# get either an AIPW, IPW, or g-computation estimator of the treatment-specific mean
# note that we can do AIPW or IPW for the missing data and AIPW, IPW, g-comp for the outcome/treatment models
# @param prob_obs the probability of being observed (used as inverse weight)
# @param prob_ps the estimated propensity score (treatment) (used as an inverse weight)
# @param outcome_regression the estimated outcome regression predictions
# @param eif_projection the projection of the full-data EIF onto the always-observed variables
# @param type whether to do aipw, ipw, or g-computation for outcome/tx models
# @param type_miss whether to do aipw or ipw for missing-data model
# @param rescale whether or not to rescale the weights for point estimation
estimate_tx_specific_mean <- function(a = NULL, a_level = 1, y = NULL, delta = NULL,
                                      outcome_regression = NULL,
                                      outcome_regression_a = NULL,
                                      prob_obs = rep(1, length(outcome_regression)),
                                      prob_ps = rep(1, length(outcome_regression)),
                                      # eif_projection = rep(1, length(outcome_regression)),
                                      eif = rep(0, length(outcome_regression)),
                                      type = "aipw", type_miss = "aipw",
                                      rescale_miss = TRUE,
                                      rescale_ps = FALSE) {
  ps_weight <- 1 / prob_ps
  if (rescale_ps) {
    # ps_weight <- 1 / prob_ps / mean(1 / prob_ps)
    ps_weight <- mean(a == a_level) / prob_ps
  }
  if (type == "g_comp") {
    outcome_treatment_contribution <- outcome_regression_a
  } else if (type == "ipw") {
    outcome_treatment_contribution <- as.numeric(a == a_level) * ps_weight * y
  } else {
    outcome_treatment_contribution <- as.numeric(a == a_level) * ps_weight * (y - outcome_regression) + outcome_regression_a
  }
  obs_weight <- 1 / prob_obs
  if (rescale_miss) {
    obs_weight <- (1 / prob_obs) / mean(1 / prob_obs)
  }
  if (length(obs_weight) == sum(delta == 1)) {
    final_weight <- obs_weight
  } else {
    final_weight <- obs_weight[delta == 1]
  }
  psi_contribution <- rep(0, length(delta))
  psi_contribution[delta == 1] <- final_weight * outcome_treatment_contribution
  if (rescale_miss) {
    est_ipw <- sum(psi_contribution) / sum(final_weight)
  } else {
    est_ipw <- mean(psi_contribution)
  }
  # est_ipw <- sum(
  #   ifelse(delta == 0, 0, final_weight * outcome_treatment_contribution)
  # ) / sum(final_weight)

  est <- est_ipw
  if (type_miss == "aipw") {
    est <- est_ipw + mean(eif)
  }
  return(est)
}

# GR for marginal parameter ----------------------------------------------------
# get nonparametric, parametric ideal-data EIFS
# @param data the dataset used for fitting the outcome regression model. complete cases
# @param Q predictions from outcome regression model (A as observed)
# @param Q_1 predictions from outcome regression model (A = 1)
# @param Q_0 predictions from outcome regression model (A = 0)
# @param g_1 propensity score for treatment (A = 1)
# @param g_0 propensity score for not treatment (A = 0)
# @param type type of outcome/tx propensity estimator
# @param psi_1 an initial estimator of outcome mean under A = 1
# @param psi_0 an initial estimator of outcome mean under A = 0
# @param coefficient_of_interest character of treatment name
# @param outcome_name the name of the outcome variable
# @param glm_fit the working glm outcome regression model fit
# @return ideal-data EIFs
get_ideal_data_eifs <- function(data = NULL,
                                Q = rep(1, nrow(data)),
                                Q_1 = rep(1, nrow(data)),
                                Q_0 = rep(1, nrow(data)),
                                g_1 = rep(1, nrow(data)),
                                g_0 = 1 - g_1,
                                type = "ipw",
                                psi_1 = 0,
                                psi_0 = 0,
                                coefficient_of_interest = "A",
                                outcome_name = "Y",
                                glm_fit = NULL, use_meat = TRUE) {
  dat_cc <- data[complete.cases(data), ]
  # estimate the nonparametric influence functions
  np_if_mu_1 <- if_tx_mean_nonparametric(
    a = dat_cc[[coefficient_of_interest]], a_level = 1, dat_cc[[outcome_name]],
    Q = Q, Q_a = Q_1, g = g_1, type = type, psi = psi_1
  )
  np_if_mu_0 <- if_tx_mean_nonparametric(
    a = dat_cc[[coefficient_of_interest]], a_level = 0, dat_cc[[outcome_name]],
    Q = Q, Q_a = Q_0, g = g_0, type = type, psi = psi_0
  )
  # estimate the parametric influence functions
  par_if_mu_1 <- eif_tx_mean_parametric(
    nonparametric_eif = np_if_mu_1, glm_fit = glm_fit,
    use_meat = use_meat
  )
  par_if_mu_0 <- eif_tx_mean_parametric(
    nonparametric_eif = np_if_mu_0, glm_fit = glm_fit,
    use_meat = use_meat
  )
  par_if_rd <- full_data_eif(mu1 = psi_1, mu0 = psi_0, eif1 = par_if_mu_1,
                             eif0 = par_if_mu_0, estimand = "RD")
  par_if_rr <- full_data_eif(mu1 = psi_1, mu0 = psi_0, eif1 = par_if_mu_1,
                             eif0 = par_if_mu_0, estimand = "RR")
  par_if_rr_nolog <- full_data_eif(mu1 = psi_1, mu0 = psi_0, eif1 = par_if_mu_1,
                             eif0 = par_if_mu_0, estimand = "nologRR")
  par_if_or <- full_data_eif(mu1 = psi_1, mu0 = psi_0, eif1 = par_if_mu_1,
                             eif0 = par_if_mu_0, estimand = "OR")
  return(list("if_1" = par_if_mu_1, "if_0" = par_if_mu_0,
              "if_rd" = par_if_rd, "if_rr" = par_if_rr, "if_or" = par_if_or,
              "np_if_1" = np_if_mu_1, "np_if_0" = np_if_mu_0,
              "if_rr_nolog" = par_if_rr_nolog))
}

# get nonparametric, parametric ideal-data EIFs and observed-data EIFs
# @param data the dataset used for fitting the outcome regression model. complete cases
# @param fully_observed_data data (complete cases) on columns that are always observed
# @param Q predictions from outcome regression model (A as observed)
# @param Q_1 predictions from outcome regression model (A = 1)
# @param Q_0 predictions from outcome regression model (A = 0)
# @param g_1 propensity score for treatment (A = 1)
# @param g_0 propensity score for not treatment (A = 0)
# @param delta indicator of being fully observed (1 = observed)
# @param prob_obs the probability of being fully observed
# @param type type of outcome/tx propensity estimator
# @param psi_1 an initial estimator of outcome mean under A = 1
# @param psi_0 an initial estimator of outcome mean under A = 0
# @param coefficient_of_interest character of treatment name
# @param outcome_name the name of the outcome variable
# @param glm_fit the working glm outcome regression model fit
# @param or_data dataset (complete cases) for the outcome regression
# @param fully_observed_cols the columns that are always observed
# @param outcome_formula the outcome regression formula
# @param tx_formula the treatment-assignment formula
# @param phase_1_formula the variables fully observed in phase 1
# @param use_meat logical, if true uses sandwich::meat (outer product of empirical estimating functions) rather than bread (empirical Hessian)
# @param use_imputation logical, if true use MI to estimate EIF projection
# @param use_hal logical, use HAL to estimate EIF projection
# @param use_basis logical, use basis expansion to estimate EIF projection
# @param use_projection_intercept logical, use an intercept in estimating EIF projection
# @return all IFs for treatment-specific means, ATE, RR, OR
get_all_eifs <- function(data = NULL,
                         fully_observed_data = NULL,
                         Q = rep(1, nrow(data)),
                         Q_1 = rep(1, nrow(data)),
                         Q_0 = rep(1, nrow(data)),
                         g_1 = rep(1, nrow(data)),
                         g_0 = 1 - g_1,
                         delta = rep(1, nrow(data)),
                         prob_obs = rep(1, nrow(data)),
                         type = "ipw",
                         psi_1 = 0,
                         psi_0 = 0,
                         coefficient_of_interest = "A",
                         missing_indicator = "Delta",
                         outcome_name = "Y",
                         family = gaussian(),
                         glm_fit = NULL,
                         or_data = NULL,
                         fully_observed_cols = NULL,
                         outcome_formula = "Y ~ X",
                         tx_formula = "X ~ Z",
                         phase_1_formula = NULL,
                         use_meat = TRUE,
                         use_imputation = FALSE,
                         use_hal = FALSE,
                         use_basis = FALSE,
                         use_projection_intercept = FALSE,
                         rescale_miss = FALSE, rescale_ps = FALSE, n_imp = 10) {
  dat_cc <- data[complete.cases(data), ]
  # get full-data EIFs
  full_data_eifs <- get_ideal_data_eifs(data = data, Q = Q, Q_1 = Q_1, Q_0 = Q_0,
                                        g_1 = g_1, g_0 = g_0, type = type,
                                        psi_1 = psi_1, psi_0 = psi_0,
                                        coefficient_of_interest = coefficient_of_interest,
                                        outcome_name = outcome_name, glm_fit = glm_fit,
                                        use_meat = use_meat)
  par_if_mu_1 <- full_data_eifs$if_1
  par_if_mu_0 <- full_data_eifs$if_0
  par_if_rd <- full_data_eifs$if_rd
  par_if_rr <- full_data_eifs$if_rr
  par_if_or <- full_data_eifs$if_or
  par_if_rr_nolog <- full_data_eifs$if_rr_nolog
  # estimate the IF projections
  if (use_imputation) {
    # estimate using MI
    mi_eifs <- estimate_eif_projection_mi(
      outcome_formula = outcome_formula, tx_formula = tx_formula,
      phase_1_formula = phase_1_formula,
      data = data, missing_indicator = missing_indicator,
      coefficient_of_interest = coefficient_of_interest,
      outcome_name = outcome_name, n_imp = n_imp, conditional = FALSE,
      use_sl = FALSE, model = "parametric", sl_args = list(),
      family = family, miss_prob = prob_obs, tx_type = type, rescale_miss = rescale_miss, rescale_ps = rescale_ps,
      use_meat = use_meat
    )
    if_projections_1 <- mi_eifs$eif_1
    if_projections_0 <- mi_eifs$eif_0
    if_projections_rd <- mi_eifs$eif_rd
    if_projections_rr <- mi_eifs$eif_rr
    if_projections_rr_nolog <- mi_eifs$eif_rr_nolog
    if_projections_or <- mi_eifs$eif_or
  } else {
    reg_eifs <- estimate_eif_projection_regression(
      eif_1 = par_if_mu_1, eif_0 = par_if_mu_0, eif_rd = par_if_rd,
      eif_rr = par_if_rr, eif_or = par_if_or, fully_observed_data = fully_observed_data,
      fully_observed_cols = fully_observed_cols, or_data = or_data, use_hal = use_hal,
      use_basis = use_basis, use_projection_intercept = use_projection_intercept
    )
    if_projections_1 <- reg_eifs$projection_1
    if_projections_0 <- reg_eifs$projection_0
    if_projections_rd <- reg_eifs$projection_rd
    if_projections_rr <- reg_eifs$projection_rr
    if_projections_rr_nolog <- rep(NA, length(reg_eifs$projection_rr))
    if_projections_or <- reg_eifs$projection_or
  }
  # estimate the observed-data influence functions
  observed_if_mu_1 <- observed_data_eif(
    delta = delta, prob_obs = prob_obs,
    full_eif = par_if_mu_1, eif_projection = if_projections_1
  )
  observed_if_mu_0 <- observed_data_eif(
    delta = delta, prob_obs = prob_obs,
    full_eif = par_if_mu_0, eif_projection = if_projections_0
  )
  observed_if_rd <- observed_data_eif(
    delta = delta, prob_obs = prob_obs,
    full_eif = par_if_rd, eif_projection = if_projections_rd
  )
  observed_if_rr <- observed_data_eif(
    delta = delta, prob_obs = prob_obs,
    full_eif = par_if_rr, eif_projection = if_projections_rr
  )
  observed_if_or <- observed_data_eif(
    delta = delta, prob_obs = prob_obs,
    full_eif = par_if_or, eif_projection = if_projections_or
  )
  observed_if_rr_nolog <- observed_data_eif(
    delta = delta, prob_obs = prob_obs,
    full_eif = par_if_rr_nolog, eif_projection = if_projections_rr_nolog
  )
  # return
  return(list("obs_if_ate" = observed_if_rd, "obs_if_rr" = observed_if_rr, "obs_if_or" = observed_if_or,
              "obs_if_mu_1" = observed_if_mu_1, "obs_if_mu_0" = observed_if_mu_0,
              "np_if_mu_1" = full_data_eifs$np_if_1, "np_if_mu_0" = full_data_eifs$np_if_0,
              "par_if_mu_1" = par_if_mu_1, "par_if_mu_0" = par_if_mu_0,
              "par_if_rd" = par_if_rd, "par_if_rr" = par_if_rr, "par_if_or" = par_if_or,
              "projection_1" = if_projections_1, "projection_0" = if_projections_0,
              "projection_rd" = if_projections_rd, "projection_rr" = if_projections_rr,
              "projection_or" = if_projections_or,
              "obs_if_rr_nolog" = observed_if_rr_nolog, "par_if_rr_nolog" = par_if_rr_nolog,
              "projection_rr_nolog" = if_projections_rr_nolog
  ))
}
# general function to run marginal estimation, either IPW or AIPW for tx propensity
# and either IPW or AIPW for missing-data model
# @param data the dataset
# @param formula the outcome regression model
# @param miss_formula the missing-data formula
# @param tx_formula the tx propensity formula
# @param phase_1_formula any additional fully-observed phase 1 variables that should be
#         included in the EIF projection model
# @param outcome_name the outcome variable
# @param family the outcome regression family (e.g., gaussian, binomial)
# @param coefficient_of_interest the treatment/exposure variable
# @param missing_indicator the missing-data indicator
# @param model whether to use a parametric model (e.g., for EIFs)
# @param use_projection_intercept logical, use an intercept in projection model
# @param use_basis logical, use basis expansion for projection model
# @param use_hal logical, use highly adaptive lasso for projection model
# @param use_imputation should we use regression (FALSE) or imputation (TRUE) to estimate the optimal raking variable?
# @param rescale logical, rescale to have mean 1
# @param use_weights logical, use IPCW weights in fitting initial tx propensity and outcome model
# @param use_meat logical, use meat function from sandwich
# @param outcome_tx_estimator "ipw" or "aipw", how to handle the treatment propensity and outcome models
# @param miss_estimator "ipw", "aipw", or "raking", how to handle the missing-data model
# @param n_imp number of imputations
# @param sl_args arguments to pass to super learner
# @param rake_on_outcome should we rake on the outcome too?
# @param calibrate_separate should we calibrate weights based on EIFs for ATE,RR, etc.?
# @param stratified is this a stratified design or Poisson sampling?
# @param strata_form the formula for stratification
# @return a tibble with the point estimate, SE, confidence interval
run_marginal <- function(data = NULL,
                         formula = "Y ~ A + W",
                         miss_formula = "Delta ~ A + Y",
                         tx_formula = "A ~ W",
                         phase_1_formula = "~ A + W",
                         outcome_name = "Y",
                         errorprone_outcome_name = "Y",
                         family = "binomial",
                         coefficient_of_interest = "A",
                         missing_indicator = "Delta",
                         model = "parametric",
                         use_projection_intercept = TRUE,
                         use_basis = FALSE,
                         use_hal = FALSE,
                         use_imputation = FALSE,
                         rescale_miss = TRUE,
                         rescale_ps = FALSE,
                         use_weights = TRUE,
                         use_meat = TRUE,
                         outcome_tx_estimator = "aipw",
                         miss_estimator = "aipw",
                         n_imp = 10,
                         sl_args = list(),
                         rake_on_outcome = TRUE,
                         calibrate_separate = FALSE,
                         stratified = FALSE,
                         strata_form = "") {
  if (is.character(family)) {
    family <- match.fun(family)()
  }
  internal_miss_estimator <- ifelse(
    grepl("MR", miss_estimator) | grepl("raking", miss_estimator) | grepl("MGR", miss_estimator), "ipw",
    miss_estimator
  )
  # get the probability of being fully observed
  if (miss_formula == "" | is.null(miss_formula)) {
    prob_obs <- rep(1, nrow(data))
    delta <- rep(1, nrow(data))
  } else {
    obs_list <- get_propensity_score(
      data = data, formula = miss_formula,
      model = model, variable = missing_indicator,
      use_weights = FALSE
    )
    prob_obs <- obs_list$ps_1
    delta <- data[[missing_indicator]]
  }
  # get complete cases and fully-observed phase 1 columns
  phase_1_terms <- terms.formula(as.formula(phase_1_formula))
  phase_1_varnames <- gsub(" ", "", rownames(attr(phase_1_terms, "factors")))
  ps_terms <- terms.formula(as.formula(tx_formula))
  ps_varnames <- gsub(" ", "", rownames(attr(ps_terms, "factors")))
  or_terms <- terms.formula(as.formula(formula))
  or_varnames <- gsub(" ", "", rownames(attr(or_terms, "factors")))
  or_plus_phase_1_data <- data[, unique(c(or_varnames, ps_varnames, phase_1_varnames))]
  cc <- complete.cases(or_plus_phase_1_data)
  fully_observed_cols <- (colSums(is.na(or_plus_phase_1_data)) == 0)
  fully_observed_data <- or_plus_phase_1_data[complete.cases(or_plus_phase_1_data), fully_observed_cols]
  dat_cc <- or_plus_phase_1_data[cc, ]
  dat0 <- dat1 <- dat_cc
  dat0[[coefficient_of_interest]] <- 0
  dat1[[coefficient_of_interest]] <- 1

  # get initial coarsening weights
  initial_weights <- 1 / prob_obs[cc]
  initial_weights_stab <- initial_weights
  if (rescale_miss) {
    initial_weights_stab <- initial_weights / mean(initial_weights)
  }
  # get the probability of being assigned treatment
  ps_list <- get_propensity_score(
    data = dat_cc[, ps_varnames], formula = tx_formula,
    model = model, variable = coefficient_of_interest, use_weights = use_weights,
    weights = initial_weights_stab
  )
  ps_1 <- ps_list$ps_1
  ps_0 <- ps_list$ps_0
  # get the outcome regression
  or_list <- get_outcome_regression(
    data = dat_cc[, or_varnames], formula = formula, coefficient_of_interest = coefficient_of_interest,
    model = model, family = family, use_weights = use_weights,
    weights = initial_weights_stab
  )
  preds <- or_list$preds
  preds_0 <- or_list$preds_0
  preds_1 <- or_list$preds_1
  # get an initial estimator of the treatment-specific means, without IF projection
  initial_mu_1 <- estimate_tx_specific_mean(
    a = dat_cc[[coefficient_of_interest]], a_level = 1, y = dat_cc[[outcome_name]],
    delta = delta, outcome_regression = preds, outcome_regression_a = preds_1,
    prob_obs = prob_obs, prob_ps = ps_1,
    eif = rep(0, length(delta)),
    type = outcome_tx_estimator, type_miss = internal_miss_estimator,
    rescale_miss = rescale_miss, rescale_ps = rescale_ps
  )
  initial_mu_0 <- estimate_tx_specific_mean(
    a = dat_cc[[coefficient_of_interest]], a_level = 0, y = dat_cc[[outcome_name]],
    delta = delta, outcome_regression = preds, outcome_regression_a = preds_0,
    prob_obs = prob_obs, prob_ps = ps_0,
    eif = rep(0, length(delta)),
    type = outcome_tx_estimator, type_miss = internal_miss_estimator,
    rescale_miss = rescale_miss, rescale_ps = rescale_ps
  )
  # estimate the ideal-data influence functions on the fully-observed data
  all_eifs <- get_all_eifs(data = data, fully_observed_data = fully_observed_data,
                           Q = preds, Q_1 = preds_1,
                           Q_0 = preds_0, g_1 = ps_1, g_0 = ps_0,
                           delta = delta, prob_obs = prob_obs,
                           type = outcome_tx_estimator,
                           psi_1 = initial_mu_1, psi_0 = initial_mu_0,
                           coefficient_of_interest = coefficient_of_interest,
                           missing_indicator = missing_indicator, family = family,
                           outcome_name = outcome_name, glm_fit = or_list$fit,
                           or_data = or_plus_phase_1_data,
                           fully_observed_cols = fully_observed_cols,
                           outcome_formula = formula, tx_formula = tx_formula,
                           phase_1_formula = phase_1_formula,
                           use_meat = use_meat, use_imputation = use_imputation,
                           use_hal = use_hal, use_basis = use_basis,
                           use_projection_intercept = use_projection_intercept,
                           rescale_miss = rescale_miss, rescale_ps = rescale_ps, n_imp = n_imp)
  observed_if_mu_1 <- all_eifs$obs_if_mu_1
  observed_if_mu_0 <- all_eifs$obs_if_mu_0

  # adjust the weights, if using MR
  if (miss_estimator == "raking" | miss_estimator == "MR" | miss_estimator == "CR" |
      miss_estimator == "MGR" | miss_estimator == "CGR") {
    all_calibrated_weights <- get_all_calibrated_weights(
      data = data, missing_indicator = missing_indicator, ip_weights = 1 / prob_obs,
      raking_variables = list("eif_1" = all_eifs$projection_1,
                              "eif_0" = all_eifs$projection_0,
                              "eif_rd" = all_eifs$projection_rd,
                              "eif_rr" = all_eifs$projection_rr,
                              "eif_or" = all_eifs$projection_or),
      rake_on_outcome = rake_on_outcome, conditional = FALSE,
      rake_on_exposure = FALSE, rake_on_covars = FALSE,
      calibrate_separate = calibrate_separate, rescale_miss = FALSE, tx_type = "ipw",
      family = family, outcome_name = outcome_name,
      errorprone_outcome_name = errorprone_outcome_name,
      exposure_name = coefficient_of_interest, stratified = stratified,
      strata_form = strata_form
    )
    final_prob_obs_1 <- 1 / all_calibrated_weights$tx_mean_1
    final_prob_obs_0 <- 1 / all_calibrated_weights$tx_mean_0
  } else {
    final_prob_obs_1 <- final_prob_obs_0 <- prob_obs
  }

  # estimate the treatment-specific means using the final weights
  mu_1 <- estimate_tx_specific_mean(
    a = dat_cc[[coefficient_of_interest]], a_level = 1, y = dat_cc[[outcome_name]],
    delta = delta, outcome_regression = preds, outcome_regression_a = preds_1,
    prob_obs = final_prob_obs_1, prob_ps = ps_1,
    eif = observed_if_mu_1,
    type = outcome_tx_estimator, type_miss = internal_miss_estimator,
    rescale_miss = rescale_miss, rescale_ps = rescale_ps
  )
  mu_0 <- estimate_tx_specific_mean(
    a = dat_cc[[coefficient_of_interest]], a_level = 0, y = dat_cc[[outcome_name]],
    delta = delta, outcome_regression = preds, outcome_regression_a = preds_0,
    prob_obs = final_prob_obs_0, prob_ps = ps_0,
    eif = observed_if_mu_0,
    type = outcome_tx_estimator, type_miss = internal_miss_estimator,
    rescale_miss = rescale_miss, rescale_ps = rescale_ps
  )
  ate <- mu_1 - mu_0
  rr <- mu_1 / mu_0
  or <- (mu_1 / (1 - mu_1)) / (mu_0 / (1 - mu_0))
  # update the EIFs using the final weights
  if (miss_estimator == "raking" | miss_estimator == "CR" | miss_estimator == "MR" |
      miss_estimator == "MGR" | miss_estimator == "CGR") {
    final_ps_list <- get_propensity_score(
      data = dat_cc[, ps_varnames], formula = tx_formula,
      model = model, variable = coefficient_of_interest, use_weights = use_weights,
      weights = final_prob_obs_1
    )
    final_or_list <-  get_outcome_regression(
      data = dat_cc[, or_varnames], formula = formula, coefficient_of_interest = coefficient_of_interest,
      model = model, family = family, use_weights = use_weights,
      weights = final_prob_obs_1
    )
    prob_obs_long <- prob_obs
    prob_obs_long[delta == 1] <- final_prob_obs_1
    final_eifs <- all_eifs
    final_eifs$obs_if_mu_1 <- observed_data_eif(
      delta = delta, prob_obs = prob_obs_long, full_eif = all_eifs$par_if_mu_1,
      eif_projection = all_eifs$projection_1
    )
    final_eifs$obs_if_mu_0 <- observed_data_eif(
      delta = delta, prob_obs = prob_obs_long, full_eif = all_eifs$par_if_mu_0,
      eif_projection = all_eifs$projection_0
    )
    final_eifs$obs_if_ate <- observed_data_eif(
      delta = delta, prob_obs = prob_obs_long, full_eif = all_eifs$par_if_rd,
      eif_projection = all_eifs$projection_rd
    )
    final_eifs$obs_if_rr <- observed_data_eif(
      delta = delta, prob_obs = prob_obs_long, full_eif = all_eifs$par_if_rr,
      eif_projection = all_eifs$projection_rr
    )
    final_eifs$obs_if_rr_nolog <- observed_data_eif(
      delta = delta, prob_obs = prob_obs_long, full_eif = all_eifs$par_if_rr_nolog,
      eif_projection = all_eifs$projection_rr_nolog
    )
    final_eifs$obs_if_or <- observed_data_eif(
      delta = delta, prob_obs = prob_obs_long, full_eif = all_eifs$par_if_or,
      eif_projection = all_eifs$projection_or
    )
  } else {
    prob_obs_long <- prob_obs
    final_eifs <- all_eifs
  }
  final_observed_if_mu_1 <- final_eifs$obs_if_mu_1
  final_observed_if_mu_0 <- final_eifs$obs_if_mu_0
  final_observed_if_ate <- final_eifs$obs_if_ate
  final_observed_if_rr <- final_eifs$obs_if_rr
  final_observed_if_rr_nolog <- final_eifs$obs_if_rr_nolog
  final_observed_if_or <- final_eifs$obs_if_or

  # return point estimates and SEs
  # se_mu_1 <- sqrt(var(final_observed_if_mu_1) / nrow(data))
  # se_mu_0 <- sqrt(var(final_observed_if_mu_0) / nrow(data))
  # se_ate <- sqrt(var(final_observed_if_ate) / nrow(data))
  # se_logrr <- sqrt(var(final_observed_if_rr) / nrow(data))
  # se_logor <- sqrt(var(final_observed_if_or) / nrow(data))
  se_mu_1 <- as.numeric(sqrt(crossprod(matrix(final_observed_if_mu_1))) / (nrow(data)))
  se_mu_0 <- as.numeric(sqrt(crossprod(matrix(final_observed_if_mu_0))) / (nrow(data)))
  se_ate <- as.numeric(sqrt(crossprod(matrix(final_observed_if_ate))) / (nrow(data)))
  se_logrr <- as.numeric(sqrt(crossprod(matrix(final_observed_if_rr))) / (nrow(data)))
  se_rr <- as.numeric(sqrt(crossprod(matrix(final_observed_if_rr_nolog))) / (nrow(data)))
  se_logor <- as.numeric(sqrt(crossprod(matrix(final_observed_if_or))) / (nrow(data)))

  result_tib <- tibble::tibble(
    "estimand" = c("mu1", "mu0", "ATE",
                   # "logRR",
                   "RR", "OR"),
    "outcome_tx_estimator" = outcome_tx_estimator,
    "miss_estimator" = miss_estimator,
    "est" = c(mu_1, mu_0, ate,
              # rr,
              rr, or),
    "SE" = c(se_mu_1, se_mu_0, se_ate, se_logrr,
             # se_rr,
             se_logor)
  )
  return(
    list(
      "result" = result_tib,
      "np_ifs" = list(
        "mu_1" = final_eifs$np_if_mu_1, "mu_0" = final_eifs$np_if_mu_0
      ),
      "par_ifs" = list(
        "mu_1" = final_eifs$par_if_mu_1, "mu_0" = final_eifs$par_if_mu_0
      ),
      "obs_ifs" = list(
        "mu_1" = final_eifs$obs_if_mu_1, "mu_0" = final_eifs$obs_if_mu_0
      ),
      "calibrated_weights" = prob_obs_long
    )
  )

}

# wrapper for glm that can handle Y as an input covariate
SL.glm.withy <- function(Y, X, newX, family, obsWeights, model = TRUE,
                         ...) {
  # check to see if Y is part of X or newX; change name to "Y.orig"
  contains_Y <- (names(X) == "Y")
  if (any(contains_Y)) {
    names(X)[contains_Y] <- "Y.orig"
    names(newX)[contains_Y] <- "Y.orig"
  }

  # fit SL.glm with the new Y names
  SL.glm(Y = Y, X = X, newX = newX, family = family, obsWeights = obsWeights,
         model = model, ...)
}
