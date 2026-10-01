# generate a dataset

# formulas for scenarios -------------------------------------------------------
# Get oracle, working formulas for outcome regression, missing-data model,
# and propensity score.
# ** THESE SHOULD BE EDITED WHEN NEW SCENARIOS ARE ADDED **
# @param scenario the scenario (e.g., 1) which determines the outcome regression and missing-data model formulas
# @param plasmode logical, whether or not this is a plasmode data analysis
get_formulas <- function(scenario = 1, include_a = FALSE) {
  if (scenario > 3) {
    include_a <- TRUE
  }
  simple_scenarios <- c(1:2, 9:10, 12, 14, 17, 18, 96, 97, 98, 99)
  misspecified_scenarios <- c(3:4, 5:8)
  complex_scenarios <- c(15, 16, 17, 18)
  census_formula <- paste0("Y ~ X + Zs + Zw + Ws + Ww", ifelse(include_a, " + As + Aw", ""))
  outcome_formula <- paste0("Y ~ X + Zs + Zw + Ws_obs + Ww_obs", ifelse(include_a, " + As + Aw", ""))
  outcome_formula_confounded <- paste0("Y ~ X + Zs + Zw", ifelse(include_a, " + As + Aw", ""))
  propensity_formula <- paste0(paste0("X ~ Zs + Zw", ifelse(!(scenario %in% c(96,97)), " + Ws_obs + Ww_obs", "")),
                                    ifelse(include_a, " + As + Aw", ""))
  phase_1_formula <- paste0("~ Y + X + Zs + Zw", ifelse(include_a, " + As + Aw", ""))
  # phase_1_formula <- paste0("~ Y + X + Zs + Zw", ifelse(include_a, " + As + Aw", ""),
                            # ifelse(scenario %in% c(5:8), "+ Us + Uw", ""))
  if (scenario %in% c(7, 8)) {
    # OR is misspecified but PS is correct
    propensity_formula <- paste0(propensity_formula, " + Us + Uw")
  } else if (scenario %in% c(5, 6)) {
    # PS is misspecified but OR is correct
    outcome_formula <- paste0(outcome_formula, " + Us + Uw")
  }
  miss_formula <- paste0("is.complete ~ X + Zs + Zw + Y", ifelse(include_a, " + As + Aw", ""))
  if (as.character(scenario) %in% as.character(simple_scenarios)) {
    oracle_formula <- paste0("Y ~ X + Zs + Zw + Ws + Ww", ifelse(include_a, " + As + Aw", ""))
  } else if (as.character(scenario) %in% as.character(misspecified_scenarios)) {
    oracle_formula <- paste0("Y ~ X + Zs + Zw + Ws + Ww + Us + Uw", ifelse(include_a, " + As + Aw", ""))
  } else if (as.character(scenario) %in% as.character(complex_scenarios)) {
    oracle_formula <- paste0("Y ~ X + Zs + I(Zs ^ 2) + Zw + Ws + Ww + Ws * Zs", ifelse(include_a, " + As + Aw", ""))
  } else {
    stop("The requested scenario has not been implemented.")
  }
  return(
    list(
      "oracle_outcome_formula" = oracle_formula, "census_outcome_formula" = census_formula,
      "working_outcome_formula" = outcome_formula, "confounded_outcome_formula" = outcome_formula_confounded,
      "miss_formula" = miss_formula, "propensity_formula" = propensity_formula,
      "phase_1_formula" = phase_1_formula
    )
  )
}


# generate an entire dataset ---------------------------------------------------
# @param N the sample size
# @param scenario the scenario number, determines covariate distribution, propensity score,
#   outcome regression model, missing-data model
# @param lowcor low correlation
# @param midcor middle correlation
# @param highcor high correlation
# @param gencor general correlation
# @return a dataset, with covariates, outcome, missing-data indicator
gen_data <- function(N = 100, scenario = 1,
                     gencor = 0.2,
                     debug_miss = FALSE) {
  # generate covariates (everything fully observed)
  covariates <- gen_x(N = N, scenario = scenario,
                      gencor = gencor)
  # generate outcome based on covariates (everything fully observed)
  Y <- gen_y(data = covariates, scenario = scenario)$Y
  df <- cbind.data.frame(Y = Y, covariates)
  # generate missing-data indicators
  full_miss <- gen_miss(data = df, scenario = scenario)
  miss <- full_miss$miss
  # create missing-data versions of confounders
  df$Ws_obs <- ifelse(miss == 1, NA, df$Ws)
  df$Ww_obs <- ifelse(miss == 1, NA, df$Ww)
  df$is.complete <- as.numeric(miss == 0)
  if (debug_miss) {
    df$true_miss_prob <- full_miss$probs
  }
  # return
  return(df)
}
# outcome regression model -----------------------------------------------------
# @param data the covariate data (expects columns X, Zw, Zs, Ww, Ws, Uw, Us)
# @param scenario the outcome regression model scenario of interest
#        1: continuous outcome, simple linear outcome model
#        2: binary outcome, simple logistic outcome model
# @return a vector of outcomes
gen_y <- function(data, scenario = 1) {
  continuous_outcome_scenarios <- c(1, 3, 5, 7, 9, 15, 17, 97, 99)
  binary_outcome_scenarios <- c(2, 4, 6, 8, 10, 12, 14, 16, 18, 96, 98)
  probit_scenarios <- c(14)
  misspecified_scenarios <- c(5:8)
  rare_scenarios <- c(3, 4, 9, 10, 12, 14, 15, 16, 17, 18)
  complex_scenarios <- c(15, 16)
  # specify regression model parameters
  # if (scenario %in% continuous_outcome_scenarios) {
  if ((scenario %% 2) == 1) { # odd numbered scenarios are continuous
    # continuous outcome, simple model
    beta0 <- -0.25
    betaX <- 1.5
    betaZs2 <- 0
    betaZw <- -1
    betaZs <- 1
    betaWw <- -1
    betaWs <- 1
    betaUw <- 0
    betaUs <- 0
    betaAs <- 0
    betaAw <- 0
    betaZsZw <- 0
    if (scenario %in% complex_scenarios) {
      betaZs2 <- 1.5
      betaZsZw <- 1
    }
    all_betas <- data.frame(
      "name" = c("Intercept", "X", "Zw", "Zs", "Ww", "Ws", "As", "Aw", "Us", "Uw", "Zs2", "ZsZw"),
      "value" = c(beta0, betaX, betaZw, betaZs, betaWw, betaWs, betaAs, betaAw, betaUs, betaUw, betaZs2, betaZsZw)
    )
  # } else if (scenario %in% binary_outcome_scenarios) {
  } else if ((scenario %% 2) == 0) {
    # binary outcome, simple models
    beta0 <- -2.6
    betaX <- 1.5
    betaZs2 <- 0
    betaZw <- -1
    betaZs <- 1
    betaWw <- -1
    betaWs <- 1
    betaUw <- 0
    betaUs <- 0
    betaAs <- 0
    betaAw <- 0
    betaZsZw <- 0
    if (scenario %in% complex_scenarios) {
      beta0 <- -5
      betaZs2 <- 1.5
      betaZsZw <- 1
    }
    all_betas <- data.frame(
      "name" = c("Intercept", "X", "Zw", "Zs", "Ww", "Ws", "As", "Aw", "Us", "Uw", "Zs2", "ZsZw"),
      "value" = c(beta0, betaX, betaZw, betaZs, betaWw, betaWs, betaAs, betaAw, betaUs, betaUw, betaZs2, betaZsZw)
    )
  } else {
    stop("This scenario has not yet been implemented.")
  }
  if (scenario %in% c(misspecified_scenarios, rare_scenarios)) {
    # rarer outcome
    beta0 <- ifelse(scenario %in% misspecified_scenarios, -8, -3)
    if (scenario %in% complex_scenarios) {
      beta0 <- -4
    }
    betaAs <- 0.75
    betaAw <- -0.75
    # non-zero coefficient for U variables
    if (scenario %in% misspecified_scenarios) {
      if (scenario %in% continuous_outcome_scenarios) {
        # only in these scenarios
        betaUw <- -1
        betaUs <- 1.5
      } else {
        betaUw <- -5
        betaUs <- 3
      }
    }
    all_betas <- data.frame(
      "name" = c("Intercept", "X", "Zw", "Zs", "Ww", "Ws", "As", "Aw", "Us", "Uw", "Zs2", "ZsZw"),
      "value" = c(beta0, betaX, betaZw, betaZs, betaWw, betaWs, betaAs, betaAw, betaUs, betaUw, betaZs2, betaZsZw)
    )
  }
  # create the true outcome regression
  if (scenario %in% c(continuous_outcome_scenarios, binary_outcome_scenarios)) {
    Xbeta <- beta0 + betaX * data$X + betaZw * data$Zw + betaZs * data$Zs +
          betaWw * data$Ww + betaWs * data$Ws + betaUw * data$Uw + betaUs * data$Us +
      betaAs * data$As + betaAw * data$Aw + betaZs2 * data$Zs ^ 2 + betaZsZw * data$Zs * data$Zw
  } else {
    stop("The current scenario has not yet been implemented.")
  }
  if (scenario %in% continuous_outcome_scenarios) {
    # continuous outcome
    Y <- as.numeric(Xbeta + rnorm(nrow(data), 0, 1))
  } else if (scenario %in% binary_outcome_scenarios) {
    # logistic outcome
    if (!(scenario %in% probit_scenarios)) {
      Y <- as.numeric(runif(nrow(data)) < stats::plogis(Xbeta))
    } else {
      # probit outcome
      Y <- as.numeric(0 < Xbeta + rnorm(nrow(data)))
    }

  }
  return(list("Y" = Y, "coefs" = all_betas))
}

# missing data model -----------------------------------------------------------
# @param data the covariate data (expects columns Y, X, Zw, Zs, Ww, Ws, Uw, Us)
# @param scenario the scenario of interest
# @return a vector of missing data indicators
gen_miss <- function(data, scenario = 1) {
  alpha0 <- 0.3
  simple_miss_25 <- c(1:2)
  simple_miss_low <- c(98:99)
  complex_scenarios <- c(17, 18)
  # if continuous outcome, need to make probabilities more reasonable
  if (scenario %in% simple_miss_25) { # approximately 20% missing data
    alpha0 <- -1.335
    alphaX <- -0.1
    alphaZw <- 0
    alphaZs <- 0.05
    alphaUs <- 0
    alphaUw <- 0
    alphaAs <- 0
    alphaAw <- 0
    alphaY <- 0.05
    alphaX2 <- 0
    alphaZsZw <- 0
  } else if (scenario %in% simple_miss_low) { # less missing data
    alpha0 <- 0.05
    alphaX <- -0.1
    alphaZw <- 0
    alphaZs <- 0.05
    alphaUs <- 0
    alphaUw <- 0
    alphaAs <- 0
    alphaAw <- 0
    alphaY <- 0.05
    alphaX2 <- 0
    alphaZsZw <- 0
  } else { # approximately 50% missing data
    alphaX <- -0.1
    alphaZw <- 0
    alphaZs <- 0.05
    alphaUs <- 0
    alphaUw <- 0
    alphaAs <- 0
    alphaAw <- 0
    alphaX2 <- 0
    alphaZsZw <- 0
    if (scenario %in% complex_scenarios) {
      alpha0 <- 0.2
      alphaX2 <- 0.1
      alphaZsZw <- 0.1
    }
    alphaY <- ifelse(scenario %% 2 == 0, 0.05, 0.025)
  }
  if (scenario > 2) {
    alphaAs <- -0.1
    alphaAw <- -0.1
  }
  if (scenario %in% c(9, 10)) {
    alphaUs <- -0.5
    alphaUw <- 0.5
  }
  XbetaM <- alpha0 + alphaX * data$X + alphaZw * data$Zw +
    alphaZs * data$Zs + alphaUs * data$Us + alphaUw * data$Uw + alphaY * data$Y +
    alphaAs * data$As + alphaAw * data$Aw + alphaX2 * data$X ^ 2 + alphaZsZw * data$Zs * data$Zw
  all_betas <- data.frame(
    "name" = c("Intercept", "X", "Zw", "Zs", "Aw", "As", "Uw", "Us", "Y", "X2", "ZsZw"),
    "value" = c(alpha0, alphaX, alphaZw, alphaZs, alphaAw, alphaAs, alphaUw, alphaUs, alphaY, alphaX2, alphaZsZw)
  )

  miss <- as.numeric(runif(nrow(data)) < stats::plogis(XbetaM))
  return(list("miss" = miss, "coefs" = all_betas, "probs" = stats::plogis(XbetaM)))
}

# generate de novo covariates --------------------------------------------------
# @param N the sample size
# @param lowcor low correlation
# @param midcor middle correlation
# @param highcor high correlation
# @param gencor general correlation
# @param XScenario the scenario for generating covariates
# @return a dataset with the covariates
gen_x <- function(N = 100, scenario = 1, gencor = 0) {
  corMat <- matrix(gencor,nrow=8,ncol=8)

  varMat <- corMat
  diag(varMat) <- 1

  ZWUA <- mvtnorm::rmvnorm(N, c(rep(0,ncol(corMat) - 1), 2), varMat)

  Zs <- ZWUA[,1]
  Zw <- ZWUA[,2]
  Ws <- ZWUA[,3]
  Ww <- ZWUA[,4]
  Us <- ZWUA[,5]
  Uw <- ZWUA[,6]
  As <- ZWUA[,7]
  Aw <- ZWUA[,8]

  # Treatment covariate of interest.
  beta <- c(rep(0.1, 5), rep(0, 4))
  no_w_tx_ps <- c(96, 97)
  # misspecified_scenarios <- c(5, 6)
  misspecified_scenarios <- c(5:8)
  include_a <- c(3:20)
  if (scenario %in% no_w_tx_ps) {
    # X does not depend on W's
    beta <- c(rep(0.1, 3), rep(0, 6))
  }
  if (scenario %in% misspecified_scenarios) {
    beta <- c(rep(0.1, 5), c(-1, -1), rep(0, 2))
  }
  if (scenario %in% include_a) {
    beta[8] <- -0.4
    beta[9] <- -0.2
  }
  pX <- stats::plogis(cbind(1, ZWUA) %*% beta)
  X <- rbinom(N, size = 1, prob = pX)

  data <- data.frame(X = X, Zs = Zs, Zw = Zw, Ws = Ws, Ww = Ww, Us = Us, Uw = Uw, As = As, Aw = Aw)
  return(data)
}
