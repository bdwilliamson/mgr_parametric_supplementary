# main simulation engine

# load required functions and packages -----------------------------------------
library("dplyr"); library("purrr")
library("future.apply") # for parallelization
library("sandwich") # for robust SEs from glms
library("mice") # for GR
library("survey") # for GR for a conditional estimand
library("marginaleffects") # for marginalizing GR for the conditional estimand
library("SuperLearner") # for flexible learners
library("splines") # for b-spline basis
library("tmle") # for TMLE
library("twoStageDesignTMLE") # for two-stage TMLE
library("xgboost") # flexible learner
library("ranger") # flexible learner
library("hal9001") # flexible learner
library("optparse") # for dynamic args
library("here") # for directory traversing
library("rprojroot") # for trajectory traversing
this_path <- normalizePath(".", mustWork = FALSE)
proj_root <- rprojroot::find_root_file(criterion = ".projectile", path = this_path)

source(paste0(proj_root, "/code/00_utils.R"))
source(paste0(proj_root, "/code/sims/01_gen_data.R"))
source(paste0(proj_root, "/code/sims/02_run_ests.R"))

# set up dynamic and static simulation args ------------------------------------
parser <- OptionParser()
parser <- add_option(parser, "--nreps-total", default = 500, help = "Number of simulation replicates in total")
parser <- add_option(parser, "--n", default = 500, help = "Sample size")
parser <- add_option(parser, "--scenario", default = 17, help = "Outcome regression model scenario")
parser <- add_option(parser, "--estimator", default = "CGR-IPTW;CGR-adjust;MGR;AIPCW-AIPTW", help = "The estimator(s) to run")
parser <- add_option(parser, "--use-cached-datasets", default = 0, help = "Should we use cached datasets or create a new dataset each time?")
parser <- add_option(parser, "--data-only", default = 0, help = "Should we only generate datasets?")
parser <- add_option(parser, "--seed", default = 1, help = "Change the seed (already depends on y, m, x scenarios)")
parser <- add_option(parser, "--testing", default = 0, help = "Are we testing?")
parser <- add_option(parser, "--num-cores", default = 24, help = "number of cores")
parser <- add_option(parser, "--eif-reg-type", default = "impute", help = "Regression estimator of EIF projection")
parser <- add_option(parser, "--link", default = "default", help = "Should we use the default (gaussian or logistic) or another family?")
parser <- add_option(parser, "--cgr-type", default = "ipw", help = "Should we use IPTW or regression adjustment (adjust) for CGR?")
args <- parse_args(parser, convert_hyphens_to_underscores = TRUE)

if (args$testing) {
  nreps_total <- 64
  n <- 500
  estimator <- "CGR;MGR;AIPCW-AIPTW"
  scenario <- 1
  use_cached_datasets <- 0
  data_only <- 0
  seed <- 1
} else {
  nreps_total <- args$nreps_total
  n <- args$n
  estimator <- args$estimator
  scenario <- args$scenario
  use_cached_datasets <- args$use_cached_datasets
  data_only <- args$data_only
  seed <- args$seed
}

est_vec <- unlist(strsplit(estimator, ";", fixed = TRUE))
est_for_save <- paste0(est_vec, collapse = "_")

output_dir <- paste0(proj_root, "/results/sims/scenario", scenario)
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}
data_dir <- paste0(proj_root, "/code/sims/data/scenario", scenario)

# run the simulation -----------------------------------------------------------
# read in cached datasets, if we're using these
if (use_cached_datasets) {
  cached_datasets <- lapply(
    as.list(1:nreps_total), function(i) {
      readRDS(file = paste0(
        data_dir, filename_prefix, i, filename_suffix
      ))
    }
  )
} else {
  cached_datasets <- NULL
}
# set up Super Learner library (only used if "TMLE" is an estimator)
xgb_tune_params <- list(max_depth = c(1, 3), shrinkage = c(1e-2, 1e-1), ntrees = 500)
num_n <- as.numeric(n)
if (is.na(num_n)) {
  num_n <- nrow(cached_datasets[[1]])
}
min_node_sizes <- round(seq.int(num_n / 100, num_n / 10, length.out = 4))
rf_tune_params <- list(num.trees = 500, min.node.size = min_node_sizes)
xgb_learners <- create.Learner("SL.xgboost", tune = xgb_tune_params,
                               detailed_names = TRUE, name_prefix = "xgb")
rf_learners <- create.Learner("SL.ranger", tune = rf_tune_params,
                              detailed_names = TRUE, name_prefix = "rf")
continuous_learner_lib <- c("SL.glm.withy", "SL.glmnet", xgb_learners$names, rf_learners$names)
binary_learner_lib <- c("SL.glm.withy", "SL.glmnet", xgb_learners$names, rf_learners$names)
simple_lib <- "SL.glm.withy"
if (args$testing) {
  tmle_args <- list("g_lib" = simple_lib,
                    "miss_lib" = simple_lib,
                    "q_lib" = simple_lib, "K" = 10,
                    "phase1_covars" = c("Y", "X", "Zs", "Zw", "As", "Aw"),
                    "phase2_covars" = c("Ws_obs", "Ww_obs"))
} else {
  tmle_args <- list("g_lib" = binary_learner_lib,
                    "miss_lib" = binary_learner_lib,
                    "q_lib" = continuous_learner_lib, "K" = 10,
                    "phase1_covars" = c("Y", "X", "Zs", "Zw", "As", "Aw"),
                    "phase2_covars" = c("Ws_obs", "Ww_obs"))
}
n_imp <- 10
raking_args <- list("n_imp" = n_imp, "rake_on_outcome" = TRUE,
                    "rake_on_exposure" = FALSE, "rake_on_covars" = FALSE,
                    "rake_on_eifs" = TRUE, "calibrate_separate" = FALSE)


# get the current seed
current_seed <- round(scenario * 100) + round(n * 1000) + round((seed - 1) * 51)

# set up parallelization
future::plan(multisession, workers = args$num_cores)
options(future.globals.maxSize = +Inf)
set.seed(current_seed)
seeds <- future_lapply(as.list(seq_len(nreps_total)), FUN = function(x) .Random.seed,
                       future.chunk.size = Inf, future.seed = current_seed)

continuous_outcome_scenarios <- c(1, 3, 5, 7, 9, 15, 17, 97, 99)
if (scenario == 12) {
  args$link <- "poisson"
}
if (scenario == 14) {
  args$link <- "probit"
}
if (scenario %in% continuous_outcome_scenarios) {
  family <- "gaussian"
  fam_txt <- "identity"
} else {
  if (args$link == "default") {
    family <- "binomial"
    fam_txt <- "logistic"
  } else if (args$link == "probit") {
    family <- binomial(link = "probit")
    fam_txt <- "probit"
  } else if (args$link == "identity") {
    family <- "gaussian"
    fam_txt <- "identity"
  } else {
    family <- "poisson"
    fam_txt <- "poisson"
  }
}
use_hal <- (args$eif_reg_type == "hal")
use_basis <- (args$eif_reg_type == "spline")
use_imputation <- (args$eif_reg_type == "impute")

# actually run the simulation
# run the simulation
cat("\nRunning simulation\n")
start <- Sys.time()
output_list <- future.apply::future_lapply(
  X = as.list(seq_len(nreps_total)), FUN = function(i) {
    tryCatch(
      investigate_performance_once(
        mc_id = i, n = n, data_args = list(
          scenario = scenario, gencor = 0
        ), outcome_name = "Y", tx_name = "X", missing_indicator = "is.complete",
        estimators = est_vec, tmle_args = tmle_args,
        raking_args = raking_args, cached_datasets = cached_datasets,
        data_only = data_only, fam = family, use_hal = use_hal,
        use_imputation = use_imputation, use_basis = use_basis,
        rescale = TRUE, rescale_ps = FALSE, projection_intercept = TRUE,
        weight_regressions = TRUE, use_meat = FALSE, parametric = TRUE,
        rake_on_outcome = TRUE, cgr_tx_type = args$cgr_type
      ), error = function(e) {
        message(conditionMessage(e))
        print(paste0("Error occurred when running Monte-Carlo iteration ", i))
        traceback()
        NULL
      }
    )
  }, future.seed = seeds, future.stdout = structure(TRUE, drop = TRUE),
  future.conditions = "message",
  future.globals = structure(TRUE, add = c(xgb_learners$names, rf_learners$names,
                                           "SL.glm.withy", "All"))
)
end <- Sys.time()
cat("\nElapsed time: ", format(end - start), "\n")
# write results ----------------------------------------------------------------
if (data_only) {
  lapply(as.list(1:length(output_list)), function(i) {
    saveRDS(output_list[[i]], file = paste0(
      data_dir, "/data_scenario_", scenario, "_n", n, "_id", i, ".rds"
    ))
  })
} else {
  output <- do.call(rbind, lapply(output_list, function(x) x))
  output <- output %>%
    mutate(scenario = scenario, seed = seed, n = n) %>%
    select(scenario, seed, n, contains("estimator"), mc_id, estimand, est, SE)
  # save each estimator separately (so that we can update one at a time if need be)
  for (i in 1:length(est_vec)) {
    this_output <- output %>%
      mutate(procedure = paste0(miss_estimator, "-", outcome_tx_estimator),
             link = fam_txt) %>%
      filter(grepl(est_vec[i], procedure)) %>%
      select(-procedure)
    # save as csvs
    readr::write_csv(this_output, file = paste0(
      output_dir, "/output_scenario_", scenario, "_link_", fam_txt,
      "_s", seed, "_n", n, "_est_", est_vec[i], ifelse(est_vec[i] == "CGR", paste0("_", this_output$outcome_tx_estimator[1]), ""),
      "_nreps_", nreps_total, ".csv"
    ))
  }
}
cat("\nSimulation complete!\n")
