# get the true values for the simulations
library("tidyr")
library("dplyr")
library("future.apply")
library("here") # for directory traversing
library("rprojroot") # for trajectory traversing
this_path <- normalizePath(".", mustWork = FALSE)
proj_root <- rprojroot::find_root_file(criterion = ".projectile", path = this_path)

source(paste0(proj_root, "/code/00_utils.R"))
source(paste0(proj_root, "/code/sims/00_sim_utils.R"))
source(paste0(proj_root, "/code/sims/01_gen_data.R"))
source(paste0(proj_root, "/code/sims/02_run_ests.R"))


# scenarios <- c(1:10, 15)
scenarios <- c(1:10, 97, 96, 99, 98)
future::plan(multisession)
set.seed(20240919)
true_values_list <- future.apply::future_lapply(
  X = as.list(seq_len(length(scenarios))), FUN = function(i) {
    get_true_values(scenario = scenarios[i], sample_size = 1e6,
                    nreps = 10, gencor = 0)
  }, future.seed = TRUE
)
true_values <- do.call(rbind.data.frame, true_values_list)
true_values_summ <- true_values %>%
  group_by(scenario, estimand) %>%
  summarize(value = mean(value), .groups = "drop")

saveRDS(true_values, file = paste0(proj_root, "/code/sims/true_values_long.rds"))
saveRDS(true_values_summ, file = paste0(proj_root, "/code/sims/true_values.rds"))

true_values_summ %>%
  print(n = Inf)
