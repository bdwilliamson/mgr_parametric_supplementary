# subsample phase 2 datasets from VCCC data

# load required packages and functions -----------------------------------------
packages <- c("dplyr", "tidyr", "here", "rprojroot")
lapply(packages, library, character.only = TRUE)

this_path <- normalizePath(".", mustWork = FALSE)
proj_root <- rprojroot::find_root_file(criterion = ".projectile", path = this_path)
# anonymized
path2 <- "path_to_VCCC_data"
data_dir <- paste0(path2, "/ade_death_mgr_analysis_dataset/")

# read in the dataset ----------------------------------------------------------
vccc_dat <- readRDS(file = paste0(data_dir, "vccc_data_clean.rds"))

# outcome-stratified sampling --------------------------------------------------
# select all cases and a subset of controls; for controls, stratify on high risk people
errorprone_cases <- vccc_dat %>%
  filter(ade_death_bin_star == 1)
errorprone_controls <- vccc_dat %>%
  filter(ade_death_bin_star == 0)

# n_total <- round(nrow(vccc_dat) * 0.75)
n_total <- round(nrow(vccc_dat) * 0.5)
n_cases <- nrow(errorprone_cases)
n_controls <- n_total - n_cases

n_datasets <- 1000
set.seed(20260130)
seeds <- sample(1e5:1e6, n_datasets)
dataset_list <- vector("list", length = n_datasets)
for (i in 1:length(dataset_list)) {
  set.seed(seeds[i])
  these_controls <- errorprone_controls[
  sample(1:nrow(errorprone_controls), n_controls, replace = FALSE),
  ]
  this_validated_subset <- dplyr::bind_rows(
    errorprone_cases, these_controls
  )

  this_dataset <- vccc_dat %>%
    mutate(validated = as.numeric(id %in% c(this_validated_subset$id)),
           ade_death_bin = ifelse(validated == 0, NA, ade_death_bin),
           contains_pi = ifelse(validated == 0, NA, contains_pi),
           cd4 = ifelse(validated == 0, NA, cd4),
           logvl = ifelse(validated == 0, NA, logvl),
           year = ifelse(validated == 0, NA, year))

  dataset_list[[i]] <- this_dataset
}

# save the datasets ------------------------------------------------------------
saveRDS(dataset_list, file = paste0(proj_root, "/data/vccc_sampled_datasets.rds"))
