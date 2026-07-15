# code to clean the VCCC data
# this code is largely based off of the code provided in the
# supplementary materials of Oh et al. (2021) Biometrical Journal

# load required packages and functions -----------------------------------------
packages <- c("dplyr", "tidyr", "here", "rprojroot", "data.table")
lapply(packages, library, character.only = TRUE)

this_path <- normalizePath(".", mustWork = FALSE)
proj_root <- rprojroot::find_root_file(criterion = ".projectile", path = this_path)
# anonymized
path2 <- "path_to_VCCC_data"
data_dir <- paste0(path2, "/ade_death_mgr_analysis_dataset/")

source(paste0(proj_root, "/code/data_analysis/00_vccc_cleaning_utils.R"))

if (!dir.exists(data_dir)) {
  dir.create(data_dir, recursive = TRUE)
}

# set the timepoint for outcome
timepoint <- 5
min_followup_time <- 3.5

# read in the data -------------------------------------------------------------

# antiretroviral therapy (ART) data
nonvalid_art_init <- read.csv(
  paste0(path2, "shaw_non_valid_art_17may2016.csv"),
  header = TRUE,
  stringsAsFactors = FALSE
)
valid_art_init <- read.csv(
  paste0(path2, "shaw_valid_art_17may2016.csv"),
  header = TRUE,
  stringsAsFactors = FALSE
)
# demographics data
demo_init <- read.csv(
  paste0(path2, "data/shaw_demo_17may2016.csv"),
  header = TRUE,
  stringsAsFactors = FALSE
)
# lab data
lab_init <- read.csv(
  paste0(path2, "data/shaw_labs_17may2016.csv"),
  header = TRUE,
  stringsAsFactors = FALSE
)
# opportunistic infections (OIs) data
nonvalid_ois_init <- read.csv(
  paste0(path2, "shaw_non_valid_ois_17may2016.csv"),
  header = TRUE,
  stringsAsFactors = FALSE
)
valid_ois_init <- read.csv(
  paste0(path2, "shaw_valid_ois_17may2016.csv"),
  header = TRUE,
  stringsAsFactors = FALSE
)

# merge demographic data into ART data and OI data
nonvalid_art_demo <- nonvalid_art_init %>%
  left_join(demo_init, by = "CFAR_PID") %>%
  filter(!is.na(AGE_AT_MED_START), AGE_AT_MED_START > 0)

valid_art_demo <- valid_art_init %>%
  left_join(demo_init, by = "CFAR_PID") %>%
  filter(!is.na(AGE_AT_MED_START), AGE_AT_MED_START > 0)

nonvalid_ois_demo <- nonvalid_ois_init %>%
  left_join(demo_init, by = "CFAR_PID")

valid_ois_demo <- valid_ois_init %>%
  left_join(demo_init, by = "CFAR_PID")

# define combination therapies: ------------------------------------------------
# if someone started multiple on the same day, or if they are explicit combinations
# pick first ART initiation
nonvalid_art_demo2 <- nonvalid_art_demo %>%
  group_by(CFAR_PID) %>%
  filter(AGE_AT_MED_START == min(AGE_AT_MED_START)) %>%
  group_by(AGE_AT_MED_START, AGE_AT_MED_STOP, CFAR_PID, SOURCE, YEAR_OF_RX_START, YEAR_OF_RX_STOP,
           AGE_AT_DEATH, AGE_AT_FIRST_VISIT, AGE_AT_LAST_VISIT, RACEETH, RISK1, RISK2,
           SEX, YEAR_OF_DEATH, YEAR_OF_ENROLLMENT, YEAR_OF_LAST_VISIT) %>%
  summarize(COMBINATION_ART_NAME = paste0(GENERIC_ART_NAME, collapse = "-"), .groups = "drop") %>%
  group_by(CFAR_PID) %>%
  arrange(AGE_AT_MED_START) %>%
  slice(1)

valid_art_demo2 <- valid_art_demo %>%
  group_by(CFAR_PID) %>%
  filter(AGE_AT_MED_START == min(AGE_AT_MED_START)) %>%
  group_by(AGE_AT_MED_START, AGE_AT_MED_STOP, CFAR_PID, YEAR_OF_RX_START, YEAR_OF_RX_STOP,
           AGE_AT_DEATH, AGE_AT_FIRST_VISIT, AGE_AT_LAST_VISIT, RACEETH, RISK1, RISK2,
           SEX, YEAR_OF_DEATH, YEAR_OF_ENROLLMENT, YEAR_OF_LAST_VISIT) %>%
  summarize(COMBINATION_ART_NAME = paste0(GENERIC_ART_NAME, collapse = "-"), .groups = "drop") %>%
  group_by(CFAR_PID) %>%
  arrange(AGE_AT_MED_START) %>%
  slice(1)

# apply exclusion criteria to the unvalidated data -----------------------------
# (1) exclude people who started on a non-combination therapy or started ART before enrollment
nonvalid_art_demo_exclude_1 <- nonvalid_art_demo2 %>%
  filter(grepl("-", COMBINATION_ART_NAME)) %>%
  filter(AGE_AT_MED_START >= AGE_AT_FIRST_VISIT)

# (2) exclude people who had an AIDS-defining event prior to enrollment or start of ART
nonvalid_merged1 <- nonvalid_ois_init %>%
  filter(!is.na(AGE_AT_DX_ONSET)) %>% # remove ADEs with undeterminate start date
  right_join(nonvalid_art_demo_exclude_1 %>% select(-SOURCE), by = "CFAR_PID") %>%
  group_by(CFAR_PID) %>%
  arrange(AGE_AT_DX_ONSET) %>% # select first ADE
  slice(1) %>%
  ungroup() %>%
  filter(AGE_AT_DX_ONSET >= AGE_AT_MED_START | is.na(AGE_AT_DX_ONSET))

# (3) require at least ~5 years of followup, unless they had an event earlier than that
exclude_3 <- lab_init %>%
  filter(CFAR_PID %in% nonvalid_merged1$CFAR_PID) %>%
  left_join(nonvalid_merged1, by = "CFAR_PID") %>%
  group_by(CFAR_PID) %>%
  mutate(time_in_study = AGE_AT_RESULT_DATE - AGE_AT_FIRST_VISIT) %>%
  slice_tail() %>%
  filter(time_in_study >= min_followup_time | !is.na(AGE_AT_DX_ONSET)) %>%
  pull(CFAR_PID)

nonvalid_merged <- nonvalid_merged1 %>%
  filter(CFAR_PID %in% exclude_3)

# apply exclusion criteria to validated data -----------------------------------
# (1) exclude people who started on a non-combination therapy or started ART before enrollment
valid_art_demo_exclude_1 <- valid_art_demo2 %>%
  filter(grepl("-", COMBINATION_ART_NAME)) %>%
  filter(AGE_AT_MED_START >= AGE_AT_FIRST_VISIT)

# (2) exclude people who had an AIDS-defining event prior to enrollment or start of ART
valid_merged1 <- valid_ois_init %>%
  filter(!is.na(AGE_AT_DX_ONSET)) %>% # remove ADEs with undeterminate start date
  right_join(valid_art_demo_exclude_1, by = "CFAR_PID") %>%
  group_by(CFAR_PID) %>%
  arrange(AGE_AT_DX_ONSET) %>% # select first ADE
  slice(1) %>%
  ungroup() %>%
  filter(AGE_AT_DX_ONSET >= AGE_AT_MED_START | is.na(AGE_AT_DX_ONSET))

# (3) require at least ~5 years of followup, unless they had an event earlier than that
valid_exclude_3 <- lab_init %>%
  filter(CFAR_PID %in% valid_merged1$CFAR_PID) %>%
  left_join(valid_merged1, by = "CFAR_PID") %>%
  group_by(CFAR_PID) %>%
  mutate(time_in_study = AGE_AT_RESULT_DATE - AGE_AT_FIRST_VISIT) %>%
  slice_tail() %>%
  filter(time_in_study >= min_followup_time | !is.na(AGE_AT_DX_ONSET)) %>%
  pull(CFAR_PID)

valid_merged <- valid_merged1 %>%
  filter(CFAR_PID %in% valid_exclude_3)

# get CD4 count and VL closest to ART start dates ------------------------------
nonvalid_cd4_vl_lst <- lapply(
  unique(nonvalid_merged$CFAR_PID),
  getClosestCD4toART,
  lab_df = lab_init,
  art_df = nonvalid_merged
)
nonvalid_cd4_vl <- data.table::rbindlist(nonvalid_cd4_vl_lst, fill = TRUE) %>%
  select(CFAR_PID, cd4, vl)

valid_cd4_vl_lst <- lapply(
  unique(valid_merged$CFAR_PID),
  getClosestCD4toART,
  lab_df = lab_init,
  art_df = valid_merged
)
valid_cd4_vl <- data.table::rbindlist(valid_cd4_vl_lst, fill = TRUE) %>%
  select(CFAR_PID, cd4, vl)

# define final dataset ---------------------------------------------------------
full_nonvalid <- nonvalid_merged %>%
  left_join(nonvalid_cd4_vl, by = "CFAR_PID") %>%
  mutate(time = AGE_AT_DX_ONSET - AGE_AT_MED_START,
         time_death = AGE_AT_DEATH - AGE_AT_MED_START,
         death_by_timepoint = ifelse(is.na(time_death), 0, as.numeric(time_death <= timepoint)),
         ade_by_timepoint = ifelse(is.na(time), 0, as.numeric(time <= timepoint)),
         ade_death_bin = as.numeric(death_by_timepoint == 1 | ade_by_timepoint == 1),
         contains_nnrti = regimen_contains_nnrti(COMBINATION_ART_NAME),
         contains_nrti = regimen_contains_nrti(COMBINATION_ART_NAME),
         contains_pi = regimen_contains_pi(COMBINATION_ART_NAME),
         contains_insti = regimen_contains_insti(COMBINATION_ART_NAME),
         race_combined = ifelse(!(RACEETH %in% c(1, 2)), 3, RACEETH),
         # race = relevel(factor(race_combined), ref = "1"),
         race = as.numeric(RACEETH == 2),
         risk_combined = ifelse(RISK1 %in% c("Other", "Unknown"), "Other or unknown", RISK1),
         # risk = relevel(factor(risk_combined), ref = "Hetero"),
         risk = as.numeric(RISK1 == "IDU"),
         sex = factor(SEX),
         year = YEAR_OF_RX_START - 1998,
         logvl = log(vl),
         cd4 = sqrt(cd4)) %>%
         # year = factor(YEAR_OF_RX_START)) %>%
  select(CFAR_PID, AGE_AT_MED_START, race, risk, sex, COMBINATION_ART_NAME,
         ade_death_bin, cd4, logvl, year, starts_with("contains_")) %>%
  rename(id = CFAR_PID, age_star = AGE_AT_MED_START, race_star = race, risk_star = risk,
         sex_star = sex, regimen_star = COMBINATION_ART_NAME, ade_death_bin_star = ade_death_bin,
         cd4_star = cd4, logvl_star = logvl, year_star = year,
         contains_nnrti_star = contains_nnrti, contains_nrti_star = contains_nrti,
         contains_insti_star = contains_insti, contains_pi_star = contains_pi)

full_valid <- valid_merged %>%
  left_join(valid_cd4_vl, by = "CFAR_PID") %>%
  mutate(time = AGE_AT_DX_ONSET - AGE_AT_MED_START,
         time_death = AGE_AT_DEATH - AGE_AT_MED_START,
         death_by_timepoint = ifelse(is.na(time_death), 0, as.numeric(time_death <= timepoint)),
         ade_by_timepoint = ifelse(is.na(time), 0, as.numeric(time <= timepoint)),
         ade_death_bin = as.numeric(death_by_timepoint == 1 | ade_by_timepoint == 1),
         contains_nnrti = regimen_contains_nnrti(COMBINATION_ART_NAME),
         contains_nrti = regimen_contains_nrti(COMBINATION_ART_NAME),
         contains_pi = regimen_contains_pi(COMBINATION_ART_NAME),
         contains_insti = regimen_contains_insti(COMBINATION_ART_NAME),
         race_combined = ifelse(!(RACEETH %in% c(1, 2)), 3, RACEETH),
         # race = relevel(factor(race_combined), ref = "1"),
         race = as.numeric(RACEETH == 2),
         risk_combined = ifelse(RISK1 %in% c("Other", "Unknown"), "Other or unknown", RISK1),
         # risk = relevel(factor(risk_combined), ref = "Hetero"),
         risk = as.numeric(RISK1 == "IDU"),
         sex = factor(SEX),
         year = YEAR_OF_RX_START - 1998,
         logvl = log(vl),
         cd4 = sqrt(cd4)) %>%
         # year = factor(YEAR_OF_RX_START)) %>%
  select(CFAR_PID, AGE_AT_MED_START, race, risk, sex, COMBINATION_ART_NAME,
         ade_death_bin, cd4, logvl, year, starts_with("contains_")) %>%
  rename(id = CFAR_PID, age = AGE_AT_MED_START, regimen = COMBINATION_ART_NAME)

summary(full_nonvalid)
summary(full_valid)

# people who have both non-validated and validated data
vccc_dat_init <- full_nonvalid %>%
  inner_join(full_valid, by = "id")

# check to see which variables have measurement error
summary(vccc_dat_init)
nrow(vccc_dat_init)
sum(vccc_dat_init$age == vccc_dat_init$age_star)
sum(vccc_dat_init$race == vccc_dat_init$race_star)
sum(vccc_dat_init$risk == vccc_dat_init$risk_star)
sum(vccc_dat_init$sex == vccc_dat_init$sex_star)
sum(vccc_dat_init$regimen == vccc_dat_init$regimen_star)
sum(vccc_dat_init$year == vccc_dat_init$year_star)

# create the final dataset; exclude people who are on a regimen with fewer than 3 medications
nonvalid_exclude_regimen_2less <- full_nonvalid %>%
  mutate(num_meds_star = nchar(regimen_star) - nchar(gsub("-", "", regimen_star))) %>%
  filter(num_meds_star >= 2)
valid_exclude_regimen_2less <- full_valid %>%
  mutate(num_meds = nchar(regimen) - nchar(gsub("-", "", regimen))) %>%
  filter(num_meds >= 2)

# subset to complete cases
nonvalid_exclude_regimen_2less %>%
  inner_join(valid_exclude_regimen_2less, by = "id")
vccc_dat <- vccc_dat_init %>%
  select(-age_star, -race_star, -risk_star, -sex_star) %>%
  mutate(num_meds = nchar(regimen) - nchar(gsub("-", "", regimen)),
         num_meds_star = nchar(regimen_star) - nchar(gsub("-", "", regimen_star))) %>%
  filter(num_meds >= 2 & num_meds_star >= 2) %>%
  select(-num_meds, -num_meds_star) %>%
  filter(complete.cases(.))

# data quality checks ----------------------------------------------------------
# check to make sure one row per participant
nrow(vccc_dat); length(unique(vccc_dat$id))

# check to make sure all are combination therapies
sum(grepl("-", vccc_dat$regimen))

# check to make sure all ART started after enrollment
nrow(nonvalid_merged); sum(nonvalid_merged$AGE_AT_MED_START >= nonvalid_merged$AGE_AT_FIRST_VISIT)
nrow(valid_merged); sum(valid_merged$AGE_AT_MED_START >= valid_merged$AGE_AT_FIRST_VISIT)

# check to make sure all ADEs occurred after enrollment and ART start
nrow(nonvalid_merged); sum(nonvalid_merged$AGE_AT_DX_ONSET >= nonvalid_merged$AGE_AT_MED_START | is.na(nonvalid_merged$AGE_AT_DX_ONSET))
nrow(valid_merged); sum(valid_merged$AGE_AT_DX_ONSET >= valid_merged$AGE_AT_MED_START | is.na(valid_merged$AGE_AT_DX_ONSET))

# check to make sure that there is a difference between non-validated and validated data
sum(vccc_dat$ade_death_bin); sum(vccc_dat$ade_death_bin_star)
summary(vccc_dat$cd4); summary(vccc_dat$cd4_star)
summary(vccc_dat$logvl); summary(vccc_dat$logvl_star)

# summary statistics -----------------------------------------------------------
exclusion_tib_nonvalidated <- tibble::tibble(
  "n_overall" = nrow(demo_init),
  "n_valid_art_data" = length(unique(nonvalid_art_demo$CFAR_PID)),
  "n_exclude_non_combo_art_preenrollment" = length(unique(nonvalid_art_demo_exclude_1$CFAR_PID)),
  "n_exclude_ade_preenrollment_pretx" = length(unique(nonvalid_merged1$CFAR_PID)),
  "n_exclude_less5yr" = length(unique(nonvalid_merged$CFAR_PID)),
  "n_with_cd4_vl" = length(unique(full_nonvalid$id)),
  "n_exclude_non_combo" = sum(complete.cases(nonvalid_exclude_regimen_2less)),
  "n_final" = nrow(vccc_dat)
)

exclusion_tib_validated <- tibble::tibble(
  "n_overall" = nrow(demo_init),
  "n_valid_art_data" = length(unique(valid_art_demo$CFAR_PID)),
  "n_exclude_non_combo_art_preenrollment" = length(unique(valid_art_demo_exclude_1$CFAR_PID)),
  "n_exclude_ade_preenrollment_pretx" = length(unique(valid_merged1$CFAR_PID)),
  "n_exclude_less5yr" = length(unique(valid_merged$CFAR_PID)),
  "n_exclude_non_combo" = sum(complete.cases(valid_exclude_regimen_2less)),
  "n_final" = nrow(vccc_dat)
)
summary_stats <- exclusion_tib_nonvalidated %>%
  mutate(dataset = "unvalidated") %>%
  bind_rows(
    exclusion_tib_validated %>% mutate(dataset = "validated")
  )
write.csv(summary_stats, paste0(data_dir, "vccc_data_cleaning_summary_statistics.csv"))

# print out the regimens
all_regimens <- unique(vccc_dat$regimen)
# get unique generic names
all_generics <- unique(unlist(strsplit(all_regimens, "-", fixed = TRUE)))
write.csv(all_regimens, paste0(data_dir, "vccc_data_unique_regimens.csv"))
write.csv(all_generics, paste0(data_dir, "vccc_data_unique_generics.csv"))

# get counts of people on regimens in each class
table(vccc_dat$contains_nnrti)
table(vccc_dat$contains_nrti)
table(vccc_dat$contains_insti)
table(vccc_dat$contains_pi)
table(!vccc_dat$contains_nnrti & !vccc_dat$contains_insti & !vccc_dat$contains_nrti & !vccc_dat$contains_pi)

# combinations
table("nnrti" = vccc_dat$contains_nnrti, "pi" = vccc_dat$contains_pi)
table("nrti" = vccc_dat$contains_nrti, "pi" = vccc_dat$contains_pi)
table("insti" = vccc_dat$contains_insti, "pi" = vccc_dat$contains_pi)
table("nnrti" = vccc_dat$contains_nnrti, "nrti" = vccc_dat$contains_nrti)
table("nnrti" = vccc_dat$contains_nnrti, "insti" = vccc_dat$contains_insti)
table("nrti" = vccc_dat$contains_nrti, "insti" = vccc_dat$contains_insti)

table(rowSums(vccc_dat %>% select(starts_with("contains_") & !contains("star"))))

three_class_regimens <- which(rowSums(vccc_dat %>% select(starts_with("contains_") & !contains("star"))) == 3)
vccc_dat[three_class_regimens, ] %>%
  select(regimen_star, starts_with("contains") & !contains("star"))
# save the dataset -------------------------------------------------------------
# clean dataset, can be used for summary statistics
saveRDS(vccc_dat, paste0(data_dir, "vccc_data_clean_all.rds"))
# analytic dataset, which doesn't contain regimen* columns and all contains_* columns but PI
saveRDS(vccc_dat %>%
          select(-contains("regimen"), -contains("contains_nnrti"),
                 -contains("contains_nrti"), -contains("contains_insti")),
        paste0(data_dir, "vccc_data_clean.rds"))

