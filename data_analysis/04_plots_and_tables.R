# compile results from data analysis, create plots and tables

# load required packages and functions -----------------------------------------
packages <- c("dplyr", "tidyr", "stringr", "here", "rprojroot", "sandwich", "survey",
              "future", "future.apply", "marginaleffects", "ggplot2",
              "kableExtra", "gtsummary")
lapply(packages, library, character.only = TRUE)

this_path <- normalizePath(".", mustWork = FALSE)
proj_root <- rprojroot::find_root_file(criterion = ".projectile", path = this_path)
plots_dir <- paste0(proj_root, "/plots/")
# anonymized
path2 <- "path_to_VCCC_data"
data_dir <- paste0(path2, "/ade_death_mgr_analysis_dataset/")

# read in analysis results -----------------------------------------------------
full_results <- readRDS(file = paste0(proj_root, "/results/data_analysis/full_data_results.rds"))
cgr_results_iptw <- readRDS(file = paste0(proj_root, "/results/data_analysis/cgr_results_iptw.rds"))
cgr_results_ra <- readRDS(file = paste0(proj_root, "/results/data_analysis/cgr_results_ra.rds"))
mgr_results <- readRDS(file = paste0(proj_root, "/results/data_analysis/mgr_results.rds"))

vccc_dat <- readRDS(file = paste0(data_dir, "vccc_data_clean.rds"))

# get summaries across the sampled datasets ------------------------------------
all_results <- bind_rows(
  data.table::rbindlist(cgr_results_iptw, use.names = TRUE),
  data.table::rbindlist(cgr_results_ra, use.names = TRUE),
  data.table::rbindlist(mgr_results, use.names = TRUE)
) %>%
  left_join(
    full_results %>% pivot_wider(names_from = dataset, values_from = value),
    by = "estimand"
  ) %>%
  filter(estimand != "OR")

ese <- all_results %>%
  group_by(estimand, miss_estimator, outcome_tx_estimator) %>%
  mutate(scaled_est = ifelse(estimand == "RR", log(est), est)) %>%
  summarize(ESE = sd(scaled_est, na.rm = TRUE), .groups = "drop")

all_results_wide <- all_results %>%
  left_join(ese, by = c("estimand", "miss_estimator", "outcome_tx_estimator")) %>%
  mutate(bias_init = est - validated,
         mse_init = (est - validated) ^ 2,
         log_est = log(est),
         ci_ll = ifelse(estimand == "RR", exp(log_est + qnorm(0.025) * SE), est + qnorm(0.025) * SE),
         ci_ul = ifelse(estimand == "RR", exp(log_est + qnorm(0.975) * SE), est + qnorm(0.975) * SE),
         ese_ci_ll = ifelse(estimand == "RR", exp(log_est + qnorm(0.025) * ESE), est + qnorm(0.025) * ESE),
         ese_ci_ul = ifelse(estimand == "RR", exp(log_est + qnorm(0.975) * ESE), est + qnorm(0.975) * ESE),
         cover_init = as.numeric(ci_ll <= validated & validated <= ci_ul),
         ese_cover_init = as.numeric(ese_ci_ll <= validated & validated <= ese_ci_ul),
         procedure_fct = factor(paste0(miss_estimator, "-", outcome_tx_estimator),
                                levels = c("CGR-ipw", "CGR-adjust", "MGR-aipw"),
                                labels = c("CGR (IPTW)", "CGR (RA)", "MGR")),
         is_na = as.numeric(is.na(est)),
         scaled_est = ifelse(estimand == "RR", log(est), est)
  ) %>%
  group_by(estimand, dataset) %>%
  mutate(problematic_id = any(is.na(est))) %>%
  ungroup()


nreps_numeric <- 1000
all_results_summ <- all_results_wide %>%
  group_by(procedure_fct, estimand) %>%
  summarize(mn_est = mean(est, na.rm = TRUE),
            mdn_est = median(est, na.rm = TRUE),
            mn_bias = mean(bias_init, na.rm = TRUE),
            med_bias = median(bias_init, na.rm = TRUE),
            mean_ASE = mean(SE, na.rm = TRUE),
            med_ASE = median(SE, na.rm = TRUE),
            ESE = mean(ESE, na.rm = TRUE),
            MAD = mad(scaled_est, na.rm = TRUE),
            cover = mean(cover_init, na.rm = TRUE),
            ese_cover = mean(ese_cover_init, na.rm = TRUE),
            num_na = sum(is_na),
            mn_est_nice = mean(est[!problematic_id]),
            mdn_est_nice = median(est[!problematic_id]),
            mn_bias_nice = mean(bias_init[!problematic_id]),
            med_bias_nice = median(bias_init[!problematic_id]),
            mean_ASE_nice = mean(SE[!problematic_id]),
            med_ASE_nice = median(SE[!problematic_id]),
            MAD_nice = mad(scaled_est[!problematic_id]),
            cover_nice = mean(cover_init[!problematic_id]),
            bias_mcse = sqrt(sum((est - mean(est)) ^ 2) / (nreps_numeric * (nreps_numeric - 1))),
            .groups = "drop") %>%
  mutate(rrmse_mn = sqrt(mn_bias ^ 2 + med_ASE ^ 2),
         rrmse_mad = sqrt(med_bias ^ 2 + MAD ^ 2),
         ese_mcse = ESE / sqrt(2 * (nreps_numeric - 1)),
         rrmse_mcse = sqrt(bias_mcse^2 + ese_mcse^2),
         rrmse_cil = rrmse_mad + rrmse_mcse * qnorm(0.025),
         rrmse_ciu = rrmse_mad + rrmse_mcse * qnorm(0.975))
all_results_summ %>%
  arrange(estimand) %>%
  filter(estimand %in% c("ATE", "RR")) %>%
  print(width = Inf)

# check MCSEs
all_results_summ %>%
  arrange(estimand) %>%
  filter(estimand %in% c("ATE", "RR")) %>%
  select(procedure_fct, estimand, med_bias, bias_mcse, ESE, ese_mcse, contains("rrmse")) %>%
  mutate(abs_bias = abs(med_bias), bias_cil = abs_bias + qnorm(0.025) * bias_mcse, bias_ciu = abs_bias + qnorm(0.975) * bias_mcse,
         ese_cil = ESE + qnorm(0.025) * ese_mcse, ese_ciu = ESE + qnorm(0.975) * ese_mcse) %>%
  select(procedure_fct, estimand, contains("bias"), contains("ese"), contains("rrmse")) %>%
  print(width = Inf)

# create table with results over the replications ------------------------------
font_size <- 9
all_results_summ %>%
  filter(estimand %in% c("ATE", "RR")) %>%
  group_by(estimand) %>%
  arrange(estimand) %>%
  select(estimand, procedure_fct, mdn_est_nice, med_bias_nice,
         # med_ASE_nice,
         ESE, MAD, rrmse_mad) %>%
  rename(Estimand = estimand, Procedure = procedure_fct, `Med. est.` = mdn_est_nice,
         `Med. bias` = med_bias_nice,
         # `Med. ASE` = med_ASE_nice,
         `rRMSE` = rrmse_mad) %>%
  group_by(Estimand) %>%
  mutate(`Rel. Eff.` = rRMSE[Procedure == "MGR"] / rRMSE) %>%
  # mutate(across(
  #   where(is.numeric),
  #   .fns = ~ as.character(ifelse(abs(.x) < 0.001, "$<$ 0.001", round(.x, digits = 3))))) %>%
  knitr::kable(digits = 3, format = "latex", escape = FALSE, booktabs = TRUE,
               # align = c("l", "l", rep("r", 7)), linesep = "",
               align = c("l", "l", rep("r", 6)), linesep = "",
               caption = paste0("Performance of estimators of the effect of ",
                                "starting a protease inhibitor (PI)-containing ",
                                "combination antiretroviral therapy (ART) regimen on ",
                                "experiencing an AIDS-defining event or death within 5 years of ART initiation ",
                                "compared to an non-PI-containing regimen, ",
                                " in the VCCC data. \\label{tab:vccc_est}")) %>%
  kableExtra::kable_styling(font_size = font_size, latex_options = "hold_position") %>%
  kableExtra::footnote(general = paste0(
    "Abbreviations: ATE: average treatment effect; RR: relative risk; ",
    "CGR: conditional generalized raking (GR); ",
    "IPTW: inverse probability of treatment weighting; ",
    "RA: regression adjustment; ",
    "MGR: marginal GR; ",
    "Med.: median; est.: estimate; ",
    # "ASE: asymptotic standard error; ",
    "ESE: empirical standard error; ",
    "MAD: median absolute deviation; ",
    "rRMSE: robust root mean squared error (square root of median bias squared + MAD squared); ",
    "Rel. Eff.: relative efficiency (rRMSE relative to that of MGR).",
    " Results are based on ", length(mgr_results), " retrospectively-sampled datasets. "
  ), general_title = "", threeparttable = TRUE, escape = FALSE) %>%
  kableExtra::collapse_rows(columns = c(1)) %>%
  kableExtra::save_kable(file = paste0(plots_dir, "vccc_analysis_results.tex"))

# plot of point estimates across replications ----------------------------------
bias_boxplots <- all_results_wide %>%
  filter(!estimand %in% c("mu0", "mu1", "OR")) %>%
  ggplot(aes(x = procedure_fct, y = bias_init)) +
  geom_boxplot() +
  scale_x_discrete(guide = guide_axis(n.dodge = 2)) +
  scale_color_viridis_d(begin = 0.5, end = 0.75) +
  labs(x = "Procedure", y = "Bias") +
  facet_grid(cols = vars(estimand), scales = "free") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red")
se_boxplots <- all_results_wide %>%
  filter(!estimand %in% c("mu0", "mu1", "OR")) %>%
  ggplot(aes(x = procedure_fct, y = SE)) +
  geom_boxplot() +
  scale_x_discrete(guide = guide_axis(n.dodge = 2)) +
  scale_color_viridis_d(begin = 0.5, end = 0.75) +
  labs(x = "Procedure", y = "SE") +
  facet_grid(cols = vars(estimand), scales = "free")
bias_boxplots
se_boxplots


# description of sample --------------------------------------------------------
table(vccc_dat$ade_death_bin)
table(vccc_dat$ade_death_bin_star)
table(vccc_dat$contains_pi)
table(vccc_dat$contains_pi_star)

table_1_data <- vccc_dat %>%
  select(-contains("logvl")) %>%
  mutate(sex = as.numeric(sex == "Female"),
         year_star = year_star + 1998,
         year = year + 1998) %>%
  pivot_longer(
    cols = -c(id, age, race, risk, sex),
    names_to = c(".value", "source"),
    names_pattern = "(ade_death_bin|cd4|year|contains_pi)(?:_(star))?"
  ) %>%
  mutate(source = ifelse(source == "star", "Unvalidated", "Validated")) %>%
  select(id, age, race, risk, sex, source, year, cd4, contains_pi, ade_death_bin) %>%
  mutate(year = factor(year))

table_1 <- table_1_data %>%
  select(-id) %>%
  gtsummary::tbl_summary(
    by = "source",
    statistic = list(
      all_continuous() ~ "{mean} ({sd})",
      all_categorical() ~ "{n} ({p}%)"
    ),
    digits = all_continuous() ~ 1,
    label = list(
      ade_death_bin ~ "ADE or death within 5 years of treatment initiation",
      contains_pi ~ "Protease inhibitor-containing regimen",
      cd4 ~ "CD4 count at treatment initiation (square root)",
      year ~ "Year of treatment initiation",
      race ~ "Self-reported Black race",
      risk ~ "Self-reported injection drug use",
      sex ~ "Female",
      age ~ "Age at treatment initiation"
    )
  ) %>%
  gtsummary::modify_header(all_stat_cols() ~ "**{level}**\nN = {n}") %>%
  gtsummary::modify_spanning_header(all_stat_cols() ~ "**Source**") %>%
  gtsummary::modify_footnote(stat_1 ~ "For categorical variables, we report n (%); for continuous variables, we report mean (SD).",
                             stat_2 = NA)

table_1 %>%
  remove_bold() %>%
  as_kable_extra(format = "latex",
                 booktabs = TRUE,
                 linesep = "",
                 strip_md_bold = TRUE,
                 caption = "Descriptive statistics for the VCCC data, by source (unvalidated vs validated)",
                 label = "vccc_descriptive") %>%
  kableExtra::save_kable(file = paste0(plots_dir, "vccc_analysis_descriptive.tex"))

# table with descriptions of measurement error ---------------------------------
# number and proportion of each variable that has an error
errors_in_variables <- vccc_dat %>%
  select(-contains("logvl")) %>%
  mutate(
    across(ends_with("_star"), ~ as.numeric(.x != get(str_remove(cur_column(), "_star"))),
           .names = "error_in_{.col}")
  ) %>%
  mutate(any_error = if_any(starts_with("error_in"), ~ .x) * 1)

num_prop_errors <- errors_in_variables %>%
  select(starts_with("error_in_"), any_error) %>%
  summarize(across(everything(), list("sum" = sum, "mean" = mean), .names = "{.col}_{.fn}")) %>%
  pivot_longer(cols = everything(),
               names_to = c("variable", "summary"), values_to = "value",
               names_pattern = "(error_in_.*|any_error)_(mean|sum)") %>%
  pivot_wider(id_cols = variable,
              names_from = summary, values_from = value) %>%
  mutate(variable = case_when(
    variable == "error_in_ade_death_bin_star" ~ "AIDS-defining event or death within 5 years",
    variable == "error_in_cd4_star" ~ "CD4 count",
    variable == "error_in_year_star" ~ "Year of ART initiation",
    variable == "error_in_contains_pi_star" ~ "Initiated PI-containing ART regimen",
    variable == "any_error" ~ "Error in any variable"
  ))

num_prop_errors %>%
  mutate(n_pct = paste0(sum, " (", round(mean * 100), "\\%)")) %>%
  select(-sum, -mean) %>%
  rename(Variable = variable, `N (percent)` = n_pct) %>%
  knitr::kable(digits = 3, format = "latex", escape = FALSE, booktabs = TRUE,
               align = c("l", "r"), linesep = "",
               caption = paste0("Number and proportion of errors in variables in the VCCC data. \\label{tab:vccc_num_prop_errors}")) %>%
  kableExtra::kable_styling(font_size = font_size, latex_options = "hold_position") %>%
  kableExtra::column_spec(1, width = "30em") %>%
  kableExtra::footnote(general = paste0(
    "Abbreviations: ART: antiretroviral treatment; PI: protease inhibitor."
  ), general_title = "", threeparttable = TRUE, escape = FALSE) %>%
  kableExtra::save_kable(file = paste0(plots_dir, "vccc_analysis_num_prop_errors.tex"))


# sensitivity, specificity, PPV, NPV for the error-prone outcome and error-prone exposure
misclassification_performance <- vccc_dat %>%
  select(ade_death_bin, ade_death_bin_star, contains_pi, contains_pi_star) %>%
  summarize(
    across(ends_with("_star"),
           list("TP" = ~ sum(.x == 1 & get(str_remove(cur_column(), "_star")) == 1),
                "TN" = ~ sum(.x == 0 & get(str_remove(cur_column(), "_star")) == 0),
                "FP" = ~ sum(.x == 1 & get(str_remove(cur_column(), "_star")) == 0),
                "FN" = ~ sum(.x == 0 & get(str_remove(cur_column(), "_star")) == 1)),
           .names = "{.fn}_{.col}")
  ) %>%
  pivot_longer(
    cols = everything(), names_to = c("summary", "variable"), values_to = "value",
    names_pattern = "(TP|TN|FP|FN)_(.*)"
  ) %>%
  pivot_wider(
    id_cols = variable, names_from = summary, values_from = value
  ) %>%
  mutate(sensitivity = TP / (TP + FN),
         specificity = TN / (TN + FP),
         ppv = TP / (TP + FP),
         npv = TN / (TN + FN))

misclassification_performance %>%
  rename(Variable = variable, Sensitivity = sensitivity,
         Specificity = specificity, PPV = ppv, NPV = npv) %>%
  mutate(Variable = case_when(
    Variable == "ade_death_bin_star" ~ "AIDS-defining event or death within 5 years",
    Variable == "contains_pi_star" ~ "PI-containing ART regiment"
  )) %>%
  knitr::kable(digits = 3, format = "latex", escape = FALSE, booktabs = TRUE,
               align = c("l", rep("r", 8)), linesep = "",
               caption = paste0("Classification performance of the unvalidated variables in the VCCC data. \\label{tab:vccc_unvalidated_class_perf}")) %>%
  kableExtra::kable_styling(font_size = font_size, latex_options = "hold_position") %>%
  kableExtra::column_spec(1, width = "15em") %>%
  kableExtra::footnote(general = paste0(
    "Abbreviations: TP: true positive; TN: true negative; FP: false positive; FN: false negative; ",
    "PPV: positive predictive value; NPV: negative predictive value; ",
    "ART: antiretroviral treatment; PI: protease inhibitor."
  ), general_title = "", threeparttable = TRUE, escape = FALSE) %>%
  kableExtra::save_kable(file = paste0(plots_dir, "vccc_analysis_unvalidated_class_perf.tex"))
