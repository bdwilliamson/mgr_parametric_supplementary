# plots of results

# load required packages -------------------------------------------------------
library("readr")
library("dplyr")
library("tibble")
library("tidyr")
library("ggplot2")
library("viridis")
library("cowplot")
theme_set(theme_cowplot())
library("rprojroot") # for trajectory traversing
this_path <- normalizePath(".", mustWork = FALSE)
proj_root <- rprojroot::find_root_file(criterion = ".projectile", path = this_path)

source(paste0(proj_root, "/code/sims/01_gen_data.R"))

format_sci <- function(x, digits = 4) {
  sci <- abs(x) < 1e-4 | abs(x) >= 1e6

  out <- format(round(x, digits))

  out[sci] <- formatC(x[sci], format = "e", digits = digits)
  out <- ifelse(grepl("e", out), paste0(gsub("e", " $\\times 10^{", out, fixed = TRUE), "}$"), out)

  out
}
# read in results --------------------------------------------------------------
output_dir <- paste0(proj_root, "/results/sims/")
plots_dir <- paste0(proj_root, "/plots/")
scenarios <- as.character(c(1,2,3,4,5,6,7,8,9,10))
nreps <- "2500"
nreps_numeric <- 2500
ns <- as.character(sort(unique(c(c(250, 500, 1000, 2000, 5000, 20000), c(1000, 2000, 4000, 8000, 20000, 40000)))))
estimators <- c("CGR-IPTW", "CGR-adjust", "MGR", "AIPCW-AIPTW")
# all_est_scenarios <- expand.grid(scenario = scenarios, estimator = estimators, n = ns)
all_est_scenarios <- list.files(output_dir, recursive = TRUE)
all_est_scenarios <- all_est_scenarios[grepl("scenario", all_est_scenarios) &
                                         !grepl("archive", all_est_scenarios)]
all_output_list <- lapply(as.list(1:length(all_est_scenarios)), function(i) {
  fam_txt <- case_when(
    grepl("logistic", all_est_scenarios[i]) ~ "logistic",
    grepl("probit", all_est_scenarios[i]) ~ "probit",
    .default = "identity"
  )
  read.csv(paste0(output_dir, all_est_scenarios[i])) %>%
    mutate(link = fam_txt)
})
all_output <- as_tibble(do.call(rbind, all_output_list))

continuous_outcome_scenarios <- c(1, 3, 5, 7, 9)
binary_outcome_scenarios <- c(2, 4, 6, 8, 10)
num_scenarios <- length(c(continuous_outcome_scenarios, binary_outcome_scenarios))
# round true values to 3 significant digits
true_values <- readRDS(paste0(proj_root, "/code/sims/true_values.rds")) %>%
  filter(scenario %in% c(continuous_outcome_scenarios, binary_outcome_scenarios)) %>%
  rename(truth = value) %>%
  mutate(truth = round(truth, 3),
         estimand = ifelse(scenario %in% continuous_outcome_scenarios & estimand == "RD", "ATE", estimand)) %>%
  filter(!(scenario %in% continuous_outcome_scenarios & estimand %in% c("OR", "RR"))) %>%
  filter(!(estimand == "OR"))

# print out a table with the true values
true_values_tab <- true_values %>%
  mutate(estimand = factor(case_when(
    estimand == "OR" ~ "$OR$",
    estimand == "RR" ~ "$RR$",
    estimand == "RD" ~ "$ATE$",
    estimand == "ATE" ~ "$ATE$",
    estimand == "tx mean 1" ~ "$\\mu_1$",
    estimand == "tx mean 0" ~ "$\\mu_0$"
  ), labels = c("$\\mu_0$", "$\\mu_1$", "$ATE$", "$RR$", "$OR$"),
  levels = c("$\\mu_0$", "$\\mu_1$", "$ATE$", "$RR$", "$OR$"), ordered = TRUE)) %>%
  mutate(scenario2 = case_when(
    scenario == 1 ~ "1a",
    scenario == 2 ~ "1b",
    scenario == 3 ~ "2a",
    scenario == 4 ~ "2b",
    scenario == 5 ~ "3a",
    scenario == 6 ~ "3b",
    scenario == 7 ~ "4a",
    scenario == 8 ~ "4b",
    # scenario == 12 ~ "5",
    scenario == 9 ~ "5a",
    scenario == 10 ~ "5b"
  )) %>%
  mutate(scenario = as.character(scenario)) %>%
  filter(!grepl("OR", estimand)) %>%
  arrange(scenario2) %>%
  pivot_wider(id_cols = "estimand", names_from = "scenario2", values_from = "truth") %>%
  arrange(estimand) %>%
  mutate(across(-estimand, .fns = ~ ifelse(is.na(.x), "", round(.x, 2)))) %>%
  rename(Estimand = estimand)
true_values_tab %>%
  knitr::kable(format = "latex", digits = 3, booktabs = TRUE, escape = FALSE,
               align = c("l", rep("r", num_scenarios)),
               caption = paste0("True values of estimands in the simulation scenarios.",
                                " Scenarios ending in ``a'' involve a continuous outcome,",
                                " so the $RR$ is not defined.",
                                "\\label{tab:true_values_all_scenarios}"),
               linesep = "") %>%
  kableExtra::add_header_above(c(" " = 1, "Scenario" = num_scenarios)) %>%
  kableExtra::kable_styling(font_size = 11) %>%
  kableExtra::save_kable(file = paste0(plots_dir, "sim_true_values_all_scenarios.tex"))
# version with only scenarios 1 and 2
true_values_tab %>%
  # select(-`3`, -`4`) %>%
  knitr::kable(format = "latex", digits = 3, booktabs = TRUE, escape = FALSE,
               align = c("l", rep("r", num_scenarios)),
               caption = paste0("True values of estimands in the simulation scenarios.",
                                " Scenarios ending in ``a'' involve a continuous outcome,",
                                " so the $RR$ is not defined.",
                                "\\label{tab:true_values}"),
               linesep = "") %>%
  kableExtra::add_header_above(c(" " = 1, "Scenario" = num_scenarios)) %>%
  kableExtra::kable_styling(font_size = 11) %>%
  kableExtra::save_kable(file = paste0(plots_dir, "sim_true_values.tex"))

# summarize over monte-carlo ID
ese <- all_output %>%
  group_by(scenario, seed, n, miss_estimator, outcome_tx_estimator, link, estimand) %>%
  mutate(scaled_est = ifelse(!(estimand %in% c("RD", "ATE", "mu1", "mu0")), log(est), est)) %>%
  summarize(ESE = sd(scaled_est, na.rm = TRUE), .groups = "drop")
expanded_output <- all_output %>%
  left_join(true_values %>% mutate(estimand = case_when(
    estimand == "RD" ~ "ATE",
    estimand == "tx mean 0" ~ "mu0",
    estimand == "tx mean 1" ~ "mu1",
    .default = estimand
  )), by = c("scenario", "estimand")) %>%
  left_join(ese, by = c("scenario", "seed", "n", "miss_estimator",
                        "outcome_tx_estimator", "link", "estimand")) %>%
  mutate(bias_init = est - truth,
         log_est = log(est),
         scaled_est = ifelse(estimand == "RR", log(est), est),
         ci_ll = ifelse(estimand == "RD" | estimand == "ATE", est + SE * qnorm(0.025),
                        exp(log_est + SE * qnorm(0.025))),
         ci_ul = ifelse(estimand == "RD" | estimand == "ATE", est + SE * qnorm(0.975),
                        exp(log_est + SE * qnorm(0.975))),
         ese_ci_ll = ifelse(estimand == "RD" | estimand == "ATE", est + ESE * qnorm(0.025),
                            exp(log_est + ESE * qnorm(0.025))),
         ese_ci_ul = ifelse(estimand == "RD" | estimand == "ATE", est + ESE * qnorm(0.975),
                            exp(log_est + ESE * qnorm(0.975))),
         cover_init = as.numeric(ci_ll <= truth & truth <= ci_ul),
         ese_cover_init = as.numeric(ese_ci_ll <= truth & truth <= ese_ci_ul))


nice_output <- expanded_output %>%
  mutate(
    n_fct = factor(n, levels = ns, labels = ns, ordered = TRUE),
    procedure = case_when(
    miss_estimator == "CGR" & outcome_tx_estimator == "adjust" ~ "CGR (RA)",
    miss_estimator == "CGR" & outcome_tx_estimator != "adjust" ~ "CGR (IPTW)",
    miss_estimator == "MGR" ~ "MGR",
    .default = paste0(miss_estimator, "-", outcome_tx_estimator)
  ),
  procedure_fct = factor(procedure,
                         levels = c("CGR", "CGR (RA)", "CGR (IPTW)", "MGR", "AIPCW-AIPTW"),
                         labels = c("CGR", "CGR (RA)", "CGR (IPTW)", "MGR", "AIPCW-AIPTW")))
summary_tib <- expanded_output %>%
  group_by(scenario, seed, n, outcome_tx_estimator, miss_estimator, link, estimand) %>%
  summarize(mn_bias = mean(bias_init, na.rm = TRUE),
            med_bias = median(bias_init, na.rm = TRUE),
            cover = mean(cover_init, na.rm = TRUE),
            ese_cover = mean(ese_cover_init, na.rm = TRUE),
            ASE = median(SE, na.rm = TRUE),
            ESE = median(ESE, na.rm = TRUE),
            MAD = mad(scaled_est, na.rm = TRUE),
            bias_to_ase = abs(mn_bias) / ASE,
            bias_to_ese = abs(mn_bias) / ESE,
            bias_mcse = sqrt(sum((est - mean(est)) ^ 2) / (nreps_numeric * (nreps_numeric - 1))),
            .groups = "drop") %>%
  mutate(rmse = sqrt(mn_bias ^ 2 + ESE ^ 2),
         rrmse = sqrt(med_bias ^ 2 + MAD ^ 2),
         ese_mcse = ESE / sqrt(2 * (nreps_numeric - 1)),
         rrmse_mcse = sqrt(bias_mcse^2 + ese_mcse^2),
         n_fct = factor(n, levels = ns, labels = ns, ordered = TRUE),
         procedure = case_when(
           miss_estimator == "CGR" & outcome_tx_estimator == "adjust" ~ "CGR (RA)",
           miss_estimator == "CGR" & outcome_tx_estimator != "adjust" ~ "CGR (IPTW)",
           miss_estimator == "MGR" ~ "MGR",
           .default = paste0(miss_estimator, "-", outcome_tx_estimator)
         ),
         procedure_fct = factor(procedure,
                                levels = c("CGR", "CGR (RA)", "CGR (IPTW)", "MGR", "AIPCW-AIPTW"),
                                labels = c("CGR", "CGR (RA)", "CGR (IPTW)", "MGR", "AIPCW-AIPTW")),
         n_eff = case_when(
           scenario %in% binary_outcome_scenarios & n == "2000" ~ 300,
           scenario %in% binary_outcome_scenarios & n == "4000" ~ 600,
           scenario %in% binary_outcome_scenarios & n == "8000" ~ 1200,
           scenario %in% binary_outcome_scenarios & n == "20000" ~ 3000,
           scenario %in% binary_outcome_scenarios & n == "40000" ~ 6000,
           scenario %in% continuous_outcome_scenarios & n == "500" ~ 300,
           scenario %in% continuous_outcome_scenarios & n == "1000" ~ 600,
           scenario %in% continuous_outcome_scenarios & n == "2000" ~ 1200,
           scenario %in% continuous_outcome_scenarios & n == "5000" ~ 3000,
           scenario %in% continuous_outcome_scenarios & n == "10000" ~ 6000
         ),
         outcome_rate = ifelse(scenario %in% c(1), 1, 0.26),
         one_over_n_eff = 1 / (n_eff * outcome_rate * (1 - 0.56)),
         n_eff_fct = factor(n_eff)) %>%
  filter(!(estimand %in% c("mu0", "mu1", "OR")))

long_summary_tib <- summary_tib %>%
  select(-contains("estimator")) %>%
  select(scenario, seed, n, n_fct, n_eff, one_over_n_eff, n_eff_fct, procedure,
         link, estimand, everything()) %>%
  pivot_longer(mn_bias:rmse, names_to = "measure", values_to = "value") %>%
  arrange(scenario, seed, n, estimand, measure, procedure, link) %>%
  group_by(scenario, seed, n, estimand, measure, link) %>%
  mutate(rel_value = value / first(value)) %>%
  ungroup() %>%
  mutate(scenario2 = case_when(
    scenario == 1 ~ "1a",
    scenario == 2 ~ "1b",
    scenario == 3 ~ "2a",
    scenario == 4 ~ "2b",
    scenario == 5 ~ "3a",
    scenario == 6 ~ "3b",
    scenario == 7 ~ "4a",
    scenario == 8 ~ "4b",
    # scenario == 12 ~ "5",
    scenario == 9 ~ "5a",
    scenario == 10 ~ "5b"
  ))

long_summary_tib2 <- long_summary_tib %>%
  filter(scenario %in% continuous_outcome_scenarios |
           ((scenario %in% binary_outcome_scenarios))) %>%
  filter(!(scenario == 1 & n == 20000))

# check bias of each estimator as a boxplot
bias_boxplots <- nice_output %>%
  filter(!estimand %in% c("mu0", "mu1", "OR"),
         !(scenario %in% continuous_outcome_scenarios & estimand == "RR")) %>%
  ggplot(aes(x = n_fct, color = procedure_fct, y = bias_init)) +
  geom_boxplot() +
  scale_x_discrete(guide = guide_axis(n.dodge = 2)) +
  scale_color_viridis_d(begin = 0.5, end = 0.75) +
  labs(x = "n", y = "Bias", color = "Estimator") +
  facet_grid(rows = vars(estimand), cols = vars(scenario), scales = "free") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red")
se_boxplots <- nice_output %>%
  filter(!estimand %in% c("mu0", "mu1", "OR")) %>%
  filter(!(scenario == 1 & estimand == "RR")) %>%
  ggplot(aes(x = n_fct, color = procedure_fct, y = SE)) +
  geom_boxplot() +
  scale_x_discrete(guide = guide_axis(n.dodge = 2)) +
  scale_color_viridis_d(begin = 0.5, end = 0.75) +
  labs(x = "n", y = "SE", color = "Estimator") +
  facet_grid(rows = vars(estimand), cols = vars(scenario), scales = "free")
bias_boxplots
se_boxplots

# create plots -----------------------------------------------------------------
bias_plot <- long_summary_tib2 %>%
  filter(measure == "med_bias") %>%
  ggplot(aes(x = n_fct, shape = procedure_fct, y = value)) +
  geom_point(position = position_dodge(width = 0.2)) +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
  scale_x_discrete(guide = guide_axis(n.dodge = 2)) +
  labs(x = "n", y = "Median bias", shape = "Estimator") +
  facet_grid(rows = vars(estimand), cols = vars(scenario), scales = "free")

bias_to_se_plot <- long_summary_tib2 %>%
  filter(grepl("bias_to_", measure)) %>%
  mutate(measure = ifelse(measure == "bias_to_ase", "ASE", "ESE")) %>%
  ggplot(aes(x = n_fct, shape = procedure_fct, y = value, color = measure)) +
  geom_point(position = position_dodge(width = 0.2)) +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
  scale_x_discrete(guide = guide_axis(n.dodge = 2)) +
  labs(x = "n", y = "Bias to SE ratio", shape = "Estimator", color = "SE type") +
  facet_grid(rows = vars(estimand), cols = vars(scenario), scales = "free")

se_plot <- long_summary_tib2 %>%
  filter(measure %in% c("ASE", "ESE")) %>%
  ggplot(aes(x = n_fct, shape = procedure_fct, y = value, color = measure)) +
  geom_point(position = position_dodge(width = 0.4)) +
  scale_x_discrete(guide = guide_axis(n.dodge = 2)) +
  labs(x = "n", y = "Standard error", shape = "Estimator", color = "SE type") +
  facet_grid(rows = vars(estimand), cols = vars(scenario), scales = "free")

# plot including link function, to check
se_plot_with_link <- long_summary_tib %>%
  filter(measure %in% c("ASE", "ESE")) %>%
  mutate(procedure2 = paste0(procedure, "_", link),
         procedure_fct2 = factor(procedure2,
                                 levels = c("CGR_identity", "MGR_identity",
                                            "CGR_logistic", "MGR_logistic",
                                            "CGR_probit", "MGR_probit", "AIPCW-AIPTW_probit"),
                                 labels = c("CGR (identity)", "MGR (identity)",
                                            "CGR (logistic)", "MGR (logistic)",
                                            "CGR (probit)", "MGR (probit)",
                                            "AIPCW-AIPTW (probit)"))) %>%
  ggplot(aes(x = n_fct, shape = procedure_fct2, y = value, color = measure)) +
  geom_point(position = position_dodge(width = 0.4)) +
  scale_shape_manual(values = c(15:19, 3, 4, 6)) +
  scale_x_discrete(guide = guide_axis(n.dodge = 2)) +
  labs(x = "n", y = "Standard error", shape = "Estimator", color = "SE type") +
  facet_grid(rows = vars(estimand), cols = vars(scenario), scales = "free")

rmse_plot <- long_summary_tib2 %>%
  filter(measure == "rmse") %>%
  ggplot(aes(x = n_fct, shape = procedure_fct, y = value)) +
  geom_point(position = position_dodge(width = 0.2)) +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
  scale_x_discrete(guide = guide_axis(n.dodge = 2)) +
  labs(x = "n", shape = "Estimation procedure", y = "RMSE") +
  facet_grid(rows = vars(estimand), cols = vars(scenario), scales = "free")

cover_plot <- long_summary_tib2 %>%
  filter(grepl("cover", measure)) %>%
  mutate(measure = ifelse(measure == "cover", "ASE", "ESE")) %>%
  ggplot(aes(x = n_fct, shape = procedure_fct, y = value, color = measure)) +
  geom_point(position = position_dodge(width = 0.4)) +
  geom_hline(yintercept = 0.95, color = "red", linetype = "dashed") +
  ylim(c(0.6, 1)) +
  scale_x_discrete(guide = guide_axis(n.dodge = 2)) +
  labs(x = "n", shape = "Estimation procedure", y = "Confidence interval coverage", color = "SE type") +
  facet_grid(rows = vars(estimand), cols = vars(scenario), scales = "free")

combined_legend <- get_plot_component(
  se_plot + theme(legend.direction = "horizontal", legend.position = "bottom"),
  "guide-box-bottom", return_all = TRUE
)

combined_plot <- plot_grid(
  plot_grid(plot_grid(
    bias_plot + theme(legend.position = "none",
                      panel.border = element_rect(fill = "transparent", color = "black", linewidth = 0.5)),
    rmse_plot + theme(legend.position = "none",
                      panel.border = element_rect(fill = "transparent", color = "black", linewidth = 0.5)),
    nrow = 2, ncol = 1
  ),
  plot_grid(
    se_plot + theme(legend.position = "none",
                    panel.border = element_rect(fill = "transparent", color = "black", linewidth = 0.5)),
    cover_plot + theme(legend.position = "none",
                       panel.border = element_rect(fill = "transparent", color = "black", linewidth = 0.5)),
    nrow = 2, ncol = 1
  ), nrow = 1, ncol = 2),
  combined_legend, rel_heights = c(1, .1), nrow = 2, ncol = 1
)

ggsave(filename = paste0(plots_dir, "/sim_combined_performance_plot.png"),
       combined_plot, dpi = 300, width = 8, height = 4.25, units = "in")
ggsave(filename = paste0(plots_dir, "/sim_bias_to_se_plot.png"),
       bias_to_se_plot, dpi = 300, width = 4, height = 2, units = "in")

# also have a series of plots and tables that are separate by scenario, for presentations
font_size <- 9
theme_set(theme_cowplot(font_size = font_size))
abbreviation_text_prefix <- "Abbreviations: "
estimand_txt <- "ATE: average treatment effect; "
abbreviation_text_suffix <- paste0("Est: Estimator; ", "Med: Median; ",
                                   "ASE: median asymptotic standard error; ",
                                   "ESE: empirical standard error; ",
                                   "MAD: mean absolute deviation; ",
                                   "rRMSE: robust root mean squared error (RMSE) (square root of med bias squared + MAD squared); ",
                                   "RMSE: RMSE (mean bias squared + ESE squared); ",
                                   "Rel. Eff.: relative efficiency (rRMSE relative to that of MGR); ",
                                   "Cover: coverage based on the ASE; ",
                                   "ESE Cover: coverage based on the ESE; ",
                                   "CGR: conditional generalized raking (GR); ",
                                   "IPTW: inverse probability of treatment weighting; ",
                                   "RA: regression adjustment; ",
                                   "MGR: marginal GR; ",
                                   "AIPCW: augmented inverse probability of coarsening weighting; ",
                                   "AIPTW: augmented IPTW.")

abbreviation_text_suffix_cgr_iptw <- paste0("Med: Median; ",
                                            "ASE: median asymptotic standard error; ",
                                            "ESE: empirical standard error; ",
                                            "MAD: mean absolute deviation; ",
                                            "rRMSE: robust root mean squared error (RMSE) (square root of med bias squared + MAD squared); ",
                                            "RMSE: RMSE (mean bias squared + ESE squared); ",
                                            "Rel. Eff.: relative efficiency (rRMSE relative to that of MGR); ",
                                            "Cover: coverage based on the ASE; ",
                                            "ESE Cover: coverage based on the ESE; ",
                                            "CGR: conditional generalized raking (GR); ",
                                            "IPTW: inverse probability of treatment weighting.")

# main manuscript tables: do not include CGR (IPTW)
# supplemental tables do include CGR (IPTW)
# also create a CGR (IPTW) table for the supplement
main_scenarios <- c(4, 6, 8)
for (i in 1:length(scenarios)) {
  this_scenario <- scenarios[i]
  this_tib <- long_summary_tib2 %>%
    filter(scenario == this_scenario, !(estimand %in% c("mu0", "mu1", "OR")))
  if (i %in% continuous_outcome_scenarios) {
    this_tib <- this_tib %>%
      filter(!(estimand %in% c("RR")))
  } else {
    estimand_txt <- paste0(estimand_txt,
                           "RR: relative risk; ")
  }
  final_scenario <- this_tib$scenario2[1]

  bias_plot <- this_tib %>%
    filter(measure == "med_bias") %>%
    ggplot(aes(x = n_fct, shape = procedure_fct, y = value)) +
    geom_point(position = position_dodge(width = 0.2)) +
    geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
    scale_x_discrete(guide = guide_axis(n.dodge = 2)) +
    labs(x = "n", y = "Median bias", shape = "Estimator") +
    facet_grid(rows = vars(estimand), scales = "free_y")

  bias_to_se_plot <- this_tib %>%
    filter(grepl("bias_to_", measure)) %>%
    mutate(measure = ifelse(measure == "bias_to_ase", "ASE", "ESE")) %>%
    ggplot(aes(x = n_fct, shape = procedure_fct, y = value, color = measure)) +
    geom_point(position = position_dodge(width = 0.2)) +
    geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
    scale_x_discrete(guide = guide_axis(n.dodge = 2)) +
    labs(x = "n", y = "Bias to SE ratio", shape = "Estimator", color = "SE type") +
    facet_grid(rows = vars(estimand), scales = "free_y")

  se_plot <- this_tib %>%
    filter(measure %in% c("ASE", "ESE")) %>%
    ggplot(aes(x = n_fct, shape = procedure_fct, y = value, color = measure)) +
    geom_point(position = position_dodge(width = 0.4)) +
    scale_x_discrete(guide = guide_axis(n.dodge = 2)) +
    labs(x = "n", y = "Standard error", shape = "Estimator", color = "SE type") +
    facet_grid(rows = vars(estimand), scales = "free_y")

  rmse_plot <- this_tib %>%
    filter(measure == "rmse") %>%
    ggplot(aes(x = n_fct, shape = procedure_fct, y = value)) +
    geom_point(position = position_dodge(width = 0.2)) +
    geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
    scale_x_discrete(guide = guide_axis(n.dodge = 2)) +
    labs(x = "n", shape = "Estimation procedure", y = "RMSE") +
    facet_grid(rows = vars(estimand), scales = "free_y")

  cover_plot <- this_tib %>%
    filter(grepl("cover", measure)) %>%
    mutate(measure = ifelse(measure == "cover", "ASE", "ESE")) %>%
    ggplot(aes(x = n_fct, shape = procedure_fct, y = value, color = measure)) +
    geom_point(position = position_dodge(width = 0.4)) +
    geom_hline(yintercept = 0.95, color = "red", linetype = "dashed") +
    ylim(c(0.6, 1)) +
    scale_x_discrete(guide = guide_axis(n.dodge = 2)) +
    labs(x = "n", shape = "Estimation procedure", y = "Confidence interval coverage", color = "SE type") +
    facet_grid(rows = vars(estimand), scales = "free_y")

  combined_legend <- get_plot_component(
    se_plot + theme(legend.direction = "horizontal", legend.position = "bottom"),
    "guide-box-bottom", return_all = TRUE
  )
  overall_title <- ggdraw() +
    draw_label(paste0("Scenario ", final_scenario),
               fontface = "bold", x = 0, hjust = 0) +
    theme(plot.margin = margin(0, 0, 0, 7))

  combined_plot <- plot_grid(
    overall_title,
    plot_grid(plot_grid(
      bias_plot + theme(legend.position = "none",
                        panel.border = element_rect(fill = "transparent", color = "black", linewidth = 0.5)),
      rmse_plot + theme(legend.position = "none",
                        panel.border = element_rect(fill = "transparent", color = "black", linewidth = 0.5)),
      nrow = 2, ncol = 1
    ),
    plot_grid(
      se_plot + theme(legend.position = "none",
                      panel.border = element_rect(fill = "transparent", color = "black", linewidth = 0.5)),
      cover_plot + theme(legend.position = "none",
                         panel.border = element_rect(fill = "transparent", color = "black", linewidth = 0.5)),
      nrow = 2, ncol = 1
    ), nrow = 1, ncol = 2),
    combined_legend, rel_heights = c(.1, 1, .1), nrow = 3, ncol = 1
  )

  ggsave(filename = paste0(plots_dir, "/sim_combined_performance_plot_s", this_scenario, ".png"),
         combined_plot, dpi = 300, width = 8, height = 4.25, units = "in")
  ggsave(filename = paste0(plots_dir, "/sim_bias_to_se_plot_s", this_scenario, ".png"),
         plot_grid(overall_title, bias_to_se_plot, rel_heights = c(0.1, 1), nrow = 2, ncol = 1),
         dpi = 300, width = 4, height = 2, units = "in")

  # also tables
  summ_tib <- this_tib %>%
    ungroup() %>%
    pivot_wider(id_cols = scenario:procedure_fct, names_from = "measure",
                values_from = c("value", "rel_value")) %>%
    rename_with(everything(), .fn = ~ gsub("value_", "", .x))
  final_summ_tib <- summ_tib %>%
    select(estimand, procedure, n,
           # one_over_n_eff,
           med_bias, ASE, ESE,
           MAD, rrmse, rmse,
           # bias_to_ase,
           cover, ese_cover) %>%
    rename(`Est.` = procedure, `Estimand` = estimand,
           # `1/$n_\\text{eff}$` = one_over_n_eff,
           `Med. bias` = med_bias,
           # `bias/SE` = bias_to_ase,
           `rRMSE` = rrmse,
           `RMSE` = rmse,
           Cover = cover,
           `ESE Cover` = ese_cover) %>%
    group_by(Estimand, n) %>%
    mutate(`Rel. Eff.` = `rRMSE`[`Est.` == "MGR"] / `rRMSE`,
           # `Rel. Eff. 2` = `RMSE`[`Est.` == "MGR"] / `RMSE`,
           .before = "Cover") %>%
    ungroup()
  estimands <- unique(final_summ_tib$Estimand)

  # describe the scenario
  scenario_descr <- paste0(" a ", ifelse(this_scenario %in% continuous_outcome_scenarios, "continuous", "binary"),
                           " outcome setting with approximately ",
                           ifelse(this_scenario %in% c(1, 2), 20, 50), "\\% missing data and ",
                           case_when(
                             this_scenario %in% 1:4 ~ "correctly specified nuisance models",
                             this_scenario %in% 5:6 ~ "a misspecified treatment propensity model",
                             this_scenario %in% 7:8 ~ "a misspecified outcome regression model",
                             this_scenario %in% 9:10 ~ "a misspecified missing-data model",
                             this_scenario %in% 12 ~ "the wrong outcome-model link function"
                           ), ".")
  # if a main manuscript table, remove CGR (IPTW) and create a separate table
  for (j in 1:length(estimands)) {
    if (this_scenario %in% main_scenarios) {
      this_final_tib <- final_summ_tib %>%
        filter(Estimand == estimands[j], `Est.` != "CGR (IPTW)") %>%
        select(-Estimand)
    } else {
      this_final_tib <- final_summ_tib %>%
        filter(Estimand == estimands[j]) %>%
        select(-Estimand)
    }
     this_final_tib %>%
      mutate(`Med. bias` = ifelse(abs(`Med. bias`) < 0.0001, "$<$ 0.0001", format_sci(`Med. bias`, digits = 4))) %>%
      # mutate(across(
        # where(is.numeric),
        # .fns = ~ ifelse(abs(.x) < 0.0001, "$<$ 0.0001", format_sci(.x, digits = 4)))) %>%
        # .fns = ~ ifelse(abs(.x) < 0.0001, "$<$ 0.0001", round(.x, digits = 4)))) %>%
      knitr::kable(digits = 3, format = "latex", escape = FALSE, booktabs = TRUE,
                   align = c("l", "l", rep("r", 9)), linesep = "",
                   caption = paste0("Performance of estimators of the ", estimands[j],
                                    " in Scenario ", final_scenario, ",",
                                    scenario_descr,
                                    "\\label{tab:sim_combined_performance_tab_s", this_scenario, "_", estimands[j], "}")) %>%
      kableExtra::kable_styling(font_size = font_size, latex_options = "hold_position") %>%
      kableExtra::footnote(general = paste0(
        abbreviation_text_prefix, estimand_txt, abbreviation_text_suffix,
        " Results are based on 2500 Monte-Carlo replications. "
      ), general_title = "", threeparttable = TRUE, escape = FALSE) %>%
      kableExtra::column_spec(column = 1, width = "3.5em") %>%
      kableExtra::collapse_rows(columns = 2) %>%
      kableExtra::save_kable(file = paste0(plots_dir, "sim_combined_performance_tab_s", this_scenario, "_", estimands[j], ".tex"))
  }
  # one table for all estimands, CGR (IPTw) in main manuscript scenarios
  if (this_scenario %in% main_scenarios) {
    cgr_iptw_tib <- final_summ_tib %>%
      filter(`Est.` == "CGR (IPTW)") %>%
      select(-`Est.`)
    cgr_iptw_tib %>%
      arrange(Estimand, n) %>%
      mutate(`Med. bias` = ifelse(abs(`Med. bias`) < 0.0001, "$<$ 0.0001", format_sci(`Med. bias`, digits = 4))) %>%
      knitr::kable(digits = 3, format = "latex", escape = FALSE, booktabs = TRUE,
                   align = c("l", "l", rep("r", 9)), linesep = "",
                   caption = paste0("Performance of CGR (IPTW) estimators of the ATE and RR",
                                    " in Scenario ", final_scenario, ",",
                                    scenario_descr,
                                    "\\label{tab:sim_combined_performance_tab_s", this_scenario, "_cgr_iptw}")) %>%
      kableExtra::kable_styling(font_size = font_size, latex_options = "hold_position") %>%
      kableExtra::footnote(general = paste0(
        abbreviation_text_prefix, abbreviation_text_suffix_cgr_iptw,
        " Results are based on 2500 Monte-Carlo replications. "
      ), general_title = "", threeparttable = TRUE, escape = FALSE) %>%
      kableExtra::column_spec(column = 1, width = "3.5em") %>%
      kableExtra::collapse_rows(columns = 2) %>%
      kableExtra::save_kable(file = paste0(plots_dir, "sim_combined_performance_tab_s", this_scenario, "_cgr_iptw.tex"))
  }
}

# table with all regression coefficients for data-generating models
all_outcome_coefs <- vector("list", length = length(scenarios))
all_missing_coefs <- vector("list", length = length(scenarios))
x <- gen_x()
for (i in 1:length(scenarios)) {
  this_scenario <- as.numeric(scenarios[i])
  these_outcome_regression_coefs <- gen_y(data = x, scenario = this_scenario)$coefs
  if (nrow(these_outcome_regression_coefs) != (length(names(x)) + 3)) {
    other_coefs <- rep(0, length(names(x)) - nrow(these_outcome_regression_coefs) + 3)
    these_nms <- names(x)[!(names(x) %in% these_outcome_regression_coefs$name)]
    these_outcome_regression_coefs <- rbind.data.frame(these_outcome_regression_coefs,
                                                       data.frame(name = these_nms, value = other_coefs))
  }
  these_missing_data_coefs <- gen_miss(data = x, scenario = this_scenario)$coefs
  # miss_nms <- c(names(x)[!(names(x) %in% c("Ww", "Ws"))], "Y")
  miss_nms <- c(names(x), "Y")
  if (nrow(these_missing_data_coefs) != (length(miss_nms) + 1)) {
    other_coefs <- rep(0, length(miss_nms) - nrow(these_missing_data_coefs) + 1)
    these_nms <- miss_nms[!(miss_nms %in% these_missing_data_coefs$name)]
    these_missing_data_coefs <- rbind.data.frame(these_missing_data_coefs,
                                                       data.frame(name = these_nms, value = other_coefs))
  }
  all_outcome_coefs[[i]] <- data.frame(scenario = this_scenario, these_outcome_regression_coefs)
  all_missing_coefs[[i]] <- data.frame(scenario = this_scenario, these_missing_data_coefs)
}

all_outcome_coefs_df <- do.call(rbind.data.frame, all_outcome_coefs) %>%
  mutate(scenario = as.character(scenario),
         scenario2 = case_when(
           scenario == 1 ~ "1a",
           scenario == 2 ~ "1b",
           scenario == 3 ~ "2a",
           scenario == 4 ~ "2b",
           scenario == 5 ~ "3a",
           scenario == 6 ~ "3b",
           scenario == 7 ~ "4a",
           scenario == 8 ~ "4b",
           # scenario == 12 ~ "5",
           scenario == 9 ~ "5a",
           scenario == 10 ~ "5b"
         )) %>%
  filter(!(name %in% c("X2", "Zs2", "ZsZw")))
all_missing_coefs_df <- do.call(rbind.data.frame, all_missing_coefs) %>%
  mutate(scenario2 = case_when(
    scenario == 1 ~ "1a",
    scenario == 2 ~ "1b",
    scenario == 3 ~ "2a",
    scenario == 4 ~ "2b",
    scenario == 5 ~ "3a",
    scenario == 6 ~ "3b",
    scenario == 7 ~ "4a",
    scenario == 8 ~ "4b",
    # scenario == 12 ~ "5",
    scenario == 9 ~ "5a",
    scenario == 10 ~ "5b"
  )) %>%
  filter(!(name %in% c("X2", "ZsZw")))

all_tx_coefs_df <- rbind.data.frame(
  data.frame(
    scenario = "1", scenario2 = "1a",
    name = c("Intercept", "Zs", "Zw", "Ws", "Ww", "Us", "Uw", "As", "Aw"),
    value = c(rep(0.1, 5), rep(0, 4))
  ),
  data.frame(
    scenario = "1", scenario2 = "1b",
    name = c("Intercept", "Zs", "Zw", "Ws", "Ww", "Us", "Uw", "As", "Aw"),
    value = c(rep(0.1, 5), rep(0, 4))
  ),
  data.frame(
    scenario = "3", scenario2 = "2a",
    name = c("Intercept", "Zs", "Zw", "Ws", "Ww", "Us", "Uw", "As", "Aw"),
    value = c(rep(0.1, 5), rep(0, 2), c(-0.4, -0.2))
  ),
  data.frame(
    scenario = "4", scenario2 = "2b",
    name = c("Intercept", "Zs", "Zw", "Ws", "Ww", "Us", "Uw", "As", "Aw"),
    value = c(rep(0.1, 5), rep(0, 2), c(-0.4, -0.2))
  ),
  data.frame(
    scenario = "5", scenario2 = "3a",
    name = c("Intercept", "Zs", "Zw", "Ws", "Ww", "Us", "Uw", "As", "Aw"),
    value = c(rep(0.1, 5), c(-1, -1, -0.4, -0.2))
  ),
  data.frame(
    scenario = "6", scenario2 = "3b",
    name = c("Intercept", "Zs", "Zw", "Ws", "Ww", "Us", "Uw", "As", "Aw"),
    value = c(rep(0.1, 5), c(-1, -1, -0.4, -0.2))
  ),
  data.frame(
    scenario = "7", scenario2 = "4a",
    name = c("Intercept", "Zs", "Zw", "Ws", "Ww", "Us", "Uw", "As", "Aw"),
    value = c(rep(0.1, 5), c(-1, -1, -0.4, -0.2))
  ),
  data.frame(
    scenario = "8", scenario2 = "4b",
    name = c("Intercept", "Zs", "Zw", "Ws", "Ww", "Us", "Uw", "As", "Aw"),
    value = c(rep(0.1, 5), c(-1, -1, -0.4, -0.2))
  ),
  data.frame(
    scenario = "9", scenario2 = "6a",
    name = c("Intercept", "Zs", "Zw", "Ws", "Ww", "Us", "Uw", "As", "Aw"),
    value = c(rep(0.1, 5), c(-1, -1, -0.4, -0.2))
  ),
  data.frame(
    scenario = "10", scenario2 = "6b",
    name = c("Intercept", "Zs", "Zw", "Ws", "Ww", "Us", "Uw", "As", "Aw"),
    value = c(rep(0.1, 5), c(-1, -1, -0.4, -0.2))
  ),
  data.frame(
    scenario = "12", scenario2 = "5",
    name = c("Intercept", "Zs", "Zw", "Ws", "Ww", "Us", "Uw", "As", "Aw"),
    value = c(rep(0.1, 5), rep(0, 2), c(-0.4, -0.2))
  )
)

outcome_tx_missing_coefs <- rbind.data.frame(
  all_outcome_coefs_df %>% mutate(model = "Outcome regression"),
  all_tx_coefs_df %>% mutate(model = "Propensity score"),
  all_missing_coefs_df %>% mutate(model = "Missing data")
) %>%
  group_by(scenario2, name, model) %>%
  slice(1) %>%
  ungroup() %>%
  arrange(scenario2, factor(name, levels = c("Intercept", "X", "Zs", "Zw", "Ws", "Ww", "Us", "Uw", "As", "Aw", "Y"))) %>%
  mutate(name2 = gsub("A", "C", gsub("s", "1", gsub("w", "2", name))))


# make a wide version
wide_coefs <- outcome_tx_missing_coefs %>%
  pivot_wider(id_cols = c(scenario2, model), names_from = name2, values_from = value) %>%
  rename(Scenario = scenario2, Model = model)

# create the table
wide_coefs %>%
  mutate(across(everything(), .fns = ~ ifelse(is.na(.x), "---", .x))) %>%
  knitr::kable(format = "latex", digits = 3, booktabs = TRUE, escape = FALSE,
               caption = paste0("Regression coefficients specifying the data-generating ",
                                "models for all simulation scenarios. ",
                                "\\label{tab:sim_model_coefs}"),
               linesep = "") %>%
  kableExtra::add_header_above(c(" " = 2, "Coefficient" = 11)) %>%
  kableExtra::collapse_rows(columns = 1) %>%
  kableExtra::kable_styling(font_size = 11) %>%
  kableExtra::save_kable(file = paste0(plots_dir, "sim_model_coefs.tex"))

