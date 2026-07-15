# simulation-specific utility functions

# get true values for a simulation
# @param scenario the data-generating scenario
# @param sample_size the sample size of the dataset to generate
# @param nreps the number of replications to generate at sample size n
# @param lowcor, midcor, highcor, gencor: correlation, to pass to DGM
get_true_values <- function(scenario = 1, sample_size = 1e6, nreps = 10,
                            lowcor = 0.2, midcor = 0.4, highcor = 0.7, gencor = 0) {
  # generate dataset, get true values for each of 1:nreps
  mu1_vec <- vector("numeric", length = nreps)
  mu0_vec <- vector("numeric", length = nreps)
  rr_vec <- vector("numeric", length = nreps)
  rd_vec <- vector("numeric", length = nreps)
  or_vec <- vector("numeric", length = nreps)
  for (i in 1:nreps) {
    # generate dataset
    df <- gen_data(N = sample_size, scenario = scenario, gencor = gencor)
    # new datasets setting X = 0, 1 for everyone
    newdf_0 <- df
    newdf_0$X <- 0
    newdf_0$Y <- gen_y(data = newdf_0 %>% select(-Y), scenario = scenario)$Y
    newdf_1 <- df
    newdf_1$X <- 1
    newdf_1$Y <- gen_y(data = newdf_1 %>% select(-Y), scenario = scenario)$Y

    # oracle values of treatment-specific means
    mu1 <- mean(newdf_1$Y)
    mu0 <- mean(newdf_0$Y)

    # oracle RD, RR, OR
    mu1_vec[i] <- mu1
    mu0_vec[i] <- mu0
    rd_vec[i] <- mu1 - mu0
    rr_vec[i] <- mu1 / mu0
    or_vec[i] <- (mu1 / (1 - mu1)) / (mu0 / (1 - mu0))
  }
  mc_id <- 1:nreps
  oracle_values <- tibble::tibble(
    scenario = scenario,
    mc_id = rep(mc_id, 5),
    estimand = rep(c("tx mean 1", "tx mean 0", "RR", "RD", "OR"), each = nreps),
    type = "oracle",
    variable = "X",
    value = c(mu1_vec, mu0_vec, rr_vec, rd_vec, or_vec)
  )
  return(oracle_values)
}
