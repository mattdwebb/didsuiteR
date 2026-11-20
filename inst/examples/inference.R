library(didsuiteR)
school_data <- read.csv(system.file("extdata", "school_policy.csv", package = "didsuiteR"))
randomized <- did_suite(school_data, "test_score", "year", "school",
                        gvar = "first_treat", covariates = "baseline_score",
                        data_mode = "rcs", estimators = "twfe",
                        inference = "full", reps = 19, seed = 2026)
replay <- did_suite(school_data, "test_score", "year", "school",
                   gvar = "first_treat", covariates = "baseline_score",
                   data_mode = "rcs", estimators = "twfe", inference = "full",
                   assignments = randomized$assignments)
p_columns <- c("native_p", "jackknife_p", "ri_beta_p", "ri_t_native_p")
print(randomized$summary[c("estimate", p_columns, "status")], row.names = FALSE)
stopifnot(randomized$summary$status == "PASS", replay$summary$status == "PASS",
          identical(randomized$assignments$hash, replay$assignments$hash),
          identical(randomized$summary[p_columns], replay$summary[p_columns]),
          randomized$summary$jackknife_df == 11,
          nrow(randomized$assignments$matrix) == 19)
