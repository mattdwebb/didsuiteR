library(didsuiteR)
school_data <- read.csv(system.file("extdata", "school_policy.csv",
                                    package = "didsuiteR"))
comparison <- did_suite(school_data, "test_score", "year", "school",
                        gvar = "first_treat", covariates = "baseline_score",
                        data_mode = "rcs", estimators = c("twfe", "did2s", "bjs"),
                        inference = "standard")
example_result <- comparison$summary[c("estimator", "estimate", "native_se",
                                      "jackknife_se", "jackknife_p", "status")]
print(example_result, digits = 7, row.names = FALSE)
