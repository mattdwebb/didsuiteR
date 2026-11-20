# Execute the installed-package example used in the README.
available <- vapply(c("fasttwfeR", "fastdid2sR", "fastbjsR"),
                    requireNamespace, logical(1), quietly = TRUE)
if (!all(available)) {
  message("README integration requires optional fasttwfeR, fastdid2sR and fastbjsR; run the release integration check after installing them.")
} else {
source(system.file("examples", "school_policy.R", package = "didsuiteR"))
expected <- read.csv(system.file("extdata", "readme_expected.csv", package = "didsuiteR"))
numeric_columns <- names(expected)[vapply(expected, is.numeric, logical(1))]
stopifnot(nrow(example_result) == nrow(expected),
          all(is.finite(as.matrix(example_result[numeric_columns]))),
          max(abs(as.matrix(example_result[numeric_columns]) - as.matrix(expected[numeric_columns]))) < 2e-6)
stopifnot(all(example_result$status == "PASS"))
}
