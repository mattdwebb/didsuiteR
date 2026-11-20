library(didsuiteR)
if(requireNamespace("fastflexR",quietly=TRUE)) {
  path <- system.file("extdata","school_policy.csv",package="didsuiteR")
  d <- read.csv(path)
  direct <- fastflexR::fastflex(d,"test_score","year","school",gvar="first_treat",covariates="baseline_score",inference="native")
  null <- direct$estimate-.5*direct$native_se
  result <- did_suite(d,"test_score","year","school",gvar="first_treat",covariates="baseline_score",estimators="flex",inference="native",null=null)$summary
  stopifnot(result$status=="PASS",result$native_df==direct$native_df,
            abs(result$native_stat-.5)<1e-10,
            abs(result$native_p-2*pt(-.5,direct$native_df))<1e-12,
            abs(result$native_p-2*pnorm(-.5))>1e-3)
  cat("FLEX_NONZERO_STUDENT_WRAPPER_PASS\n")
} else cat("SKIP: optional FLEX backend not installed\n")
