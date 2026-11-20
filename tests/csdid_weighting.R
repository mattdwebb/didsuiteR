library(didsuiteR)
if(requireNamespace("fastrcsdid",quietly=TRUE)) {
  base <- if(dir.exists("tests")) "tests" else "."
  d <- read.csv(file.path(base,"csdid_weighting_data.csv"))
  expected <- read.csv(file.path(base,"csdid_weighting_reference.csv"))
  for(method in c("reg","dr")) for(control in c("never","not_yet")) for(weighting in c("cohort","comparison")) {
    z <- did_suite(d,"y","time","cluster",gvar="g",covariates=c("x","z"),weights="w",estimators="csdid",control_group=control,estimator_options=list(csdid=list(method=method,aggregation_weights=weighting)))
    row <- z$summary
    ref <- expected[expected$method==method & expected$control==control & expected$weighting==weighting,]
    stopifnot(row$status=="PASS",abs(row$estimate-ref$estimate[1])<1e-8,abs(row$native_se-ref$native_se[1])<1e-8,abs(row$jackknife_se-ref$jk_se[1])<1e-8,row$jackknife_df==11)
    vals <- z$details$csdid$methods$jackknife$value$deletions
    stopifnot(max(abs(vals-ref$estimate[match(as.integer(names(vals)),ref$cluster)]))<1e-8)
  }
  cat("CSDID_ALL_WEIGHTS_CONTROLS_SURVIVING_COHORTS_PARITY_PASS\n")
} else cat("SKIP: optional CSDID backend not installed\n")
