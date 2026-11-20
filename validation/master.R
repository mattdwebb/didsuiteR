# Run from didsuiteR with R_LIBS_USER pointing at installed native dependencies.
library(didsuiteR)
spec <- if (length(commandArgs(TRUE))) commandArgs(TRUE)[1] else "../didsuite_spec"
d <- read.csv(file.path(spec,"fixtures","stochastic_rcs.csv"),stringsAsFactors=FALSE)
a <- as.matrix(read.csv(file.path(spec,"fixtures","shared_assignments.csv"),check.names=FALSE))
meta <- jsonlite::read_json(file.path(spec,"fixtures","shared_assignments.json"),simplifyVector=TRUE)
r <- did_suite(d,"y","year","state",gvar="g",data_mode="rcs",inference="full",
  assignments=a,estimator_options=list(stacked=list(event_window=c(-1,1),control_group="both"),
                                      lpdid=list(event_window=c(0,1))))
stopifnot(r$assignments$hash=="bea08f08f138e2343a8e1c618a7542169a9c56f966dd174dbac3d240abc3b4ec")
export_results(r,"validation/master_r.csv")
export_results(r,"validation/master_r.json","json")
py <- read.csv(file.path(spec,"validation","results","master_python.csv"),stringsAsFactors=FALSE)
fields <- intersect(c("estimate","native_se","jackknife_se","ri_beta_p","ri_t_native_p"),
                    names(py))
diffs <- do.call(rbind,lapply(fields,function(v) {
  rhs <- py[[v]][match(r$summary$estimator,py$estimator)]
  lhs <- r$summary[[v]]
  data.frame(estimator=r$summary$estimator,field=v,R=lhs,Python=rhs,
    absolute_gap=abs(lhs-rhs),status=ifelse(is.na(lhs)&is.na(rhs),"BOTH_NA",
      ifelse(is.na(lhs)|is.na(rhs),"AVAILABILITY_DIFF",ifelse(abs(lhs-rhs)<1e-8,"PASS","DIFFERENCE"))))
}))
write.csv(diffs,"validation/master_python_comparison.csv",row.names=FALSE,na="")
print(r)
print(diffs[diffs$status %in% c("DIFFERENCE","AVAILABILITY_DIFF"),])
stopifnot(!any(diffs$status=="DIFFERENCE"))
cat("DIDSUITE_R_MASTER_COMPARISON_PASS\n")
