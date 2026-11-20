library(didsuiteR)
base <- if(dir.exists("tests")) "tests" else "."
if(requireNamespace("fastrcsdid",quietly=TRUE)) {
bad <- read.csv(file.path(base,"unidentified_cell.csv"))
partial <- read.csv(file.path(base,"partial_deletion.csv"))
run <- function(d,method,weighting,inference="standard",policy="changing") {
 did_suite(d,"y","time","cluster",gvar="g",weights="w",covariates=c("x","z"),
  data_mode="rcs",estimators="csdid",inference=inference,
  estimator_options=list(csdid=list(method=method,aggregation_weights=weighting,jackknife_support=policy)))
}
records <- list()
for(method in c("reg","dr")) for(weighting in c("cohort","comparison")) {
 for(inference in c("point","standard")) {
  z <- run(bad,method,weighting,inference)
  stopifnot(z$summary$status=="IDENTIFICATION_FAILURE",!is.finite(z$summary$estimate))
 }
 z <- run(partial,method,weighting)
 stopifnot(z$summary$status=="PASS",z$summary$jackknife_df==12)
 jk <- z$details$csdid$methods$jackknife$value
 stopifnot(length(jk$deletions)==13,all(is.finite(jk$deletions)),jk$support_changed['13'])
 stopifnot(abs(jk$se-sqrt(12/13*sum((jk$deletions-mean(jk$deletions))^2)))<1e-12)
 physical <- run(partial[partial$cluster!=13,],method,weighting,"point")
 stopifnot(physical$summary$status=="IDENTIFICATION_FAILURE")
 fixed <- run(partial,method,weighting,policy="fixed")
 stopifnot(fixed$summary$status!="PASS",!is.finite(fixed$summary$jackknife_se))
 records[[length(records)+1L]] <- data.frame(method,weighting,cluster=as.numeric(names(jk$deletions)),
  deletion=as.numeric(jk$deletions),att=z$summary$estimate,se=jk$se)
}
out <- Sys.getenv("SUPPORT_RESULTS")
if(nzchar(out)) write.csv(do.call(rbind,records),out,row.names=FALSE)
cat("FULL_SAMPLE_FAILURE_AND_SURVIVING_DELETIONS_PASS\n")
} else cat("CSDID integration checks skipped: optional backend unavailable\n")
p <- list(data=data.frame(cluster=1:2,.didsuite_g=c(2,0),y=c(1,3)),cluster="cluster")
jk <- didsuiteR:::.jack(function(pp,options,native) list(estimate=mean(pp$data$y)),p,list(),list(estimate=2),0)
stopifnot(jk$df==1,is.finite(jk$se),length(jk$deletions)==2)
