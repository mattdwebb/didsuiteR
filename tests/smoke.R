library(didsuiteR)
if (!requireNamespace("fastrcsdid",quietly=TRUE)) {
  cat("SKIP canonical integration: fastrcsdid is not installed\n")
  quit(status=0)
}
n <- 0L
check <- function(name,expr) { force(expr); n <<- n+1L; cat("PASS",name,"\n") }
eq <- function(a,b,tol=1e-9) stopifnot(isTRUE(all.equal(a,b,tolerance=tol,check.attributes=FALSE)))
fails <- function(expr) stopifnot(inherits(tryCatch({force(expr); NULL},error=identity),"error"))
d <- expand.grid(person=1:5,time=1:5,state=1:9)
d$state <- sprintf("s%02d",d$state)
s <- as.integer(sub("s","",d$state))
d$g <- c(3,3,3,4,4,4,0,0,0)[s]
d$x <- sin(d$person*2+s*.3+d$time*.7)
d$w <- 1+(d$person%%3)/4
d$D <- as.integer(d$g>0 & d$time>=d$g)
d$y <- s*.3+d$time*.2+.4*d$x+2*d$D+sin(s*d$time)*.3+cos(d$person+s)*.2
d$unit <- paste(d$state,d$person,sep=":")
run <- function(data=d,...) did_suite(data,"y","time","state",gvar="g",estimators="csdid",...)
check("nine discovery rows",stopifnot(nrow(discover_estimators())==9))
check("canonical installed",stopifnot(requireNamespace("fastrcsdid",quietly=TRUE)))
for (method in c("none","reg","ipw","dr")) {
  cov <- if (method=="none") character() else "x"
  got <- run(inference="point",covariates=cov,weights="w",estimator_options=list(csdid=list(method=method)))
  ref <- fastrcsdid::fast_csdid(d,"y","time","g","state",weights="w",inference="none",
    est_method=method,xformla=if (method=="none") NULL else ~x,aggregation="simple")
  check(paste("canonical point",method),{
    stopifnot(got$summary$status=="PASS")
    eq(got$summary$estimate,unname(ref$aggregation$estimates["simple"]))
  })
}
a <- c(3,3,3,4,4,4,0,0,0); ids <- sprintf("s%02d",1:9)
draws <- rbind(a[c(4:9,1:3)],a[c(7:9,1:6)],a[c(2:9,1)])
colnames(draws) <- ids
aa <- assignment_set(a,ids,supplied=draws)
check("shared assignment shape",stopifnot(nrow(aa$matrix)==3,nchar(aa$hash)==64))
check("RNG preserved",{set.seed(331); old <- .Random.seed; assignment_set(a,ids,reps=3,seed=17); stopifnot(identical(old,.Random.seed))})
check("seed deterministic",eq(assignment_set(a,ids,reps=3,seed=17)$matrix,
                             assignment_set(a,ids,reps=3,seed=17)$matrix,0))
check("exact support",stopifnot(nrow(assignment_set(c(1,1,0,0),exact=TRUE)$matrix)==5))
check("duplicate assignment rejected",fails(assignment_set(a,ids,supplied=rbind(draws,draws[1,]))))
check("observed excluded",fails(assignment_set(a,ids,supplied=rbind(draws,a))))
check("changed dates rejected",{bad<-draws;bad[1,1]<-7;fails(assignment_set(a,ids,supplied=bad))})
check("column order rejected",fails(assignment_set(a,ids,supplied=draws[,9:1])))
check("treatment alternative",eq(run(inference="point")$summary$estimate,
  did_suite(d,"y","time","state",treatment="D",estimators="csdid",inference="point")$summary$estimate))
check("inconsistent treatment rejected",{bad<-d;bad$D[1]<-1;fails(did_suite(bad,"y","time","state",gvar="g",treatment="D"))})
check("nonabsorbing rejected",{bad<-d;bad$D[bad$state=="s01"&bad$time==5]<-0;fails(did_suite(bad,"y","time","state",treatment="D"))})
check("missing required rejected",{bad<-d;bad$x[1]<-NA;fails(run(bad,covariates="x"))})
check("none cannot ignore X",stopifnot(run(covariates="x",estimator_options=list(csdid=list(method="none")))$summary$status=="UNSUPPORTED"))
check("panel CSDID explicit",stopifnot(run(unit="unit",data_mode="panel",inference="point")$summary$status=="UNSUPPORTED"))
check("duplicate panel rows rejected",fails(run(rbind(d,d[1,]),unit="unit",data_mode="panel")))
check("harmonized explicit",stopifnot(run(target_mode="harmonized",event_window=c(0,1))$summary$status=="UNSUPPORTED"))
check("sample common explicit",stopifnot(run(sample_mode="common")$summary$status=="UNSUPPORTED"))
check("different assignment units explicit",{dd<-d;dd$assignment<-dd$state;
  stopifnot(run(dd,assignment_unit="assignment")$summary$status=="UNSUPPORTED")})
full <- run(inference="full",assignments=aa,estimator_options=list(csdid=list(method="none")))
print(full)
check("CSDID all methods pass",stopifnot(full$summary$status=="PASS"))
check("CV3J physical canonical",{
  b <- vapply(ids,function(id) {
    f<-fastrcsdid::fast_csdid(d[d$state!=id,],"y","time","g","state",inference="none",aggregation="simple")
    unname(f$aggregation$estimates["simple"])
  },numeric(1))
  eq(full$details$csdid$methods$jackknife$value$deletions,b)
  eq(full$summary$jackknife_se,sqrt(8/9*sum((b-mean(b))^2)))
  eq(full$summary$jackknife_p,2*pt(-abs(full$summary$jackknife_t),8))
})
check("RI coefficient physical canonical",{
  vals<-vapply(1:3,function(i) {
    dd<-d;dd$g<-draws[i,match(dd$state,ids)]
    unname(fastrcsdid::fast_csdid(dd,"y","time","g","state",inference="none",aggregation="simple")$aggregation$estimates["simple"])
  },numeric(1))
  eq(full$details$csdid$methods$ri_beta$value$statistics,vals)
})
check("native direct canonical",{
  f<-fastrcsdid::fast_csdid(d,"y","time","g","state",inference="none",aggregation="simple")
  eq(full$summary$native_se,sqrt(f$vcov$analytical["simple","simple"]))
})
check("native RI recomputed denominators",{
  expected<-vapply(1:3,function(i) {
    dd<-d;dd$g<-draws[i,match(dd$state,ids)]
    run(dd,inference="native",estimator_options=list(csdid=list(method="none")))$summary$native_stat
  },numeric(1))
  eq(full$details$csdid$methods$ri_t_native$value$statistics,expected)
})
check("row order invariant",eq(full$summary$estimate,
 run(d[nrow(d):1,],inference="point",estimator_options=list(csdid=list(method="none")))$summary$estimate))
check("singleton deletion retains point under explicit fixed support",{
  dd<-d[d$state %in% c("s01","s04","s05","s06","s07","s08","s09"),]
  ans<-run(dd,inference="standard",estimator_options=list(csdid=list(method="none",jackknife_support="fixed")))
  stopifnot(is.finite(ans$summary$estimate),is.na(ans$summary$jackknife_p),ans$summary$status=="PARTIAL")
})
check("nonzero RI explicitly unsupported",{
  ans<-run(inference=c("native","ri_beta"),assignments=aa,null=1)
  stopifnot(is.finite(ans$summary$native_p),is.na(ans$summary$ri_beta_p),ans$summary$status=="PARTIAL")
})
check("failure isolation",{
  broken<-function(...) stop("UNSUPPORTED: deliberately unavailable")
  ans<-did_suite(d,"y","time","state",gvar="g",estimators=c("csdid","twfe"),
    inference="point",adapters=list(twfe=broken))
  stopifnot(ans$summary["csdid","status"]=="PASS",ans$summary["twfe","status"]=="UNSUPPORTED")
})
check("method failure isolation",{
  native_fail<-function(p,options,native) {
    if(native) stop("UNSUPPORTED: native missing")
    b<-mean(p$data[[p$outcome]])
    list(estimate=b,target_label="test only",covariate_mode="test",weighting_scheme="test")
  }
  ans<-did_suite(d,"y","time","state",gvar="g",estimators="twfe",inference="standard",adapters=list(twfe=native_fail))
  stopifnot(is.finite(ans$summary$jackknife_p),is.na(ans$summary$native_p))
})
check("all rows retained",stopifnot(nrow(did_suite(d,"y","time","state",gvar="g",inference="point")$summary)==9))
check("CSV and JSON export",{
  paths<-c(tempfile(fileext=".csv"),tempfile(fileext=".json"))
  export_results(full,paths[1]);export_results(full,paths[2],"json")
  stopifnot(nrow(read.csv(paths[1]))==1,jsonlite::validate(paste(readLines(paths[2]),collapse="")))
  unlink(paths)
})
# Native installed sibling APIs are independently compared to direct same-language calls.
for (id in c("twfe","did2s","bjs","flex","undid","didint")) {
  catalog<-discover_estimators();r<-catalog[catalog$estimator==id,]
  if (!r$available) {cat("SKIP optional",id,"\n");next}
  ans<-did_suite(d,"y","time","state",gvar="g",estimators=id,inference="point")
  ref<-do.call(getExportedValue(r$package,r$entrypoint),list(data=d,outcome="y",time="time",
    cluster="state",gvar="g",inference="point",data_mode="rcs"))
  check(paste("native routing",id),{stopifnot(ans$summary$status=="PASS");eq(ans$summary$estimate,ref$estimate)})
}
check("canonical assignment SHA256",stopifnot(assignment_set(c(1,1,0,0),c("a","b","c","d"),exact=TRUE)$hash == "4e306ced29240f402f52b029e9253bb9cbab1cf173a35316e0e4f92dbb49f174"))
for (method in c("reg","ipw","dr")) {
  ans <- run(inference="full",covariates="x",weights="w",assignments=draws[1,,drop=FALSE],
    estimator_options=list(csdid=list(method=method)))
  check(paste("weighted covariate full inference",method),stopifnot(ans$summary$status=="PASS"))
}
if (requireNamespace("faststackedR",quietly=TRUE)) {
  ans <- did_suite(d,"y","time","state",gvar="g",estimators="stacked",inference="full",
    event_window=c(-1,1),assignments=aa)
  check("Stacked translated never cohort and full inference",stopifnot(ans$summary$status=="PASS"))
  dd <- d; dd$g[dd$g==0] <- NA_real_
  ref <- faststackedR::fit(dd,y="y",time="time",cluster="state",gvar="g",data_mode="rcs",
    event_window=c(-1,1),control_group="both",inference="native")
  check("Stacked direct routing",eq(ans$summary$estimate,ref$estimate))
}
check("zero JK variance retains point",{
  constant<-function(...) list(estimate=2,target_label="test",covariate_mode="none",weighting_scheme="test")
  ans<-did_suite(d,"y","time","state",gvar="g",estimators="twfe",inference="jackknife",adapters=list(twfe=constant))
  stopifnot(ans$summary$estimate==2,is.na(ans$summary$jackknife_p),ans$summary$status=="PARTIAL")
})
check("invalid RI draws cannot create p value",{
  oracle<-function(p,options,native) {
    if(p$data$.didsuite_g[1]!=3) stop("IDENTIFICATION_FAILURE: deliberate draw failure")
    list(estimate=1,target_label="test",covariate_mode="none",weighting_scheme="test")
  }
  ans<-did_suite(d,"y","time","state",gvar="g",estimators="twfe",inference="ri_beta",assignments=aa,adapters=list(twfe=oracle))
  stopifnot(is.na(ans$summary$ri_beta_p),ans$summary$status=="PARTIAL")
})
check("two independent R processes",{
  cl<-parallel::makePSOCKcluster(2)
  parallel::clusterCall(cl,function(lib) .libPaths(lib),.libPaths())
  got<-parallel::parLapply(cl,list(d,d),function(dd) didsuiteR::did_suite(dd,"y","time","state",gvar="g",
    estimators="csdid",inference="point")$summary$estimate)
  parallel::stopCluster(cl)
  eq(got[[1]],got[[2]],0)
})
cat("TOTAL PASSED:",n,"\n")

