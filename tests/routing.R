library(didsuiteR)
d <- expand.grid(person=1:5,time=1:5,state=1:9)
d$state <- sprintf("s%02d",d$state)
s <- as.integer(sub("s","",d$state))
d$g <- c(3,3,3,4,4,4,0,0,0)[s]
d$x <- sin(d$person*2+s*.3+d$time*.7)
d$D <- as.integer(d$g>0 & d$time>=d$g)
d$y <- s*.3+d$time*.2+.4*d$x+2*d$D+sin(s*d$time)*.3+cos(d$person+s)*.2
d$unit <- paste(d$state,d$person,sep=":")
available <- discover_estimators()
has <- function(id) available$available[match(id,available$estimator)]
run <- function(...) did_suite(d,"y","time","state",gvar="g",...)
if (has("stacked") && has("did2s")) {
  r <- run(estimators=c("did2s","stacked"),inference="point",
           estimator_options=list(stacked=list(event_window=c(-1,1))))
  stopifnot(all(r$summary$status=="PASS"),
    r$summary$control_group[1]=="untreated observations: never and not-yet",
    r$summary$control_group[2]=="both",r$summary$event_window[1]=="",
    r$summary$event_window[2]=="-1,1")
  h <- run(estimators="stacked",inference="point",target_mode="harmonized",event_window=c(0,1))
  stopifnot(h$summary$status=="PASS",!h$summary$target_compatible,
            abs(h$summary$estimate-r$summary$estimate[2])<1e-12,
            h$summary$event_window=="-1,1")
  bad <- run(estimators="did2s",inference="point",control_group="never")
  stopifnot(bad$summary$status=="UNSUPPORTED")
  cat("PASS native controls, per-estimator windows and honest harmonized translation\n")
}
if (has("lpdid")) {
  a <- c(3,3,3,4,4,4,0,0,0)
  assignments <- rbind(a[c(4:9,1:3)])
  colnames(assignments) <- sprintf("s%02d",1:9)
  r <- run(estimators="lpdid",unit="unit",data_mode="panel",inference="full",
    assignments=assignments,estimator_options=list(lpdid=list(event_window=c(0,1))))
  ref <- fastlpdidR::fit(d,"y","time",unit="unit",gvar="g",cluster="state",
    event_window=c(0,1),data_mode="panel",inference="native")
  stopifnot(r$summary$status=="PASS",abs(r$summary$estimate-ref$estimate)<1e-10,
            abs(r$summary$native_se-ref$native_se)<1e-10,
            r$summary$control_group=="native clean through pooled horizon")
  ra <- run(estimators="lpdid",unit="unit",data_mode="panel",covariates="x",
    inference="standard",event_window=c(0,1))
  stopifnot(ra$summary$status=="PASS",is.infinite(ra$summary$native_df))
  rcs <- run(estimators="lpdid",data_mode="rcs",inference="point",event_window=c(0,1))
  stopifnot(rcs$summary$status=="UNSUPPORTED")
  cat("PASS LP panel full inference, RA normal reference and RCS exclusion\n")
}
cat("DIDSUITE_ROUTING_PASS\n")
