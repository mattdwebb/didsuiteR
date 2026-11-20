library(didsuiteR)
d <- expand.grid(person=1:5,time=1:5,state=1:9)
d$g <- c(3,3,3,4,4,4,0,0,0)[d$state]
d$x <- sin(d$person*2+d$state*.3+d$time*.7)
d$w <- 1+(d$person%%3)/4
d$y <- d$state*.3+d$time*.2+.4*d$x+2*(d$g>0 & d$time>=d$g)+
  sin(d$state*d$time)*.3+cos(d$person+d$state)*.2
d <- d[!d$state %in% 2:3,]
run <- function(method,policy=NULL,...) {
  opts <- list(method=method)
  if(!is.null(policy)) opts$jackknife_support <- policy
  did_suite(d,'y','time','state',gvar='g',estimators='csdid',weights='w',
    covariates=if(method=='none') character() else 'x',
    estimator_options=list(csdid=opts),...)
}
if (requireNamespace('fastrcsdid', quietly = TRUE)) {
for(method in c('none','reg','ipw','dr')) {
  got <- run(method,inference='standard')
  explicit <- run(method,'changing',inference='jackknife')
  fixed <- run(method,'fixed',inference='standard')
  stopifnot(got$summary$status=='PASS',got$summary$jackknife_support=='changing',
    is.finite(got$summary$jackknife_p),fixed$summary$status=='PARTIAL',
    fixed$summary$jackknife_support=='fixed',is.na(fixed$summary$jackknife_p),
    identical(got$summary$estimate,fixed$summary$estimate),
    identical(got$summary$native_p,fixed$summary$native_p),
    identical(got$summary$jackknife_p,explicit$summary$jackknife_p))
  jk <- got$details$csdid$methods$jackknife$value
  stopifnot(jk$support_policy=='changing',jk$support_changed['1'],length(jk$deletions)==7)
  physical <- vapply(names(jk$deletions),function(id) {
    f <- fastrcsdid::fast_csdid(d[as.character(d$state)!=id,], 'y','time','g','state',
      weights='w',est_method=method,xformla=if(method=='none') NULL else ~x,
      aggregation='simple',inference='none')
    unname(f$aggregation$estimates['simple'])
  },numeric(1))
  stopifnot(max(abs(jk$deletions-physical))<1e-9,
    abs(jk$se-sqrt(6/7*sum((physical-mean(physical))^2)))<1e-9)
}
stopifnot(run('none','silent',inference='point')$summary$status=='UNSUPPORTED')
cat('CSDID_SUPPORT_POLICY_INTEGRATION_PASS: four methods; physical deletions; defaults; fixed opt-in\n')
} else {
  unavailable <- run('none', inference = 'standard')
  stopifnot(unavailable$summary$status == 'UNSUPPORTED',
            is.na(unavailable$summary$estimate),
            is.na(unavailable$summary$jackknife_p))
  cat('CSDID_BACKEND_UNAVAILABLE_PASS: optional backend reported explicitly; physical-deletion integration requires fastrcsdid\n')
}
# Changing support never authorizes discarding an undefined deletion.
p <- list(data=data.frame(cluster=1:6,.didsuite_g=c(3,3,3,0,0,0)),cluster='cluster')
bad <- function(pp,options,native) list(estimate=if(1 %in% pp$data$cluster) 2 else NA_real_)
err <- tryCatch(didsuiteR:::.jack(bad,p,list(),list(estimate=2,jackknife_support='changing'),0),error=identity)
stopifnot(err$status=='IDENTIFICATION_FAILURE',is.na(err$deletions['1']),
          all(err$deletions[as.character(2:6)]==2),grepl('Undefined deletion',err$failures['1']))
cat('CSDID_SUPPORT_POLICY_FAILURE_GUARD_PASS: undefined deletion is never discarded\n')
