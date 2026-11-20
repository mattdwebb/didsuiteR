# Dependency-independent harness tests use an explicit mathematical adapter.
library(didsuiteR)
d <- expand.grid(time=1:4,cluster=letters[1:6])
d$g <- c(3,3,4,4,0,0)[match(d$cluster,letters[1:6])]
d$y <- d$time + match(d$cluster,letters[1:6])/10
adapter <- function(p,options,native=FALSE) list(estimate=mean(p$data$y),
  native_se=.25,native_df=5,native_method="test fixture",target_label="test mean")
r <- did_suite(d,"y","time","cluster",gvar="g",estimators="twfe",inference="full",
               reps=2,adapters=list(twfe=adapter))
stopifnot(r$summary$status=="PASS",r$summary$jackknife_df==5,
          r$summary$ri_beta_p==1,nrow(r$assignments$matrix)==2)
z <- assignment_set(c(1,1,0,0),letters[1:4],exact=TRUE)
stopifnot(z$hash=="4e306ced29240f402f52b029e9253bb9cbab1cf173a35316e0e4f92dbb49f174")
cat("DIDSUITE_CONTRACTS_PASS\n")
