# No private estimator dependency: validate the harness with a fixed adapter.
if (dir.exists("R")) {
  env <- new.env(parent = globalenv())
  for (file in list.files("R", pattern = "\\.[rR]$", full.names = TRUE))
    sys.source(file, envir = env)
} else env <- asNamespace("didsuiteR")

reference_draws <- function(a, reps, seed) {
  set.seed(seed)
  key <- function(x) paste(sprintf("%.17g", x), collapse = ",")
  seen <- key(a)
  rows <- list()
  while (length(rows) < reps) {
    draw <- a[sample.int(length(a))]
    candidate <- key(draw)
    if (!candidate %in% seen) {
      seen <- c(seen, candidate)
      rows[[length(rows) + 1L]] <- draw
    }
  }
  do.call(rbind, rows)
}

a <- c(3, 3, rep(0, 5998))
long_key <- paste(a, collapse = ",")
stopifnot(nchar(long_key, type = "bytes") > 10000L)
expected <- reference_draws(a, 4L, 19L)
units <- sprintf("s%04d", seq_along(a))
units[1L] <- paste0(units[1L], strrep("\u00e9", 6000L))
set.seed(291)
before <- .Random.seed
actual <- env$assignment_set(a, units, reps = 4L, seed = 19L)
supplied <- env$assignment_set(a, units, supplied = expected)
stopifnot(identical(actual, supplied), identical(before, .Random.seed))
rm(".Random.seed", envir = .GlobalEnv)
stopifnot(identical(env$assignment_set(a, units, reps = 4L, seed = 19L), actual),
          !exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))

seen <- new.env(parent = emptyenv())
stopifnot(identical(env$.ri_key_bucket("abc"), env$.ri_key_bucket("b`d")))
env$.ri_key_add(seen, "abc")
stopifnot(!env$.ri_key_has(seen, "b`d"))
env$.ri_key_add(seen, "b`d")
env$.ri_key_add(seen, "abc")
stopifnot(length(seen[[env$.ri_key_bucket("abc")]]) == 2L)
for (key in c("", "\u00e9\u4e2d", long_key, paste0(long_key, ",3"))) {
  env$.ri_key_add(seen, key)
  stopifnot(env$.ri_key_has(seen, key))
}
stopifnot(!env$.ri_key_has(seen, paste0(long_key, ",missing")),
          all(nchar(ls(seen), type = "bytes") < 100L))
colliding <- new.env(parent = env)
colliding$.ri_key_bucket <- function(key) "one-bucket"
for (name in c(".ri_key_has", ".ri_key_add", "assignment_set")) {
  fun <- get(name, envir = env)
  environment(fun) <- colliding
  assign(name, fun, envir = colliding)
}
stopifnot(identical(colliding$assignment_set(a, units, reps = 4L, seed = 19L), actual))

# This deterministic adapter exercises all harness inference paths; its simple
# mean contrast and fixed native SE are test fixtures, not a claimed DID model.
adapter <- function(p, options, native) {
  d <- p$data
  treated <- d$.didsuite_g > 0 & d[[p$time]] >= d$.didsuite_g
  list(estimate = mean(d[[p$outcome]][treated]) - mean(d[[p$outcome]][!treated]),
       native_se = 0.25, native_df = Inf, native_method = "fixture",
       target_label = "fixture mean contrast", covariate_mode = "none",
       weighting_scheme = "uniform")
}
d <- expand.grid(s = 1:20, t = 1:4)
cohorts <- rep(c(0, 3), each = 10L)
d$g <- cohorts[d$s]
d$y <- d$s/10 + d$t/3 + 2*(d$g > 0 & d$t >= d$g) + sin(d$s*d$t)/5
ids <- sort(as.character(unique(d$s)))
ordered_cohorts <- cohorts[as.integer(ids)]
draws <- reference_draws(ordered_cohorts, 3L, 29L)
args <- list(data = d, outcome = "y", time = "t", cluster = "s", gvar = "g",
             estimators = "twfe", adapters = list(twfe = adapter), inference = "full")
set.seed(192)
before <- .Random.seed
generated <- do.call(env$did_suite, c(args, list(reps = 3L, seed = 29L)))
stopifnot(identical(before, .Random.seed), generated$summary$status == "PASS")
supplied <- do.call(env$did_suite, c(args, list(assignments = draws)))
fields <- setdiff(names(generated$summary), grep("^runtime", names(generated$summary), value = TRUE))
stopifnot(identical(generated$summary[fields], supplied$summary[fields]),
          identical(generated$assignments, supplied$assignments))
for (method in c("native", "jackknife", "ri_beta", "ri_t_native")) {
  stopifnot(identical(generated$details$twfe$methods[[method]]$value,
                      supplied$details$twfe$methods[[method]]$value))
}
stopifnot(env$assignment_set(c(1,1,0,0), c("a","b","c","d"), exact = TRUE)$hash ==
            "4e306ced29240f402f52b029e9253bb9cbab1cf173a35316e0e4f92dbb49f174")
cat("PASS: didsuiteR long keys/labels, collisions, RNG, canonical hashes, seeded harness inference\n")
