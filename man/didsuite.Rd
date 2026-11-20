\name{did_suite}
\alias{did_suite}
\alias{discover_estimators}
\alias{assignment_set}
\alias{export_results}
\title{Native Difference-in-Differences Comparison Harness}
\description{Discover native R estimators, validate inputs, and isolate estimator
and inference failures. Development interface, not universal certification.}
\usage{
did_suite(data, outcome, time, cluster, estimators = "all",
  gvar = NULL, treatment = NULL, assignment_unit = cluster, unit = NULL,
  covariates = character(), weights = NULL, data_mode = "auto",
  inference = "standard", assignments = NULL, reps = 999L, seed = 0L,
  exact = FALSE, null = 0, target_mode = "native", sample_mode = "native",
  event_window = NULL, control_group = "native", estimator_options = list(),
  adapters = list())
discover_estimators()
assignment_set(observed, units = names(observed), reps = 999L, seed = 0L,
  exact = FALSE, supplied = NULL)
export_results(x, path, format = c("csv", "json"))
}
\arguments{
\item{data}{Input data frame with complete required variables.}
\item{outcome, time, cluster}{Column names for outcome, time and inference cluster.}
\item{estimators}{Family names or all.}
\item{gvar, treatment}{Cohort or absorbing binary treatment column; zero cohort means never.}
\item{assignment_unit, unit}{Assignment cluster and optional panel unit columns.}
\item{covariates, weights}{Numeric covariate names and optional positive weight column.}
\item{data_mode}{auto, panel or rcs.}
\item{inference}{point, standard, full, or individual method names.}
\item{assignments, supplied}{Supplied alternative cohort matrix, rows by assignment and columns by sorted units.}
\item{reps, seed, exact}{Bounded alternative count, seed and exact enumeration switch.}
\item{null}{Scalar null. RI currently requires zero.}
\item{target_mode, sample_mode}{Native or explicitly unsupported alternatives; see README.}
\item{event_window}{Two endpoints; supported natively for Stacked and LP-DiD.}
\item{control_group}{Native by default. Unsupported alternate rules fail.}
\item{estimator_options}{Named per-estimator lists; window and control overrides do not leak.}
\item{adapters}{Optional named explicit adapter functions, primarily for contract tests.}
\item{observed, units}{Observed cohort labels and unique unit labels.}
\item{x}{Harness result.}
\item{path, format}{Explicit output path and CSV or JSON format.}
}
\value{did_suite returns summary, method details, event studies and assignments.
discover_estimators returns installed-entrypoint metadata, not certification.
assignment_set returns validated draws, design and SHA256 hash. export_results
invisibly returns the output path.}
\details{Jackknife uses all physical cluster deletions with mean centering and
t(G-1). RI recomputes each assigned denominator; invalid draws invalidate its
p-value. Native controls and targets remain estimator specific.}
\examples{
discover_estimators()
assignment_set(c(1,1,0,0), letters[1:4], exact=TRUE)
}
