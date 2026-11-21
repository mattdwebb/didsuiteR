# DiD Suite for R

**Compare difference-in-differences estimates and their uncertainty in one table.**

Suppose schools adopt a policy in different years. You want to see how the
estimated effect and standard error change when you use TWFE, two-stage DiD,
imputation, or another supported estimator. `did_suite()` calls the installed
R implementations with a common set of input names and returns a table you can
inspect or export. It keeps each estimator's own definition of the effect and
reports unavailable methods explicitly.

## Install and run a first comparison

Use R 4.1 or later. Install the three estimators used in the first example and
the suite; `remotes` installs its public R dependencies:

```r
install.packages("remotes")  # once, if needed
remotes::install_github(c("mattdwebb/fasttwfeR", "mattdwebb/fastdid2sR",
                          "mattdwebb/fastbjsR", "mattdwebb/didsuiteR"))
```

To install from downloaded source directories instead, run `R CMD INSTALL` on each estimator
and then on `didsuiteR`; the suite needs `digest` and `jsonlite`. R runs the
estimators directly and does not require Python or Stata.

The included teaching data contain 960 pupil observations in 12 schools over
2010–2017. Four schools adopt in 2013, four in 2015, and four never adopt.
`test_score` is measured in score points and `baseline_score` is a predetermined
covariate. These are fictional repeated cross sections, not a panel of pupils.
The artificial policy effect is two points; sampling disturbances make the fitted
estimates differ from that value. [Variable definitions and generation code](inst/extdata/DATA_GUIDE.md)
are included.

```r
library(didsuiteR)
school_data <- read.csv(system.file("extdata", "school_policy.csv",
                                    package = "didsuiteR"))
comparison <- did_suite(school_data, "test_score", "year", "school",
                        gvar = "first_treat", covariates = "baseline_score",
                        data_mode = "rcs", estimators = c("twfe", "did2s", "bjs"),
                        inference = "standard")
example_result <- comparison$summary[c("estimator", "estimate", "native_se",
                                      "jackknife_se", "jackknife_p", "status")]
print(example_result, digits = 7, row.names = FALSE)
```

Actual output, rounded to six decimal places:

| Estimator | Estimate | Native SE | Jackknife SE | Jackknife p-value | Status |
| --- | --- | --- | --- | --- | --- |
| twfe | 2.012599 | 0.367339 | 0.401294 | 0.000393 | PASS |
| did2s | 2.316869 | 0.285375 | 0.336814 | 0.000027 | PASS |
| bjs | 2.316869 | 0.290170 | 0.336814 | 0.000027 | PASS |

TWFE estimates an adjusted difference of about **2.01 score points**, compared
with **2.32 points** for two-stage DiD and imputation. The last two estimates
coincide here because they use the same untreated model and static aggregate;
their native standard errors differ. `jackknife_se` deletes schools one at a time.
Its p-value uses Student t(11), reflecting the 12 schools, not the 960 pupil rows.
`PASS` means the requested computations succeeded, not that the causal assumptions
have been established. Save the complete comparison with:

```r
export_results(comparison, "school_comparison.csv")
```

## Add estimators

Install only the packages you need. `discover_estimators()` shows which are
available in your R library.

| Suite name | R package / repository | Important scope |
|---|---|---|
| `twfe` | `fasttwfeR` | Pooled fixed-effect coefficient |
| `did2s` | `fastdid2sR` | Untreated first stage, static aggregate |
| `bjs` | `fastbjsR` | Untreated-outcome imputation |
| `flex` | `fastflexR` | School/state-group flexible effects |
| `csdid` | `fastrcsdid` | Repeated cross sections; none, RA, IPW or DR |
| `undid` | `suiteundidR` from `suiteundidR` | Suite adapter; no separate native analytic SE |
| `didint` | `suitedidintR` from `suitedidintR` | Covariate standardization; no separate native analytic SE |
| `stacked` | `faststackedR` | Explicit event window containing −1 |
| `lpdid` | `fastlpdidR` | Genuine panel and explicit pooled window |

Repositories are under [mattdwebb](https://github.com/mattdwebb).
The UN-DID adapter installs as `suiteundidR`, separately from the original CRAN `undidR` package. CSDID uses the canonical `fastrcsdid` package from the `R` subdirectory of `mattdwebb/fastcsdid`.

For stacked DiD, for example, use `estimators = c("twfe", "stacked")` and
`estimator_options = list(stacked = list(event_window = c(-1, 1)))` in the call
above after installing `faststackedR`. For LP-DiD load the included
`school_panel.csv`, pass `unit = "school", data_mode = "panel"`, select
`estimators = "lpdid"`, and set `event_window = c(0, 1)`. That file tracks the
same schools over time; pupil row numbers must not be invented as panel IDs.

`estimators = "all"` asks for every adapter. It cannot make LP-DiD support
repeated cross sections or make CSDID use panel equations. Expect an explicit
`UNSUPPORTED` row when a package or requested feature is unavailable.

## Choose covariates and inference

Omit `covariates` for no additional controls. Otherwise supply names of numeric
columns, as above. Unsupported covariate requests fail rather than being ignored.
Never-treated cohorts are zero in the common interface; the suite translates
that convention when a standalone estimator uses a different one.

Every jackknife requires an estimate for each input cluster deletion. When one
cohort disappears, an estimable surviving aggregate is reweighted and retained.
An undefined deletion invalidates that inference method; it does not erase an
otherwise valid original point estimate. CSDID additionally permits an explicit
fixed-support policy through its estimator options.

An unavailable native standard error is reported as missing, not replaced with a
jackknife standard error under a native label. Inspect `comparison$details` and
the status/failure columns when a row is partial. In observational work,
randomization inference still needs a defensible reference assignment mechanism.

## Interpret comparisons carefully

Different estimators can weight cohorts and periods differently. A common input
table does not ensure a common ATT under heterogeneous effects. The suite keeps
those native targets and labels them; automatic common-sample intersections and
general harmonized targets are not implemented. Missing required inputs are
rejected rather than silently selecting a different sample for each estimator.

Use `?did_suite`, `?assignment_set`, and `?export_results` for the API. The
[method specification](METHOD_SPECIFICATION.md), [inference guide](INFERENCE.md),
and [validation record](VALIDATION.md) describe the supported calculations.
[Developer details](docs/DEVELOPER_REFERENCE.md) are separate from this walkthrough.
Cite the estimators used and record package versions and source commits.
The wrapper's original code is [MIT licensed](LICENSE); estimator licenses remain
separate.

After installation, `Rscript tests/readme.R` checks the worked comparison when
its three optional estimator packages are installed. The [executable example](inst/examples/school_policy.R)
and its expected numerical results ship with the package. Release integration
checks install all supported dependencies; a wrapper-only CI run does not certify
every optional estimator.

Current candidate checks: [release validation](RELEASE_VALIDATION.md).

### CSDID dependency

Install the canonical public CSDID package with
`remotes::install_github("mattdwebb/fastcsdid", ref = "v1.0.1", subdir = "R")`.
The separate research repository named `fastrcsdid` is not a suite dependency.


## Inference scope and deletion support

The suite provides coefficient RI (`ri_beta`) and native-studentized RI
(`ri_t_native`) in addition to native and standalone jackknife inference.
These tests are implemented by the suite. Standalone Fast CSDID intentionally
has a narrower public inference API, so its API checker rejects those option
names; that restriction does not remove the suite's separately tested methods.
See [jackknife support and identification](JACKKNIFE_SUPPORT.md) for surviving
ATT aggregation, full-sample failures, overlap limits, and interface differences.
