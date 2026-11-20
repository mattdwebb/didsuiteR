# Earlier development notes

These are the previous landing-page notes, preserved on 19 November 2026 for
provenance. They may describe superseded restrictions. Use the [current README](../README.md),
[inference guide](../INFERENCE.md), and [change history](../CHANGELOG.md) for current behavior.

# didsuiteR

Native R comparison harness, version 0.1.0. Development milestone, not full
nine-estimator certification. Original code is MIT, copyright 2026 Matt Webb.
Estimator packages remain separately licensed and independently installed.

## Install And Run

Install digest and jsonlite, then run `R CMD INSTALL didsuiteR`.
Use `discover_estimators()` to inspect optional native package entrypoints.
No Python, Stata, or Julia runtime is invoked by the harness.

```r
library(didsuiteR)
result <- did_suite(data, outcome="y", time="year", cluster="state",
  gvar="g", estimators="all", data_mode="rcs", inference="full",
  assignments=common_matrix,
  estimator_options=list(stacked=list(event_window=c(-1,1),control_group="both"),
                         lpdid=list(event_window=c(0,1))))
result$summary
export_results(result,"results.csv")
```

Use treatment="D" instead of gvar for absorbing binary treatment; both may be
supplied consistently. Cohort 0 means never-treated in the common interface.
Assignment columns follow sorted character cluster labels; named columns are
validated. A unit argument with repeated units enables panel auto-detection.
Each row remains available when another estimator or inference method fails.

## Estimator Contracts

| Family | Native R entrypoint | Main limitation |
|---|---|---|
| TWFE | fasttwfeR::fasttwfe | Native coefficient, not common ATT |
| DID2S | fastdid2sR::fastdid2s | Untreated never/not-yet first stage |
| BJS | fastbjsR::fastbjs | Native imputation target |
| FLEX | fastflexR::fastflex | Historical state-group target |
| CSDID | fastrcsdid::fast_csdid | RCS only; none/reg/ipw/dr; dr default |
| UN-DID | suiteundidR::suite_undid | Native VCE unavailable |
| DID-INT | suitedidintR::suite_didint | Native VCE unavailable |
| Stacked | faststackedR::fit | Explicit full window containing -1 |
| LP-DiD | fastlpdidR::fit | Panel only; explicit pooled window |

control_group="native" is the default. Actual control metadata comes from the
adapter, not from an ignored common rule. Explicit unsupported rules fail.
Stacked defaults to both clean future-treated and never-treated controls;
the adapter translates common cohort 0 to missing. Its point calls use
inference="point". CSDID/UN-DID/DID-INT native defaults are never controls.

Per-estimator options override event_window and control_group for that estimator
only. CSDID additionally accepts method, ps_fit, aggregation_weights and jackknife_support; UN-DID aggregation_weights; DID-INT ccc; LP-DiD weighting
and covariate_lag. Other unknown options fail. Observation weights are native
when supported; LP/UN-DID/DID-INT reject unavailable weighted variants.

## Inference

The jackknife refits the estimator after each cluster deletion. If a cohort
disappears, the estimator recomputes the aggregate and weights over its
estimable effects. All cluster deletions must yield estimates for the
jackknife standard error. CSDID also accepts an explicit fixed-support
policy through its estimator option `jackknife_support="fixed"`.

For UN-DID, `estimator_options=list(undid=list(aggregation_weights="comparison"))`
uses the two comparison periods' treated sample sizes within a cohort-year.
The default uses treated post-period counts. For CSDID, `aggregation_weights="comparison"`
uses treated cohort mass in the comparison periods in place of full-window mass.
For adjusted DID-INT, the control's fitted mean at the treated covariate
mean must be identified at both dates. Each state-year effect uses its
eligible controls, and surviving state-year effects are reweighted after deletion.

RI preserves the cohort-count multiset, excludes observed/duplicate draws and
uses inclusive absolute tails with plus-one correction. Every reassignment
recomputes its own estimate and native or complete-JK denominator. Invalid
draws produce NA, never a p-value based on survivors. RI supports zero null only.
assignment_set supports supplied matrices and bounded exact/seeded generation,
with RNG restoration and SHA256 canonical serialization shared with Python.
Large assignment-space counts are floating approximations; no arbitrary precision.

## Scope And Evidence

target_mode="native" does not imply equal targets. target_compatible remains
FALSE. Harmonized Stacked (0,H) maps to its (-1,H) window but is not a claim of
common target equivalence. Other harmonized adapters, common sample intersection,
and distinct assignment/inference units are unsupported.

Strict complete common inputs are required; this is not automatic complete-case
selection. Canonical CSDID formula support is narrowed to named numeric covariates.
Production algorithms are ordinary native R; no large Monte Carlo run was made.

tests/contracts.R runs without optional estimators. tests/smoke.R contains 53
local checks with canonical CSDID and installed native ports. tests/routing.R
checks control metadata, window isolation, Stacked translation and LP panel
inference. validation/master.R compares the shared RCS fixture to supplied Python
outputs; it requires the external sibling spec fixtures. Check CSVs for actual
availability and discrepancies. Frozen comparison success is not external
certification of every design or package option.
