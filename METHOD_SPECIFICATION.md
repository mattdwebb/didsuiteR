# What the R suite compares

`did_suite()` prepares a common input description and calls independent R
estimators. It returns one row for each requested estimator, including its point
estimate, available inference, target label, support diagnostics and failure
reason. It does not average estimates across methods.

## Data and treatment

Give column names for the outcome, period, inference cluster and first-treatment
date. First-treatment dates are constant within a cluster; zero means never
treated. Alternatively supply an absorbing binary treatment history. If both
are provided they must agree. Assignment and inference clusters currently must
coincide. Required inputs must be complete; automatic estimator-specific missing-
value selection is not performed.

Repeated cross sections and genuine panels are distinct. Panel mode requires a
unit identifier with unique unit-period records. CSDID uses repeated-cross-section
equations in this implementation; LP-DiD requires a panel. The suite reports an
unsupported combination instead of manufacturing longitudinal observations.

## Estimator targets

| Estimator | Calculation and reported scalar |
|---|---|
| TWFE | Pooled treatment coefficient with fixed effects and optional covariates. |
| did2s | Untreated first-stage regression followed by the static treated-outcome adjustment. |
| BJS | Imputed untreated outcomes, averaged over supported treated observations. |
| FLEX | Historical school/state-group flexible treatment effects and its native aggregate. |
| CSDID | Group-time comparisons and the supported simple aggregate; method can be none, RA, IPW or DR. |
| UN-DID | School/state changes and comparisons, aggregated with the selected count convention. |
| DID-INT | Target-standardized fitted school-year means, contrasts and surviving-effect aggregation. |
| Stacked | Corrective-weighted stacked event regression over the explicitly requested window. |
| LP-DiD | Supported pooled local-projection effect over the explicitly requested horizon. |

`target_mode = "native"` preserves these definitions. With heterogeneous effects,
the rows need not target the same weighted ATT. `target_compatible` therefore does
not assert common-target equivalence. General target harmonization and automatic
common-sample intersections are unsupported. The limited stacked window
translation in harmonized mode is not a proof of equal estimands.

## Controls, covariates and weights

The default control rule is each estimator's native rule. TWFE uses the pooled
sample; did2s/BJS use untreated observations; CSDID/UN-DID/DID-INT default to
never-treated controls; stacked DiD uses eligible clean controls; LP-DiD requires
controls that remain untreated through the horizon. Explicit unsupported rules
raise a classified failure rather than being ignored.

Numeric covariates are passed to the selected implementation. The suite does not
silently freeze contemporaneous repeated-cross-section controls at baseline.
LP-DiD's default is lag-one covariates. `estimator_options` selects supported
method-specific settings, such as CSDID's nuisance method, DID-INT's slope model,
and separate stacked/LP windows. Unknown options are rejected.

Observation weights are used only where the native implementation supports them.
UN-DID, DID-INT and LP-DiD reject unavailable observation-weight specifications.
No unsupported weighted estimator is substituted.

## Failures and reproducibility

A failed inference method need not erase a valid point estimate. Inspect the
summary's status and failure reason, and the per-method objects in `details`.
Record source commits and package versions together with the data definition,
options and assignment object. The [README](README.md) contains an executable
first comparison; [INFERENCE.md](INFERENCE.md) defines its uncertainty calculations.
