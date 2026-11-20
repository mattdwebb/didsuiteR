# Inference in the R suite

## Native inference

The suite uses each estimator's own clustered standard error and reported
reference degrees of freedom. An infinite degree of freedom means a normal
reference; a finite value means Student t. These conventions can differ across
estimators. Covariances that require first-stage estimation uncertainty retain
that adjustment in the corresponding native implementation.

## Delete-cluster jackknife

Refit the estimator after deleting each of the G input inference clusters. Let
theta[-g] be its deleted estimate and theta_bar their mean. The reported variance is

    V_J = (G - 1) / G * sum_g (theta[-g] - theta_bar)^2.

The standard error is sqrt(V_J). The standalone two-sided p-value compares
(theta - null)/sqrt(V_J) with Student t(G-1). This is the mean-centered CV3J
variance, not a variance centered on the full-sample estimate.

Each deletion re-estimates the nuisance models and aggregation weights. If a
cohort disappears but supported effects remain estimable, the aggregate is
recomputed over those effects. CSDID also accepts an explicit
`estimator_options = list(csdid = list(jackknife_support = "fixed"))` policy
requiring its original comparison support. All scalar deletion estimates are
required under either policy. Failed deletions are retained in diagnostics;
they are never omitted to produce a smaller-G standard error. Numerical-zero
variance is an explicit failure.

The suite physically refits its point adapter for every deletion. In particular,
it does not substitute CSDID's optional linearized nonlinear jackknife for these
complete refits.

## Randomization inference

Assignments preserve the observed number of clusters in each adoption cohort,
including never treated. The observed assignment and duplicate alternatives are
excluded. `reps` requests the number of alternatives, capped by available support;
`exact = TRUE` requests the bounded supported enumeration. Seeded sampling restores
the caller's random-number state. Full serialized assignments remain values in
short hash buckets, with full-key comparison, so long identifiers and hash
collisions cannot silently merge different assignments.

    (1 + number of alternative absolute statistics >= observed) / (B + 1).

Inclusive comparisons allow only a roundoff tolerance of 64 machine eps times
the maximum of one and the two statistic magnitudes. An invalid required draw
makes the affected p-value unavailable; inference is not based on surviving draws.
Randomization currently supports a zero sharp null. A defensible assignment
mechanism remains a substantive research assumption.

## Generate and replay the same assignments

After installing `didsuiteR` and `fasttwfeR`, this small example runs all five
procedures and reuses the exact assignments:



Nineteen alternatives keep this demonstration short and give a coarse 0.05
p-value grid; use an appropriately larger count for research. Save the complete
assignment object with `saveRDS(randomized$assignments, "assignments.rds")` and
reload it with `readRDS()` to preserve its identifiers and hash. Supplied matrix
columns must match the suite's returned unit order, which uses sorted character
cluster identifiers. Seeds alone do not imply matching draws across languages;
share the actual assignments for cross-language comparisons.
