# Validation Milestone

2026-09-29, R 4.5.3 Windows. OMP/BLAS threads capped at one; two-process smoke
test only. R CMD build and R CMD check --no-manual: **Status OK**, zero
errors/warnings/notes. tests/contracts.R, tests/routing.R and all 53 named
tests/smoke.R checks passed.

The 549-row shared stochastic RCS fixture has 14 clusters and nine times.
The three supplied alternatives preserve observed cohort counts, exclude observed
assignment and are unique. Canonical matrix SHA256:
bea08f08f138e2343a8e1c618a7542169a9c56f966dd174dbac3d240abc3b4ec

validation/master.R runs all nine families in one call with only Stacked
window(-1,1), both controls, and LP window(0,1). Eight point estimates and all
available native/JK SE and three RI p-values agree with parent Python CSV:
maximum absolute gap 1.95e-14. All valid RI p-values are .25.
UN-DID/DID-INT rows are PARTIAL due to unsupported native VCE/RI-t, matching
Python numeric missingness but deliberately stricter status. LP is unsupported
on RCS, not missing implementation; panel routing and full inference are tested.

The CSV/JSON exports in validation/ retain actual results. Fixtures are external
sibling didsuite_spec inputs, never modified or redistributed here. Source
packages are independently installed. Canonical CSDID checks compare native
package outputs, not a fresh external official-oracle certification.

Remaining: broader adversarial/external-oracle coverage, harmonized common target
and sample intersection, distinct assignment/inference units, full package-option
parity, optimized production-size performance. GitHub CI only runs public harness
dependencies unless estimator backends are separately installed.
