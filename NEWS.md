# hdrm 1.0.0

## Methods

- Added one-group inference for high-dimensional repeated-measures data.
- Added multiple-group inference for heterogeneous covariance matrices.
- Added a multiple-group procedure under equal covariance matrices.
- Added predefined and user-supplied projection-matrix hypotheses.
- Added exact and subsampling-based trace estimators.
- Uses `B` as a base subsampling budget in grouped procedures; both grouped
  third-trace estimators use `a * B` draws, with the equal-covariance
  budget distributed across groups.

## Reliability and interface

- Standardized wide matrix input to the usual R orientation: subjects are
  represented by rows and repeated-measurement dimensions by columns. Returned
  processed data matrices use the same orientation; internal calculations
  continue to use the transposed representation.
- Added explicit validation of data, hypotheses, group and subject identifiers,
  subsampling budgets, and random seeds.
- Added subject-wise handling of incomplete repeated-measures data.
- Added reproducible stochastic calculations through the `seed` argument.
- Added numerically stable upper-tail p-value calculations.
- Added safeguards for degenerate trace estimates and undefined test
  statistics.
- Added independent R reference tests for the principal R and C++ estimators.
