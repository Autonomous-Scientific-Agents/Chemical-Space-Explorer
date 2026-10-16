# ADR 0006: Hybrid exact plus sampled estimation, with stratified sampling (Chen) as the default

- Status: accepted. It supersedes the original plan's Knuth-only estimator.
- Date: 2026-10-14

## Context

Exact enumeration is practical only to about 7 atoms. The approved plan used Knuth's random-probe estimator below the exact frontier. In testing it was heavy-tailed: at 10 atoms, 300 probes gave 4.96×10⁷ and 1000 probes gave 1.59×10⁸ ± 1.3×10⁸, so the error bars were unreliable.

## Decision

- Count exactly up to `k_exact` (default 6).
- Below that, use **Chen's stratified sampling** by default.
  - Per level, keep one weighted representative per stratum. A stratum is (number of atoms with each free-valence count, ring count).
  - The representative is replaced with probability w/total, which keeps the estimate unbiased.
  - Each pass is one independent estimate, and the reported SE is computed across passes.
- Knuth's estimator stays available as `--method knuth`.

## Alternatives

- Knuth with many more probes. The heavy tail persists, and the error bars stay unreliable.
- Regression or extrapolation from exact counts only. It can't be checked.

## Consequences

- CNOF up to 11 atoms gives 6.29×10⁸ ± 4.2×10⁷ (seed 0) and 6.23×10⁸ ± 3.9×10⁷ (seed 1), which agree.
- A pass costs more than a probe, and the cost per level about doubles at each level (ADR 0009).
