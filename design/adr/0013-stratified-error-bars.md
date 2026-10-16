# ADR 0013: Treat small-pass stratified error bars as optimistic; use hundreds of passes for reported numbers

- Status: accepted (finding); a site-data rebuild with more passes is proposed
- Date: 2026-10-16

## Context

The first exact counts beyond 7 atoms came from Crux (ADR 0012): 634,889,614 C·N·O·F molecules at 11 atoms. Earlier stratified estimates at 11 atoms were all below it, by 1.5–3 reported standard errors:
- 5.73×10⁸ ± 4.2×10⁷ (10 passes);
- 5.69×10⁸ ± 3.3×10⁷;
- 5.53×10⁸ ± 5% (site data, 12 passes).

## Evidence

400 passes from the 6-atom frontier, compared with exact counts:

| n | mean / exact | z | median pass / exact | max pass / mean |
|---|---|---|---|---|
| 7 | 0.996 | −0.6 | 0.989 | 1.5 |
| 9 | 0.988 | −1.4 | 0.976 | 1.7 |
| 11 | 0.990 | −0.8 | 0.960 | 1.9 |

The estimator is unbiased within error, as expected from Chen's construction. But the distribution of single passes is right-skewed: the median pass is about 4% low at 10–11 atoms, and rare passes reach about 2× the mean. With about 10 passes, a run usually misses those rare high passes, so it comes out low and its sample standard error is too small.

## Decision

- Numbers we report or publish use hundreds of passes, which is cheap on Crux with `shard plan --kind passes`.
- Error bars from fewer than about 50 passes are described as optimistic.
- Tests keep comparing samplers with exact counts. Exact counts up to 11 atoms (C·N·O·F) are now available as references.

## Consequences

- The explorer's 12-pass data is probably several percent low at 11–17 atoms. That doesn't change any order of magnitude.
- Proposed: rebuild the site data with about 200–500 passes per element set on Crux. That's about 40 sets × 300 passes × ≤ 50 core-seconds, roughly 170 core-hours.
