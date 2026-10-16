# ADR 0011: Compute the explorer data up to 17 heavy atoms

- Status: accepted. It amends the parameters in ADR 0007.
- Date: 2026-10-16

## Context

ADR 0007 computed the site data to 12 atoms with 8 passes, which took about 15 minutes on 10 cores in pure Python. The Rust core (ADR 0010) makes the same sampling about 90× faster on 10 cores, so the computed range could be extended. The user approved rebuilding the data on 2026-10-16.

## Decision

- `build_site_data.py --max-atoms 17 --samples 12`, with exact counts to 6 atoms. That is 40 element sets and about 14,500 core-seconds, about 25 minutes on 10 cores.
- 17 atoms is chosen so that the GDB-17 element set can be compared directly at a computed size, not an extrapolated one.

## Results and checks

- Relative standard error of the total ranges from 3.5% to 32%. The largest element set (C, N, O, S and 4 halogens) is the hardest; per size at 17 atoms the error ranges from 4% to 34%.
- C, N, O, F up to 11 atoms: 6.11×10⁸, which agrees with the reference value of about 6.3×10⁸.
- One subset/superset inversion between totals (CNO_h0 against CNO_h1) at z = 0.1, which is noise.
- The GDB-17 element set gives 1.29×10¹⁷ molecules up to 17 atoms. GDB-17 lists 1.66×10¹¹, so our looser rules give about 10⁶× more.

## Consequences

- The growth factor per atom keeps rising through 17 atoms: carbon goes from 5.4 to 13, and the largest set from 10 to 18.5.
- The page's "upper edge" extrapolation, which holds the last factor constant, is therefore not a true bound. At 20 atoms the GDB-like preset moved from about 4×10¹⁸ (12-atom data) to about 3×10²⁰ (17-atom data), above the earlier envelope.
- The method page now says this. A growth model in which the factor itself keeps rising would be a separate decision.
