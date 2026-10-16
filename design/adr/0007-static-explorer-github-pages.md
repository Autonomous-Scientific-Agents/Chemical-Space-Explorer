# ADR 0007: Static web explorer on GitHub Pages with precomputed data

- Status: accepted. The computed range was extended to 17 atoms by ADR 0011.
- Date: 2026-10-15

## Context

The user supplied an `index.html` visualization of chemical space and asked to integrate it and publish it with GitHub Pages. Its original growth model was a heuristic interpolation of GDB numbers.

## Decision

- Serve Pages from `main` `/docs`, with no backend.
- `scripts/build_site_data.py` precomputes 40 runs: C plus any subset of N, O and S, times 0–4 halogens. Halogens are interchangeable under the rules, so 128 element choices reduce to 40. It writes `docs/data.js`. Parameters: max 12 atoms, exact to 6, 8 passes.
- The page shows computed rows, with a ±2 SE band in the sampled range. Beyond 12 atoms, `series()` extrapolates in JavaScript:
  - the high edge keeps the last growth factor per atom;
  - the low edge decays toward the alkane rate (2.9529), floored at the alkane counts;
  - the central value is the midpoint in log space.
  The extrapolation is labeled as such.
- The GDB interpolation model was removed; GDB-13/17/20 remain as reference points.
- `docs/method/` (added 2026-10-15) is a teaching guide whose example numbers are hard-coded from `enumerate_exact` output. `docs/mode/` redirects to it.

## Alternatives

- Running Python in the browser (Pyodide). Too slow beyond tiny sizes.
- A server backend. It needs hosting and adds nothing for precomputable data.

## Consequences

Updating the site means re-running the build script, which takes about 15 minutes on 10 cores. Data beyond 12 atoms is extrapolated, not computed.
