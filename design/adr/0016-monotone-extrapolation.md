# ADR 0016: Explorer extrapolation is monotone in the element set

- Status: accepted. Amends the extrapolation in ADR 0007.
- Date: 2026-10-30

## Context

`series()` in `docs/index.html` extrapolates each data key from the growth factor between the last two computed sizes (16 → 17 atoms). At 17 atoms the sampled counts have 4–34% relative standard error, so that factor carries noise of a few hundredths in log₁₀. Over ~100 extrapolated atoms this becomes several orders of magnitude, enough to reverse the order of element sets. With the GDB-like preset at 120 atoms, CNOS plus 4 halogens gave 10¹⁰⁸, plus 3 halogens 10¹⁰⁴·⁵ and plus 2 halogens 10¹⁰⁵·⁶: removing an element *raised* the estimate. This cannot be true, because every molecule over a subset of elements is also a molecule over the larger set.

Averaging the factor over the last 3 or 5 sizes reduces the noise but does not remove the reversals.

## Decision

- `rawSeries(key)` keeps the ADR 0007 extrapolation unchanged.
- `series()` returns, at each size, the row with the largest central value among the selected key and all of its subset keys: any subset of the selected N, O, S, and any halogen count up to the selected one.
- This guarantees that adding an element never lowers the shown count, at any size.

## Consequences

- Within the computed range the effect is small (it removes sampling-noise reversals). In the extrapolated range a set can show the curve of one of its subsets.
- The noise itself remains. The large-size numbers are still order-of-magnitude illustrations, and the ordering is now consistent.
