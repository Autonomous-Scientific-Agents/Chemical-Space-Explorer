# ADR 0004: McKay canonical augmentation with hand-written canonical labeling

- Status: accepted
- Date: 2026-10-14

## Context

Each skeleton must be generated exactly once without a global hash set of everything seen. That rules out global memory, and the tree must be splittable so it can be sampled and parallelized.

## Decision

- Grow skeletons by adding one atom bonded to a non-empty subset of the existing atoms.
- Accept a child only if the new atom lies in the automorphism orbit of the child's *canonical deletion vertex*. That is the non-cut vertex that minimizes (element, degree, refined color), with ties broken by canonical position.
- Deduplicate siblings with a Weisfeiler–Lehman hash, canonicalizing only on collisions.
- Canonical labeling (`canon.py`) works by color refinement followed by an individualization–refinement search, with pruning by the automorphisms found along the way.
- The cheap keys are checked first in `_accept`. This ordering was the main 3× speedup (commit `2ca9f96`).

## Alternatives

- Orderly generation (Read/Faradzev). It needs a canonicity test at every node, which is slower.
- A global canonical-form set. Memory grows with the count, and the tree can't be sampled.

## Consequences

- The tree is exact and contains no duplicates; this is tested against brute force for ≤ 4 atoms and against known graph counts.
- Subtrees are independent, so exact enumeration parallelizes without coordination.
