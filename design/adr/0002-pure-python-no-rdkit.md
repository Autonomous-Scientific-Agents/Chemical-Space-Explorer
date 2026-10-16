# ADR 0002: Pure-Python implementation without RDKit

- Status: accepted. ADR 0009 adds an optional Rust core; the pure-Python code stays the reference.
- Date: 2026-10-14

## Context

RDKit is not available in the target environment. In the design session the user said "numpy is ok, use networkx if necessary".

## Decision

Write everything in pure Python: canonical labeling, ring finding and rules. networkx was first used for cycle finding. It was dropped on 2026-10-15 (commit `2ca9f96`) in favour of a bitmask DFS (`graph.simple_cycles`), and numpy turned out to be unnecessary. There are no runtime dependencies.

## Alternatives

- RDKit for canonical SMILES and sanitization. Not available, and too slow per molecule for enumeration anyway.
- nauty for canonical labeling. It is a C dependency, and molecules this small don't need it.

## Consequences

The package is easy to install and read, but it runs at about 130 µs per molecule. That limits exact enumeration to about 7–8 heavy atoms (see ADR 0009).
