# Chemical-Universe

Rough estimates of the number of chemically stable molecules made of a given set of elements, for every size up to a maximum number of heavy atoms. Counting follows GDB conventions:
- hydrogens are implicit;
- molecules are neutral, closed-shell and connected;
- only constitutional isomers count (no stereoisomers).

```bash
pip install -e ".[test]"
python -m chemuniverse --max-atoms 11 --elements C N O F --exact-up-to 6 --samples 10 --seed 0
```

```
  n    molecules       ± SE     mode   cumulative
  1            4               exact            4
  ...
  6         8738               exact    1.027e+04
  7    6.859e+04    2.7e+03 estimate    7.886e+04
 ...
 11    6.601e+08    6.9e+07 estimate    7.274e+08
total: 7.274e+08 ± 7e+07   (298.1 s)
```

From Python: `chemuniverse.estimate(max_atoms, elements, RuleConfig(...), k_exact, samples, seed, method)` returns one result per size.

## Method

1. **Skeleton tree.** The nodes are connected graphs whose vertices are labeled by element.
   - Each child adds one atom, using McKay canonical augmentation, so every skeleton appears exactly once.
   - Skeleton-level rules prune the tree as it is built: degree ≤ valence, banned element pairs, and optionally no rings.
2. **Payload.** For each skeleton, the program:
   - enumerates bond orders 1–3 by backtracking;
   - derives implicit H counts;
   - applies the rules;
   - counts the distinct `(element, H)` labelings up to symmetry.

   Resonance/Kekulé structures (e.g. benzene) therefore count once.
3. **Exact counts** up to `--exact-up-to` heavy atoms come from a full traversal of the tree.
4. **Estimates** for larger sizes come from sampling the tree below the exact frontier.
   - **Stratified sampling** (default, Chen 1992): each pass works level by level.
     - Children are grouped by a coarse type: the number of atoms with each free valence, plus the ring count.
     - One weighted representative is kept per group. Its weight is the total of everything merged into it, and it is swapped in with probability w/total.
     - Each pass is an unbiased estimate; the standard error comes from the spread across `--samples` passes.
     - With C, N, O, F the default 10 passes give about 10% error at 11 atoms.
   - **Knuth probes** (`--method knuth`): random root-to-leaf walks, each weighted by the product of branching factors along the way.
     - Unbiased and much cheaper per sample, but heavy-tailed: rare huge-weight probes dominate.
     - At equal cost its error bars are unreliable. A 300-probe and a 1000-probe run to 10 atoms disagreed by 3×.

## Rules

Each group can be turned off from the command line.

| Group | Rules | Flag |
|---|---|---|
| Valence | C4, N3, O2, S2, P3, Si4, B3, halogens 1; optional S 4/6 and P 5 | `--hypervalent` |
| Bond orders | ≤ 3; no triple bonds on O, S or halogens | `--no-bond-order-limits` |
| Strain | no triple bond in a ring < 8 atoms; no double bond in a 3- or 4-membered ring; no cumulated bonds at an atom in a ring < 9; Bredt (no multiple bond at a bridgehead of rings ≤ 7) | `--no-strain` |
| Heteroatoms | no O–O, N–O or halogen–(N, O, S, halogen) bonds; no enols; no hemiacetals, hemiaminals or gem-diols | `--no-heteroatom`, `--allow-n-o` |
| Validation | acyclic only / single bonds only | `--no-rings`, `--no-multiple-bonds` |

## Checks

`pytest` checks:
- alkane counts against OEIS A000602;
- carbon skeleton counts against the known numbers of connected graphs with max degree 4;
- the generation tree against brute-force enumeration of every labeled graph with ≤ 4 atoms;
- hand counts for small molecules and for each rule;
- the estimator against exact counts.

With C, N, O, F up to 11 atoms the total is about 7×10⁸. This is about 30× GDB-11 (2.6×10⁷), which uses much stricter functional-group and ring filters.

## Limitations

- Tautomers count as separate molecules (apart from the enol ban).
- There are no charged species or radicals.
- The strain and Bredt rules are heuristics.
- Exact enumeration is pure Python and practical up to about 7 heavy atoms (6 takes about 2 s for CNOF).
- Stratified passes get slower with depth, because the number of strata grows: about 30 s per pass from 6 to 11 atoms for CNOF.
