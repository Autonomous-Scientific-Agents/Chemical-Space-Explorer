# ADR 0005: Rule set for chemical stability

- Status: accepted
- Date: 2026-10-14

## Context

The user asked for "basic chemistry rules" only. In the design session they chose four groups: valence and connectivity, bond-order limits, GDB-style strain filters and heteroatom filters.

## Decision

Implement `RuleConfig` (`rules.py`), with each group toggleable from the CLI:
- **Valence:** C4, N3, O2, S2, P3, Si4, B3, halogens 1. `--hypervalent` adds S 4/6 and P 5.
- **Bond orders:** at most 3; O and S at most 2; halogens 1.
- **Strain:**
  - no double bond in a 3- or 4-membered ring;
  - no triple bond in a ring of fewer than 8 atoms;
  - no cumulated bonds at an atom in a ring of fewer than 9 atoms;
  - approximate Bredt rule: no multiple bond at a bridgehead of rings ≤ 7.
- **Heteroatoms:**
  - no O–O or N–O bonds (`--allow-n-o` lifts the N–O ban);
  - no bond from a halogen to N, O, S or another halogen;
  - no enols, hemiacetals, hemiaminals or gem-diols.
- **Validation switches:** `--no-rings` and `--no-multiple-bonds`.

Element-pair bans are hereditary, so they are enforced in the tree. They are also re-checked in `skeleton_ok` inside the payload, so hand-built skeletons are filtered too; this was a bug fix on 2026-10-14.

## Alternatives

- GDB-17's full functional-group filter list. Out of scope for "basic rules"; it may come back as an optional rule set.
- Energy-based strain from force fields. Needs 3D coordinates and is far too slow.

## Consequences

Counts are about 25× GDB-11 at 11 atoms, which is expected given GDB's stricter filters. The strain and Bredt rules are heuristics.
