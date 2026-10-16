# ADR 0003: Two-layer model of skeleton tree plus bond-order payload; a molecule is the skeleton plus H counts

- Status: accepted
- Date: 2026-10-14

## Context

Enumerating molecules directly, as graphs with bond orders, mixes two kinds of symmetry and makes rules that apply to fragments hard to use for pruning.

## Decision

- **Layer 1, skeletons:** connected graphs with element labels and no bond orders. These are the tree nodes. Only *hereditary* rules prune here: degree ≤ max valence, banned element pairs, and optionally no rings.
- **Layer 2, payload:** for each skeleton, enumerate bond orders 1–3 by backtracking, derive implicit H counts and apply the full rules.
- **Molecule identity:** the canonical form of the skeleton relabeled with (element, H count). Kekulé and resonance structures, which have the same H pattern, collapse into one molecule. A molecule passes if *any* of its bond assignments passes.

## Alternatives

- One layer with bond orders in the tree. That gives a bigger tree and duplicate resonance forms.
- Aromaticity perception. It needs a model of aromaticity that we would have to invent.

## Consequences

- Benzene counts once.
- The payload is cached per skeleton.
- Exporting a molecule needs one stored passing bond assignment (PRD 0001).
