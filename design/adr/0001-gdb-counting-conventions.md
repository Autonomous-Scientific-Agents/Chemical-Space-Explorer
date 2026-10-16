# ADR 0001: Count molecules with GDB conventions

- Status: accepted
- Date: 2026-10-14

## Context

"The number of stable molecules" is meaningless until we say what a molecule is and how size is measured. The user asked for rough estimates from a maximum atom count and an element list. In the design session they chose heavy atoms as the size and constitutional isomers as the identity.

## Decision

Use the conventions of the GDB databases (GDB-11/13/17):
- size is the number of **heavy (non-hydrogen) atoms**;
- hydrogens are implicit;
- only **constitutional isomers** count, with no stereoisomers;
- molecules are neutral, closed-shell and connected.

## Alternatives

- Count all atoms, including H. This makes sizes incomparable with GDB and inflates N.
- Count stereoisomers. That needs 3D or parity perception and roughly multiplies counts by 2^(stereocenters). It is out of scope for a rough estimate.
- Include ions and radicals. This would make the count enormous and the rules ill-defined.

## Consequences

Numbers can be compared in order of magnitude with GDB-11 (2.6×10⁷ for C, N, O, F up to 11 atoms). Tautomers count as separate molecules.
