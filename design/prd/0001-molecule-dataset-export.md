# PRD 0001: Molecule dataset export (GDB-style SMILES files)

- Status: implemented 2026-10-16 (`chemuniverse/export.py`); open questions remain
- Date: 2026-10-16
- Decisions: ADR 0003, ADR 0004, ADR 0009

## Problem

Exact enumeration already builds every molecule exactly once, but keeps only a count. Users want the molecules themselves, like the GDB datasets, for machine learning, screening or inspection.

## Goals

- Write every molecule from an exact enumeration as SMILES, exactly once.
- Support any element set and rule configuration the counter supports.
- Stream to disk: memory must stay flat however many molecules there are.
- Run in parallel across cores, with output independent of `--jobs`.

## Non-goals

- Canonical SMILES. Each molecule is emitted once, so canonical strings aren't needed for deduplication.
- Aromatic SMILES. We write one Kekulé form.
- GDB's exact filter set. That would be a separate ADR if wanted.
- Exporting from the sampler (weighted samples). Possible later.

## Requirements

1. `python -m chemuniverse export --max-atoms N --elements … --out DIR [--jobs J] [rule flags]`.
2. Output is `DIR/n{size}_part{k}.smi.gz`, one SMILES per line, plus `DIR/manifest.json` (elements, rules, per-size counts, package version, date).
3. SMILES:
   - organic-subset atoms are written without H, since implicit H follows the same lowest-allowed-valence rule we use;
   - other elements are written in brackets with an explicit H count;
   - bond orders are written explicitly (`=`, `#`), with ring-closure digits.
4. The number of lines per size equals `enumerate_exact` counts.

## Success metrics

- The test round-trip SMILES parser rebuilds each molecule with an identical canonical form, and every form is unique.
- Line counts match the exact counts (C up to 8 atoms: 1, 3, 5, 17, 55, 229, 1035, 5599).
- Throughput is at least the counting throughput times the number of cores, minus write overhead.

## Open questions

- Optional GDB-17 filter set?
- Other output formats, such as Parquet with per-molecule descriptors?
