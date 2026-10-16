# Design records

This folder records why the project is built the way it is.

- **ADRs** (`adr/`, Architecture Decision Records) cover one decision each: its context, the decision itself, the alternatives considered and the consequences. Never rewrite an accepted ADR. To change a decision, add a new ADR that supersedes it and update the old one's status line.
- **PRDs** (`prd/`, Product Requirements Documents) describe a feature before or while it is built: the problem, goals and non-goals, requirements, success metrics and open questions.

To add a record, copy the structure of an existing one and use the next number. Then add a line to the index below.

## Decisions

| # | Decision | Status |
|---|---|---|
| [0001](adr/0001-gdb-counting-conventions.md) | Count with GDB conventions (heavy atoms, constitutional isomers) | accepted |
| [0002](adr/0002-pure-python-no-rdkit.md) | Pure Python, no RDKit | accepted |
| [0003](adr/0003-skeleton-tree-and-payload.md) | Skeleton tree plus bond-order payload; molecule = skeleton + H counts | accepted |
| [0004](adr/0004-canonical-augmentation.md) | McKay canonical augmentation, hand-written canonical labeling | accepted |
| [0005](adr/0005-stability-rules.md) | Stability rule set | accepted |
| [0006](adr/0006-stratified-sampling.md) | Hybrid exact + stratified sampling (Knuth demoted) | accepted |
| [0007](adr/0007-static-explorer-github-pages.md) | Static explorer on GitHub Pages | accepted |
| [0008](adr/0008-agent-traces.md) | Sanitized agent-trace replays, HTML-only publishing | accepted |
| [0009](adr/0009-performance-strategy.md) | Processes + Rust core; no GPU, no threads | accepted, refined by 0010 |
| [0010](adr/0010-rust-core-design.md) | Rust core: exact port, sampling loops in Python, scoped threads, mimalloc | accepted |

## Features

| # | Feature | Status |
|---|---|---|
| [0001](prd/0001-molecule-dataset-export.md) | Molecule dataset export (SMILES) | implemented |
| [0002](prd/0002-performance.md) | Performance | stages 1 and 2 done |
