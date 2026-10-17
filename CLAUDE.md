# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Status

The design comes from the plan in `~/.claude/plans/goal-is-to-write-vivid-walrus.md`. The code follows it, with one addition: `skeleton_ok` re-checks the hereditary rules inside the payload, so hand-built skeletons are also filtered.

## Purpose

Rough estimator of the number of chemically stable molecules, given a maximum **heavy-atom** count N and an element list:
- Counts follow GDB-17 conventions: hydrogens are implicit, and only constitutional isomers count (no stereoisomers).
- Output is a per-size count for every n ≤ N, plus the cumulative total.

Stack: pure Python with no runtime dependencies (numpy/networkx are allowed but unused; rings come from `graph.simple_cycles`). **No RDKit**, which is unavailable and must not be added. There is an optional Rust core in `rust/` (PyO3 + maturin), and the Python code remains the reference implementation.

## Commands

```bash
pip install -e ".[test]"              # no runtime deps; pytest for tests
pytest                                # all tests
pytest tests/test_counts.py::test_alkanes   # single test
ruff check --line-length 120 .
python -m chemuniverse --max-atoms 12 --elements C N O F --exact-up-to 6 --samples 10 --seed 0 --jobs 10   # --method knuth for probes
python -m chemuniverse export --max-atoms 8 --elements C N O F --out DIR --jobs 10   # SMILES dataset
python scripts/bench.py --jobs 1 10                      # performance numbers for design/prd/0002
pip install ./rust                                        # optional Rust core (needs rustup + maturin); rebuild after editing rust/
CHEMUNIVERSE_BACKEND=python pytest                        # force pure Python; =rust requires the core
```

## Design records

Decisions live in `design/adr/` and features in `design/prd/`, with an index in `design/README.md`. Whenever a significant decision is made (method, rules, dependencies, publishing, performance), add an ADR, or supersede an old one by adding a new ADR and updating the old one's status line, in the same change.

## Architecture (package `chemuniverse/`)

There is one generation tree, and it is used in two ways:

- **Tree nodes** are labeled skeletons: connected simple graphs with element labels.
  - They are stored as a tuple of labels plus adjacency bitmasks (`graph.py`).
  - Only *hereditary* constraints prune the tree: degree ≤ max valence, and no banned heteroatom pairs.
- **Children** (`tree.py`) come from McKay canonical augmentation:
  - Add one atom bonded to a non-empty subset of the existing atoms.
  - Deduplicate children by canonical form.
  - Accept a child only if the new vertex is in the same automorphism orbit as the canonically chosen non-cut vertex.
  - The canonical vertex minimizes (element, degree, refined color), then canonical position. Cheap keys are checked first in `_accept`; this ordering is the main speed lever.
  - Siblings are deduplicated by `invariant_hash`, with canonicalization only on hash collisions.
  - The samplers use `cached_children` / `cached_count` (LRU caches).
  - This yields each isomorphism class exactly once.
- **Canonical labeling** (`canon.py`) is written by hand in pure Python: color refinement, then an individualization–refinement search. Its API is `canonical_form()` and `same_orbit()`.
- **Payload per node** (`bonds.py`, `count_molecules`):
  - Enumerate bond orders 1–3 per edge by backtracking, then apply the full rules (`rules.py`).
  - Each molecule is keyed by the canonical form of the skeleton labeled with `(element, H count)`, so Kekulé/resonance structures collapse into one molecule.
- **Exact mode** is a breadth-first traversal of the tree up to `k_exact`, which defaults to 6.
- **Parallelism** (`jobs`, through `tree._pmap` with a `fork` context): exact subtrees and sampling passes or probe batches run in worker processes.
  - Every pass or probe gets its own seed drawn from the caller's RNG, so results never depend on `jobs`.
  - Keep that property, because `tests/test_parallel.py` checks it.
- **Rust core** (`rust/`, optional package `chemuniverse_core`, selected in `backend.py`):
  - It is a line-by-line port of `graph`, `canon`, `rules`, `bonds`, the `tree` children and `export.to_smiles`. It must give the *same* children in the same order, the same counts and the same SMILES as Python, not merely equal totals.
  - `tests/test_backend_parity.py` enforces this. Any algorithm change goes into both backends.
  - Rust provides `children`, `count_molecules`, `molecule_smiles` and `expand` (exact mode, scoped threads).
  - The sampler loops stay in Python, so estimates are bit-identical across backends.
  - It must not leave threads running, because Python forks afterwards.
  - Python functions keep `_py`/`_rust` suffixed variants (`children_py`, `count_molecules_py`, `molecule_smiles_py`).
- **Shards** (`shards.py`, `python -m chemuniverse shard plan|run|reduce|finalize`) are the multi-node layer.
  - Exact/export inputs are stored in blocks of 1000 shards (`inputs/<block>.json.gz`), so a worker never loads the whole frontier.
  - `finalize` (export only) merges the SMILES parts into ~2 GB checksummed chunks and writes `dataset/manifest.json` (ADR 0014).
  - `plan` writes `plan.json` plus a gzip frontier and assignment. Each shard writes its result file atomically; `reduce` must equal the single-machine functions, which `tests/test_shards.py` checks.
  - Crux: `scripts/hpc/run_ensemble.py` (ensemble_launcher), `scripts/hpc/crux/{setup_env.sh,run_shards.sh,iri_job.py}`.
  - Submit through IRI with `~/git/genesis-skills/skills/amsc-skills/iri-api/scripts/iri_api_call.py`. The Crux resource ID is `8b9b42f7-572a-4909-8472-a0453436304c`. Crux jobs clone the repo from GitHub, so push before submitting.
- **Export** (`export.py`) does a depth-first walk of the same tree and writes `molecules(..., with_orders=True)` as Kekulé SMILES (`to_smiles`).
  - It splits into a fixed 128 subtree tasks, so the files are identical for any `jobs`.
  - `tests/test_export.py` parses the SMILES back and checks that every canonical form is unique and that the counts equal the exact counts.
- **Estimate mode** covers levels deeper than `k_exact`. `estimate(..., method=)` chooses the sampler.
  - `stratified_estimate` (the default) is Chen's stratified sampling. Per level, it keeps one weighted representative per `stratum()`, where a stratum is (atoms per free-valence count, ring count).
  - `knuth_estimate` uses random probes. It is heavy-tailed, so its error bars are unreliable.
  - Both return independent samples, which `_summarize` turns into mean ± SE.
- **Rules** (`RuleConfig`, each group toggleable):
  - valence and connectivity;
  - bond-order limits;
  - GDB-style strain filters: small-ring triple bonds and allenes, double bonds in 3- or 4-membered rings, an approximate Bredt rule;
  - heteroatom filters: banned O–O, N–O and halogen–heteroatom bonds, no enols, no hemiacetals or hemiaminals.
  - There are validation switches `allow_rings` and `allow_multiple_bonds`.

## Web explorer (`docs/`, GitHub Pages)

- `docs/index.html` is a static page that reads `window.CHEMUNIVERSE_DATA` from `docs/data.js`.
- `docs/data.js` is generated by `scripts/build_site_data.py`; never edit it by hand.
- Data keys look like `CNO_h2`: C plus the non-halogens present, plus the number of halogens. Halogens are interchangeable under the rules, so a key's counts don't depend on which halogens are chosen.
- Beyond the computed sizes, the page extrapolates in JS: `series()` in `index.html`.
- `docs/method/index.html` is a hand-written teaching page about the method. Its example counts and rule-effect table are hard-coded from `enumerate_exact` output; keep them in sync if the rules change. `docs/mode/` redirects to it.
- `docs/agent-traces/` holds published session replays (HTML only).

## Correctness anchors for tests

- Alkanes (C only, no rings, no multiple bonds) for n = 1..10 must give 1,1,1,2,3,5,9,18,35,75 (OEIS A000602).
- C n=2 → 3, and C,O n=2 → 5.
- Benzene counts once.
- For n ≤ 6, both samplers must fall within 4 SE of the exact count.
- A C,N,O,F run up to 11 atoms gives ≈6×10⁸, about 25× GDB-11 (2.6×10⁷).
- Exact C,N,O,F counts computed on Crux (ADR 0012):
  - n = 7–12: 67,403; 583,112; 5,526,145; 56,997,371; 634,889,614; 7,605,193,456.
  - C,N,O,S exact: n = 1–9: 4, 19, 80, 463, 3,155, 24,994, 219,566, 2,129,403, 22,382,806.
  - Use them to check samplers. ADR 0013 notes that small-pass stratified estimates tend to come out a few percent low.
- `tests/test_counts.py::test_tree_matches_brute_force` checks the canonical augmentation against brute-force enumeration.
