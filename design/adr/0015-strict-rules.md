# ADR 0015: Optional "strict" rule tier (S1–S6)

- Status: accepted. The basic tier of ADR 0005 is unchanged and remains the default.
- Date: 2026-10-17

## Context

The first public dataset (`chemuniverse-CNOS-le12`, ADR 0014) uses the basic rules. Its samples include chemistry a public user would not read as "plausible molecules":
- heteroatom chains such as `S(N(N(N(NC(S)=O)C)C)C)C`;
- unstable tautomers such as `C(=S)O` next to `C(S)=O`;
- thioaldehydes, enethiols and ynamines.

The user asked for stricter filters. Each candidate was measured on all 248,281 basic C·N·O·S molecules up to 7 atoms before being chosen.

## Decision

Add six individually switchable `RuleConfig` fields, all off by default. `RuleConfig.strict()` and `--rules strict` turn them all on. Heteroatoms here means N, O and S.

| Rule | Field | Removes | Kind |
|---|---|---|---|
| S1 heteroatom chains: no heteroatom with ≥ 2 heteroatom neighbours | `hetero_chains` | X–X–X (N–N–N, S–N–N, …) | hereditary: pruned in the tree |
| S2 no S–O bond | `ban_s_o` | sulfenic acids and esters | hereditary (element pair) |
| S3 no 3- or 4-membered ring with ≥ 2 heteroatoms | `small_ring_heteroatoms` | diaziridines, thiaziridines, 1,2-diazetidines, 1,3-dithietanes, … | hereditary: pruned in the tree |
| S4 unstable tautomers | `tautomers` | enethiol HS–C=C, NH-enamine HN–C=C, ynol/ynamine/ynethiol X–C≡C, imidic acid HO/HS–C=N, thiono acid HO–C=S, ketene C=C=O | molecule |
| S5 hemithioacetals | `thio_gem` | an sp³ C with ≥ 2 heteroatoms, one of them S, one of them with an H | molecule |
| S6 thiocarbonyls whose only heteroatom is the S | `thiocarbonyl` | thioaldehydes, thioketones, thioketenes | molecule |

**Deliberately kept:**
- isolated N–N, S–S and S–N bonds (hydrazines, hydrazones, azo compounds, disulfides, sulfenamides);
- cyanamides, cyanates and thiocyanates (X–C≡N);
- amides, thioamides, thioesters, thioureas, CS₂ and COS;
- orthoesters, and acetals without OH/NH/SH.

**Why pruning is exact:** deleting a vertex can neither create an X–X–X chain nor a small ring. So a valid skeleton's canonical parent is valid, and an invalid skeleton has no valid descendants. `rules.hetero_ok` therefore checks only the new vertex of each child, in both `tree.children` and the Rust core.

## Effect (exact)

| Elements | n = 4 | 6 | 7 | 8 | skeletons at 8, strict vs basic |
|---|---|---|---|---|---|
| C,N,O,S | 195 (42%) | 7,517 (30%) | 58,431 (27%) | 512,270 (24%) | 318,195 vs 1,055,896 |
| C,N,O,F | 160 (72%) | 5,432 (62%) | 39,534 (59%) | 324,498 (56%) | 155,598 vs 316,644 |

## Verification

- `tests/test_rules_strict.py` checks one kept and one removed case per rule. Its brute-force test checks that the counts from the pruned tree equal the basic molecules filtered by an independent re-implementation of S1–S6, for every size up to 6 atoms (C,N,O,S).
- The Rust core implements the same rules; strict configurations are in `tests/test_backend_parity.py`.

## Alternatives

- **Banning all heteroatom–heteroatom bonds,** closer to GDB. It would remove real hydrazines and disulfides; kept as an open option.
- **Filtering the published basic dataset by parsing SMILES.** Re-running the export with `--rules strict` is deterministic and cheaper, because the tree is about 3× smaller.

## Consequences

- Datasets and counts must state their tier. Strict datasets are named `…-strict`.
- These are still heuristics, not stability calculations. S3 also removes the rare but known 1,3-dithietane.
