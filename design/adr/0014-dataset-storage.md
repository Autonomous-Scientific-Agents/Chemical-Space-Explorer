# ADR 0014: Dataset storage on ALCF Eagle; public release; C, N, O, S first

- Status: accepted (sizes and costs are estimates until the calibration runs)
- Date: 2026-10-16

## Context

The user wants to publish exported SMILES datasets for everyone, starting with **C, N, O, S**, and pointed to two storage areas:
- `/lus/eagle/projects/datascience/keceli`, with more than 100 TB available;
- IQC's `/lus/grand/projects/IQC/keceli`, with more than 5 TB.

A probe job (362501, Crux, `filesystems=home:grand:eagle`) showed:
- Crux compute nodes mount `/lus/eagle` read-write: 500 TB total, 188 TB free.
- On Crux, `/lus/grand` is a 10 TB clone hosted on Eagle (`/eagle/clone/g2`), with 5.6 TB free.
- The default stripe count is 1.

## Decision

- **Runs and datasets live on Eagle:** `/lus/eagle/projects/datascience/keceli/chemuniverse/<run>`. That's the default in `iri_job.py`, whose JobSpecs request `filesystems=home:grand:eagle`. Small summaries stay on `/home` (readable through IRI).
- **Format:** gzip Kekulé SMILES, one molecule per line, about 2.9 bytes per molecule (measured for CNOS and CNOF to 9 atoms).
  - After a run, the shard parts of each size are merged by concatenating the gzip files into chunks of about 2 GB, which are then striped.
  - A `manifest.json` records the elements, rules, per-size counts, code commit, job IDs and a SHA-256 per file.
  - Manifests (small) are committed to the repo as the provenance record.
- **Release:** a Globus guest collection on Eagle for everyone. A Hugging Face mirror only for sets up to about 100 GB (≤ 12 atoms).
- **Depth:** store full datasets up to **13 atoms** (about 1.5 TB) and **14 atoms** (about 21 TB) if the calibration confirms the costs. At 15 atoms and beyond (≥ 300 TB) keep counts only. Export is deterministic, so any shard can be regenerated on request.

## First dataset (2026-10-17)

`chemuniverse-CNOS-le12` holds 42,655,220,742 molecules: 132.4 GB in 74 chunks, **3.10 bytes per molecule**.
- **Export:** job 362511, 4 nodes, 40,000 shards. It took 1,441,389 core-seconds (400 core-hours; estimate 508) in 49 minutes, with 95% efficiency and a slowest shard of 531 s.
- **Serial finalize:** opening 160,008 parts one by one on Lustre ran at about 2 MB/s, so the job was cancelled. Replaced by a parallel finalize: chunks are planned from the shard counts, then 64 writers run at once. Job 362533 took 21 minutes on 1 node, about 100 MB/s.
- **Verification:** `shard verify` (job 362542) re-read every chunk: 74 files, 0 bad, line counts summing to exactly 42,655,220,742 molecules, in 72 s on one node.
- **Scaling concern:** at 100 MB/s per node, a 14-atom dataset (about 21 TB) would need about 2.5 node-days just to merge. Before such a run, the merge must be distributed across nodes (as ensemble tasks), or shards must write larger parts.

## Consequences

- No transfer is needed after a run, since Crux writes Eagle directly.
- Grand (IQC) is not used for datasets; on Crux it is a small clone anyway.
- Many-small-file pressure on Lustre is removed by the merge step, which still needs to be implemented.
