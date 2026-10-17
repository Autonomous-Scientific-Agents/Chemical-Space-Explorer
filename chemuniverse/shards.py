"""Split a run into independent shards for multi-node execution, run them, and reduce the results.

Launcher-agnostic: any task farm (ensemble_launcher, job arrays, mpi4py, --jobs)
runs `run_shard(run_dir, i)` for every pending shard, then `reduce` combines them.

    plan    exact / export: count the tree serially up to a split level, then
            assign the split-level nodes round-robin (sorted by a cost proxy) to
            shards. passes: exact counts up to k_exact plus one seed per pass.
    run     each shard writes results/<i // 1000>/<i>.json atomically when done,
            so finished shards are skipped on restart (preemptable queues).
    reduce  sums shard rows (exact, export) or combines passes (mean ± SE);
            the result is identical to the single-machine functions.
"""
import gzip
import json
import os
import random
import socket
import time
from dataclasses import asdict
from functools import lru_cache
from pathlib import Path

from .backend import cfg_tuple, core
from .bonds import count_molecules
from .graph import Skeleton
from .rules import RuleConfig
from .tree import LevelResult, _expand, _stratified_pass, _summarize, cached_count, children, roots

KINDS = ("exact", "export", "passes")
VERSION = 2
BLOCK = 1000  # shards per input block / results directory / SMILES part directory
BYTES_PER_MOLECULE = 3  # gzip Kekule SMILES, measured 2.8-2.9 for CNOF/CNOS up to 9 atoms


def _write_json(path, data):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    tmp.write_text(json.dumps(data))
    os.replace(tmp, path)


def result_path(run_dir, i):
    return Path(run_dir) / "results" / f"{i // BLOCK:04d}" / f"{i:07d}.json"


def _block_path(run_dir, b):
    return Path(run_dir) / "inputs" / f"{b:05d}.json.gz"


def plan(run_dir, kind, max_atoms, elements, cfg=None, shards=1000, k_exact=6, passes=10, seed=0, jobs=1,
         log=print, split_level=None, nodes_per_shard=16):
    """Write run_dir/plan.json plus the shard inputs; returns the plan.

    passes: frontier.json.gz (every pass starts from the whole k_exact level).
    exact/export: inputs/<block>.json.gz with the nodes of BLOCK shards each, so a worker
    loads only its block (a full frontier of ~10^6 nodes per worker would exhaust memory).

    exact/export split at `split_level`, or else at the first level with at least
    nodes_per_shard * shards nodes. A deeper split makes the largest single subtree, which
    bounds the run time, smaller: on Crux CNOF<=12 split at level 7 with ~4 nodes per shard
    had one 735 s shard in a 761 s run (job 362424).
    """
    if kind not in KINDS:
        raise ValueError(f"kind must be one of {KINDS}")
    cfg = cfg or RuleConfig()
    elements = list(dict.fromkeys(elements))
    run_dir = Path(run_dir)
    if (run_dir / "plan.json").exists():
        raise FileExistsError(f"{run_dir}/plan.json exists; use a new directory")
    run_dir.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    level = roots(elements)
    rows = []  # exact (molecules, skeletons) for levels 1..split
    split = min(k_exact, max_atoms) if kind == "passes" else None
    for n in range(1, max_atoms + 1):
        if n > 1:
            level = [c for s in level for c in children(s, elements, cfg)]
        rows.append((sum(count_molecules(s, cfg) for s in level), len(level)))
        if kind != "passes" and (n == max_atoms or n == split_level
                                 or (split_level is None and len(level) >= nodes_per_shard * shards)):
            split = n
        if n == split:
            break
    if kind == "passes":
        rng = random.Random(seed)
        seeds = [rng.getrandbits(64) for _ in range(passes)]  # same draws as tree.stratified_estimate
        frontier, n_shards = level, passes
    else:
        seeds = None
        # cost proxy: number of children; deal nodes round-robin from the most expensive down
        cost = [len(children(s, elements, cfg)) for s in level] if split < max_atoms else [0] * len(level)
        order = sorted(range(len(level)), key=lambda i: -cost[i])
        n_shards = max(1, min(shards, len(level))) if split < max_atoms else 0
        frontier = level
        for b in range(0, n_shards, BLOCK):
            block = [[[list(level[k].labels), list(level[k].adj)] for k in order[i::n_shards]]
                     for i in range(b, min(b + BLOCK, n_shards))]
            _block_path(run_dir, b // BLOCK).parent.mkdir(exist_ok=True)
            with gzip.open(_block_path(run_dir, b // BLOCK), "wt") as f:
                json.dump(block, f)
    if kind == "passes":
        with gzip.open(run_dir / "frontier.json.gz", "wt") as f:
            json.dump([[list(s.labels), list(s.adj)] for s in frontier], f)
    p = {
        "version": VERSION,
        "kind": kind,
        "elements": elements,
        "rules": asdict(cfg),
        "max_atoms": max_atoms,
        "split_level": split,
        "n_shards": n_shards,
        "exact_rows": rows,
        "seeds": seeds,
        "seed": seed,
        "frontier_size": len(frontier),
        "created": time.strftime("%Y-%m-%d %H:%M:%S"),
        "plan_seconds": round(time.time() - t0, 1),
    }
    _write_json(run_dir / "plan.json", p)
    log(f"plan: {kind}, split at level {split} ({len(frontier)} nodes), {n_shards} shards "
        f"({p['plan_seconds']} s)")
    return p


@lru_cache(maxsize=4)
def _load(run_dir):
    """Plan, rule config and (passes only) the frontier; cached per worker process."""
    run_dir = Path(run_dir)
    p = json.loads((run_dir / "plan.json").read_text())
    frontier = None
    if p["kind"] == "passes":
        with gzip.open(run_dir / "frontier.json.gz", "rt") as f:
            frontier = [Skeleton(tuple(lab), tuple(adj)) for lab, adj in json.load(f)]
    return p, RuleConfig(**p["rules"]), frontier


@lru_cache(maxsize=2)
def _block(run_dir, b):
    with gzip.open(_block_path(run_dir, b), "rt") as f:
        return [[Skeleton(tuple(lab), tuple(adj)) for lab, adj in nodes] for nodes in json.load(f)]


def pending(run_dir):
    p = json.loads((Path(run_dir) / "plan.json").read_text())
    return [i for i in range(p["n_shards"]) if not result_path(run_dir, i).exists()]


def run_shard(run_dir, i):
    """Run shard i unless its result exists; returns a short summary (picklable, small)."""
    run_dir = str(run_dir)
    out = result_path(run_dir, i)
    if out.exists():
        return {"shard": i, "skipped": True}
    p, cfg, frontier = _load(run_dir)
    t0 = time.time()
    elements, split, max_atoms = p["elements"], p["split_level"], p["max_atoms"]
    if p["kind"] == "passes":
        rows = _stratified_pass(p["seeds"][i], frontier, split, max_atoms, elements, cfg, cached_count)
    else:
        nodes = _block(run_dir, i // BLOCK)[i % BLOCK]
        if p["kind"] == "exact":
            if core:
                rows, _ = core.expand(nodes, max_atoms - split, elements, cfg_tuple(cfg), False, 1)
            else:
                rows, _ = _expand(nodes, max_atoms - split, elements, cfg, count_molecules)
        else:
            from .export import _export_subtrees

            smiles_dir = Path(run_dir) / "smiles" / f"{i // BLOCK:04d}"
            smiles_dir.mkdir(parents=True, exist_ok=True)
            counts = _export_subtrees(i + 1, nodes, split, max_atoms, elements, cfg, str(smiles_dir))
            rows = [(c, None) for c in counts]
    data = {"shard": i, "rows": [list(r) for r in rows], "seconds": round(time.time() - t0, 3),
            "host": socket.gethostname(), "job": os.environ.get("PBS_JOBID")}
    _write_json(out, data)
    return {"shard": i, "seconds": data["seconds"]}


def reduce(run_dir):
    """Combine all shard results; returns (list[LevelResult], stats). Raises if shards are missing."""
    run_dir = Path(run_dir)
    p = json.loads((run_dir / "plan.json").read_text())
    missing = pending(run_dir)
    if missing:
        raise RuntimeError(f"{len(missing)} of {p['n_shards']} shards missing, e.g. {missing[:5]}")
    shards = [json.loads(result_path(run_dir, i).read_text()) for i in range(p["n_shards"])]
    split = p["split_level"]
    results = [LevelResult(n + 1, mols, skeletons=sk) for n, (mols, sk) in enumerate(p["exact_rows"])]
    depth = p["max_atoms"] - split
    if p["kind"] == "passes":
        if depth:
            results += _summarize([[tuple(r) for r in s["rows"]] for s in shards], split)
    else:
        for d in range(depth):
            sk = None if p["kind"] == "export" else sum(s["rows"][d][1] for s in shards)
            results.append(LevelResult(split + d + 1, sum(s["rows"][d][0] for s in shards), skeletons=sk))
    secs = [s["seconds"] for s in shards]
    stats = {"shards": len(shards), "core_seconds": round(sum(secs), 1), "max_shard_seconds": max(secs, default=0),
             "hosts": len({s["host"] for s in shards})}
    out = {"plan": {k: v for k, v in p.items() if k not in ("seeds", "exact_rows")}, "stats": stats,
           "rows": [[r.n, r.count, r.se, r.exact, r.skeletons] for r in results]}
    _write_json(run_dir / "result.json", out)
    if p["kind"] == "export":
        _write_levels(run_dir, p, split)
    return results, stats


def _write_levels(run_dir, p, split):
    """Export: the molecules of levels <= split come from the plan's serial part."""
    from .export import molecule_smiles

    cfg = RuleConfig(**p["rules"])
    smiles_dir = Path(run_dir) / "smiles" / "serial"
    smiles_dir.mkdir(parents=True, exist_ok=True)
    level = roots(p["elements"])
    for n in range(1, split + 1):
        if n > 1:
            level = [c for s in level for c in children(s, p["elements"], cfg)]
        with gzip.open(smiles_dir / f"n{n:02d}_part0000.smi.gz", "wt") as f:
            f.write("".join(smi + "\n" for s in level for smi in molecule_smiles(s, cfg)))


def _code_version():
    import subprocess

    try:
        return subprocess.run(["git", "-C", str(Path(__file__).resolve().parent), "rev-parse", "HEAD"],
                              capture_output=True, text=True, timeout=10).stdout.strip() or None
    except Exception:
        return None


def _write_chunk(out_path, parts):
    """Concatenate gzip parts into one chunk; returns (bytes, sha256). Runs in a worker process."""
    import hashlib

    digest, written = hashlib.sha256(), 0
    with open(out_path, "wb") as chunk:
        for f in parts:
            with open(f, "rb") as src:
                while True:
                    buf = src.read(1 << 24)
                    if not buf:
                        break
                    chunk.write(buf)
                    digest.update(buf)
                    written += len(buf)
    return written, digest.hexdigest()


def finalize(run_dir, chunk_bytes=2 * 10**9, delete_parts=False, jobs=1, log=print):
    """Export runs: merge each size's SMILES parts into ~chunk_bytes gzip chunks (gzip members
    concatenate into a valid gzip file), checksum them and write dataset/manifest.json.

    Chunks are planned from the per-shard molecule counts (about BYTES_PER_MOLECULE each), so
    planning needs no stat calls, and built in parallel: on Lustre, opening ~10^4-10^5 small
    part files one by one is metadata-bound (job 362511 merged only ~2 MB/s serially). The plan
    depends only on the counts, so the files are identical for any number of jobs.
    """
    from .tree import _pmap

    run_dir = Path(run_dir)
    p = json.loads((run_dir / "plan.json").read_text())
    if p["kind"] != "export":
        raise ValueError("finalize applies to export runs")
    if not (run_dir / "result.json").exists():
        raise RuntimeError("run reduce first")
    split, n_shards = p["split_level"], p["n_shards"]
    shard_rows = [json.loads(result_path(run_dir, i).read_text()) for i in range(n_shards)]
    out_dir = run_dir / "dataset"
    out_dir.mkdir(exist_ok=True)
    t0 = time.time()
    target = max(1, chunk_bytes // BYTES_PER_MOLECULE)
    chunks, counts = [], {}  # (name, size, molecules, [part paths])
    for n in range(1, p["max_atoms"] + 1):
        if n <= split:
            parts = [(run_dir / "smiles" / "serial" / f"n{n:02d}_part0000.smi.gz", p["exact_rows"][n - 1][0])]
        else:
            parts = [(run_dir / "smiles" / f"{i // BLOCK:04d}" / f"n{n:02d}_part{i + 1:04d}.smi.gz",
                      shard_rows[i]["rows"][n - split - 1][0]) for i in range(n_shards)]
        parts = [(f, c) for f, c in parts if c]
        counts[n] = sum(c for _, c in parts)
        group, mols = [], 0
        for f, c in parts:
            group.append(str(f))
            mols += c
            if mols >= target:
                chunks.append((f"n{n:02d}_{sum(ch[1] == n for ch in chunks):05d}.smi.gz", n, mols, group))
                group, mols = [], 0
        if group:
            chunks.append((f"n{n:02d}_{sum(ch[1] == n for ch in chunks):05d}.smi.gz", n, mols, group))
    log(f"finalize: {len(chunks)} chunks from {sum(len(ch[3]) for ch in chunks)} parts, {jobs} jobs")
    written = _pmap(_write_chunk, [(str(out_dir / name), group) for name, _, _, group in chunks], jobs)
    files = [{"file": f"dataset/{name}", "size": n, "molecules": mols, "bytes": nbytes, "sha256": sha}
             for (name, n, mols, _), (nbytes, sha) in zip(chunks, written)]
    for n in counts:
        log(f"n={n}: {counts[n]} molecules in {sum(f['size'] == n for f in files)} chunk(s)")
    result = json.loads((run_dir / "result.json").read_text())
    manifest = {
        "name": f"chemuniverse-{''.join(p['elements'])}-le{p['max_atoms']}",
        "description": "All molecules up to max_atoms heavy atoms that pass the chemuniverse stability rules, "
                       "each exactly once, as Kekule SMILES (implicit hydrogens, constitutional isomers only).",
        "elements": p["elements"],
        "max_atoms": p["max_atoms"],
        "rules": p["rules"],
        "counts": {str(n): c for n, c in counts.items()},
        "total": sum(counts.values()),
        "files": files,
        "code": {"package_version": VERSION, "git_commit": _code_version()},
        "jobs": sorted({s["job"] for s in shard_rows if s.get("job")}),
        "stats": result["stats"],
        "created": time.strftime("%Y-%m-%d"),
        "format": "gzip text, one SMILES per line; files are concatenated gzip members",
    }
    _write_json(out_dir / "manifest.json", manifest)
    if delete_parts:
        import shutil

        shutil.rmtree(run_dir / "smiles")
    log(f"dataset: {manifest['total']} molecules in {len(files)} files ({time.time() - t0:.0f} s)")
    return manifest
