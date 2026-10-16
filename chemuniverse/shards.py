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
VERSION = 1


def _write_json(path, data):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    tmp.write_text(json.dumps(data))
    os.replace(tmp, path)


def result_path(run_dir, i):
    return Path(run_dir) / "results" / f"{i // 1000:04d}" / f"{i:07d}.json"


def plan(run_dir, kind, max_atoms, elements, cfg=None, shards=1000, k_exact=6, passes=10, seed=0, jobs=1,
         log=print):
    """Write run_dir/plan.json and run_dir/frontier.json.gz; returns the plan."""
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
        if kind != "passes" and (n == max_atoms or len(level) >= 4 * shards):
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
        assignment = [order[k::n_shards] for k in range(n_shards)]
        with gzip.open(run_dir / "assignment.json.gz", "wt") as f:
            json.dump(assignment, f)
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
    """Plan, rule config, frontier and assignment; cached per worker process."""
    run_dir = Path(run_dir)
    p = json.loads((run_dir / "plan.json").read_text())
    with gzip.open(run_dir / "frontier.json.gz", "rt") as f:
        frontier = [Skeleton(tuple(lab), tuple(adj)) for lab, adj in json.load(f)]
    assignment = None
    if p["kind"] != "passes":
        with gzip.open(run_dir / "assignment.json.gz", "rt") as f:
            assignment = json.load(f)
    return p, RuleConfig(**p["rules"]), frontier, assignment


def pending(run_dir):
    p = json.loads((Path(run_dir) / "plan.json").read_text())
    return [i for i in range(p["n_shards"]) if not result_path(run_dir, i).exists()]


def run_shard(run_dir, i):
    """Run shard i unless its result exists; returns a short summary (picklable, small)."""
    run_dir = str(run_dir)
    out = result_path(run_dir, i)
    if out.exists():
        return {"shard": i, "skipped": True}
    p, cfg, frontier, assignment = _load(run_dir)
    t0 = time.time()
    elements, split, max_atoms = p["elements"], p["split_level"], p["max_atoms"]
    if p["kind"] == "passes":
        rows = _stratified_pass(p["seeds"][i], frontier, split, max_atoms, elements, cfg, cached_count)
    else:
        nodes = [frontier[k] for k in assignment[i]]
        if p["kind"] == "exact":
            if core:
                rows, _ = core.expand(nodes, max_atoms - split, elements, cfg_tuple(cfg), False, 1)
            else:
                rows, _ = _expand(nodes, max_atoms - split, elements, cfg, count_molecules)
        else:
            from .export import _export_subtrees

            smiles_dir = Path(run_dir) / "smiles"
            smiles_dir.mkdir(exist_ok=True)
            counts = _export_subtrees(i + 1, nodes, split, max_atoms, elements, cfg, str(smiles_dir))
            rows = [(c, None) for c in counts]
    data = {"shard": i, "rows": [list(r) for r in rows], "seconds": round(time.time() - t0, 3),
            "host": socket.gethostname()}
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
        manifest = {"elements": p["elements"], "max_atoms": p["max_atoms"], "rules": p["rules"],
                    "counts": {str(r.n): r.count for r in results}, "total": sum(r.count for r in results),
                    "format": "smiles/n{size}_part{k}.smi.gz, one Kekule SMILES per line; each molecule once"}
        _write_json(run_dir / "manifest.json", manifest)
        _write_levels(run_dir, p, split)
    return results, stats


def _write_levels(run_dir, p, split):
    """Export: the molecules of levels <= split come from the plan's serial part."""
    from .export import molecule_smiles

    cfg = RuleConfig(**p["rules"])
    smiles_dir = Path(run_dir) / "smiles"
    smiles_dir.mkdir(exist_ok=True)
    level = roots(p["elements"])
    for n in range(1, split + 1):
        if n > 1:
            level = [c for s in level for c in children(s, p["elements"], cfg)]
        with gzip.open(smiles_dir / f"n{n:02d}_part0000.smi.gz", "wt") as f:
            f.write("".join(smi + "\n" for s in level for smi in molecule_smiles(s, cfg)))
