import argparse
import json
import sys
import time
from pathlib import Path

from .rules import RuleConfig
from .tree import estimate


def add_rule_args(p):
    p.add_argument("--elements", nargs="+", default=["C", "N", "O", "F"])
    p.add_argument("--jobs", type=int, default=1, help="worker processes (results do not depend on it)")
    p.add_argument("--hypervalent", action="store_true", help="allow S(IV/VI) and P(V)")
    p.add_argument("--allow-n-o", action="store_true", help="allow N-O bonds")
    p.add_argument("--no-bond-order-limits", action="store_true")
    p.add_argument("--no-strain", action="store_true")
    p.add_argument("--no-heteroatom", action="store_true")
    p.add_argument("--no-rings", action="store_true")
    p.add_argument("--no-multiple-bonds", action="store_true")


def rule_config(a):
    return RuleConfig(
        hypervalent=a.hypervalent,
        bond_order_limits=not a.no_bond_order_limits,
        strain=not a.no_strain,
        heteroatom=not a.no_heteroatom,
        ban_n_o=not a.allow_n_o,
        allow_rings=not a.no_rings,
        allow_multiple_bonds=not a.no_multiple_bonds,
    )


def export_main(argv):
    from .export import export

    p = argparse.ArgumentParser(prog="chemuniverse export",
                                description="Write every molecule up to --max-atoms as gzipped SMILES files.")
    p.add_argument("--max-atoms", type=int, required=True, help="maximum number of heavy atoms")
    p.add_argument("--out", required=True, help="output directory (must not contain .smi.gz files)")
    add_rule_args(p)
    a = p.parse_args(argv)
    counts = export(a.max_atoms, a.elements, rule_config(a), a.out, a.jobs)
    print(f"total: {sum(counts)} molecules written to {a.out}")


def print_results(results, elements, max_atoms, seconds):
    print(f"elements: {' '.join(elements)}   max heavy atoms: {max_atoms}")
    print(f"{'n':>3} {'molecules':>12} {'± SE':>10} {'mode':>8} {'cumulative':>12}")
    total = total_var = 0.0
    for r in results:
        total += r.count
        total_var += r.se ** 2
        mode = "exact" if r.exact else "estimate"
        se = "" if r.exact else f"{r.se:.2g}"
        print(f"{r.n:>3} {r.count:>12.4g} {se:>10} {mode:>8} {total:>12.4g}")
    print(f"total: {total:.4g} ± {total_var ** 0.5:.2g}   ({seconds:.1f} s)")


def shard_main(argv):
    from . import shards

    p = argparse.ArgumentParser(prog="chemuniverse shard",
                                description="Split a run into shards for a task farm (see chemuniverse/shards.py).")
    sub = p.add_subparsers(dest="cmd", required=True)
    pl = sub.add_parser("plan", help="create RUN_DIR with the shard plan")
    pl.add_argument("run_dir")
    pl.add_argument("--kind", choices=shards.KINDS, required=True)
    pl.add_argument("--max-atoms", type=int, required=True)
    pl.add_argument("--shards", type=int, default=1000, help="exact/export: number of subtree shards")
    pl.add_argument("--exact-up-to", type=int, default=6, help="passes: exact frontier level")
    pl.add_argument("--passes", type=int, default=10, help="passes: number of stratified passes (= shards)")
    pl.add_argument("--seed", type=int, default=0)
    add_rule_args(pl)
    ru = sub.add_parser("run", help="run pending shards locally")
    ru.add_argument("run_dir")
    ru.add_argument("--index", type=int, nargs="*", help="shard indices (default: all pending)")
    ru.add_argument("--jobs", type=int, default=1)
    rd = sub.add_parser("reduce", help="combine shard results into RUN_DIR/result.json")
    rd.add_argument("run_dir")
    a = p.parse_args(argv)
    t0 = time.time()
    if a.cmd == "plan":
        shards.plan(a.run_dir, a.kind, a.max_atoms, a.elements, rule_config(a), a.shards, a.exact_up_to, a.passes,
                    a.seed)
    elif a.cmd == "run":
        from .tree import _pmap

        todo = a.index if a.index else shards.pending(a.run_dir)
        _pmap(shards.run_shard, [(a.run_dir, i) for i in todo], a.jobs)
        print(f"ran {len(todo)} shards ({time.time() - t0:.1f} s); {len(shards.pending(a.run_dir))} pending")
    else:
        results, stats = shards.reduce(a.run_dir)
        plan = json.loads((Path(a.run_dir) / "plan.json").read_text())
        print_results(results, plan["elements"], plan["max_atoms"], time.time() - t0)
        print(f"{stats['shards']} shards on {stats['hosts']} hosts, {stats['core_seconds']} core-seconds, "
              f"slowest shard {stats['max_shard_seconds']} s")


def main(argv=None):
    argv = sys.argv[1:] if argv is None else argv
    if argv[:1] == ["export"]:
        return export_main(argv[1:])
    if argv[:1] == ["shard"]:
        return shard_main(argv[1:])
    p = argparse.ArgumentParser(prog="chemuniverse", description="Estimate the number of stable molecules.",
                                epilog="Subcommands: 'chemuniverse export --help' writes the molecules as SMILES; "
                                       "'chemuniverse shard --help' splits a run for a multi-node task farm.")
    p.add_argument("--max-atoms", type=int, required=True, help="maximum number of heavy atoms")
    p.add_argument("--exact-up-to", type=int, default=6, help="enumerate exactly up to this size")
    p.add_argument("--method", choices=["stratified", "knuth"], default="stratified",
                   help="sampler for sizes above --exact-up-to")
    p.add_argument("--samples", type=int, default=None,
                   help="independent samples: stratified passes (default 10) or Knuth probes (default 1000)")
    p.add_argument("--seed", type=int, default=0)
    add_rule_args(p)
    a = p.parse_args(argv)
    cfg = rule_config(a)
    t0 = time.time()
    results = estimate(a.max_atoms, a.elements, cfg, a.exact_up_to, a.samples, a.seed, a.method, a.jobs)
    print_results(results, a.elements, a.max_atoms, time.time() - t0)


if __name__ == "__main__":
    main()
