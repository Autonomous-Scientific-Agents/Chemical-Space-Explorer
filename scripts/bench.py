"""Benchmark exact enumeration, stratified sampling and export (see design/prd/0002-performance.md).

    python scripts/bench.py --jobs 1 10
    CHEMUNIVERSE_BACKEND=python python scripts/bench.py --jobs 1 10
"""
import argparse
import random
import sys
import tempfile
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from chemuniverse import RuleConfig  # noqa: E402
from chemuniverse.backend import core  # noqa: E402
from chemuniverse.export import export  # noqa: E402
from chemuniverse.tree import cached_children, cached_count, enumerate_exact, stratified_estimate  # noqa: E402

EL = ["C", "N", "O", "F"]


def timed(fn):
    # cold caches, so forked workers do not inherit results from an earlier run
    cached_children.cache_clear()
    cached_count.cache_clear()
    t0 = time.perf_counter()
    out = fn()
    return out, time.perf_counter() - t0


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--jobs", type=int, nargs="+", default=[1])
    p.add_argument("--exact", type=int, default=7, help="exact CNOF enumeration up to this size")
    p.add_argument("--sample-to", type=int, default=12, help="stratified passes from the 6-atom frontier to this size")
    p.add_argument("--passes", type=int, default=10)
    p.add_argument("--export", type=int, default=7, help="export CNOF molecules up to this size")
    a = p.parse_args()
    cfg = RuleConfig()
    _, frontier = enumerate_exact(6, EL, cfg)
    print(f"backend: {'rust' if core else 'python'}   (set CHEMUNIVERSE_BACKEND=python|rust)")
    print(f"{'jobs':>4} {'task':<34} {'seconds':>8} {'molecules/s':>12}")
    for jobs in a.jobs:
        res, t = timed(lambda: enumerate_exact(a.exact, EL, cfg, jobs=jobs)[0])
        total = sum(r.count for r in res)
        print(f"{jobs:>4} {f'exact CNOF <= {a.exact}':<34} {t:>8.1f} {total / t:>12.0f}")
        _, t = timed(lambda: stratified_estimate(frontier, 6, a.sample_to, EL, cfg, a.passes, random.Random(0),
                                                 jobs=jobs))
        print(f"{jobs:>4} {f'{a.passes} stratified passes 6 -> {a.sample_to}':<34} {t:>8.1f} "
              f"{'(' + format(t * jobs / a.passes, '.1f') + ' s/pass/core)':>12}")
        with tempfile.TemporaryDirectory() as d:
            counts, t = timed(lambda: export(a.export, EL, cfg, d, jobs=jobs, log=lambda *_: None))
        print(f"{jobs:>4} {f'export CNOF <= {a.export}':<34} {t:>8.1f} {sum(counts) / t:>12.0f}")


if __name__ == "__main__":
    main()
