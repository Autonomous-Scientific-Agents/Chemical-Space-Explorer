import argparse
import time

from .rules import RuleConfig
from .tree import estimate


def main(argv=None):
    p = argparse.ArgumentParser(prog="chemuniverse", description="Estimate the number of stable molecules.")
    p.add_argument("--max-atoms", type=int, required=True, help="maximum number of heavy atoms")
    p.add_argument("--elements", nargs="+", default=["C", "N", "O", "F"])
    p.add_argument("--exact-up-to", type=int, default=6, help="enumerate exactly up to this size")
    p.add_argument("--method", choices=["stratified", "knuth"], default="stratified",
                   help="sampler for sizes above --exact-up-to")
    p.add_argument("--samples", type=int, default=None,
                   help="independent samples: stratified passes (default 10) or Knuth probes (default 1000)")
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--hypervalent", action="store_true", help="allow S(IV/VI) and P(V)")
    p.add_argument("--allow-n-o", action="store_true", help="allow N-O bonds")
    p.add_argument("--no-bond-order-limits", action="store_true")
    p.add_argument("--no-strain", action="store_true")
    p.add_argument("--no-heteroatom", action="store_true")
    p.add_argument("--no-rings", action="store_true")
    p.add_argument("--no-multiple-bonds", action="store_true")
    a = p.parse_args(argv)
    cfg = RuleConfig(
        hypervalent=a.hypervalent,
        bond_order_limits=not a.no_bond_order_limits,
        strain=not a.no_strain,
        heteroatom=not a.no_heteroatom,
        ban_n_o=not a.allow_n_o,
        allow_rings=not a.no_rings,
        allow_multiple_bonds=not a.no_multiple_bonds,
    )
    t0 = time.time()
    results = estimate(a.max_atoms, a.elements, cfg, a.exact_up_to, a.samples, a.seed, a.method)
    print(f"elements: {' '.join(a.elements)}   max heavy atoms: {a.max_atoms}")
    print(f"{'n':>3} {'molecules':>12} {'± SE':>10} {'mode':>8} {'cumulative':>12}")
    total = total_var = 0.0
    for r in results:
        total += r.count
        total_var += r.se ** 2
        mode = "exact" if r.exact else "estimate"
        se = "" if r.exact else f"{r.se:.2g}"
        print(f"{r.n:>3} {r.count:>12.4g} {se:>10} {mode:>8} {total:>12.4g}")
    print(f"total: {total:.4g} ± {total_var ** 0.5:.2g}   ({time.time() - t0:.1f} s)")


if __name__ == "__main__":
    main()
