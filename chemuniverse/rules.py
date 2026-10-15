"""Stability rules.

Hereditary rules (degree caps, banned element pairs, no rings) hold for every
connected subgraph, so they prune the generation tree. The remaining rules are
checked on full molecules with bond orders and implicit hydrogens.
"""
from dataclasses import dataclass

from . import elements
from .graph import bits, simple_cycles

MAX_RING = 8  # largest ring any strain rule looks at


@dataclass(frozen=True)
class RuleConfig:
    hypervalent: bool = False
    bond_order_limits: bool = True
    strain: bool = True
    heteroatom: bool = True
    ban_n_o: bool = True
    allow_rings: bool = True
    allow_multiple_bonds: bool = True

    def valences(self, symbol):
        if self.hypervalent and symbol in elements.HYPERVALENT:
            return elements.HYPERVALENT[symbol]
        return elements.get(symbol).valences

    def max_valence(self, symbol):
        return self.valences(symbol)[-1]

    def max_bond_order(self, symbol):
        if not self.allow_multiple_bonds:
            return 1
        return elements.get(symbol).max_bond_order if self.bond_order_limits else 3

    def pair_allowed(self, a, b):
        if not self.heteroatom:
            return True
        pair = {a, b}
        if pair == {"O"}:
            return False
        if self.ban_n_o and pair == {"N", "O"}:
            return False
        ea, eb = elements.get(a), elements.get(b)
        for x, y in ((ea, eb), (eb, ea)):
            if x.halogen and (y.halogen or y.symbol in ("N", "O", "S")):
                return False
        return True


def skeleton_ok(cfg, skel):
    """Hereditary rules; the generation tree enforces these while building skeletons."""
    edges = skel.edges()
    if not cfg.allow_rings and len(edges) >= skel.n:
        return False
    if any(bin(a).count("1") > cfg.max_valence(lab) for lab, a in zip(skel.labels, skel.adj)):
        return False
    return all(cfg.pair_allowed(skel.labels[u], skel.labels[v]) for u, v in edges)


class RingInfo:
    """Ring facts about a skeleton needed by the strain rules."""

    def __init__(self, skel):
        n = skel.n
        self.edge_ring = {}  # (u, v) -> smallest ring size (<= MAX_RING) containing the edge
        self.atom_ring = [None] * n  # smallest ring size containing the atom
        self.bridgeheads = set()
        if len(skel.edges()) < n:  # tree
            return
        cycles = simple_cycles(skel.adj, MAX_RING)
        for c in cycles:
            k = len(c)
            for i, u in enumerate(c):
                v = c[(i + 1) % k]
                e = (min(u, v), max(u, v))
                if self.edge_ring.get(e, 99) > k:
                    self.edge_ring[e] = k
                if (self.atom_ring[u] or 99) > k:
                    self.atom_ring[u] = k
        # Bredt: two rings <= 7 sharing a path of >= 3 atoms form a bridged
        # system; the ends of the shared path are bridgeheads.
        small = [c for c in cycles if len(c) <= 7]
        for i, a in enumerate(small):
            sa = set(a)
            for b in small[i + 1:]:
                shared = sa.intersection(b)
                if len(shared) < 3:
                    continue
                for ring in (a, b):
                    k = len(ring)
                    for j, x in enumerate(ring):
                        if x in shared and (ring[j - 1] not in shared or ring[(j + 1) % k] not in shared):
                            self.bridgeheads.add(x)


def edge_cap(cfg, skel, rings, u, v):
    """Largest bond order allowed on edge (u, v) by element and strain rules."""
    cap = min(cfg.max_bond_order(skel.labels[u]), cfg.max_bond_order(skel.labels[v]), 3)
    if cfg.strain:
        size = rings.edge_ring.get((min(u, v), max(u, v)))
        if size is not None:
            cap = min(cap, 1 if size <= 4 else 2 if size < 8 else 3)
        if u in rings.bridgeheads or v in rings.bridgeheads:
            cap = 1
    return cap


def molecule_ok(cfg, skel, rings, orders, hcounts):
    """Rules that need the full bond-order assignment.

    orders: dict (u, v) with u < v -> bond order.
    """
    n = skel.n
    labels = skel.labels
    if cfg.strain:
        multi = [0] * n
        for (u, v), o in orders.items():
            if o > 1:
                multi[u] += 1
                multi[v] += 1
        for x in range(n):
            if multi[x] >= 2 and rings.atom_ring[x] is not None and rings.atom_ring[x] < 9:
                return False  # cumulated bonds in a small ring
    if cfg.heteroatom:
        def order(u, v):
            return orders[(u, v) if u < v else (v, u)]

        for x in range(n):
            if labels[x] != "C":
                continue
            nbrs = list(bits(skel.adj[x]))
            # enol: OH on a carbon that is double bonded to carbon
            if any(labels[y] == "C" and order(x, y) == 2 for y in nbrs):
                if any(labels[y] == "O" and hcounts[y] > 0 for y in nbrs):
                    return False
            # hemiacetal / hemiaminal / gem-diol on an sp3 carbon
            if all(order(x, y) == 1 for y in nbrs):
                het = [y for y in nbrs if labels[y] in ("N", "O")]
                if len(het) >= 2 and any(hcounts[y] > 0 for y in het):
                    return False
    return True
