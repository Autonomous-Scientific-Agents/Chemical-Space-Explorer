"""Stability rules.

Hereditary rules (degree caps, banned element pairs, no rings) hold for every
connected subgraph, so they prune the generation tree. The remaining rules are
checked on full molecules with bond orders and implicit hydrogens.
"""
from dataclasses import dataclass

from . import elements
from .graph import bits, simple_cycles

MAX_RING = 8  # largest ring any strain rule looks at
HETERO = frozenset("NOS")  # heteroatoms for the strict-tier rules (ADR 0015)


@dataclass(frozen=True)
class RuleConfig:
    hypervalent: bool = False
    bond_order_limits: bool = True
    strain: bool = True
    heteroatom: bool = True
    ban_n_o: bool = True
    allow_rings: bool = True
    allow_multiple_bonds: bool = True
    # strict tier (ADR 0015); off by default so the basic tier is unchanged
    hetero_chains: bool = False  # S1: no N/O/S atom with >= 2 N/O/S neighbours (X-X-X chains)
    ban_s_o: bool = False  # S2: no S-O bond
    small_ring_heteroatoms: bool = False  # S3: no 3- or 4-membered ring with >= 2 N/O/S
    tautomers: bool = False  # S4: no enethiol, NH-enamine, ynol/ynamine/ynethiol, imidic/thiono acid, ketene
    thio_gem: bool = False  # S5: no hemithioacetal / gem-thiol on an sp3 carbon
    thiocarbonyl: bool = False  # S6: no C=S on a carbon without a heteroatom substituent

    @classmethod
    def strict(cls, **overrides):
        """All strict-tier rules on (ADR 0015)."""
        on = dict(hetero_chains=True, ban_s_o=True, small_ring_heteroatoms=True, tautomers=True, thio_gem=True,
                  thiocarbonyl=True)
        on.update(overrides)
        return cls(**on)

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
        pair = {a, b}
        if self.ban_s_o and pair == {"S", "O"}:
            return False
        if not self.heteroatom:
            return True
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
    if not all(cfg.pair_allowed(skel.labels[u], skel.labels[v]) for u, v in edges):
        return False
    return all(hetero_ok(cfg, skel.labels, skel.adj, v) for v in range(skel.n))


def hetero_ok(cfg, labels, adj, v):
    """Strict hereditary rules S1 and S3 around vertex v. Deleting a vertex can neither create an
    X-X-X chain nor a small ring, so checking only the new vertex keeps the tree exact."""
    het_v = labels[v] in HETERO
    if cfg.hetero_chains:
        het_nbrs = [u for u in bits(adj[v]) if labels[u] in HETERO]
        if het_v and len(het_nbrs) >= 2:
            return False
        if het_v and any(sum(labels[w] in HETERO for w in bits(adj[u])) >= 2 for u in het_nbrs):
            return False
    if cfg.small_ring_heteroatoms:
        nbrs = list(bits(adj[v]))
        for i, a in enumerate(nbrs):
            for b in nbrs[i + 1:]:
                rings = [()] if adj[a] >> b & 1 else []  # triangle v-a-b
                rings += [(c,) for c in bits(adj[a] & adj[b] & ~(1 << v))]  # 4-ring v-a-c-b
                for extra in rings:
                    if sum(labels[x] in HETERO for x in (v, a, b, *extra)) >= 2:
                        return False
    return True


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
    def order(u, v):
        return orders[(u, v) if u < v else (v, u)]

    if cfg.tautomers or cfg.thio_gem or cfg.thiocarbonyl:
        for x in range(n):
            if labels[x] == "C" and not strict_carbon_ok(cfg, labels, list(bits(skel.adj[x])), x, order, hcounts):
                return False
    if cfg.heteroatom:
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


def strict_carbon_ok(cfg, labels, nbrs, x, order, hcounts):
    """Strict molecule rules S4-S6 at carbon x (ADR 0015)."""
    double = [y for y in nbrs if order(x, y) == 2]
    single_het = [y for y in nbrs if labels[y] in HETERO and order(x, y) == 1]
    single_het_h = [y for y in single_het if hcounts[y] > 0]
    if cfg.tautomers:
        c_double = any(labels[y] == "C" for y in double)
        if c_double and any(labels[y] in ("S", "N") for y in single_het_h):
            return False  # enethiol HS-C=C, NH-enamine HN-C=C (enols are basic-tier)
        if single_het and any(labels[y] == "C" and order(x, y) == 3 for y in nbrs):
            return False  # ynol / ynamine / ynethiol X-C#C (cyanamides, cyanates X-C#N stay)
        if any(labels[y] == "N" for y in double) and any(labels[y] in ("O", "S") for y in single_het_h):
            return False  # imidic acid HO-C=N / HS-C=N (the amide / thioamide is kept)
        if any(labels[y] == "S" for y in double) and any(labels[y] == "O" for y in single_het_h):
            return False  # thiono acid HO-C=S (the thiol acid O=C-SH is kept)
        if c_double and any(labels[y] == "O" for y in double):
            return False  # ketene C=C=O
    if cfg.thio_gem and len(double) == 0 and all(order(x, y) == 1 for y in nbrs):
        if len(single_het) >= 2 and single_het_h and any(labels[y] == "S" for y in single_het):
            return False  # hemithioacetal / gem-thiol
    if cfg.thiocarbonyl and any(labels[y] == "S" for y in double) and sum(labels[y] in HETERO for y in nbrs) == 1:
        return False  # thioketone / thioaldehyde / thioketene: S is the carbon's only heteroatom
        # (thioamides, thioesters, thioureas, CS2 and COS stay)
    return True
