"""Count distinct stable molecules on one skeleton (the tree-node payload)."""
from .canon import canonical_form, refined_colors
from .rules import RingInfo, edge_cap, molecule_ok, skeleton_ok


def molecules(skel, cfg, with_orders=False):
    """Yield the distinct (element, H) labelings of `skel` that pass all rules.

    Bond assignments that give the same hydrogen counts (Kekule/resonance
    structures) are one molecule; it counts if any of its assignments passes.
    With `with_orders`, yield (key, orders, h) instead of key, where `orders`
    is one passing assignment {(u, v): order} and `h` the hydrogen counts.
    """
    if not skeleton_ok(cfg, skel):
        return
    labels = skel.labels
    n = skel.n
    edges = skel.edges()
    rings = RingInfo(skel)
    deg = [bin(a).count("1") for a in skel.adj]
    spare = [cfg.max_valence(labels[v]) - deg[v] for v in range(n)]
    caps = [edge_cap(cfg, skel, rings, u, v) for u, v in edges]
    valences = [cfg.valences(lab) for lab in labels]
    accepted = {}  # hydrogen counts -> first passing bond orders
    extra = [0] * len(edges)  # bond order minus one, per edge

    def finish():
        bond_sum = [deg[v] for v in range(n)]
        for (u, v), x in zip(edges, extra):
            bond_sum[u] += x
            bond_sum[v] += x
        h = tuple(next(val for val in valences[v] if val >= bond_sum[v]) - bond_sum[v] for v in range(n))
        if h in accepted:
            return
        orders = {e: x + 1 for e, x in zip(edges, extra)}
        if molecule_ok(cfg, skel, rings, orders, h):
            accepted[h] = orders

    def assign(i):
        if i == len(edges):
            finish()
            return
        u, v = edges[i]
        for x in range(min(caps[i] - 1, spare[u], spare[v]) + 1):
            extra[i] = x
            spare[u] -= x
            spare[v] -= x
            assign(i + 1)
            spare[u] += x
            spare[v] += x
        extra[i] = 0

    assign(0)
    if len(set(refined_colors(labels, skel.adj))) == n:
        # no symmetry: every hydrogen pattern is a distinct molecule
        for h, orders in accepted.items():
            key = tuple(zip(labels, h))
            yield (key, orders, h) if with_orders else key
        return
    seen = set()
    for h, orders in accepted.items():
        form = canonical_form(list(zip(labels, h)), skel.adj)
        if form not in seen:
            seen.add(form)
            yield (form, orders, h) if with_orders else form


def count_molecules(skel, cfg):
    return sum(1 for _ in molecules(skel, cfg))
