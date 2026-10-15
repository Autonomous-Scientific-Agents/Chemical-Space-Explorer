"""Count distinct stable molecules on one skeleton (the tree-node payload)."""
from .canon import canonical_form
from .rules import RingInfo, edge_cap, molecule_ok, skeleton_ok


def molecules(skel, cfg):
    """Yield the distinct (element, H) labelings of `skel` that pass all rules.

    Bond assignments that give the same hydrogen counts (Kekule/resonance
    structures) are one molecule; it counts if any of its assignments passes.
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
    accepted = set()
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
            accepted.add(h)

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
    seen = set()
    for h in accepted:
        form = canonical_form(list(zip(labels, h)), skel.adj)
        if form not in seen:
            seen.add(form)
            yield form


def count_molecules(skel, cfg):
    return sum(1 for _ in molecules(skel, cfg))
