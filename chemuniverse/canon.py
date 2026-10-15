"""Pure-Python canonical labeling of vertex-labeled simple graphs.

Color refinement followed by an individualization-refinement search that keeps
the lexicographically smallest relabeled adjacency. Subtrees that are images
of already explored ones under discovered automorphisms are pruned.
"""
from .graph import bits


def _nbr_lists(adj):
    return [list(bits(a)) for a in adj]


def initial_colors(labels):
    order = {lab: i for i, lab in enumerate(sorted(set(labels)))}
    return [order[lab] for lab in labels]


def refine(nbrs, colors):
    """Coarsest equitable refinement of an ordered partition (colors are ranks)."""
    k = len(set(colors))
    while True:
        sigs = [(colors[v], tuple(sorted(colors[u] for u in nb))) for v, nb in enumerate(nbrs)]
        uniq = sorted(set(sigs))
        rank = {s: i for i, s in enumerate(uniq)}
        colors = [rank[s] for s in sigs]
        if len(uniq) == k:
            return colors
        k = len(uniq)


def refined_colors(labels, adj):
    """Isomorphism-invariant vertex colors (stable color refinement)."""
    return refine(_nbr_lists(adj), initial_colors(labels))


def _orbit_of(v, gens):
    orbit, stack = {v}, [v]
    while stack:
        x = stack.pop()
        for g in gens:
            y = g[x]
            if y not in orbit:
                orbit.add(y)
                stack.append(y)
    return orbit


def canonicalize(labels, adj):
    """Return (form, perm): a hashable canonical form and perm[v] = canonical position."""
    n = len(labels)
    nbrs = _nbr_lists(adj)
    best = None  # (cert, perm)
    first = None
    autos = []

    def leaf(colors):
        nonlocal best, first
        perm = colors
        inv = [0] * n
        for v, p in enumerate(perm):
            inv[p] = v
        cert = tuple(sum(1 << perm[u] for u in nbrs[inv[p]]) for p in range(n))
        for ref in (first, best):
            if ref is not None and ref[0] == cert:
                rinv = [0] * n
                for v, p in enumerate(ref[1]):
                    rinv[p] = v
                autos.append([rinv[perm[v]] for v in range(n)])
                break
        if first is None:
            first = (cert, perm)
        if best is None or cert < best[0]:
            best = (cert, perm)

    def search(colors, path):
        counts = {}
        for c in colors:
            counts[c] = counts.get(c, 0) + 1
        target = min((c for c, k in counts.items() if k > 1), default=None)
        if target is None:
            leaf(colors)
            return
        tried = set()
        for v in (v for v in range(n) if colors[v] == target):
            if v in tried:
                continue
            child = [2 * c + (c == target and u != v) for u, c in enumerate(colors)]
            search(refine(nbrs, child), path + (v,))
            # skip vertices equivalent to v under automorphisms fixing the path
            gens = [g for g in autos if all(g[p] == p for p in path)]
            tried |= _orbit_of(v, gens)

    search(refine(nbrs, initial_colors(labels)), ())
    cert, perm = best
    return (tuple(sorted(labels)), cert), perm


def canonical_form(labels, adj):
    return canonicalize(labels, adj)[0]


def same_orbit(labels, adj, u, v):
    """True if some automorphism maps u to v."""
    if u == v:
        return True
    mu = [(lab, w == u) for w, lab in enumerate(labels)]
    mv = [(lab, w == v) for w, lab in enumerate(labels)]
    return canonical_form(mu, adj) == canonical_form(mv, adj)
