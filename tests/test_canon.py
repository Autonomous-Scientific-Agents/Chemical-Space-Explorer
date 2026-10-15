import random

from chemuniverse.canon import canonical_form, same_orbit
from chemuniverse.graph import Skeleton


def permuted(skel, perm):
    """Relabel vertex v as perm[v]."""
    n = skel.n
    labels = [None] * n
    adj = [0] * n
    for v in range(n):
        labels[perm[v]] = skel.labels[v]
        adj[perm[v]] = sum(1 << perm[u] for u in range(n) if skel.adj[v] >> u & 1)
    return Skeleton(tuple(labels), tuple(adj))


def from_edges(labels, edges):
    adj = [0] * len(labels)
    for u, v in edges:
        adj[u] |= 1 << v
        adj[v] |= 1 << u
    return Skeleton(tuple(labels), tuple(adj))


CUBANE = from_edges("C" * 8, [(0, 1), (1, 2), (2, 3), (3, 0), (4, 5), (5, 6), (6, 7), (7, 4),
                              (0, 4), (1, 5), (2, 6), (3, 7)])
PYRIDINE_RING = from_edges("CCCCCN", [(i, (i + 1) % 6) for i in range(6)])
NORBORNANE = from_edges("C" * 7, [(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 0), (0, 6), (3, 6)])


def test_permutation_invariance():
    rng = random.Random(0)
    for skel in (CUBANE, PYRIDINE_RING, NORBORNANE):
        form = canonical_form(skel.labels, skel.adj)
        for _ in range(20):
            perm = list(range(skel.n))
            rng.shuffle(perm)
            p = permuted(skel, perm)
            assert canonical_form(p.labels, p.adj) == form


def test_distinguishes_non_isomorphic():
    # 3-regular graphs on 6 vertices: prism vs K3,3
    prism = from_edges("C" * 6, [(0, 1), (1, 2), (2, 0), (3, 4), (4, 5), (5, 3), (0, 3), (1, 4), (2, 5)])
    k33 = from_edges("C" * 6, [(a, b) for a in range(3) for b in range(3, 6)])
    assert canonical_form(prism.labels, prism.adj) != canonical_form(k33.labels, k33.adj)
    # same graph, different labels
    other = from_edges("CCCNCC", [(i, (i + 1) % 6) for i in range(6)])
    assert canonical_form(other.labels, other.adj) == canonical_form(PYRIDINE_RING.labels, PYRIDINE_RING.adj)
    oxa = from_edges("CCCCCO", [(i, (i + 1) % 6) for i in range(6)])
    assert canonical_form(oxa.labels, oxa.adj) != canonical_form(PYRIDINE_RING.labels, PYRIDINE_RING.adj)


def test_orbits():
    s = PYRIDINE_RING  # N at 5; orbits {5}, {0,4}, {1,3}, {2}
    assert same_orbit(s.labels, s.adj, 0, 4)
    assert same_orbit(s.labels, s.adj, 1, 3)
    assert not same_orbit(s.labels, s.adj, 0, 1)
    assert not same_orbit(s.labels, s.adj, 2, 5)
    assert all(same_orbit(CUBANE.labels, CUBANE.adj, 0, v) for v in range(8))
    n = NORBORNANE  # bridgeheads 0 and 3, bridge 6
    assert same_orbit(n.labels, n.adj, 0, 3)
    assert not same_orbit(n.labels, n.adj, 0, 6)
