"""Labeled skeleton graphs stored as a labels tuple plus adjacency bitmasks."""
from typing import NamedTuple


class Skeleton(NamedTuple):
    labels: tuple  # element symbol per vertex
    adj: tuple  # adj[v] is a bitmask of v's neighbors

    @property
    def n(self):
        return len(self.labels)

    def edges(self):
        return [(u, v) for u in range(self.n) for v in bits(self.adj[u]) if u < v]

    def add_vertex(self, label, nbr_mask):
        n = self.n
        adj = tuple(a | (1 << n) if nbr_mask >> v & 1 else a for v, a in enumerate(self.adj))
        return Skeleton(self.labels + (label,), adj + (nbr_mask,))


def bits(mask):
    while mask:
        low = mask & -mask
        yield low.bit_length() - 1
        mask ^= low


def is_connected(adj, mask=None):
    """True if the subgraph induced by `mask` (default: all vertices) is connected."""
    if mask is None:
        mask = (1 << len(adj)) - 1
    if not mask:
        return True
    seen = frontier = mask & -mask
    while frontier:
        nxt = 0
        for v in bits(frontier):
            nxt |= adj[v]
        frontier = nxt & mask & ~seen
        seen |= frontier
    return seen == mask


def non_cut_vertices(adj):
    full = (1 << len(adj)) - 1
    if len(adj) == 1:
        return [0]
    return [v for v in range(len(adj)) if is_connected(adj, full & ~(1 << v))]
