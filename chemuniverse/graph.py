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
    """Vertices whose removal leaves a connected graph (complement of Tarjan's articulation points)."""
    n = len(adj)
    if n <= 2:
        return list(range(n))
    nbrs = [list(bits(a)) for a in adj]
    disc = [-1] * n
    low = [0] * n
    cut = [False] * n
    disc[0] = 0
    t = 1
    root_children = 0
    stack = [(0, -1, iter(nbrs[0]))]
    while stack:
        v, parent, it = stack[-1]
        for w in it:
            if disc[w] < 0:
                disc[w] = low[w] = t
                t += 1
                stack.append((w, v, iter(nbrs[w])))
                break
            if w != parent and disc[w] < low[v]:
                low[v] = disc[w]
        else:
            stack.pop()
            if parent >= 0:
                if low[v] < low[parent]:
                    low[parent] = low[v]
                if parent == 0:
                    root_children += 1
                elif low[v] >= disc[parent]:
                    cut[parent] = True
    cut[0] = root_children > 1
    return [v for v in range(n) if not cut[v]]


def simple_cycles(adj, max_len):
    """All simple cycles with at most `max_len` vertices, each listed once."""
    out = []
    for s in range(len(adj)):
        higher = -1 << (s + 1)  # cycles are rooted at their smallest vertex
        path = [s]

        def extend(used):
            v = path[-1]
            if len(path) >= 3 and adj[v] >> s & 1 and path[1] < v:
                out.append(list(path))
            if len(path) == max_len:
                return
            for w in bits(adj[v] & higher & ~used):
                path.append(w)
                extend(used | 1 << w)
                path.pop()

        extend(1 << s)
    return out
