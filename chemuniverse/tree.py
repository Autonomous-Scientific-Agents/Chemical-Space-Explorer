"""Generation tree of labeled skeletons and the two ways of using it.

Children are built by McKay-style canonical augmentation: add one atom, keep
the child only if the new atom is in the orbit of the child's canonically
chosen non-cut vertex. Every connected skeleton then appears exactly once.
"""
import math
import random
from dataclasses import dataclass
from itertools import combinations

from .bonds import count_molecules
from .canon import canonicalize, initial_colors, refine, same_orbit, _nbr_lists
from .graph import Skeleton, non_cut_vertices
from .rules import RuleConfig


def roots(elements):
    return [Skeleton((e,), (0,)) for e in elements]


def _accept(child, new):
    """Canonical-deletion test: is `new` in the orbit of the canonical non-cut vertex?"""
    colors = refine(_nbr_lists(child.adj), initial_colors(child.labels))
    nc = non_cut_vertices(child.adj)
    cmin = min(colors[v] for v in nc)
    if colors[new] != cmin:
        return False, None
    cands = [v for v in nc if colors[v] == cmin]
    if len(cands) == 1:
        return True, None
    form, perm = canonicalize(child.labels, child.adj)
    m = min(cands, key=perm.__getitem__)
    return same_orbit(child.labels, child.adj, m, new), form


def children(skel, elements, cfg):
    n = skel.n
    spare = [cfg.max_valence(lab) - bin(a).count("1") for lab, a in zip(skel.labels, skel.adj)]
    out = {}
    for e in elements:
        eligible = [v for v in range(n) if spare[v] > 0 and cfg.pair_allowed(e, skel.labels[v])]
        max_size = min(cfg.max_valence(e), len(eligible), 1 if not cfg.allow_rings else n)
        for size in range(1, max_size + 1):
            for subset in combinations(eligible, size):
                child = skel.add_vertex(e, sum(1 << v for v in subset))
                ok, form = _accept(child, n)
                if not ok:
                    continue
                if form is None:
                    form = canonicalize(child.labels, child.adj)[0]
                out.setdefault(form, child)
    return list(out.values())


@dataclass
class LevelResult:
    n: int
    count: float
    se: float = 0.0
    exact: bool = True
    skeletons: float = 0.0


def enumerate_exact(max_atoms, elements, cfg, payload=count_molecules):
    """Breadth-first traversal; returns (results per level, nodes of the last level)."""
    level = roots(elements)
    results = []
    for n in range(1, max_atoms + 1):
        if n > 1:
            level = [c for s in level for c in children(s, elements, cfg)]
        results.append(LevelResult(n, sum(payload(s, cfg) for s in level), skeletons=len(level)))
    return results, level


def _summarize(samples, start_level):
    """samples: one (molecules, skeletons) vector per independent sample, per level."""
    results = []
    m = len(samples)
    for d in range(len(samples[0])):
        xs = [smp[d][0] for smp in samples]
        mean = sum(xs) / m
        var = sum((x - mean) ** 2 for x in xs) / (m - 1) if m > 1 else 0.0
        skel = sum(smp[d][1] for smp in samples) / m
        results.append(LevelResult(start_level + d + 1, mean, math.sqrt(var / m), False, skel))
    return results


def knuth_estimate(start_nodes, start_level, max_atoms, elements, cfg, probes, rng, payload=count_molecules):
    """Knuth's random-probe estimate of per-level payload sums below `start_level`.

    Unbiased but heavy-tailed: rare probes with huge weights dominate the mean.
    """
    depth = max_atoms - start_level
    samples = []
    for _ in range(probes):
        smp = [(0.0, 0.0)] * depth
        node = rng.choice(start_nodes)
        weight = float(len(start_nodes))
        for d in range(depth):
            ch = children(node, elements, cfg)
            if not ch:
                break
            weight *= len(ch)
            node = rng.choice(ch)
            smp[d] = (weight * payload(node, cfg), weight)
        samples.append(smp)
    return _summarize(samples, start_level)


def stratum(skel, cfg):
    """Coarse node type for stratified sampling: atoms per free-valence count, ring count."""
    free = [0] * 7
    for lab, a in zip(skel.labels, skel.adj):
        free[cfg.max_valence(lab) - bin(a).count("1")] += 1
    return tuple(free), len(skel.edges()) - skel.n


def stratified_estimate(start_nodes, start_level, max_atoms, elements, cfg, passes, rng, payload=count_molecules):
    """Chen's stratified sampling of the tree; each pass is one unbiased estimate.

    Per level, each stratum keeps one representative node whose weight is the
    total weight of all nodes merged into it (replaced with probability
    w / total, which keeps the estimate unbiased).
    """
    def merge(level, node, w):
        key = stratum(node, cfg)
        slot = level.get(key)
        if slot is None:
            level[key] = [node, w]
        else:
            slot[1] += w
            if rng.random() < w / slot[1]:
                slot[0] = node

    samples = []
    for _ in range(passes):
        level = {}
        for node in start_nodes:
            merge(level, node, 1.0)
        smp = []
        for _n in range(start_level + 1, max_atoms + 1):
            nxt = {}
            for node, w in level.values():
                for c in children(node, elements, cfg):
                    merge(nxt, c, w)
            level = nxt
            smp.append((sum(w * payload(c, cfg) for c, w in level.values()), sum(w for _, w in level.values())))
        samples.append(smp)
    return _summarize(samples, start_level)


METHODS = {"stratified": (stratified_estimate, 10), "knuth": (knuth_estimate, 1000)}


def estimate(max_atoms, elements, cfg=None, k_exact=6, samples=None, seed=0, method="stratified"):
    """Exact counts up to `k_exact` heavy atoms, sampled estimates beyond.

    samples: stratified passes or Knuth probes (default 10 or 1000).
    """
    cfg = cfg or RuleConfig()
    elements = list(dict.fromkeys(elements))
    sampler, default = METHODS[method]
    k = min(k_exact, max_atoms)
    results, frontier = enumerate_exact(k, elements, cfg)
    if max_atoms > k and frontier:
        results += sampler(frontier, k, max_atoms, elements, cfg, samples or default, random.Random(seed))
    return results
