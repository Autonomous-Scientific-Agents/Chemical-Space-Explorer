"""Generation tree of labeled skeletons and the two ways of using it.

Children are built by McKay-style canonical augmentation: add one atom, keep
the child only if the new atom is in the orbit of the child's canonically
chosen non-cut vertex. Every connected skeleton then appears exactly once.
"""
import math
import multiprocessing
import random
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass
from functools import lru_cache
from itertools import combinations

from .bonds import count_molecules
from .canon import _nbr_lists, canonicalize, initial_colors, invariant_hash, refine, same_orbit
from .graph import Skeleton, non_cut_vertices
from .rules import RuleConfig


def roots(elements):
    return [Skeleton((e,), (0,)) for e in elements]


def _accept(child, new):
    """Canonical-deletion test: is `new` in the orbit of the canonical non-cut vertex?

    The canonical vertex minimizes (element, degree, refined color), then
    canonical position. Cheap keys come first so most children are rejected
    before color refinement or canonical labeling runs.
    """
    labels, adj = child.labels, child.adj
    deg = [bin(a).count("1") for a in adj]
    key_new = (labels[new], deg[new])
    # leaves are never cut vertices
    if any(d == 1 and (labels[v], 1) < key_new for v, d in enumerate(deg)):
        return False, None
    nc = non_cut_vertices(adj)
    if min((labels[v], deg[v]) for v in nc) < key_new:
        return False, None
    cands = [v for v in nc if (labels[v], deg[v]) == key_new]
    if len(cands) == 1:
        return True, None
    colors = refine(_nbr_lists(adj), initial_colors(labels))
    cmin = min(colors[v] for v in cands)
    if colors[new] != cmin:
        return False, None
    cands = [v for v in cands if colors[v] == cmin]
    if len(cands) == 1:
        return True, None
    form, perm = canonicalize(labels, adj)
    m = min(cands, key=perm.__getitem__)
    return same_orbit(labels, adj, m, new), form


def children(skel, elements, cfg):
    n = skel.n
    spare = [cfg.max_valence(lab) - bin(a).count("1") for lab, a in zip(skel.labels, skel.adj)]
    buckets = {}  # invariant hash -> [[child, canonical form or None], ...]
    for e in elements:
        eligible = [v for v in range(n) if spare[v] > 0 and cfg.pair_allowed(e, skel.labels[v])]
        max_size = min(cfg.max_valence(e), len(eligible), 1 if not cfg.allow_rings else n)
        for size in range(1, max_size + 1):
            for subset in combinations(eligible, size):
                child = skel.add_vertex(e, sum(1 << v for v in subset))
                ok, form = _accept(child, n)
                if not ok:
                    continue
                # siblings can be isomorphic; canonicalize only on hash collisions
                bucket = buckets.setdefault(invariant_hash(child.labels, child.adj), [])
                if bucket:
                    if form is None:
                        form = canonicalize(child.labels, child.adj)[0]
                    for entry in bucket:
                        if entry[1] is None:
                            entry[1] = canonicalize(entry[0].labels, entry[0].adj)[0]
                    if any(entry[1] == form for entry in bucket):
                        continue
                bucket.append([child, form])
    return [child for bucket in buckets.values() for child, _ in bucket]


@lru_cache(maxsize=1 << 16)
def cached_children(skel, elements, cfg):
    return tuple(children(skel, elements, cfg))


@lru_cache(maxsize=1 << 18)
def cached_count(skel, cfg):
    return count_molecules(skel, cfg)


@dataclass
class LevelResult:
    n: int
    count: float
    se: float = 0.0
    exact: bool = True
    skeletons: float = 0.0


def _pmap(fn, tasks, jobs):
    """map(fn, tasks) in order, in `jobs` worker processes when jobs > 1."""
    if jobs <= 1 or len(tasks) <= 1:
        return [fn(*t) for t in tasks]
    # fork (where available) starts faster and, unlike spawn, works from notebooks and stdin scripts
    ctx = multiprocessing.get_context("fork" if "fork" in multiprocessing.get_all_start_methods() else None)
    with ProcessPoolExecutor(min(jobs, len(tasks)), mp_context=ctx) as pool:
        return list(pool.map(fn, *zip(*tasks)))


def _chunks(items, k):
    """Split `items` into at most k contiguous, nearly equal chunks."""
    size = -(-len(items) // k)
    return [items[i:i + size] for i in range(0, len(items), size)]


def _expand(nodes, depth, elements, cfg, payload):
    """Breadth-first expansion of `nodes` by `depth` levels; per-level (payload, skeletons) and last level."""
    rows = []
    for _ in range(depth):
        nodes = [c for s in nodes for c in children(s, elements, cfg)]
        rows.append((sum(payload(s, cfg) for s in nodes), len(nodes)))
    return rows, nodes


def enumerate_exact(max_atoms, elements, cfg, payload=count_molecules, jobs=1):
    """Breadth-first traversal; returns (results per level, nodes of the last level).

    With jobs > 1 the tree is grown serially until it is wide enough, then the
    subtrees are expanded in worker processes. Canonical augmentation makes the
    subtrees disjoint, so their counts simply add up; level order is unchanged.
    """
    level = roots(elements)
    results = []
    for n in range(1, max_atoms + 1):
        if n > 1:
            if jobs > 1 and len(level) >= 16 * jobs:
                tasks = [(c, max_atoms - n + 1, elements, cfg, payload) for c in _chunks(level, 4 * jobs)]
                parts = _pmap(_expand, tasks, jobs)
                level = [s for _, nodes in parts for s in nodes]
                for d in range(max_atoms - n + 1):
                    results.append(LevelResult(n + d, sum(p[0][d][0] for p in parts),
                                               skeletons=sum(p[0][d][1] for p in parts)))
                return results, level
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


def _knuth_probes(seeds, start_nodes, depth, elements, cfg, payload):
    samples = []
    for seed in seeds:
        rng = random.Random(seed)
        smp = [(0.0, 0.0)] * depth
        node = rng.choice(start_nodes)
        weight = float(len(start_nodes))
        for d in range(depth):
            ch = cached_children(node, tuple(elements), cfg)
            if not ch:
                break
            weight *= len(ch)
            node = rng.choice(ch)
            smp[d] = (weight * payload(node, cfg), weight)
        samples.append(smp)
    return samples


def knuth_estimate(start_nodes, start_level, max_atoms, elements, cfg, probes, rng, payload=cached_count, jobs=1):
    """Knuth's random-probe estimate of per-level payload sums below `start_level`.

    Unbiased but heavy-tailed: rare probes with huge weights dominate the mean.
    Each probe gets its own seed drawn from `rng`, so results do not depend on `jobs`.
    """
    seeds = [rng.getrandbits(64) for _ in range(probes)]
    tasks = [(c, start_nodes, max_atoms - start_level, elements, cfg, payload) for c in _chunks(seeds, 4 * jobs)]
    samples = [smp for part in _pmap(_knuth_probes, tasks, jobs) for smp in part]
    return _summarize(samples, start_level)


def stratum(skel, cfg):
    """Coarse node type for stratified sampling: atoms per free-valence count, ring count."""
    free = [0] * 7
    for lab, a in zip(skel.labels, skel.adj):
        free[cfg.max_valence(lab) - bin(a).count("1")] += 1
    return tuple(free), len(skel.edges()) - skel.n


def _stratified_pass(seed, start_nodes, start_level, max_atoms, elements, cfg, payload):
    """One pass of Chen's stratified sampling; per level (payload estimate, skeleton estimate).

    Per level, each stratum keeps one representative node whose weight is the
    total weight of all nodes merged into it (replaced with probability
    w / total, which keeps the estimate unbiased).
    """
    rng = random.Random(seed)

    def merge(level, node, w):
        key = stratum(node, cfg)
        slot = level.get(key)
        if slot is None:
            level[key] = [node, w]
        else:
            slot[1] += w
            if rng.random() < w / slot[1]:
                slot[0] = node

    level = {}
    for node in start_nodes:
        merge(level, node, 1.0)
    smp = []
    for _n in range(start_level + 1, max_atoms + 1):
        nxt = {}
        for node, w in level.values():
            for c in cached_children(node, tuple(elements), cfg):
                merge(nxt, c, w)
        level = nxt
        smp.append((sum(w * payload(c, cfg) for c, w in level.values()), sum(w for _, w in level.values())))
    return smp


def stratified_estimate(start_nodes, start_level, max_atoms, elements, cfg, passes, rng, payload=cached_count,
                        jobs=1):
    """Chen's stratified sampling of the tree; each pass is one unbiased estimate.

    Each pass gets its own seed drawn from `rng`, so results do not depend on `jobs`.
    """
    tasks = [(rng.getrandbits(64), start_nodes, start_level, max_atoms, elements, cfg, payload)
             for _ in range(passes)]
    return _summarize(_pmap(_stratified_pass, tasks, jobs), start_level)


METHODS = {"stratified": (stratified_estimate, 10), "knuth": (knuth_estimate, 1000)}


def estimate(max_atoms, elements, cfg=None, k_exact=6, samples=None, seed=0, method="stratified", jobs=1):
    """Exact counts up to `k_exact` heavy atoms, sampled estimates beyond.

    samples: stratified passes or Knuth probes (default 10 or 1000).
    jobs: worker processes; results are identical for any value.
    """
    cfg = cfg or RuleConfig()
    elements = list(dict.fromkeys(elements))
    sampler, default = METHODS[method]
    k = min(k_exact, max_atoms)
    results, frontier = enumerate_exact(k, elements, cfg, jobs=jobs)
    if max_atoms > k and frontier:
        results += sampler(frontier, k, max_atoms, elements, cfg, samples or default, random.Random(seed), jobs=jobs)
    return results
