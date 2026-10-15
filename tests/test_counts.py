import random

import pytest

from chemuniverse import RuleConfig, count_molecules, estimate
from chemuniverse.tree import enumerate_exact, knuth_estimate, stratified_estimate

from test_canon import NORBORNANE, from_edges


def counts(max_atoms, elements, cfg=RuleConfig()):
    return [r.count for r in enumerate_exact(max_atoms, elements, cfg)[0]]


def test_alkanes():  # OEIS A000602
    cfg = RuleConfig(allow_rings=False, allow_multiple_bonds=False)
    assert counts(10, ["C"], cfg) == [1, 1, 1, 2, 3, 5, 9, 18, 35, 75]


def test_carbon_skeletons():  # connected graphs with max degree <= 4
    skel = [r.skeletons for r in enumerate_exact(8, ["C"], RuleConfig())[0]]
    assert skel == [1, 1, 2, 6, 21, 78, 353, 1929]


def test_small_hand_counts():
    assert counts(2, ["C"]) == [1, 3]  # CC, C=C, C#C
    assert counts(2, ["C", "O"]) == [2, 5]  # + CO, C=O; O-O banned
    assert counts(3, ["C"])[2] == 5  # CCC, C=CC, C#CC, C=C=C, cyclopropane


def test_six_ring_with_resonance_merging():
    ring = from_edges("C" * 6, [(i, (i + 1) % 6) for i in range(6)])
    # cyclohexane, cyclohexene, 1,3- and 1,4-cyclohexadiene, benzene (one Kekule class)
    assert count_molecules(ring, RuleConfig()) == 5


def test_strain_rules():
    tri = from_edges("CCC", [(0, 1), (1, 2), (2, 0)])
    assert count_molecules(tri, RuleConfig()) == 1  # no cyclopropene
    # without strain rules: + cyclopropene, -yne, ring allene, C3 with three double bonds
    assert count_molecules(tri, RuleConfig(strain=False)) == 5
    # norbornane, norbornene, norbornadiene; no double bond at a bridgehead
    assert count_molecules(NORBORNANE, RuleConfig()) == 3


def test_heteroatom_rules():
    co2 = from_edges("CO", [(0, 1)])
    ok = RuleConfig()
    assert count_molecules(co2, ok) == 2  # methanol, formaldehyde
    enol = from_edges("CCO", [(0, 1), (1, 2)])
    # only the enol C=C-OH is removed by the heteroatom rules
    with_enol = count_molecules(enol, RuleConfig(heteroatom=False))
    assert count_molecules(enol, ok) == with_enol - 1
    gem = from_edges("OCO", [(0, 1), (1, 2)])
    assert count_molecules(gem, ok) == 2  # formic acid, CO2; the gem-diol is banned
    assert count_molecules(gem, RuleConfig(heteroatom=False)) == 3
    no = from_edges("NO", [(0, 1)])
    assert count_molecules(no, ok) == 0
    assert count_molecules(no, RuleConfig(ban_n_o=False)) == 2  # NH2OH, HN=O


@pytest.mark.parametrize("sampler,samples", [(knuth_estimate, 2000), (stratified_estimate, 20)])
def test_sampler_matches_exact(sampler, samples):
    cfg = RuleConfig()
    el = ["C", "N", "O"]
    _, frontier = enumerate_exact(3, el, cfg)
    exact = enumerate_exact(6, el, cfg)[0]
    for r in sampler(frontier, 3, 6, el, cfg, samples, random.Random(0)):
        assert r.se > 0
        assert abs(r.count - exact[r.n - 1].count) < 4 * r.se
        assert abs(r.skeletons - exact[r.n - 1].skeletons) < 0.5 * exact[r.n - 1].skeletons


def test_estimate_api():
    for method in ("stratified", "knuth"):
        res = estimate(5, ["C", "N"], k_exact=3, samples=5, seed=1, method=method)
        assert [r.exact for r in res] == [True, True, True, False, False]
        assert res == estimate(5, ["C", "N"], k_exact=3, samples=5, seed=1, method=method)


def test_tree_matches_brute_force():
    """Every valid connected skeleton appears exactly once in the tree."""
    from itertools import combinations, product

    from chemuniverse.canon import canonical_form
    from chemuniverse.graph import is_connected
    from chemuniverse.rules import skeleton_ok

    el = ["C", "N", "O", "F"]
    cfg = RuleConfig()
    tree = [r.skeletons for r in enumerate_exact(4, el, cfg)[0]]
    for n in range(1, 5):
        pairs = list(combinations(range(n), 2))
        forms = set()
        for labels in product(el, repeat=n):
            for k in range(len(pairs) + 1):
                for es in combinations(pairs, k):
                    s = from_edges(labels, es)
                    if is_connected(s.adj) and skeleton_ok(cfg, s):
                        forms.add(canonical_form(s.labels, s.adj))
        assert tree[n - 1] == len(forms)
