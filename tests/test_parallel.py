import random

from chemuniverse import RuleConfig, estimate
from chemuniverse.tree import enumerate_exact, knuth_estimate, stratified_estimate

EL = ["C", "N", "O", "F"]


def test_exact_parallel_matches_serial():
    cfg = RuleConfig()
    serial, last = enumerate_exact(6, EL, cfg)
    par, last_par = enumerate_exact(6, EL, cfg, jobs=3)
    assert par == serial
    assert last_par == last  # same nodes in the same order


def test_samplers_do_not_depend_on_jobs():
    cfg = RuleConfig()
    _, frontier = enumerate_exact(4, ["C", "N", "O"], cfg)
    for sampler, n in ((stratified_estimate, 6), (knuth_estimate, 50)):
        one = sampler(frontier, 4, 7, ["C", "N", "O"], cfg, n, random.Random(3))
        three = sampler(frontier, 4, 7, ["C", "N", "O"], cfg, n, random.Random(3), jobs=3)
        assert one == three


def test_estimate_jobs():
    assert estimate(7, EL, k_exact=5, samples=4, seed=2) == estimate(7, EL, k_exact=5, samples=4, seed=2, jobs=2)
