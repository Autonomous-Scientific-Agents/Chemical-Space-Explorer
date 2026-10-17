"""Strict rule tier (ADR 0015): per-rule cases and an independent brute-force check of the pruned tree."""
import pytest

from chemuniverse import RuleConfig
from chemuniverse.bonds import molecules
from chemuniverse.export import molecule_smiles
from chemuniverse.graph import bits, simple_cycles
from chemuniverse.tree import children, enumerate_exact, roots

from test_canon import from_edges

BASIC, STRICT = RuleConfig(), RuleConfig.strict()


def removed(labels, edges, cfg=STRICT):
    skel = from_edges(labels, edges)
    basic = set(molecule_smiles(skel, BASIC))
    return basic, basic - set(molecule_smiles(skel, cfg))


@pytest.mark.parametrize("labels,edges,gone,kept", [
    ("NNN", [(0, 1), (1, 2)], {"NNN", "NN=N"}, set()),  # S1 heteroatom chain
    ("CNN", [(0, 1), (1, 2)], set(), {"CNN", "CN=N", "C=NN"}),  # isolated N-N kept
    ("CSSC", [(0, 1), (1, 2), (2, 3)], set(), {"CSSC"}),  # disulfide kept
    ("CSO", [(0, 1), (1, 2)], {"CSO"}, set()),  # S2 S-O
    ("CNS", [(0, 1), (1, 2), (2, 0)], {"C1NS1"}, set()),  # S3 3-ring, 2 heteroatoms
    ("CCNN", [(0, 1), (1, 2), (2, 3), (3, 0)], {"C1CNN1"}, set()),  # S3 4-ring
    ("CCN", [(0, 1), (1, 2), (2, 0)], set(), {"C1CN1"}),  # aziridine kept
    ("CCN", [(0, 1), (1, 2)], {"C=CN", "C#CN"}, {"CCN", "CC=N", "CC#N"}),  # S4 enamine, ynamine
    ("CNO", [(0, 1), (0, 2)], {"C(=N)O"}, {"C(N)=O", "C(#N)O"}),  # S4 imidic acid; amide kept
    ("CSO", [(0, 1), (0, 2)], {"C(=S)O", "C(S)O"}, {"C(S)=O", "C(=S)=O"}),  # S4 thiono acid, S5; COS kept
    ("CCO", [(0, 1), (1, 2)], {"C=C=O", "C#CO"}, {"CCO", "CC=O"}),  # S4 ketene, ynol
    ("SCS", [(0, 1), (1, 2)], {"SCS"}, {"S=C=S", "SC=S"}),  # S5 gem-dithiol; CS2 kept
    ("CCS", [(0, 1), (1, 2)], {"CC=S", "C=CS", "C#CS", "C=C=S"}, {"CCS"}),  # S6, S4; thiol kept
    ("CCNS", [(0, 1), (0, 2), (0, 3)], {"C(C)(N)S", "C(C)(=N)S"}, {"C(C)(N)=S"}),  # thioamide kept
])
def test_rule_cases(labels, edges, gone, kept):
    basic, out = removed(labels, edges)
    assert gone <= basic and kept <= basic
    assert gone <= out and not (kept & out)


def test_single_rule_switches():
    _, out = removed("CSO", [(0, 1), (1, 2)], RuleConfig.strict(ban_s_o=False))
    assert out == set()
    _, out = removed("CCS", [(0, 1), (1, 2)], RuleConfig.strict(thiocarbonyl=False, tautomers=False))
    assert out == set()


HET = set("NOS")


def strict_violation(skel, orders, h):
    """Independent re-implementation of S1-S6 on a finished molecule."""
    lab, n = skel.labels, skel.n
    nb = [list(bits(a)) for a in skel.adj]

    def o(u, v):
        return orders[(u, v) if u < v else (v, u)]

    for x in range(n):
        if lab[x] in HET and sum(lab[y] in HET for y in nb[x]) >= 2:
            return "S1"
        if any(x < y and {lab[x], lab[y]} == {"S", "O"} for y in nb[x]):
            return "S2"
    for c in simple_cycles(skel.adj, 4):
        if sum(lab[v] in HET for v in c) >= 2:
            return "S3"
    for x in range(n):
        if lab[x] != "C":
            continue
        dbl = [y for y in nb[x] if o(x, y) == 2]
        het1 = [y for y in nb[x] if lab[y] in HET and o(x, y) == 1]
        het1h = [y for y in het1 if h[y]]
        cdbl = any(lab[y] == "C" for y in dbl)
        if (cdbl and any(lab[y] in "SN" for y in het1h)) or (het1 and any(lab[y] == "C" and o(x, y) == 3 for y in nb[x])) \
                or (any(lab[y] == "N" for y in dbl) and any(lab[y] in "OS" for y in het1h)) \
                or (any(lab[y] == "S" for y in dbl) and any(lab[y] == "O" for y in het1h)) \
                or (cdbl and any(lab[y] == "O" for y in dbl)):
            return "S4"
        if all(o(x, y) == 1 for y in nb[x]) and len(het1) >= 2 and het1h and any(lab[y] == "S" for y in het1):
            return "S5"
        if any(lab[y] == "S" for y in dbl) and sum(lab[y] in HET for y in nb[x]) == 1:
            return "S6"
    return None


def test_pruned_tree_equals_filtered_basic():
    """Strict counts from the pruned tree == basic molecules that pass an independent strict check."""
    el = ["C", "N", "O", "S"]
    strict_counts = [r.count for r in enumerate_exact(6, el, STRICT)[0]]
    level, filtered = roots(el), []
    for n in range(1, 7):
        if n > 1:
            level = [c for s in level for c in children(s, el, BASIC)]
        filtered.append(sum(1 for s in level for _, orders, h in molecules(s, BASIC, with_orders=True)
                            if strict_violation(s, orders, h) is None))
    assert strict_counts == filtered
