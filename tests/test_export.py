import gzip
import json
import re

import pytest

from chemuniverse import RuleConfig
from chemuniverse.canon import canonical_form
from chemuniverse.export import export, molecule_smiles, to_smiles
from chemuniverse.graph import Skeleton
from chemuniverse.tree import enumerate_exact

from test_canon import from_edges

TOKEN = re.compile(r"\[[^\]]+\]|Cl|Br|[BCNOPSFI]|%\d\d|\d|[=#()]")
ORDER = {"": 1, "=": 2, "#": 3}


def parse(smiles, cfg=RuleConfig()):
    """Minimal parser for the SMILES we write: returns (canonical form of (element, H) graph, orders)."""
    labels, explicit_h, orders = [], [], {}
    stack, prev, bond, rings = [], None, "", {}
    pos = 0
    for tok in TOKEN.findall(smiles):
        pos += len(tok)
        if tok in "=#" and tok:
            bond = tok
        elif tok == "(":
            stack.append(prev)
        elif tok == ")":
            prev = stack.pop()
        elif tok[0].isdigit() or tok[0] == "%":
            d = int(tok.lstrip("%"))
            if d in rings:
                other, b = rings.pop(d)
                orders[(min(other, prev), max(other, prev))] = ORDER[bond or b]
            else:
                rings[d] = (prev, bond)
            bond = ""
        else:
            m = re.fullmatch(r"\[([A-Z][a-z]?)(H(\d*))?\]", tok)
            labels.append(m.group(1) if m else tok)
            explicit_h.append((int(m.group(3) or 1) if m.group(2) else 0) if m else None)
            v = len(labels) - 1
            if prev is not None:
                orders[(prev, v)] = ORDER[bond]
            prev, bond = v, ""
    assert pos == len(smiles) and not rings and not stack, smiles
    n = len(labels)
    adj = [0] * n
    bond_sum = [0] * n
    for (u, v), o in orders.items():
        adj[u] |= 1 << v
        adj[v] |= 1 << u
        bond_sum[u] += o
        bond_sum[v] += o
    h = [eh if eh is not None else next(x for x in cfg.valences(lab) if x >= s) - s
         for lab, eh, s in zip(labels, explicit_h, bond_sum)]
    return canonical_form(list(zip(labels, h)), tuple(adj))


def test_simple_smiles():
    assert sorted(molecule_smiles(from_edges("CC", [(0, 1)]), RuleConfig())) == ["C#C", "C=C", "CC"]
    benzene_ring = from_edges("C" * 6, [(i, (i + 1) % 6) for i in range(6)])
    smiles = molecule_smiles(benzene_ring, RuleConfig())
    assert len(smiles) == 5 and sum(s.count("=") == 3 for s in smiles) == 1  # benzene once
    si = Skeleton(("Si", "C"), (0b10, 0b01))
    assert to_smiles(si.labels, si.adj, {(0, 1): 1}, (3, 3)) == "[SiH3]C"


@pytest.mark.parametrize("elements,max_atoms", [(["C"], 8), (["C", "N", "O", "F"], 6), (["C", "S", "Cl"], 5)])
def test_export_round_trip(tmp_path, elements, max_atoms):
    cfg = RuleConfig()
    counts = export(max_atoms, elements, cfg, tmp_path, jobs=2, log=lambda *_: None)
    exact = [int(r.count) for r in enumerate_exact(max_atoms, elements, cfg)[0]]
    assert counts == exact
    manifest = json.loads((tmp_path / "manifest.json").read_text())
    assert manifest["total"] == sum(exact)
    forms = set()
    for n in range(1, max_atoms + 1):
        lines = [s for f in sorted(tmp_path.glob(f"n{n:02d}_part*.smi.gz"))
                 for s in gzip.open(f, "rt").read().split()]
        assert len(lines) == exact[n - 1]
        for s in lines:
            form = parse(s)
            assert form not in forms, s
            forms.add(form)


def test_export_does_not_depend_on_jobs(tmp_path):
    a, b = tmp_path / "a", tmp_path / "b"
    export(7, ["C", "O"], RuleConfig(), a, jobs=1, log=lambda *_: None)
    export(7, ["C", "O"], RuleConfig(), b, jobs=3, log=lambda *_: None)
    files = sorted(p.name for p in a.glob("*.smi.gz"))
    assert files == sorted(p.name for p in b.glob("*.smi.gz"))
    for name in files:
        assert gzip.open(a / name).read() == gzip.open(b / name).read()
