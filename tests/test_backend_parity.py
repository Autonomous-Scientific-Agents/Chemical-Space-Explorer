"""The Rust core must reproduce the Python reference exactly: same children in the same order,
same counts, same SMILES. Skipped when chemuniverse-core is not installed."""
import pytest

from chemuniverse import RuleConfig
from chemuniverse.backend import core
from chemuniverse.bonds import count_molecules_py, count_molecules_rust
from chemuniverse.export import molecule_smiles_py, molecule_smiles_rust
from chemuniverse.tree import children_py, children_rust, roots

pytestmark = pytest.mark.skipif(core is None, reason="chemuniverse-core (Rust) not installed")

CONFIGS = [
    RuleConfig(),
    RuleConfig(strain=False, heteroatom=False),
    RuleConfig(hypervalent=True, ban_n_o=False),
    RuleConfig(allow_rings=False, allow_multiple_bonds=False),
    RuleConfig(bond_order_limits=False),
]


@pytest.mark.parametrize("elements,max_atoms", [
    (["C", "N", "O", "F"], 6), (["C"], 8), (["C", "S", "P", "Cl", "Si", "B"], 5), (["N", "O", "Br", "I"], 6)])
@pytest.mark.parametrize("cfg", CONFIGS, ids=range(len(CONFIGS)))
def test_tree_payload_and_smiles_identical(elements, max_atoms, cfg):
    level = roots(elements)
    for n in range(1, max_atoms + 1):
        for s in level:
            assert count_molecules_rust(s, cfg) == count_molecules_py(s, cfg), s
            assert molecule_smiles_rust(s, cfg) == molecule_smiles_py(s, cfg), s
        if n < max_atoms:
            nxt = []
            for s in level:
                py = children_py(s, elements, cfg)
                assert children_rust(s, elements, cfg) == py, s
                nxt += py
            level = nxt
