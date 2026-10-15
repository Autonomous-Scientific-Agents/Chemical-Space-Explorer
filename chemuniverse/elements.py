"""Element table: allowed valences and bond-order limits."""
from dataclasses import dataclass


@dataclass(frozen=True)
class Element:
    symbol: str
    valences: tuple  # allowed neutral valences, ascending
    max_bond_order: int = 3
    halogen: bool = False


ELEMENTS = {
    "C": Element("C", (4,)),
    "N": Element("N", (3,)),
    "O": Element("O", (2,), max_bond_order=2),
    "S": Element("S", (2,), max_bond_order=2),
    "P": Element("P", (3,)),
    "Si": Element("Si", (4,)),
    "B": Element("B", (3,)),
    "F": Element("F", (1,), max_bond_order=1, halogen=True),
    "Cl": Element("Cl", (1,), max_bond_order=1, halogen=True),
    "Br": Element("Br", (1,), max_bond_order=1, halogen=True),
    "I": Element("I", (1,), max_bond_order=1, halogen=True),
}

HYPERVALENT = {"S": (2, 4, 6), "P": (3, 5)}


def get(symbol):
    try:
        return ELEMENTS[symbol]
    except KeyError:
        raise ValueError(f"unknown element {symbol!r}; known: {', '.join(ELEMENTS)}") from None
