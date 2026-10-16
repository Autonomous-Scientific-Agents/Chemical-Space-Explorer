"""Optional Rust core (package chemuniverse-core, built from rust/).

CHEMUNIVERSE_BACKEND=auto (default) uses the Rust core when it is installed,
=python forces the pure-Python reference code, =rust requires the core.
Both backends give identical results: same tree, same child order, same counts.
"""
import os
from dataclasses import astuple
from functools import lru_cache

NAME = os.environ.get("CHEMUNIVERSE_BACKEND", "auto")
if NAME not in ("auto", "python", "rust"):
    raise ValueError(f"CHEMUNIVERSE_BACKEND must be auto, python or rust, not {NAME!r}")

core = None
if NAME != "python":
    try:
        import chemuniverse_core as core  # noqa: F401  (re-exported for bonds, tree, export)
    except ImportError:
        if NAME == "rust":
            raise


@lru_cache(maxsize=None)
def cfg_tuple(cfg):
    """RuleConfig as the tuple of booleans the Rust core expects (field order)."""
    return astuple(cfg)
