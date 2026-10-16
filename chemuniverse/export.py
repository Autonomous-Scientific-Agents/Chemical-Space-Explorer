"""Write every molecule of an exact enumeration as SMILES (a GDB-style dataset).

Each molecule is generated exactly once (canonical augmentation, one key per
hydrogen pattern), so no canonical SMILES or global deduplication is needed.
Output: DIR/n{size}_part{k}.smi.gz (one SMILES per line) and DIR/manifest.json.
"""
import gzip
import json
import time
from dataclasses import asdict
from pathlib import Path

from .backend import cfg_tuple, core
from .bonds import molecules
from .graph import bits
from .tree import _chunks, _pmap, children, roots

ORGANIC = {"B", "C", "N", "O", "P", "S", "F", "Cl", "Br", "I"}
BOND = {1: "", 2: "=", 3: "#"}
TASKS = 128  # subtree chunks; fixed so the files do not depend on the number of jobs


def to_smiles(labels, adj, orders, h):
    """Kekule SMILES for one molecule.

    orders: {(u, v): order} with u < v. Organic-subset atoms are written bare:
    SMILES gives them the lowest normal valence >= their bond sum, the same rule
    used for implicit hydrogens here (including hypervalent S and P). Other
    elements are bracketed with an explicit H count.
    """
    n = len(labels)
    nbrs = [list(bits(a)) for a in adj]
    order = lambda u, v: orders[(u, v) if u < v else (v, u)]  # noqa: E731
    state = [0] * n  # 0 unseen, 1 on the DFS stack, 2 done
    kids = [[] for _ in range(n)]
    opens = [[] for _ in range(n)]  # ring bonds opened at this atom (it is the ancestor)
    closes = [[] for _ in range(n)]  # ring bonds closed at this atom

    def dfs(v, parent):
        state[v] = 1
        for w in nbrs[v]:
            if w == parent:
                continue
            if state[w] == 0:
                kids[v].append(w)
                dfs(w, v)
            elif state[w] == 1:  # back edge to an ancestor
                ring = (w, v)
                opens[w].append(ring)
                closes[v].append(ring)
        state[v] = 2

    dfs(0, -1)
    digit = {}
    free = []
    next_digit = [1]

    def take():
        if free:
            free.sort()
            return free.pop(0)
        next_digit[0] += 1
        return next_digit[0] - 1

    def label(d):
        return str(d) if d < 10 else f"%{d}"

    def emit(v):
        lab = labels[v]
        out = [lab if lab in ORGANIC else f"[{lab}{'H' + (str(h[v]) if h[v] > 1 else '') if h[v] else ''}]"]
        closed = [digit.pop(ring) for ring in closes[v]]
        out += [BOND[order(*ring)] + label(d) for ring, d in zip(closes[v], closed)]
        for ring in opens[v]:
            digit[ring] = take()
            out.append(label(digit[ring]))
        free.extend(closed)  # reusable only after this atom, so it never prints "11"
        for i, w in enumerate(kids[v]):
            branch = BOND[order(v, w)] + emit(w)
            out.append(branch if i == len(kids[v]) - 1 else f"({branch})")
        return "".join(out)

    return emit(0)


def molecule_smiles_py(skel, cfg):
    return [to_smiles(skel.labels, skel.adj, orders, h) for _, orders, h in molecules(skel, cfg, with_orders=True)]


def molecule_smiles_rust(skel, cfg):
    return core.molecule_smiles(skel.labels, skel.adj, cfg_tuple(cfg))


molecule_smiles = molecule_smiles_rust if core else molecule_smiles_py


def _export_subtrees(part, nodes, start, max_atoms, elements, cfg, out):
    """Depth-first export below `nodes` (at level `start`); returns molecule counts per level."""
    counts = [0] * (max_atoms - start)
    files = {}

    def write(n, smiles):
        if not smiles:
            return
        f = files.get(n)
        if f is None:
            f = files[n] = gzip.open(Path(out) / f"n{n:02d}_part{part:04d}.smi.gz", "wt")
        f.write("\n".join(smiles) + "\n")
        counts[n - start - 1] += len(smiles)

    def visit(skel):
        n = skel.n
        if n > start:
            write(n, molecule_smiles(skel, cfg))
        if n < max_atoms:
            for c in children(skel, elements, cfg):
                visit(c)

    for s in nodes:
        visit(s)
    for f in files.values():
        f.close()
    return counts


def export(max_atoms, elements, cfg, out, jobs=1, log=print):
    """Write all molecules up to `max_atoms` heavy atoms; returns the per-size counts."""
    elements = list(dict.fromkeys(elements))
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    if any(out.glob("*.smi.gz")):
        raise FileExistsError(f"{out} already contains .smi.gz files")
    t0 = time.time()
    counts = []
    level = roots(elements)
    for n in range(1, max_atoms + 1):
        if n > 1:
            level = [c for s in level for c in children(s, elements, cfg)]
        smiles = [smi for s in level for smi in molecule_smiles(s, cfg)]
        with gzip.open(out / f"n{n:02d}_part0000.smi.gz", "wt") as f:
            f.write("".join(smi + "\n" for smi in smiles))
        counts.append(len(smiles))
        log(f"n={n}: {len(smiles)} molecules ({time.time() - t0:.1f} s)")
        # grow serially until the tree is wide enough to split into subtree tasks
        if n < max_atoms and len(level) >= TASKS:
            tasks = [(k + 1, c, n, max_atoms, elements, cfg, str(out)) for k, c in enumerate(_chunks(level, TASKS))]
            parts = _pmap(_export_subtrees, tasks, jobs)
            for d in range(max_atoms - n):
                counts.append(sum(p[d] for p in parts))
                log(f"n={n + d + 1}: {counts[-1]} molecules")
            break
    manifest = {
        "elements": elements,
        "max_atoms": max_atoms,
        "rules": asdict(cfg),
        "counts": {str(i + 1): c for i, c in enumerate(counts)},
        "total": sum(counts),
        "format": "n{size}_part{k}.smi.gz, one Kekule SMILES per line; each molecule appears once",
        "generated": time.strftime("%Y-%m-%d"),
        "seconds": round(time.time() - t0, 1),
    }
    (out / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    return counts
