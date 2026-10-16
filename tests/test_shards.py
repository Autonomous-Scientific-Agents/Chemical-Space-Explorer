import gzip

from chemuniverse import RuleConfig, estimate
from chemuniverse import shards
from chemuniverse.export import export
from chemuniverse.tree import enumerate_exact

EL = ["C", "N", "O", "F"]


def run_all(run_dir):
    for i in shards.pending(run_dir):
        shards.run_shard(run_dir, i)


def test_exact_shards_match_enumeration(tmp_path):
    shards.plan(tmp_path, "exact", 6, EL, shards=50, nodes_per_shard=4, log=lambda *_: None)
    run_all(tmp_path)
    results, stats = shards.reduce(tmp_path)
    assert results == enumerate_exact(6, EL, RuleConfig())[0]
    assert stats["shards"] == 50


def test_restart_runs_only_missing_shards(tmp_path):
    shards.plan(tmp_path, "exact", 7, ["C", "O"], shards=20, log=lambda *_: None)
    for i in range(0, 20, 2):
        shards.run_shard(tmp_path, i)
    assert shards.pending(tmp_path) == list(range(1, 20, 2))
    assert shards.run_shard(tmp_path, 0) == {"shard": 0, "skipped": True}
    run_all(tmp_path)
    assert shards.reduce(tmp_path)[0] == enumerate_exact(7, ["C", "O"], RuleConfig())[0]


def test_pass_shards_match_estimate(tmp_path):
    shards.plan(tmp_path, "passes", 7, EL, k_exact=4, passes=6, seed=4, log=lambda *_: None)
    run_all(tmp_path)
    assert shards.reduce(tmp_path)[0] == estimate(7, EL, k_exact=4, samples=6, seed=4)


def test_export_shards_match_export(tmp_path):
    shards.plan(tmp_path / "run", "export", 6, EL, shards=30, log=lambda *_: None)
    run_all(tmp_path / "run")
    results, _ = shards.reduce(tmp_path / "run")
    counts = export(6, EL, RuleConfig(), tmp_path / "ref", log=lambda *_: None)
    assert [r.count for r in results] == counts

    def smiles(d):
        return sorted(s for f in d.glob("*.smi.gz") for s in gzip.open(f, "rt").read().split())

    assert smiles(tmp_path / "run" / "smiles") == smiles(tmp_path / "ref")
