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


def test_export_shards_finalize_match_export(tmp_path):
    import hashlib
    import json

    run = tmp_path / "run"
    shards.plan(run, "export", 6, EL, shards=30, nodes_per_shard=4, log=lambda *_: None)
    run_all(run)
    results, _ = shards.reduce(run)
    counts = export(6, EL, RuleConfig(), tmp_path / "ref", log=lambda *_: None)
    assert [r.count for r in results] == counts

    one = shards.finalize(run, chunk_bytes=2000, log=lambda *_: None)
    manifest = shards.finalize(run, chunk_bytes=2000, delete_parts=True, jobs=3, log=lambda *_: None)
    assert [f["sha256"] for f in one["files"]] == [f["sha256"] for f in manifest["files"]]  # independent of jobs
    assert not (run / "smiles").exists()
    assert manifest["name"] == "chemuniverse-CNOF-le6" and manifest["tier"] == "basic"
    assert shards.refresh_manifest(run, log=lambda *_: None)["name"] == "chemuniverse-CNOF-le6"
    assert shards.dataset_name({"elements": ["C", "N", "O", "S"], "max_atoms": 12,
                                "rules": vars(RuleConfig.strict())}) == "chemuniverse-CNOS-le12-strict"
    assert manifest["total"] == sum(counts) and [manifest["counts"][str(n)] for n in range(1, 7)] == counts
    assert json.loads((run / "dataset" / "manifest.json").read_text()) == manifest
    assert any(f["size"] == 6 for f in manifest["files"][1:]) and len(manifest["files"]) > 6  # chunked

    def smiles(files):
        return sorted(s for f in files for s in gzip.open(f, "rt").read().split())

    for f in manifest["files"]:
        data = (run / f["file"]).read_bytes()
        assert hashlib.sha256(data).hexdigest() == f["sha256"] and len(data) == f["bytes"]
        assert gzip.decompress(data).count(b"\n") == f["molecules"]
    assert smiles(run.glob("dataset/*.smi.gz")) == smiles((tmp_path / "ref").glob("*.smi.gz"))
    assert shards.verify(run, jobs=2, log=lambda *_: None)["bad"] == []
    first = run / manifest["files"][0]["file"]
    first.write_bytes(first.read_bytes() + gzip.compress(b"CC\n"))  # corrupt: extra molecule
    bad = shards.verify(run, log=lambda *_: None)["bad"]
    assert len(bad) == 1 and "sha256 mismatch" in bad[0]["errors"]
