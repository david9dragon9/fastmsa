"""End-to-end wiring of the MMseqs2 pipeline, driven by a stub binary.

These tests do not check alignment quality -- they check that each module is
handed the database it is supposed to read, that the flat and clustered code
paths differ only where they should, and that the results land in files named
after the queries.
"""

import os
import sys
from pathlib import Path

import pytest

from fastmsa.a3m import TAXONOMY_M8
from fastmsa.search import SearchConfig, search

STUB = Path(__file__).parent / "stub_mmseqs.py"
QUERIES = {"chainA": "MKKD", "chain B": "QQQ"}


@pytest.fixture
def calls(tmp_path, monkeypatch):
    """Run searches against the stub and expose the recorded invocations."""
    log = tmp_path / "calls.log"
    monkeypatch.setenv("STUB_MMSEQS_LOG", str(log))

    def recorded():
        return [line.split("\t") for line in log.read_text().splitlines()]

    return recorded


def make_db(directory: Path, name: str, *, clustered: bool, indexed: bool) -> Path:
    directory.mkdir(parents=True, exist_ok=True)
    Path(directory / f"{name}.dbtype").write_text("stub")
    if indexed:
        Path(directory / f"{name}.idx").write_text("stub")
    if clustered:
        for suffix in ("_seq", "_aln"):
            Path(directory / f"{name}{suffix}.dbtype").write_text("stub")
    return directory


def run_search(tmp_path, *, clustered=False, indexed=False, taxonomy=False):
    db_dir = make_db(tmp_path / "dbs", "testDB", clustered=clustered, indexed=indexed)
    config = SearchConfig(
        clustered=clustered,
        taxonomy=taxonomy,
        mmseqs=Path(sys.executable),
        threads=4,
    )
    # Every module runs as `python stub_mmseqs.py <module> ...`, so the logged
    # arguments line up one-for-one with what the pipeline passed to mmseqs.
    from fastmsa import search as search_module

    original_check_call = search_module.subprocess.check_call

    def check_call(command):
        return original_check_call([command[0], str(STUB), *command[1:]])

    search_module.subprocess.check_call = check_call
    try:
        return search(QUERIES, db_dir, "testDB", tmp_path / "out", config)
    finally:
        search_module.subprocess.check_call = original_check_call


def modules(calls):
    return [call[0] for call in calls]


def call_for(calls, module):
    return next(call for call in calls if call[0] == module)


def test_flat_search_writes_one_a3m_per_query(tmp_path, calls):
    written = run_search(tmp_path)

    assert [path.name for path in written] == ["chainA.a3m", "chain_B.a3m"]
    assert all(path.is_file() for path in written)


def test_flat_search_realigns_against_the_database_itself(tmp_path, calls):
    run_search(tmp_path)
    recorded = calls()

    assert "expandaln" not in modules(recorded)
    for module in ("align", "filterresult", "result2msa"):
        assert call_for(recorded, module)[2].endswith("testDB")


def test_clustered_search_expands_cluster_members(tmp_path, calls):
    run_search(tmp_path, clustered=True)
    recorded = calls()

    expand = call_for(recorded, "expandaln")
    assert expand[2].endswith("testDB_seq")
    assert expand[4].endswith("testDB_aln")
    assert call_for(recorded, "align")[2].endswith("testDB_seq")


def test_an_indexed_clustered_database_is_read_through_its_index(tmp_path, calls):
    run_search(tmp_path, clustered=True, indexed=True)
    recorded = calls()

    expand = call_for(recorded, "expandaln")
    assert expand[2].endswith("testDB.idx")
    assert expand[4].endswith("testDB.idx")


def test_an_unindexed_database_is_not_preloaded(tmp_path, calls):
    run_search(tmp_path)
    search_call = call_for(calls(), "search")

    assert search_call[search_call.index("--db-load-mode") + 1] == "0"


def test_taxonomy_report_is_optional(tmp_path, calls):
    run_search(tmp_path)
    assert "convertalis" not in modules(calls())
    assert not (tmp_path / "out" / TAXONOMY_M8).exists()


def test_taxonomy_report_is_written_next_to_the_alignments(tmp_path, calls):
    run_search(tmp_path, taxonomy=True)
    assert (tmp_path / "out" / TAXONOMY_M8).is_file()


def test_intermediates_are_cleaned_up(tmp_path, calls):
    run_search(tmp_path)
    leftovers = sorted(p.name for p in (tmp_path / "out").iterdir())

    assert leftovers == ["chainA.a3m", "chain_B.a3m"]


def test_a_missing_database_is_reported_before_any_work(tmp_path, calls):
    with pytest.raises(FileNotFoundError, match="absentDB"):
        search(
            QUERIES,
            make_db(tmp_path / "dbs", "testDB", clustered=False, indexed=False),
            "absentDB",
            tmp_path / "out",
            SearchConfig(mmseqs=Path(sys.executable)),
        )
