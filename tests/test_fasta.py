"""Queries may arrive as FASTA or as a JSON mapping; both must agree."""

import json

import pytest

from fastmsa.fasta import load_queries, parse_fasta, safe_filename


def test_multiline_and_described_fasta_entries(tmp_path):
    path = tmp_path / "queries.fasta"
    path.write_text(">chainA description here\nMKKD\nVRIL\n\n>chainB\nQQQ\n")
    assert load_queries(path) == {"chainA": "MKKDVRIL", "chainB": "QQQ"}


def test_json_queries_match_fasta_queries(tmp_path):
    fasta = tmp_path / "queries.fasta"
    fasta.write_text(">chainA\nmkkd\n")
    as_json = tmp_path / "queries.json"
    as_json.write_text(json.dumps({"chainA": "MKKD"}))
    assert load_queries(fasta) == load_queries(as_json)


@pytest.mark.parametrize(
    "contents, message",
    [
        (">chainA\nMKKD\n>chainA\nQQQ\n", "duplicate"),
        (">chainA\n", "empty sequence"),
        ("", "no queries"),
    ],
)
def test_malformed_queries_are_rejected(tmp_path, contents, message):
    path = tmp_path / "queries.fasta"
    path.write_text(contents)
    with pytest.raises(ValueError, match=message):
        load_queries(path)


def test_parse_fasta_ignores_comments_and_blank_lines():
    assert parse_fasta("# a comment\n\n>a\nMK\n") == [("a", "MK")]


def test_safe_filename_keeps_names_usable_on_disk():
    assert safe_filename("sp|P12345|A B/C") == "sp_P12345_A_B_C"
