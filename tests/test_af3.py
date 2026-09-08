"""The generated input JSON has to be a valid AlphaFold 3 external-MSA record."""

import pytest

from fastmsa.af3 import DIALECT, VERSION, af3_input


def test_inline_msas_are_embedded():
    content = af3_input("q", "MKKD", unpaired_msa=">q\nMKKD\n", paired_msa=">q\nMKKD\n")
    assert content["dialect"] == DIALECT
    assert content["version"] == VERSION
    protein = content["sequences"][0]["protein"]
    assert protein["sequence"] == "MKKD"
    assert protein["unpairedMsa"] == ">q\nMKKD\n"
    assert protein["templates"] == []


def test_missing_msas_are_empty_rather_than_null():
    # An empty string tells AF3 to fold without the MSA; null would make it run
    # its own data pipeline, which is the thing we are avoiding.
    protein = af3_input("q", "MKKD")["sequences"][0]["protein"]
    assert protein["unpairedMsa"] == ""
    assert protein["pairedMsa"] == ""


def test_paths_replace_the_inline_fields():
    protein = af3_input(
        "q", "MKKD", unpaired_msa_path="/msas/q.a3m", paired_msa_path="/paired/q.a3m"
    )["sequences"][0]["protein"]
    assert protein["unpairedMsaPath"] == "/msas/q.a3m"
    assert protein["pairedMsaPath"] == "/paired/q.a3m"
    assert "unpairedMsa" not in protein and "pairedMsa" not in protein


def test_an_msa_cannot_be_both_inline_and_a_path():
    with pytest.raises(ValueError, match="not both"):
        af3_input("q", "MKKD", unpaired_msa="x", unpaired_msa_path="/msas/q.a3m")
