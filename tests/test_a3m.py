"""The taxonomy rewrite has to produce headers AlphaFold 3 can actually parse."""

import re

import pytest

from fastmsa.a3m import add_taxonomy_ids, iter_a3m, read_taxonomy_m8, taxonomy_header

# AlphaFold 3's own pattern for pulling a species ID out of an A3M header,
# copied verbatim from its MSA parser.
UNIPROT_ENTRY_NAME_REGEX = re.compile(
    r"(?:cb|tr|sp)\|"
    r"(?:[A-Z0-9]{6,10})"
    r"(?:_\d+)?\|"
    r"(?:[A-Z0-9]{1,10}_)"
    r"(?P<SpeciesId>[A-Z0-9]{1,5})"
)

# An MMseqs2 alignment as written by `fastmsa search`: the query first, then
# hits headed by a bare accession followed by tab-separated alignment stats.
SAMPLE_A3M = """\
>query
MKKDVRILLVGEPRVGKTSLIMSLVSEEFPEEVPPR
>Q6NVC5\t917\t0.859\t4.922E-295\t0\t617\t618\t0\t618\t619
MRKDVRILLVGEPKVGKTSLIMSLVSEEFPDEVPPR
>A0A6P6JVL9\t910\t0.857\t1.939E-292\t0\t617\t618\t0\t618\t619
MRKDVRILLVGEPKVGKTSLIMSLVSEEFPDEVPlR
>UniRef90_W5N8K8\t891\t0.813\t5.096E-286\t0\t617\t618\t0\t659\t660
MRKDVRILLVGEPKVGKTSLIMSLVSEEFPDEVPLR
>P00000\t880\t0.800\t1.000E-280\t0\t617\t618\t0\t659\t660
MRKDVRILLVGEPKVGKTSLIMSLVSEEFPDEVPLR
"""

TAXONOMY = {
    "Q6NVC5": "9606",  # 6-character accession, 4-digit taxid
    "A0A6P6JVL9": "10090",  # 10-character accession, 5-digit taxid
    "UniRef90_W5N8K8": "9031",  # cluster accession, prefix must be stripped
    # P00000 is deliberately absent: not every hit has a known taxid.
}


@pytest.fixture
def annotated(tmp_path):
    a3m = tmp_path / "query.a3m"
    a3m.write_text(SAMPLE_A3M)
    output = tmp_path / "annotated.a3m"
    count = add_taxonomy_ids(a3m, TAXONOMY, output)
    return count, list(iter_a3m(output))


def test_annotated_headers_parse_as_alphafold3_species_ids(annotated):
    count, entries = annotated
    assert count == 3

    species = {}
    for header, _sequence in entries:
        match = UNIPROT_ENTRY_NAME_REGEX.search(header)
        if match:
            species[header] = match.group("SpeciesId")

    assert sorted(species.values()) == ["10090", "9031", "9606"]


def test_hits_without_a_taxid_are_left_alone(annotated):
    _count, entries = annotated
    headers = [header for header, _ in entries]

    assert headers[0] == "query"
    assert headers[-1].startswith("P00000\t")
    assert UNIPROT_ENTRY_NAME_REGEX.search(headers[-1]) is None


def test_insertion_columns_survive_the_rewrite(annotated):
    _count, entries = annotated
    assert entries[2][1].endswith("EVPlR")


def test_uniref_prefix_is_stripped():
    assert taxonomy_header("UniRef90_W5N8K8", "9031") == "cb|W5N8K8|W5N8K8_9031"
    assert taxonomy_header("Q6NVC5", "9606") == "cb|Q6NVC5|Q6NVC5_9606"


def test_read_taxonomy_m8_keys_on_the_target_accession(tmp_path):
    m8 = tmp_path / "taxonomy.m8"
    m8.write_text(
        "query\tQ6NVC5\t9606\tHomo sapiens\tEukaryota;Metazoa\t0.859\n"
        "query\tP00000\t\tunclassified\t\t0.800\n"
    )
    assert read_taxonomy_m8(m8) == {"Q6NVC5": "9606"}
