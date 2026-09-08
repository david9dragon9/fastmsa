"""Reading A3M alignments and annotating them with NCBI taxonomy IDs.

AlphaFold 3 pairs chains by species, and it recovers the species of a hit by
matching its A3M header against a UniProt-style pattern. MMseqs2 emits bare
accessions instead, so a paired MSA has to be rewritten into the form AF3 can
parse. :func:`add_taxonomy_ids` does that rewrite using the taxonomy report
that :mod:`fastmsa.search` writes alongside the alignments.
"""

from __future__ import annotations

import re
from pathlib import Path
from typing import Iterator

#: Name of the taxonomy report written by a search with ``--taxonomy``.
TAXONOMY_M8 = "taxonomy.m8"

#: Columns requested from ``mmseqs convertalis`` for the taxonomy report.
#: ``read_taxonomy_m8`` relies on ``target`` and ``taxid`` being the second and
#: third columns.
TAXONOMY_COLUMNS = (
    "query,target,taxid,taxname,taxlineage,fident,alnlen,mismatch,"
    "gapopen,qstart,qend,tstart,tend,evalue,bits,cigar"
)

#: Strips the ``UniRef100_``/``UniRef90_``/... prefix off a cluster accession.
_UNIREF_PREFIX = re.compile(r"UniRef\d+_(?P<accession>.+)")


def iter_a3m(path: str | Path) -> Iterator[tuple[str, str]]:
    """Yield ``(header, sequence)`` pairs from an A3M file.

    The header is returned without its leading ``>``. Sequences split over
    several lines are concatenated; case is preserved, so lower-case insertion
    columns survive a round trip.
    """
    header: str | None = None
    chunks: list[str] = []
    with open(path) as handle:
        for line in handle:
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(chunks)
                header = line[1:].rstrip("\n")
                chunks = []
            elif header is not None:
                chunks.append(line.strip())
    if header is not None:
        yield header, "".join(chunks)


def read_taxonomy_m8(path: str | Path) -> dict[str, str]:
    """Map hit accession to NCBI taxonomy ID from an MMseqs2 taxonomy report."""
    taxonomy: dict[str, str] = {}
    with open(path) as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 3:
                continue
            _query, target, taxid = fields[0], fields[1], fields[2]
            if taxid:
                taxonomy[target] = taxid
    return taxonomy


def taxonomy_header(hit_id: str, taxid: str) -> str:
    """Build a header AlphaFold 3 can extract a species ID from.

    AF3 expects ``cb|ACCESSION|ACCESSION_SPECIESID``; any UniRef cluster prefix
    on the accession is dropped so that the accession matches AF3's pattern.
    """
    match = _UNIREF_PREFIX.match(hit_id)
    accession = match.group("accession") if match else hit_id
    return f"cb|{accession}|{accession}_{taxid}"


def add_taxonomy_ids(
    a3m_path: str | Path,
    taxonomy: dict[str, str],
    output_path: str | Path,
) -> int:
    """Rewrite an A3M so AF3 can read species IDs off the hit headers.

    Hits with no known taxonomy ID -- including the query itself -- keep their
    original header and are left for the unpaired MSA. Returns the number of
    annotated hits.
    """
    annotated = 0
    with open(output_path, "w") as out:
        for header, sequence in iter_a3m(a3m_path):
            hit_id = header.split("\t")[0]
            taxid = taxonomy.get(hit_id)
            if taxid is not None:
                header = taxonomy_header(hit_id, taxid)
                annotated += 1
            out.write(f">{header}\n{sequence}\n")
    return annotated
