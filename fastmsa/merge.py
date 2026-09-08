"""Merging the per-database alignments of one query into a single MSA."""

from __future__ import annotations

from pathlib import Path
from typing import Sequence

#: mmCIF polymer types accepted by AlphaFold 3, plus short aliases.
CHAIN_POLY_TYPES = {
    "protein": "polypeptide(L)",
    "rna": "polyribonucleotide",
    "dna": "polydeoxyribonucleotide",
    "polypeptide(L)": "polypeptide(L)",
    "polyribonucleotide": "polyribonucleotide",
    "polydeoxyribonucleotide": "polydeoxyribonucleotide",
}


def resolve_chain_poly_type(name: str) -> str:
    """Map a polymer type name or alias onto its mmCIF spelling."""
    key = name.strip().lower()
    for candidate, resolved in CHAIN_POLY_TYPES.items():
        if key == candidate.lower():
            return resolved
    raise ValueError(
        f"Unknown chain polymer type {name!r}; choose from: {', '.join(CHAIN_POLY_TYPES)}"
    )


def _load_msa_class():
    """Import AlphaFold 3's MSA class, which is an optional dependency."""
    try:
        from alphafold3.data.msa import Msa
    except ImportError as error:  # pragma: no cover - depends on the environment
        raise ImportError(
            "Merging uses AlphaFold 3's own MSA de-duplication, so the alphafold3 "
            "package must be importable. Run this step inside the AlphaFold 3 "
            "container or environment."
        ) from error
    return Msa


def merge_a3ms(
    a3m_paths: Sequence[str | Path],
    chain_poly_type: str,
    *,
    deduplicate: bool = True,
) -> tuple[str, int]:
    """Concatenate several A3Ms for one query into a single alignment.

    All inputs must be alignments of the same query sequence. Returns the merged
    A3M text and its depth.
    """
    if not a3m_paths:
        raise ValueError("No A3M files to merge")

    msa_class = _load_msa_class()
    a3m_strings = [Path(path).read_text() for path in a3m_paths]
    msa = msa_class.from_multiple_a3ms(
        a3ms=a3m_strings,
        chain_poly_type=resolve_chain_poly_type(chain_poly_type),
        deduplicate=deduplicate,
    )
    return msa.to_a3m(), msa.depth
