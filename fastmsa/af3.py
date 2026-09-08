"""Building AlphaFold 3 input JSON around a precomputed MSA."""

from __future__ import annotations

from pathlib import Path
from typing import Sequence

#: AF3 input dialect. Version 2 is the first that accepts external MSAs.
DIALECT = "alphafold3"
VERSION = 2


def af3_input(
    name: str,
    sequence: str,
    *,
    unpaired_msa: str | None = None,
    paired_msa: str | None = None,
    unpaired_msa_path: str | Path | None = None,
    paired_msa_path: str | Path | None = None,
    model_seeds: Sequence[int] = (1,),
) -> dict:
    """Return an AF3 input record for a single protein chain.

    An MSA may be supplied inline or as a path, but not both. Omitting one
    entirely writes an empty string, which tells AF3 to fold without that MSA
    rather than to compute it itself.
    """
    if unpaired_msa is not None and unpaired_msa_path is not None:
        raise ValueError("Give the unpaired MSA either inline or as a path, not both")
    if paired_msa is not None and paired_msa_path is not None:
        raise ValueError("Give the paired MSA either inline or as a path, not both")

    protein: dict = {
        "id": "A",
        "sequence": sequence,
        "modifications": [],
        "templates": [],
    }
    if unpaired_msa_path is not None:
        protein["unpairedMsaPath"] = str(unpaired_msa_path)
    else:
        protein["unpairedMsa"] = unpaired_msa or ""
    if paired_msa_path is not None:
        protein["pairedMsaPath"] = str(paired_msa_path)
    else:
        protein["pairedMsa"] = paired_msa or ""

    return {
        "dialect": DIALECT,
        "version": VERSION,
        "name": name,
        "sequences": [{"protein": protein}],
        "modelSeeds": list(model_seeds),
        "bondedAtomPairs": None,
        "userCCD": None,
    }
