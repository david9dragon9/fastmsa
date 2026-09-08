"""Loading query sequences and deriving safe output filenames."""

from __future__ import annotations

import json
from pathlib import Path

FASTA_SUFFIXES = (".fasta", ".fa", ".faa", ".fas")
JSON_SUFFIXES = (".json",)

QUERY_SUFFIXES = FASTA_SUFFIXES + JSON_SUFFIXES


def parse_fasta(fasta_string: str) -> list[tuple[str, str]]:
    """Parse a FASTA string into ``(name, sequence)`` pairs.

    The name is the first whitespace-delimited token of the description line.
    Sequences spanning several lines are concatenated.
    """
    entries: list[tuple[str, list[str]]] = []
    for line in fasta_string.splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        if line.startswith(">"):
            description = line[1:].strip()
            name = description.split()[0] if description else ""
            entries.append((name, []))
        elif entries:
            entries[-1][1].append(line)
    return [(name, "".join(chunks)) for name, chunks in entries]


def load_queries(path: str | Path) -> dict[str, str]:
    """Load queries from a FASTA file or a JSON ``{name: sequence}`` mapping.

    Both formats are accepted everywhere queries are needed, so the same file
    can drive the search, the merge and the JSON-writing steps.
    """
    path = Path(path)
    if not path.is_file():
        raise FileNotFoundError(f"Query file not found: {path}")

    suffix = path.suffix.lower()
    if suffix in JSON_SUFFIXES:
        loaded = json.loads(path.read_text())
        if not isinstance(loaded, dict):
            raise ValueError(f"{path} must contain a JSON object mapping name to sequence")
        entries = [(str(name), str(sequence)) for name, sequence in loaded.items()]
    elif suffix in FASTA_SUFFIXES:
        entries = parse_fasta(path.read_text())
    else:
        raise ValueError(
            f"Unsupported query file format {path.suffix!r}; expected one of {', '.join(QUERY_SUFFIXES)}"
        )

    queries: dict[str, str] = {}
    for name, sequence in entries:
        if not name:
            raise ValueError(f"{path} contains an entry without a name")
        if not sequence:
            raise ValueError(f"{path} contains an empty sequence for query {name!r}")
        if name in queries:
            raise ValueError(f"{path} contains duplicate query name {name!r}")
        queries[name] = sequence.upper()

    if not queries:
        raise ValueError(f"{path} contains no queries")
    return queries


def safe_filename(name: str) -> str:
    """Replace characters that are awkward in filenames with underscores."""
    return "".join(c if c.isalnum() or c in ("_", ".", "-") else "_" for c in name)
