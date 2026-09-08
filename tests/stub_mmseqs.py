#!/usr/bin/env python3
"""A stand-in for the mmseqs binary, for exercising the search pipeline.

Records every invocation to ``$STUB_MMSEQS_LOG`` and fabricates just enough
output for the next module in the pipeline to be satisfied.
"""

import os
import sys
from pathlib import Path

MODULE_OUTPUT_POS = {
    "align": 4,
    "expandaln": 5,
    "filterresult": 4,
    "lndb": 2,
    "mvdb": 2,
    "result2msa": 4,
    "search": 3,
}


def touch_db(db: Path) -> None:
    db.parent.mkdir(parents=True, exist_ok=True)
    Path(f"{db}.dbtype").write_text("stub")


def main(argv: list[str]) -> int:
    module, params = argv[0], argv[1:]
    with open(os.environ["STUB_MMSEQS_LOG"], "a") as log:
        log.write("\t".join(argv) + "\n")

    if module == "createdb":
        query_db = Path(params[1])
        names = [
            line.split(">")[1].strip()
            for line in Path(params[0]).read_text().splitlines()
            if line.startswith(">")
        ]
        touch_db(query_db)
        touch_db(Path(f"{query_db}_h"))
        query_db.with_suffix(".lookup").write_text(
            "".join(f"{i}\t{name}\t{i}\n" for i, name in enumerate(names))
        )
    elif module == "search":
        touch_db(Path(params[1]))
        touch_db(Path(params[3]) / "latest" / "profile_1")
    elif module == "convertalis":
        Path(params[3]).write_text("query\tQ6NVC5\t9606\n")
    elif module == "unpackdb":
        packed, destination = Path(params[0]), Path(params[1])
        lookup = packed.parent / "qdb.lookup"
        for line in lookup.read_text().splitlines():
            key = line.split("\t")[0]
            (destination / f"{key}.a3m").write_text(f">{key}\nMKKD\n")
    elif module == "rmdb":
        Path(f"{params[0]}.dbtype").unlink(missing_ok=True)
    elif module in MODULE_OUTPUT_POS:
        touch_db(Path(params[MODULE_OUTPUT_POS[module] - 1]))
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
