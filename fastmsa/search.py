"""Bulk MMseqs2 search producing one A3M alignment per query sequence.

The same pipeline serves both shapes of MMseqs2 database:

* **flat** -- a database created straight from a FASTA file, which is what the
  AlphaFold 3 database set gives you. Results are realigned against the
  database itself.
* **clustered** -- a database that ships ``_seq``/``_aln`` sub-databases (the
  ColabFold-style databases). Cluster members are recovered with ``expandaln``
  before realignment, which deepens the alignment at no extra search cost.

Derived from ColabFold's ``colabfold_search`` (MIT licensed).
"""

from __future__ import annotations

import logging
import math
import os
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Sequence

from fastmsa.a3m import TAXONOMY_COLUMNS, TAXONOMY_M8
from fastmsa.fasta import safe_filename

logger = logging.getLogger(__name__)

#: Position of the output database in each MMseqs2 module's argument list, used
#: to skip work that a previous run already completed.
MODULE_OUTPUT_POS = {
    "align": 4,
    "convertalis": 4,
    "expandaln": 5,
    "filterresult": 4,
    "lndb": 2,
    "mvdb": 2,
    "result2msa": 4,
    "search": 3,
}

QUERY_FASTA = "query.fas"
QUERY_DB = "qdb"
UNPACK_DIR = "unpacked"


def dbtype_path(db: Path) -> Path:
    """Path of the ``.dbtype`` sidecar MMseqs2 writes next to every database."""
    return Path(f"{db}.dbtype")


def run_mmseqs(mmseqs: Path, params: Sequence[str | Path]) -> None:
    """Run one MMseqs2 module, skipping it if its output already exists."""
    module = str(params[0])
    if module in MODULE_OUTPUT_POS:
        output_path = dbtype_path(Path(params[MODULE_OUTPUT_POS[module]]))
        if output_path.exists():
            logger.info("Skipping %s because %s already exists", module, output_path)
            return

    logger.info("Running %s %s", mmseqs, " ".join(str(param) for param in params))
    # Hide the verbose MMseqs2 parameter list that otherwise clogs up the log.
    os.environ["MMSEQS_CALL_DEPTH"] = "1"
    subprocess.check_call([str(mmseqs)] + [str(param) for param in params])


def remove_db(mmseqs: Path, db: Path) -> None:
    """Delete an intermediate database, tolerating one that was never created."""
    if dbtype_path(db).exists():
        run_mmseqs(mmseqs, ["rmdb", db])


@dataclass
class SearchConfig:
    """Tunables for :func:`search`, mirroring the MMseqs2 flags they feed."""

    clustered: bool = False
    taxonomy: bool = False
    filter_msa: bool = True
    align_eval: float = 10.0
    expand_eval: float = math.inf
    diff: int = 3000
    qsc: float | None = None
    max_accept: int | None = None
    prefilter_mode: int = 0
    sensitivity: float | None = None
    db_load_mode: int = 0
    threads: int = 32
    gpu: bool = False
    gpu_server: bool = False
    mmseqs: Path = Path("mmseqs")

    def resolved_qsc(self) -> float:
        """Minimum score threshold, defaulting to the active filter mode."""
        if self.qsc is not None:
            return self.qsc
        return 0.8 if self.filter_msa else -20.0

    def resolved_max_accept(self) -> int:
        """Alignment cap, defaulting to the active filter mode."""
        if self.max_accept is not None:
            return self.max_accept
        return 100_000 if self.filter_msa else 1_000_000


@dataclass
class _DatabaseLayout:
    """Which database MMseqs2 should read at each stage of the pipeline."""

    search: Path
    align: Path
    cluster_alignments: Path | None
    db_load_mode: int


def _resolve_database(
    dbbase: Path, db_name: str, *, clustered: bool, db_load_mode: int
) -> _DatabaseLayout:
    """Locate a database and work out whether its index can be used."""
    if not dbbase.joinpath(f"{db_name}.dbtype").is_file():
        raise FileNotFoundError(f"Database {db_name} does not exist in {dbbase}")

    indexed = (
        dbbase.joinpath(f"{db_name}.idx").is_file()
        or dbbase.joinpath(f"{db_name}.idx.index").is_file()
    ) and not os.environ.get("MMSEQS_IGNORE_INDEX")
    if not indexed:
        logger.info("Database %s has no usable index; searching without one", db_name)
        db_load_mode = 0

    if not clustered:
        db = dbbase.joinpath(db_name)
        return _DatabaseLayout(db, db, None, db_load_mode)

    suffix = ".idx" if indexed else "_seq"
    return _DatabaseLayout(
        search=dbbase.joinpath(db_name),
        align=dbbase.joinpath(f"{db_name}{suffix}"),
        cluster_alignments=dbbase.joinpath(f"{db_name}{'.idx' if indexed else '_aln'}"),
        db_load_mode=db_load_mode,
    )


def _search_params(config: SearchConfig, db_load_mode: int) -> list[str]:
    params = [
        "--num-iterations", "3",
        "--db-load-mode", str(db_load_mode),
        "-a",
        "-e", "0.1",
        "--max-seqs", "10000",
    ]
    if config.gpu:
        # The GPU code path only supports the ungapped prefilter and always
        # runs at maximum sensitivity.
        params += ["--gpu", "1", "--prefilter-mode", "1"]
    else:
        params += ["--prefilter-mode", str(config.prefilter_mode)]
        if config.sensitivity is not None:
            params += ["-s", f"{config.sensitivity:.1f}"]
        else:
            params += ["--k-score", "seq:96,prof:80"]
    if config.gpu_server:
        params += ["--gpu-server", "1"]
    return params


def write_query_db(queries: dict[str, str], base: Path, mmseqs: Path) -> None:
    """Write the queries to a FASTA file and turn it into an MMseqs2 database."""
    base.mkdir(parents=True, exist_ok=True)
    query_fasta = base.joinpath(QUERY_FASTA)
    with query_fasta.open("w") as handle:
        for name, sequence in queries.items():
            handle.write(f">{name}\n{sequence}\n")

    run_mmseqs(
        mmseqs,
        ["createdb", query_fasta, base.joinpath(QUERY_DB), "--shuffle", "0", "--dbtype", "1"],
    )
    # Rewrite the lookup so intermediate databases carry the query names.
    with base.joinpath(f"{QUERY_DB}.lookup").open("w") as handle:
        for index, name in enumerate(queries):
            handle.write(f"{index}\t{name}\t{index}\n")
    query_fasta.unlink()


def search(
    queries: dict[str, str],
    dbbase: Path,
    db_name: str,
    output_dir: Path,
    config: SearchConfig | None = None,
) -> list[Path]:
    """Search every query against one database and write one A3M per query.

    Intermediate MMseqs2 databases live in ``output_dir`` and are removed on
    success, so an interrupted run can be restarted with the same arguments and
    will pick up where it stopped. Returns the written A3M paths, in query
    order.
    """
    config = config or SearchConfig()
    base = Path(output_dir)
    layout = _resolve_database(
        Path(dbbase),
        db_name,
        clustered=config.clustered,
        db_load_mode=config.db_load_mode,
    )
    mmseqs = config.mmseqs

    write_query_db(queries, base, mmseqs)

    query_db = base.joinpath(QUERY_DB)
    msa_db = base.joinpath("msa.a3m")
    threads = str(config.threads)
    db_load_mode = str(layout.db_load_mode)

    if not dbtype_path(msa_db).exists():
        run_mmseqs(
            mmseqs,
            ["search", query_db, layout.search, base.joinpath("res"), base.joinpath("tmp"),
             "--threads", threads]
            + _search_params(config, layout.db_load_mode),
        )
        # The last search iteration leaves behind the query profile; realigning
        # against it is what makes the alignment sensitive.
        run_mmseqs(mmseqs, ["mvdb", base.joinpath("tmp/latest/profile_1"), base.joinpath("prof_res")])
        run_mmseqs(mmseqs, ["lndb", base.joinpath(f"{QUERY_DB}_h"), base.joinpath("prof_res_h")])

        hits = base.joinpath("res")
        if layout.cluster_alignments is not None:
            run_mmseqs(
                mmseqs,
                ["expandaln", query_db, layout.align, hits, layout.cluster_alignments,
                 base.joinpath("res_exp"),
                 "--db-load-mode", db_load_mode,
                 "--threads", threads,
                 "--expansion-mode", "0",
                 "-e", str(config.expand_eval),
                 "--expand-filter-clusters", "1" if config.filter_msa else "0",
                 "--max-seq-id", "0.95"],
            )
            hits = base.joinpath("res_exp")

        run_mmseqs(
            mmseqs,
            ["align", base.joinpath("prof_res"), layout.align, hits, base.joinpath("res_realign"),
             "--db-load-mode", db_load_mode,
             "-e", str(config.align_eval),
             "--max-accept", str(config.resolved_max_accept()),
             "--threads", threads,
             "--alt-ali", "10",
             "-a"],
        )
        run_mmseqs(
            mmseqs,
            ["filterresult", query_db, layout.align, base.joinpath("res_realign"),
             base.joinpath("res_filter"),
             "--db-load-mode", db_load_mode,
             "--qid", "0",
             "--qsc", str(config.resolved_qsc()),
             "--diff", "0",
             "--threads", threads,
             "--max-seq-id", "1.0",
             "--filter-min-enable", "100"],
        )
        if config.taxonomy:
            run_mmseqs(
                mmseqs,
                ["convertalis", query_db, layout.align, base.joinpath("res_filter"),
                 base.joinpath(TAXONOMY_M8),
                 "--format-output", TAXONOMY_COLUMNS,
                 "--db-load-mode", db_load_mode,
                 "--threads", threads],
            )
        run_mmseqs(
            mmseqs,
            ["result2msa", query_db, layout.align, base.joinpath("res_filter"), msa_db,
             "--msa-format-mode", "6",
             "--db-load-mode", db_load_mode,
             "--threads", threads,
             "--filter-msa", "1" if config.filter_msa else "0",
             "--filter-min-enable", "1000",
             "--diff", str(config.diff),
             "--qid", "0.0,0.2,0.4,0.6,0.8,1.0",
             "--qsc", "0",
             "--max-seq-id", "0.95"],
        )
        for intermediate in ("res_filter", "res_realign", "res_exp", "res"):
            remove_db(mmseqs, base.joinpath(intermediate))
    else:
        logger.info("Skipping the %s search because %s already exists", db_name, msa_db)

    written = _unpack(queries, base, msa_db, mmseqs)

    for intermediate in ("prof_res", "prof_res_h", QUERY_DB, f"{QUERY_DB}_h"):
        remove_db(mmseqs, base.joinpath(intermediate))
    base.joinpath(f"{QUERY_DB}.lookup").unlink(missing_ok=True)
    shutil.rmtree(base.joinpath("tmp"), ignore_errors=True)
    return written


def _unpack(queries: dict[str, str], base: Path, msa_db: Path, mmseqs: Path) -> list[Path]:
    """Split the packed MSA database into one named A3M file per query."""
    unpack_dir = base.joinpath(UNPACK_DIR)
    unpack_dir.mkdir(parents=True, exist_ok=True)
    run_mmseqs(
        mmseqs,
        ["unpackdb", msa_db, unpack_dir, "--unpack-name-mode", "0", "--unpack-suffix", ".a3m"],
    )

    written = []
    for index, name in enumerate(queries):
        unpacked = unpack_dir.joinpath(f"{index}.a3m")
        if not unpacked.is_file():
            raise FileNotFoundError(
                f"No alignment came back for query {name!r}. If {base} holds results "
                "for a different set of queries, search into a fresh directory."
            )
        destination = base.joinpath(f"{safe_filename(name)}.a3m")
        unpacked.replace(destination)
        written.append(destination)

    shutil.rmtree(unpack_dir, ignore_errors=True)
    remove_db(mmseqs, msa_db)
    return written
