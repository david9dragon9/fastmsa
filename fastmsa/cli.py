"""Command line interface: ``fastmsa <search|taxid|merge|json>``."""

from __future__ import annotations

import argparse
import json
import logging
import math
from pathlib import Path

from fastmsa import __version__
from fastmsa.a3m import TAXONOMY_M8, add_taxonomy_ids, read_taxonomy_m8
from fastmsa.af3 import af3_input
from fastmsa.fasta import load_queries, safe_filename
from fastmsa.merge import CHAIN_POLY_TYPES, merge_a3ms
from fastmsa.search import SearchConfig, search

logger = logging.getLogger(__name__)


def _add_search_parser(subparsers: argparse._SubParsersAction) -> None:
    parser = subparsers.add_parser(
        "search",
        help="search queries against one MMseqs2 database, one A3M per query",
        description=(
            "Search every query sequence against a single MMseqs2 database and "
            "write one A3M alignment per query into the output directory."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("queries", type=Path, help="FASTA file, or JSON mapping name to sequence.")
    parser.add_argument("db_dir", type=Path, help="Directory holding the MMseqs2 databases.")
    parser.add_argument("output_dir", type=Path, help="Directory for the A3M files and intermediates.")
    parser.add_argument("--db", required=True, help="Name of the database to search within DB_DIR.")
    parser.add_argument(
        "--db-type",
        choices=["flat", "clustered"],
        default="flat",
        help="'flat' for a database built from a FASTA file; 'clustered' for one "
        "with _seq/_aln sub-databases, whose cluster members are expanded.",
    )
    parser.add_argument(
        "--taxonomy",
        action="store_true",
        help=f"Also write {TAXONOMY_M8}, the taxonomy report needed by 'fastmsa taxid'. "
        "Requires a database with taxonomy information.",
    )
    parser.add_argument("--threads", type=int, default=32, help="Number of threads to use.")
    parser.add_argument(
        "--gpu", action="store_true", help="Search on GPU (select devices with CUDA_VISIBLE_DEVICES)."
    )
    parser.add_argument("--gpu-server", action="store_true", help="Use a running MMseqs2 GPU server.")

    tuning = parser.add_argument_group("search tuning")
    tuning.add_argument(
        "--no-filter", dest="filter_msa", action="store_false",
        help="Keep every hit instead of filtering the alignment for diversity.",
    )
    tuning.add_argument(
        "-s", "--sensitivity", type=float, default=None,
        help="MMseqs2 sensitivity. Lowering it is much faster but yields sparser "
        "alignments. Left unset, the k-mer thresholds are pinned to the ones the "
        "ColabFold server uses, equivalent to a sensitivity of ~8. Ignored with --gpu.",
    )
    tuning.add_argument(
        "--prefilter-mode", type=int, default=0, choices=[0, 1, 2],
        help="Prefilter: 0 k-mer (memory-hungry), 1 ungapped (CPU-hungry), 2 exhaustive.",
    )
    tuning.add_argument("--align-eval", type=float, default=10.0, help="E-value threshold for 'align'.")
    tuning.add_argument(
        "--expand-eval", type=float, default=math.inf,
        help="E-value threshold for 'expandaln'. Clustered databases only.",
    )
    tuning.add_argument(
        "--diff", type=int, default=3000, help="Keep at least this many sequences per alignment block."
    )
    tuning.add_argument(
        "--qsc", type=float, default=None,
        help="Minimum alignment score for 'filterresult'. Default: 0.8 with --filter, -20.0 without.",
    )
    tuning.add_argument(
        "--max-accept", type=int, default=None,
        help="Alignments accepted per query before 'align' stops. "
        "Default: 100000 with --filter, 1000000 without.",
    )
    tuning.add_argument(
        "--db-load-mode", type=int, default=0, choices=[0, 1, 2, 3],
        help="Database preload mode: 0 auto, 1 fread, 2 mmap, 3 mmap+touch.",
    )
    tuning.add_argument("--mmseqs", type=Path, default=Path("mmseqs"), help="Path to the mmseqs binary.")


def _add_taxid_parser(subparsers: argparse._SubParsersAction) -> None:
    parser = subparsers.add_parser(
        "taxid",
        help="rewrite A3M headers with NCBI taxonomy IDs for AlphaFold 3 pairing",
        description=(
            "Rewrite the hit headers of each A3M into the cb|ACCESSION|ACCESSION_TAXID "
            "form AlphaFold 3 parses species IDs from, so the alignment can be used as "
            "a paired MSA. Hits with no known taxonomy ID pass through unchanged."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "msa_dir", type=Path, help=f"Search output directory containing *.a3m and {TAXONOMY_M8}."
    )
    parser.add_argument("output_dir", type=Path, help="Directory for the rewritten A3M files.")
    parser.add_argument(
        "--taxonomy-m8", type=Path, default=None,
        help=f"Taxonomy report to read. Default: {TAXONOMY_M8} inside MSA_DIR.",
    )


def _add_merge_parser(subparsers: argparse._SubParsersAction) -> None:
    parser = subparsers.add_parser(
        "merge",
        help="merge the per-database alignments of each query into one MSA",
        description=(
            "Concatenate and de-duplicate the A3Ms a query has in each of several "
            "search output directories, producing one unpaired MSA per query. Uses "
            "AlphaFold 3's own MSA handling, so it must run where alphafold3 is importable."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("queries", type=Path, help="FASTA file, or JSON mapping name to sequence.")
    parser.add_argument("output_dir", type=Path, help="Directory for the merged A3M files.")
    parser.add_argument(
        "--msa-dirs", type=Path, nargs="+", required=True,
        help="Search output directories, each containing {query}.a3m.",
    )
    parser.add_argument(
        "--chain-poly-type", default="protein", choices=sorted(CHAIN_POLY_TYPES),
        help="Polymer type of the queries.",
    )
    parser.add_argument(
        "--no-deduplicate", dest="deduplicate", action="store_false",
        help="Keep sequences that appear in more than one input alignment.",
    )


def _add_json_parser(subparsers: argparse._SubParsersAction) -> None:
    parser = subparsers.add_parser(
        "json",
        help="write AlphaFold 3 input JSON referencing the computed MSAs",
        description=(
            "Write one AlphaFold 3 input JSON per query, wiring in the unpaired and "
            "paired MSAs so that AF3 skips its own data pipeline."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("queries", type=Path, help="FASTA file, or JSON mapping name to sequence.")
    parser.add_argument("output_dir", type=Path, help="Directory for the input JSON files.")
    parser.add_argument(
        "--unpaired-msa-dir", type=Path, default=None, help="Directory containing {query}.a3m."
    )
    parser.add_argument(
        "--paired-msa-dir", type=Path, default=None,
        help="Directory containing {query}.a3m with taxonomy-annotated headers.",
    )
    parser.add_argument(
        "--msa-as-path", action="store_true",
        help="Reference the A3M files by path instead of inlining them, keeping the JSON small.",
    )
    parser.add_argument(
        "--model-seeds", type=int, nargs="+", default=[1], help="Model seeds to fold with."
    )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="fastmsa",
        description="Bulk MSAs from the AlphaFold 3 database set, at MMseqs2 speed.",
    )
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    verbosity = parser.add_mutually_exclusive_group()
    verbosity.add_argument(
        "-v", "--verbose", action="store_true", help="Report progress for every query."
    )
    verbosity.add_argument(
        "-q", "--quiet", action="store_true", help="Only report warnings and errors."
    )
    subparsers = parser.add_subparsers(dest="command", required=True)
    _add_search_parser(subparsers)
    _add_taxid_parser(subparsers)
    _add_merge_parser(subparsers)
    _add_json_parser(subparsers)
    return parser


def _run_search(args: argparse.Namespace) -> None:
    queries = load_queries(args.queries)
    config = SearchConfig(
        clustered=args.db_type == "clustered",
        taxonomy=args.taxonomy,
        filter_msa=args.filter_msa,
        align_eval=args.align_eval,
        expand_eval=args.expand_eval,
        diff=args.diff,
        qsc=args.qsc,
        max_accept=args.max_accept,
        prefilter_mode=args.prefilter_mode,
        sensitivity=args.sensitivity,
        db_load_mode=args.db_load_mode,
        threads=args.threads,
        gpu=args.gpu,
        gpu_server=args.gpu_server,
        mmseqs=args.mmseqs,
    )
    written = search(queries, args.db_dir, args.db, args.output_dir, config)
    logger.info("Wrote %d alignments to %s", len(written), args.output_dir)


def _run_taxid(args: argparse.Namespace) -> None:
    taxonomy_m8 = args.taxonomy_m8 or args.msa_dir.joinpath(TAXONOMY_M8)
    if not taxonomy_m8.is_file():
        raise SystemExit(
            f"No taxonomy report at {taxonomy_m8}. Re-run 'fastmsa search --taxonomy' "
            "against a database that carries taxonomy information."
        )
    taxonomy = read_taxonomy_m8(taxonomy_m8)
    logger.info("Read %d taxonomy assignments from %s", len(taxonomy), taxonomy_m8)

    args.output_dir.mkdir(parents=True, exist_ok=True)
    a3m_paths = sorted(args.msa_dir.glob("*.a3m"))
    if not a3m_paths:
        raise SystemExit(f"No A3M files found in {args.msa_dir}")
    for a3m_path in a3m_paths:
        annotated = add_taxonomy_ids(
            a3m_path, taxonomy, args.output_dir.joinpath(a3m_path.name)
        )
        logger.debug("%s: annotated %d hits", a3m_path.name, annotated)
    logger.info("Wrote %d alignments to %s", len(a3m_paths), args.output_dir)


def _run_merge(args: argparse.Namespace) -> None:
    queries = load_queries(args.queries)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    for name in queries:
        filename = f"{safe_filename(name)}.a3m"
        paths = [d.joinpath(filename) for d in args.msa_dirs if d.joinpath(filename).is_file()]
        if not paths:
            raise SystemExit(
                f"No alignment named {filename} in any of: "
                + ", ".join(str(d) for d in args.msa_dirs)
            )
        merged, depth = merge_a3ms(paths, args.chain_poly_type, deduplicate=args.deduplicate)
        args.output_dir.joinpath(filename).write_text(merged)
        logger.debug("%s: merged %d alignments into depth %d", name, len(paths), depth)
    logger.info("Wrote %d alignments to %s", len(queries), args.output_dir)


def _msa_for(msa_dir: Path | None, name: str) -> Path | None:
    if msa_dir is None:
        return None
    path = msa_dir.joinpath(f"{safe_filename(name)}.a3m")
    if not path.is_file():
        raise SystemExit(f"No alignment for query {name!r} at {path}")
    return path


def _run_json(args: argparse.Namespace) -> None:
    queries = load_queries(args.queries)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    for name, sequence in queries.items():
        unpaired = _msa_for(args.unpaired_msa_dir, name)
        paired = _msa_for(args.paired_msa_dir, name)
        if args.msa_as_path:
            content = af3_input(
                name, sequence,
                unpaired_msa_path=unpaired.resolve() if unpaired else None,
                paired_msa_path=paired.resolve() if paired else None,
                model_seeds=args.model_seeds,
            )
        else:
            content = af3_input(
                name, sequence,
                unpaired_msa=unpaired.read_text() if unpaired else None,
                paired_msa=paired.read_text() if paired else None,
                model_seeds=args.model_seeds,
            )
        destination = args.output_dir.joinpath(f"{safe_filename(name)}.json")
        destination.write_text(json.dumps(content, indent=4) + "\n")
    logger.info("Wrote %d input files to %s", len(queries), args.output_dir)


COMMANDS = {
    "search": _run_search,
    "taxid": _run_taxid,
    "merge": _run_merge,
    "json": _run_json,
}


def main(argv: list[str] | None = None) -> None:
    args = build_parser().parse_args(argv)
    if args.quiet:
        level = logging.WARNING
    elif args.verbose:
        level = logging.DEBUG
    else:
        level = logging.INFO
    logging.basicConfig(level=level, format="%(asctime)s %(levelname)s %(message)s")
    COMMANDS[args.command](args)


if __name__ == "__main__":
    main()
