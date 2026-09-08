# fastmsa

**TL;DR: MMseqs2-speed MSAs using AlphaFold 3 databases, including paired MSAs for AF3 complex prediction.**

AlphaFold 3 uses `jackhmmer` and
`nhmmer` against UniRef90, MGnify, BFD and UniProt, which are slow and may cost hours of CPU time per
target.

MMseqs2 is much faster, and `colabfold_search` already applies it to monomer MSAs, but it uses ColabFold's *clustered* databases, and it does not produce alignments AlphaFold 3 can pair by species. 

`fastmsa` searches flat databases straight on the AF3 FASTA files, with no `_seq`/`_aln` cluster sub-databases, rewrites hit headers into the form AF3 parses species IDs from, and writes the
input JSON for AF3 folding.

## Install

```bash
pip install -e .
```

The only external requirement is [MMseqs2](https://github.com/soedinglab/MMseqs2)
on your `PATH` (or pass `--mmseqs /path/to/mmseqs`). The `merge` step additionally
needs the [`alphafold3`](https://github.com/google-deepmind/alphafold3) package
importable, since it reuses AF3's own MSA de-duplication; run that one step inside
the AF3 container or environment.

## The pipeline

| Command          | Does                                                                     |
| ---------------- | ------------------------------------------------------------------------ |
| `fastmsa search` | Searches every query against one MMseqs2 database; writes `{query}.a3m`. |
| `fastmsa taxid`  | Rewrites hit headers so AF3 can read species IDs, makes a paired MSA.   |
| `fastmsa merge`  | Merges and de-duplicates each query's per-database A3Ms into one MSA.    |
| `fastmsa json`   | Writes AF3 input JSON wiring in the unpaired and paired MSAs.            |

Run `fastmsa <command> --help` for the full flag list.

## Preparing databases

Any MMseqs2 database works; you select one per search with `--db`. To reproduce
AlphaFold 3's own database set, build one database per FASTA file in the AF3
bundle:

| AF3 FASTA                                   | Used for                    |
| ------------------------------------------- | --------------------------- |
| `uniref90_*.fa`                             | unpaired MSA                |
| `mgy_clusters_*.fa`                         | unpaired MSA                |
| `bfd-first_non_consensus_sequences.fasta`   | unpaired MSA                |
| `uniprot_all_*.fa`                          | paired MSA (has taxonomy)   |

```bash
mmseqs createdb uniref90_2022_05.fa                      unirefDB
mmseqs createdb mgy_clusters_2022_05.fa                  mgyDB
mmseqs createdb bfd-first_non_consensus_sequences.fasta  bfdDB
mmseqs createdb uniprot_all_2021_04.fa                   uniprotDB

# Taxonomy IDs, needed only for the paired MSA.
mmseqs createtaxdb uniprotDB tmp

# Optional but a large speedup for repeated searches: prebuilt indexes.
for db in unirefDB mgyDB bfdDB uniprotDB; do
  mmseqs createindex "$db" tmp --threads 24
done
```

## System requirements
- The full AF3 set with indexes needs a few TB of disk, and index-backed searching is fastest when the index fits in RAM (order of 1 TB for
the whole set). Without an index, `fastmsa` falls back to reading the database
from disk automatically which is slower per query, but it runs on far less memory.

## Quickstart

Queries can be a FASTA file or a JSON `{"name": "SEQUENCE"}` mapping.

```fasta
>chainA
MKKDVRILLVGEPRVGKTSLIMSLVSEEF...
>chainB
MSEQLTNAVKVIDVAKRNGVSVAAIAKEL...
```

**1. Unpaired MSA**: search the three unpaired databases, then merge.

```bash
for db in unirefDB mgyDB bfdDB; do
  fastmsa search queries.fasta /path/to/databases msas_$db --db $db --threads 24
done

fastmsa merge queries.fasta msas_unpaired \
    --msa-dirs msas_unirefDB msas_mgyDB msas_bfdDB
```

**2. Paired MSA**: search UniProt with taxonomy, then annotate the headers.

```bash
fastmsa search queries.fasta /path/to/databases msas_uniprot \
    --db uniprotDB --taxonomy --threads 24

fastmsa taxid msas_uniprot msas_paired
```

**3. AlphaFold 3 input**: one JSON per query, MSAs already filled in via filename or expanded MSA.

```bash
fastmsa json queries.fasta af3_input \
    --unpaired-msa-dir msas_unpaired \
    --paired-msa-dir msas_paired
```

**4. Fold**, with AF3's data pipeline switched off:

```bash
python run_alphafold.py --input_dir=af3_input --output_dir=predictions \
    --norun_data_pipeline
```

Use `--msa-as-path` in step 3 to reference the A3M files instead of inlining
them, which keeps the JSON small for deep alignments.

## Notes

- `--db-type flat` (the default) realigns hits
against the database itself. `--db-type clustered` expects `_seq`/`_aln`
sub-databases (ColabFold style).
- An interrupted run can be restarted with the same arguments.
- AF3 recovers a hit's species by matching its A3M
header against a UniProt entry-name pattern. MMseqs2 emits bare accessions, so `fastmsa taxid` uses `cb|ACCESSION|ACCESSION_TAXID` as header for pairing. Hits with no known taxonomy ID are left alone.
- `-s/--sensitivity` trades depth for time (left unset, the k-mer
thresholds match the ColabFold server, ≈ sensitivity 8). `--gpu` runs the search
on GPU. `--no-filter` keeps every hit instead of filtering for diversity.


## References

MMseqs2 pipeline is derived from
[ColabFold](https://github.com/sokrypton/ColabFold)'s `colabfold_search`

AlphaFold 3: [AlphaFold 3 repository](https://github.com/google-deepmind/alphafold3)
