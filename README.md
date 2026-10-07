# nanoTimeSort

[![CI](https://github.com/duceppemo/nanoTimeSort/actions/workflows/ci.yml/badge.svg)](https://github.com/duceppemo/nanoTimeSort/actions/workflows/ci.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Python 3.8+](https://img.shields.io/badge/python-3.8%2B-blue.svg)](https://www.python.org/downloads/)
[![No dependencies](https://img.shields.io/badge/dependencies-none-brightgreen.svg)](pyproject.toml)

Bin Oxford Nanopore reads by **cumulative sequencing time intervals**.

Every Nanopore read carries its acquisition start time in the FASTQ header. nanoTimeSort uses it
to split a sequencing run into cumulative time slices — e.g. with `-i 1h`, you get one FASTQ with
the reads of the first hour, one with the first two hours, and so on. This is useful to answer
questions like *"How long did I actually need to sequence to close this genome / detect this
pathogen?"* or to benchmark real-time analysis pipelines.

```text
fastq_pass/  ──►  sample_0-1h_125437reads_1103093674bp.fastq.gz
                  sample_0-2h_225891reads_2087456221bp.fastq.gz
                  sample_0-3h_301255reads_2812345678bp.fastq.gz
                  ...
```

## Features

- **Guppy/MinKNOW *and* Dorado compatible** — reads the legacy `start_time=<ISO 8601>` header
  field as well as the Dorado SAM-style `st:Z:<ISO 8601>` tag, in any timestamp flavor
  (`...22Z`, `...45.731+00:00`, `...28.681306+00:00`).
- **Fast** — each read is compressed exactly once; cumulative files are assembled by raw gzip
  member concatenation instead of recompressing. Roughly **30–100× faster** than v0.1
  (75 s → 0.8 s on a 20k-read test set).
- **Low memory** — reads are streamed, never held in RAM, so full PromethION runs are fine.
- **Parallel** — FASTQ files are scanned and binned across multiple processes.
- **Zero dependencies** — pure Python standard library.

## Installation

```bash
git clone https://github.com/duceppemo/nanoTimeSort.git
cd nanoTimeSort
pip install .
```

Requires Python ≥ 3.8. No other dependencies.

## Usage

```bash
nanotimesort -f /path/to/fastq_pass/ -o /path/to/output/ -i 1h -p my_sample
```

```text
usage: nanotimesort [-h] -f /basecalled/folder/ -o /output/folder/ -i 1h
                    [-p my_sample] [-t N] [-c 4] [-v]

options:
  -f, --fastq              Input folder (searched recursively) or single FASTQ file.
                           Accepts .fastq, .fq, .fastq.gz and .fq.gz.
  -o, --output             Output folder. Created if it does not exist.
  -i, --interval           Time interval for the bins, e.g. '1h', '30m' or '90s'.
                           Bins are cumulative.
  -p, --prefix             Output file prefix. Default: interval
  -t, --threads            Number of FASTQ files to process in parallel. Default: all CPUs
  -c, --compression-level  Gzip level for outputs (1=fastest, 9=smallest). Default: 4
  -v, --version            Show version and exit.
```

Output files are named `{prefix}_0-{end}{unit}_{reads}reads_{bp}bp.fastq.gz`, so read and
base-pair counts are visible at a glance.

> **Note** — bins are *cumulative*: the `0-2h` file also contains the reads from the `0-1h`
> file. The last file contains the complete run.

📖 **Full documentation, tutorial and worked example:** see the
[**wiki**](https://github.com/duceppemo/nanoTimeSort/wiki).

## Supported FASTQ headers

| Basecaller | Header style | Example time field |
|---|---|---|
| Guppy / Albacore | key=value | `start_time=2019-07-16T19:51:22Z` |
| MinKNOW (live basecalling, incl. Dorado) | key=value | `start_time=2025-01-13T10:45:28.681306+00:00` |
| Dorado (standalone) | SAM tags | `st:Z:2023-09-01T11:13:45.731+00:00` |

Reads without a recognizable start time are skipped with a warning.

## How it works (and why it's fast now)

1. **Scan** — all FASTQ files are streamed once (in parallel) to find the first and last read
   start time of the run.
2. **Chunk** — files are streamed a second time; each read is gzip-compressed *once* into the
   chunk of the single interval it belongs to.
3. **Assemble** — cumulative output *k* = output *k-1* + chunks of interval *k*, concatenated at
   the byte level. Concatenated gzip members form a valid gzip stream, so nothing is ever
   recompressed, and read/bp counts are tallied along the way.

The previous version recompressed every read into every cumulative bin at gzip level 9 and held
the entire run in memory — the slowness was CPU-bound compression, not disk I/O.

> Output files contain multiple gzip members. This is fully standard (`gzip`, `zcat`, `seqkit`,
> `minimap2`, `samtools`, etc. all handle it); it's the same trick used by bgzip.

## Development

```bash
pip install -e .[dev]
pytest        # run the test suite
ruff check nanotimesort tests
```

## License

[MIT](LICENSE) © Marc-Olivier Duceppe
