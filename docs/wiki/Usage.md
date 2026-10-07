# Usage

```bash
nanotimesort -f /path/to/fastq_pass/ -o /path/to/output/ -i 1h -p my_sample
```

## Command-line options

| Flag | Required | Description |
|---|---|---|
| `-f, --fastq` | yes | Input folder (searched recursively) or a single FASTQ file. Accepts `.fastq`, `.fq`, `.fastq.gz` and `.fq.gz`. |
| `-o, --output` | yes | Output folder. Created if it does not exist. |
| `-i, --interval` | yes | Time interval for the bins, e.g. `1h`, `30m` or `90s`. Fractional values like `0.5h` work too. Bins are cumulative. |
| `-p, --prefix` | no | Output file prefix. Default: `interval`. |
| `-t, --threads` | no | Number of FASTQ files to process in parallel. Default: all CPUs. |
| `-c, --compression-level` | no | Gzip level for output files, 1 (fastest) to 9 (smallest). Default: 4. |
| `-v, --version` | no | Show version and exit. |

## Output

Files are named `{prefix}_0-{end}{unit}_{reads}reads_{bp}bp.fastq.gz`, so read and base-pair
counts are visible straight from a directory listing:

```text
my_sample_0-1h_125437reads_1103093674bp.fastq.gz
my_sample_0-2h_225891reads_2087456221bp.fastq.gz
my_sample_0-3h_301255reads_2812345678bp.fastq.gz
```

Bins are **cumulative**: the `0-2h` file also contains the reads of the `0-1h` file. The last
file holds the complete run. Elapsed time is measured from the earliest read start time found
across all input files.

## Supported FASTQ headers

| Basecaller | Header style | Example time field |
|---|---|---|
| Albacore / Guppy | key=value | `start_time=2019-07-16T19:51:22Z` |
| MinKNOW live basecalling (incl. built-in Dorado) | key=value | `start_time=2025-01-13T10:45:28.681306+00:00` |
| Dorado standalone | SAM tags | `st:Z:2023-09-01T11:13:45.731+00:00` |

Any ISO 8601 / RFC 3339 timestamp flavor is parsed (with `Z` suffix or numeric offset, with or
without fractional seconds). Reads without a recognizable start time are skipped, with a single
warning giving the total count.

Basecalled to BAM with Dorado? Convert while keeping the tags:

```bash
samtools fastq -T '*' calls.bam | gzip > calls.fastq.gz
```

## Tips

- Parallelism is per *file*, so `-t` helps most with the many small files MinKNOW produces
  (4000 reads per file by default). A single huge FASTQ is processed on one core.
- Use `-c 1` when the bins are intermediate files you will delete anyway; use `-c 9` for
  long-term archiving.
- The deprecated `nanopore_reads_binner.py` wrapper still accepts the same flags, so old
  pipelines keep working unchanged.
