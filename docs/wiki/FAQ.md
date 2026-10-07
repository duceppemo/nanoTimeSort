# FAQ

## Which basecallers are supported?

All of them, as far as read start times go:

| Basecaller | Header style | Time field |
|---|---|---|
| Albacore / Guppy | key=value | `start_time=2019-07-16T19:51:22Z` |
| MinKNOW live basecalling (incl. built-in Dorado) | key=value | `start_time=2025-01-13T10:45:28.681306+00:00` |
| Dorado standalone (`dorado basecaller --emit-fastq`, `dorado demux`, BAM converted to FASTQ with `samtools fastq -T '*'`) | SAM tags | `st:Z:2023-09-01T11:13:45.731+00:00` |

nanoTimeSort looks for either a `start_time=` field or an `st:Z:` tag in each header and parses
any ISO 8601 / RFC 3339 timestamp (with `Z` suffix or numeric offset, with or without fractional
seconds). Naive timestamps are assumed to be UTC.

## My reads were basecalled by Dorado into BAM. How do I get FASTQ with the time tag?

Convert with samtools, keeping the tags in the header:

```bash
samtools fastq -T '*' calls.bam | gzip > calls.fastq.gz
```

(`-T st` keeps just the start-time tag.)

## Why are the bins cumulative?

The tool's main use case is "what would my analysis have looked like after N hours of
sequencing?" — each file is a self-contained snapshot of the run at that point, ready to feed to
an assembler or classifier without concatenating anything.

## What happens to reads without a start time?

They are skipped, and a single warning with the total count is printed at the end of the scan
step. If *no* read has a recognizable start time, the run aborts with an explanatory error.

## The output files contain multiple gzip "members". Is that OK?

Yes. Cumulative files are assembled by concatenating independently compressed gzip blocks, which
the gzip format explicitly allows (it's how `bgzip` works too). `zcat`, `seqkit`, `minimap2`,
`samtools`, BioPython, etc. all read them transparently. The only tools that can trip on it are
hand-rolled scripts that call `zlib.decompress()` once instead of using a proper gzip reader.

## How much memory does it need?

Almost none — reads are streamed and never held in memory. Peak RSS is a few tens of MB
regardless of run size. (v0.1 loaded the entire run into RAM; that is gone.)

## Why is it so much faster than v0.1?

v0.1 recompressed every read into every cumulative bin at gzip level 9 — with 10 bins, the first
hour's reads were compressed 10 times each. v1.0 compresses each read once and builds the
cumulative files by raw byte concatenation, parallelized across input files. On a 20k-read test
set this went from 75 s to 0.8 s. The job was never I/O bound; it was CPU-bound on gzip.

## Does the read order in the output matter?

Reads within an interval keep their input-file order, and intervals are concatenated in
chronological order, so files are roughly time-sorted. There is no strict per-read sort — no
downstream tool needs FASTQ sorted by time, and skipping the sort keeps memory flat.

## Can I still use `nanopore_reads_binner.py`?

Yes, it remains as a deprecated wrapper with the same flags, but prints a notice. Switch to the
`nanotimesort` command when convenient.
