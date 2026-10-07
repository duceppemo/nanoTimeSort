# How it works

nanoTimeSort streams the data in three stages; reads are never held in memory.

```text
            Stage 1: SCAN                Stage 2: CHUNK               Stage 3: ASSEMBLE
   ┌──────────────────────────┐  ┌──────────────────────────────┐  ┌──────────────────────┐
   │ stream all FASTQ files   │  │ stream files again; compress │  │ output k =           │
   │ (parallel) to find the   │─►│ each read ONCE into the gzip │─►│   output k-1         │
   │ first & last start time  │  │ chunk of its interval        │  │   + chunks of bin k  │
   └──────────────────────────┘  └──────────────────────────────┘  │ (byte concatenation) │
                                                                   └──────────────────────┘
```

1. **Scan** — every FASTQ file is streamed once (in parallel across files) to find the earliest
   and latest read start time, which define the run's time range and the number of bins.
2. **Chunk** — files are streamed a second time. Each read's elapsed time places it in exactly
   one interval, and it is gzip-compressed exactly once into that interval's chunk file.
   Read and base-pair counts are tallied here, so outputs never need to be re-read.
3. **Assemble** — cumulative output *k* is built by concatenating output *k−1* and the chunks of
   interval *k* at the raw byte level. Concatenated gzip members form a valid gzip stream
   (the same property bgzip relies on), so nothing is ever decompressed or recompressed.

## Why v1.0 is ~100× faster than v0.1

v0.1's bottleneck was never disk I/O — it was CPU-bound gzip compression:

| | v0.1 | v1.0 |
|---|---|---|
| Compressions per read | once per cumulative bin (×10 for 10 bins) | exactly 1 |
| Compression level | 9 (hardcoded) | 4 (configurable 1–9) |
| Timestamp parsing | `dateutil.parser.parse` per read | `datetime.fromisoformat` |
| Memory | entire run in RAM | streaming, flat ~tens of MB |
| Parallelism | none | per input file (scan + chunk) |
| Read/bp counts | re-read every output file | tallied while writing |
| Dependencies | numpy, python-dateutil | none |

Benchmark (20,000 reads × 1 kb, 5 h run, 30 min bins, 4 files):

| Version | Wall time |
|---|---|
| v0.1 | 75 s |
| v1.0, 1 thread | 2.5 s |
| v1.0, 4 threads | **0.8 s** |

Output read sets per bin are identical between versions, and file naming is unchanged.

## Multi-member gzip outputs

Because cumulative files are assembled by concatenation, each output contains multiple gzip
members. This is fully standard: `zcat`, `seqkit`, `minimap2`, `samtools`, BioPython and every
proper gzip reader handle it transparently. Only hand-rolled code that calls
`zlib.decompress()` a single time instead of using a gzip reader would stop at the first member.
