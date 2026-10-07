# Tutorial

This tutorial walks through a complete run of nanoTimeSort. If you don't have Nanopore data at
hand, step 1 generates a realistic mock dataset so you can follow along on any machine.

## 1. Get some input data

### Option A — your own run

Point nanoTimeSort at the `fastq_pass/` folder of a basecalled run (MinKNOW live basecalling,
Guppy, or Dorado). The folder is searched recursively, and `.fastq`, `.fq`, `.fastq.gz` and
`.fq.gz` files are all accepted — barcoded subfolders included.

### Option B — generate a mock run

Save this as `make_mock_run.py` and run `python make_mock_run.py mock_run/`. It creates four
gzipped FASTQ files simulating a 5-hour run with 20,000 reads, mixing Guppy-style and
Dorado-style headers:

```python
import gzip, os, random, sys

out_dir = sys.argv[1]
os.makedirs(out_dir, exist_ok=True)
random.seed(42)
run_seconds = 5 * 3600  # 5 h run

for f in range(4):
    with gzip.open(f"{out_dir}/reads_{f}.fastq.gz", "wt") as fh:
        for r in range(5000):
            t = random.uniform(0, run_seconds)
            hh, rem = divmod(int(t), 3600)
            mm, ss = divmod(rem, 60)
            stamp = f"2024-05-01T{10 + hh:02d}:{mm:02d}:{ss:02d}"
            seq = "".join(random.choices("ACGT", k=1000))
            qual = "".join(random.choices("FGHI5:?J", k=1000))
            if f % 2 == 0:  # Guppy/MinKNOW style
                header = f"@read_{f}_{r} runid=abc ch={r % 512} start_time={stamp}Z"
            else:           # Dorado style
                header = f"@read_{f}_{r} qs:f:20.1 ch:i:{r % 512} st:Z:{stamp}.000+00:00"
            fh.write(f"{header}\n{seq}\n+\n{qual}\n")

print("Mock run written to", out_dir)
```

## 2. Bin the reads

Bin the run into cumulative 30-minute intervals:

```bash
nanotimesort -f mock_run/ -o binned/ -i 30m -p mock
```

Expected console output:

```text
Scanning 4 FASTQ file(s)... 20000 reads in 0.09s
Run spans 4h59m57s -> 10 intervals of 30m
Binning reads... done in 0.50s
Writing cumulative interval files... done in 0.12s
Total run time: 0.72s
```

## 3. Inspect the output

```bash
ls -1 binned/
```

```text
mock_0-30m_1932reads_1932000bp.fastq.gz
mock_0-60m_4002reads_4002000bp.fastq.gz
mock_0-90m_6031reads_6031000bp.fastq.gz
mock_0-120m_8114reads_8114000bp.fastq.gz
mock_0-150m_10082reads_10082000bp.fastq.gz
mock_0-180m_12091reads_12091000bp.fastq.gz
mock_0-210m_14110reads_14110000bp.fastq.gz
mock_0-240m_16090reads_16090000bp.fastq.gz
mock_0-270m_18011reads_18011000bp.fastq.gz
mock_0-300m_20000reads_20000000bp.fastq.gz
```

Three things to notice:

1. **Bins are cumulative** — `mock_0-60m_...` contains everything in `mock_0-30m_...` plus the
   reads from minutes 30–60. The last file is the complete run.
2. **Counts are in the names** — `..._4002reads_4002000bp...` means 4,002 reads totalling
   4.002 Mbp. You can plot yield-over-time straight from a directory listing.
3. **Elapsed time starts at the first read** — interval 0 begins at the earliest `start_time`
   found across all input files, not at a clock hour.

## 4. Use the bins

Each output is a standard gzipped FASTQ. For example, assemble each time slice to find the
point of diminishing returns:

```bash
for fq in binned/mock_0-*.fastq.gz; do
    name=$(basename "$fq" .fastq.gz)
    flye --nano-hq "$fq" -o "asm_${name}" -t 16
done
```

Or map each slice and track reference coverage over sequencing time:

```bash
for fq in binned/mock_0-*.fastq.gz; do
    minimap2 -ax map-ont ref.fasta "$fq" | samtools sort -o "${fq%.fastq.gz}.bam" -
    samtools coverage "${fq%.fastq.gz}.bam"
done
```

## 5. Options worth knowing

| Flag | What it does |
|---|---|
| `-t 8` | Process up to 8 FASTQ files in parallel. Defaults to all CPUs. Parallelism is per *file*, so it helps most with the many small files MinKNOW produces (4000 reads per file by default). |
| `-c 1` | Fastest output compression (bigger files). `-c 9` for smallest files. Default `4` is a good balance. |
| `-i 90s` | Intervals also accept minutes (`m`) and seconds (`s`), and fractional values like `0.5h`. |
| `-f run/reads.fastq.gz` | A single file works too, not just folders. |
