<p align="center">
  <img src="https://raw.githubusercontent.com/duceppemo/nanoTimeSort/master/docs/images/logo.svg" alt="nanoTimeSort" width="520">
</p>

**nanoTimeSort** bins Oxford Nanopore reads by cumulative sequencing time intervals, using the
read start time embedded in every FASTQ header by MinKNOW, Guppy and Dorado.

## Pages

- [[Installation]] — requirements and install options
- [[Usage]] — CLI reference, output naming, supported FASTQ headers
- [[Tutorial]] — a complete worked example, from input data to binned output
- [[How it works|How-it-works]] — the three-stage pipeline and why v1.0 is ~100× faster
- [[FAQ]] — common questions (header compatibility, cumulative bins, multi-member gzip)

## What is it for?

A Nanopore run produces reads continuously over hours or days. Splitting the run into cumulative
time slices lets you re-analyze "what I would have had after 1 h, 2 h, 3 h…" without re-running
the sequencer. Typical uses:

- **Time-to-answer studies** — how much sequencing time is needed to detect a pathogen, type a
  strain, or reach a target depth?
- **Assembly saturation** — at what point do more reads stop improving the assembly?
- **Real-time pipeline benchmarking** — replay a run as realistic incremental datasets.
- **Run QC** — spot a pore-death or air-bubble event by comparing per-interval yields (the read
  and bp counts are in the file names).

## Quick start

```bash
pip install nanotimesort
nanotimesort -f fastq_pass/ -o binned/ -i 1h -p my_sample
```

See the [[Tutorial]] for a full walk-through.
