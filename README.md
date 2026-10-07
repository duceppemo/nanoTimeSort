<p align="center">
  <img src="docs/images/logo.svg" alt="nanoTimeSort" width="520">
</p>

<p align="center">
  <a href="https://github.com/duceppemo/nanoTimeSort/actions/workflows/ci.yml"><img src="https://github.com/duceppemo/nanoTimeSort/actions/workflows/ci.yml/badge.svg" alt="CI"></a>
  <a href="LICENSE"><img src="https://img.shields.io/badge/License-MIT-yellow.svg" alt="License: MIT"></a>
  <a href="https://www.python.org/downloads/"><img src="https://img.shields.io/badge/python-3.8%2B-blue.svg" alt="Python 3.8+"></a>
  <a href="pyproject.toml"><img src="https://img.shields.io/badge/dependencies-none-brightgreen.svg" alt="No dependencies"></a>
</p>

Bin Oxford Nanopore reads by **cumulative sequencing time intervals**, using the read start time
that MinKNOW, Guppy and Dorado embed in every FASTQ header. Useful to answer *"how long did I
actually need to sequence?"* — for time-to-detection studies, assembly saturation curves, or
benchmarking real-time pipelines.

```text
fastq_pass/  ──►  sample_0-1h_125437reads_1103093674bp.fastq.gz
                  sample_0-2h_225891reads_2087456221bp.fastq.gz   (bins are cumulative)
                  sample_0-3h_301255reads_2812345678bp.fastq.gz
```

Compatible with **Guppy/MinKNOW** (`start_time=`) and **Dorado** (`st:Z:`) headers. Fast
(each read compressed once, ~100× faster than v0.1), low-memory (streaming), parallel, and
pure standard library.

## Install

```bash
pip install git+https://github.com/duceppemo/nanoTimeSort.git
```

Requires Python ≥ 3.8. No other dependencies.

## Usage

```bash
nanotimesort -f /path/to/fastq_pass/ -o /path/to/output/ -i 1h -p my_sample
```

📖 **Full documentation** — CLI reference, tutorial with example data, header compatibility,
design notes and FAQ — is in the [**wiki**](https://github.com/duceppemo/nanoTimeSort/wiki).

## License

[MIT](LICENSE) © Marc-Olivier Duceppe
