# Installation

## Requirements

- Python ≥ 3.8
- That's it — nanoTimeSort has **no third-party dependencies** (pure standard library).

## From GitHub (recommended)

```bash
pip install git+https://github.com/duceppemo/nanoTimeSort.git
```

or clone first:

```bash
git clone https://github.com/duceppemo/nanoTimeSort.git
cd nanoTimeSort
pip install .
```

Both install the `nanotimesort` command.

## In a conda environment

```bash
conda create -n nanotimesort python=3.12 pip
conda activate nanotimesort
pip install git+https://github.com/duceppemo/nanoTimeSort.git
```

## For development

```bash
git clone https://github.com/duceppemo/nanoTimeSort.git
cd nanoTimeSort
pip install -e .[dev]   # adds pytest and ruff
pytest                  # run the test suite
```

## Verify the installation

```bash
nanotimesort --version
# nanoTimeSort v1.0.0
```

## Upgrading from v0.1 (`nanopore_reads_binner.py`)

The old script is deprecated but still works as a thin wrapper around the new CLI, with the same
flags (`-f`, `-o`, `-i`, `-p`, `-t`). Just replace:

```bash
python nanopore_reads_binner.py -f fastq/ -o out/ -i 1h -p sample
```

with:

```bash
nanotimesort -f fastq/ -o out/ -i 1h -p sample
```

Output file naming is unchanged (`{prefix}_0-{end}{unit}_{reads}reads_{bp}bp.fastq.gz`).
The `numpy` and `python-dateutil` dependencies of v0.1 are no longer needed.
