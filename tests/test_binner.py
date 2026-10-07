import gzip
import os

import pytest

from nanotimesort.binner import NanoTimeSort, find_fastq_files

BASE = "2024-05-01T{:02d}:{:02d}:00Z"


def make_read(read_id, minutes, style="guppy", seq="ACGT" * 10):
    """Build one FASTQ record with a start time `minutes` after 10:00 UTC."""
    hour, minute = divmod(600 + minutes, 60)
    stamp = BASE.format(hour, minute)
    if style == "guppy":
        header = "@{} runid=abc ch=1 start_time={}".format(read_id, stamp)
    else:  # dorado
        stamp = stamp.replace("Z", ".000+00:00")
        header = "@{} qs:f:20.0 ch:i:1 st:Z:{}".format(read_id, stamp)
    return "{}\n{}\n+\n{}\n".format(header, seq, "I" * len(seq))


@pytest.fixture
def run_folder(tmp_path):
    """Two FASTQ files (one gzipped, one plain) spanning ~2.5 hours."""
    fastq_dir = tmp_path / "fastq"
    fastq_dir.mkdir()
    # File 1 (gzipped, Guppy headers): reads at 0, 30 and 70 minutes.
    content1 = (
        make_read("read1", 0)
        + make_read("read2", 30)
        + make_read("read3", 70)
    )
    with gzip.open(fastq_dir / "part1.fastq.gz", "wt") as handle:
        handle.write(content1)
    # File 2 (plain, Dorado headers): reads at 90 and 150 minutes.
    content2 = make_read("read4", 90, style="dorado") + make_read("read5", 150, style="dorado")
    (fastq_dir / "part2.fastq").write_text(content2)
    return fastq_dir


def read_ids(path):
    with gzip.open(path, "rt") as handle:
        return [line.split()[0][1:] for i, line in enumerate(handle) if i % 4 == 0]


def test_find_fastq_files(run_folder):
    files = find_fastq_files(str(run_folder))
    assert [os.path.basename(f) for f in files] == ["part1.fastq.gz", "part2.fastq"]


def test_cumulative_binning(run_folder, tmp_path):
    out_dir = tmp_path / "out"
    binner = NanoTimeSort(
        input_path=str(run_folder),
        output_folder=str(out_dir),
        interval="1h",
        prefix="test",
        threads=1,
    )
    outputs = binner.run()

    # Run spans 150 minutes -> three cumulative 1 h bins.
    names = [os.path.basename(p) for p in outputs]
    bp = 40  # each read is 40 bp
    assert names == [
        "test_0-1h_2reads_{}bp.fastq.gz".format(2 * bp),
        "test_0-2h_4reads_{}bp.fastq.gz".format(4 * bp),
        "test_0-3h_5reads_{}bp.fastq.gz".format(5 * bp),
    ]

    # Bin contents are cumulative and mix Guppy + Dorado reads.
    assert sorted(read_ids(outputs[0])) == ["read1", "read2"]
    assert sorted(read_ids(outputs[1])) == ["read1", "read2", "read3", "read4"]
    assert sorted(read_ids(outputs[2])) == ["read1", "read2", "read3", "read4", "read5"]


def test_parallel_matches_serial(run_folder, tmp_path):
    results = {}
    for threads in (1, 2):
        out_dir = tmp_path / "out_t{}".format(threads)
        binner = NanoTimeSort(
            input_path=str(run_folder),
            output_folder=str(out_dir),
            interval="30m",
            prefix="par",
            threads=threads,
        )
        outputs = binner.run()
        results[threads] = {
            os.path.basename(p): sorted(read_ids(p)) for p in outputs
        }
    assert results[1] == results[2]


def test_single_file_input(run_folder, tmp_path):
    out_dir = tmp_path / "out_single"
    binner = NanoTimeSort(
        input_path=str(run_folder / "part2.fastq"),
        output_folder=str(out_dir),
        interval="2h",
        prefix="single",
        threads=1,
    )
    outputs = binner.run()
    assert len(outputs) == 1
    assert sorted(read_ids(outputs[0])) == ["read4", "read5"]


def test_invalid_interval():
    with pytest.raises(ValueError):
        NanoTimeSort("in", "out", interval="1x")
    with pytest.raises(ValueError):
        NanoTimeSort("in", "out", interval="-5m")


def test_reads_without_timestamp_are_skipped(tmp_path, capsys):
    fastq_dir = tmp_path / "fastq"
    fastq_dir.mkdir()
    content = make_read("good1", 0) + "@orphan length=4\nACGT\n+\nIIII\n" + make_read("good2", 10)
    (fastq_dir / "mixed.fastq").write_text(content)
    out_dir = tmp_path / "out"
    binner = NanoTimeSort(str(fastq_dir), str(out_dir), interval="1h", threads=1)
    outputs = binner.run()
    assert sorted(read_ids(outputs[-1])) == ["good1", "good2"]


def test_handle_eviction_with_many_bins(run_folder, tmp_path, monkeypatch):
    """With more bins than allowed open handles, chunks are reopened in
    append mode and no read is lost."""
    import nanotimesort.binner as binner_module

    monkeypatch.setattr(binner_module, "MAX_OPEN_CHUNKS", 1)
    out_dir = tmp_path / "out_evict"
    binner = NanoTimeSort(
        input_path=str(run_folder),
        output_folder=str(out_dir),
        interval="10m",  # 150 min run -> 16 bins, far above the cap of 1
        prefix="evict",
        threads=1,
    )
    outputs = binner.run()
    assert len(outputs) == 16
    assert sorted(read_ids(outputs[-1])) == ["read1", "read2", "read3", "read4", "read5"]
    # Every cumulative file must contain its predecessor's reads.
    previous = set()
    for path in outputs:
        current = set(read_ids(path))
        assert previous <= current
        previous = current
