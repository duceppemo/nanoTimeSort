import gzip

import pytest

from nanotimesort import __version__
from nanotimesort.cli import main


@pytest.fixture
def fastq_folder(tmp_path):
    folder = tmp_path / "fastq"
    folder.mkdir()
    reads = "".join(
        "@read{} runid=abc start_time=2024-05-01T1{}:00:00Z\nACGT\n+\nIIII\n".format(i, i)
        for i in range(3)
    )
    (folder / "run.fastq").write_text(reads)
    return folder


def test_main_success(fastq_folder, tmp_path):
    out_dir = tmp_path / "out"
    rc = main(["-f", str(fastq_folder), "-o", str(out_dir), "-i", "1h", "-p", "cli"])
    assert rc == 0
    outputs = sorted(out_dir.glob("cli_0-*.fastq.gz"))
    assert len(outputs) == 3
    with gzip.open(outputs[-1], "rt") as handle:
        assert sum(1 for line in handle if line.startswith("@read")) >= 1


def test_main_no_fastq(tmp_path):
    empty = tmp_path / "empty"
    empty.mkdir()
    rc = main(["-f", str(empty), "-o", str(tmp_path / "out"), "-i", "1h"])
    assert rc == 1


def test_main_bad_interval(fastq_folder, tmp_path):
    rc = main(["-f", str(fastq_folder), "-o", str(tmp_path / "out"), "-i", "1x"])
    assert rc == 1


def test_main_no_timestamps(tmp_path):
    folder = tmp_path / "fastq"
    folder.mkdir()
    (folder / "run.fastq").write_text("@read1 length=4\nACGT\n+\nIIII\n")
    rc = main(["-f", str(folder), "-o", str(tmp_path / "out"), "-i", "1h"])
    assert rc == 1


def test_version(capsys):
    with pytest.raises(SystemExit) as exc:
        main(["--version"])
    assert exc.value.code == 0
    assert __version__ in capsys.readouterr().out


def test_main_bad_max_time(fastq_folder, tmp_path):
    rc = main(["-f", str(fastq_folder), "-o", str(tmp_path / "out"), "-i", "1h", "-m", "2x"])
    assert rc == 1


def test_main_zero_max_time(fastq_folder, tmp_path):
    rc = main(["-f", str(fastq_folder), "-o", str(tmp_path / "out"), "-i", "1h", "-m", "0h"])
    assert rc == 1


def test_main_prefix_with_path_separator(fastq_folder, tmp_path):
    rc = main(["-f", str(fastq_folder), "-o", str(tmp_path / "out"), "-i", "1h", "-p", "a/b"])
    assert rc == 1
