"""Core binning engine.

The pipeline runs in three stages, none of which holds reads in memory:

1. **Scan** (parallel): stream every FASTQ once to find the earliest and
   latest read start time and the total read count.
2. **Chunk** (parallel): stream every FASTQ again; each read is gzip-
   compressed exactly once, into the chunk file of the single interval it
   belongs to.
3. **Assemble** (sequential): cumulative output *k* is built by raw byte
   concatenation of output *k-1* and the chunks of interval *k*. Gzip members
   concatenate into a valid gzip stream, so no data is ever recompressed.

Compared to the original implementation (which recompressed every read into
every cumulative bin at gzip level 9), this is typically one to two orders of
magnitude faster.
"""

from __future__ import annotations

import gzip
import os
import shutil
import sys
import tempfile
import zlib
from collections import OrderedDict
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass
from datetime import datetime
from math import ceil, floor
from time import time
from typing import IO, Iterator, List, Optional, Tuple

from .timestamps import extract_start_time, find_time_field, parse_timestamp

UNIT_SECONDS = {"h": 3600.0, "m": 60.0, "s": 1.0}
UNIT_NAMES = {"h": "hour", "m": "minute", "s": "second"}

_COPY_BUFFER = 4 * 1024 * 1024

# Cap on simultaneously open chunk files per worker. Small intervals over a
# long run can produce tens of thousands of bins, which would exhaust the
# process file-descriptor limit if every chunk stayed open. Evicted chunks
# are reopened in append mode, which adds a new gzip member -- exactly what
# the assembly stage concatenates anyway.
MAX_OPEN_CHUNKS = 128


@dataclass
class ScanResult:
    """Summary of one scan pass over a FASTQ file."""

    reads: int = 0
    missing: int = 0  # reads with no start-time field in the header
    malformed: int = 0  # reads whose start-time value failed to parse
    t_min: Optional[datetime] = None
    t_max: Optional[datetime] = None


def open_fastq(path: str) -> IO[bytes]:
    """Open a plain or gzipped FASTQ for binary reading."""
    if path.endswith(".gz"):
        return gzip.open(path, "rb")
    return open(path, "rb", buffering=1024 * 1024)


def fastq_records(handle: IO[bytes]) -> Iterator[Tuple[bytes, bytes, bytes]]:
    """Yield (header, sequence, quality) byte lines, newline-stripped."""
    while True:
        header = handle.readline()
        if not header:
            return
        seq = handle.readline()
        handle.readline()  # '+' separator line
        qual = handle.readline()
        if not qual:
            return  # truncated final record: skip it
        yield header.rstrip(), seq.rstrip(), qual.rstrip()


def find_fastq_files(input_path: str) -> List[str]:
    """Collect FASTQ files from a folder (recursively) or a single file path."""
    extensions = (".fastq", ".fastq.gz", ".fq", ".fq.gz")
    if os.path.isfile(input_path):
        return [input_path] if input_path.endswith(extensions) else []
    fastq_list = []
    for root, _dirs, filenames in os.walk(input_path):
        for filename in filenames:
            if filename.endswith(extensions):
                fastq_list.append(os.path.join(root, filename))
    return sorted(fastq_list)


def scan_file(path: str) -> ScanResult:
    """Pass 1: stream one file and record read count and time extremes."""
    result = ScanResult()
    if os.path.getsize(path) == 0:
        return result
    try:
        with open_fastq(path) as handle:
            for header, _seq, _qual in fastq_records(handle):
                raw = find_time_field(header)
                if raw is None:
                    result.missing += 1
                    continue
                try:
                    start_time = parse_timestamp(raw)
                except (ValueError, UnicodeDecodeError):
                    result.malformed += 1
                    continue
                result.reads += 1
                if result.t_min is None or start_time < result.t_min:
                    result.t_min = start_time
                if result.t_max is None or start_time > result.t_max:
                    result.t_max = start_time
    except (EOFError, zlib.error, gzip.BadGzipFile) as err:
        raise ValueError("Corrupt or truncated input file {}: {}".format(path, err)) from err
    return result


def chunk_file(
    path: str,
    file_index: int,
    t_min: datetime,
    bin_seconds: float,
    num_bins: int,
    cutoff_seconds: Optional[float],
    chunk_dir: str,
    compresslevel: int,
) -> Tuple[List[int], List[int]]:
    """Pass 2: write each read of one file into its interval chunk.

    Chunk files are opened lazily, so intervals with no reads in this file
    cost nothing. Reads at or past ``cutoff_seconds`` (when set) are
    discarded. Returns per-interval read and base-pair counts.
    """
    reads_per_bin = [0] * num_bins
    bp_per_bin = [0] * num_bins
    handles: "OrderedDict[int, IO[bytes]]" = OrderedDict()

    if os.path.getsize(path) == 0:
        return reads_per_bin, bp_per_bin

    try:
        with open_fastq(path) as handle:
            for header, seq, qual in fastq_records(handle):
                start_time = extract_start_time(header)
                if start_time is None:
                    continue
                elapsed = (start_time - t_min).total_seconds()
                if cutoff_seconds is not None and elapsed >= cutoff_seconds:
                    continue
                bin_index = min(floor(elapsed / bin_seconds), num_bins - 1)
                out = handles.get(bin_index)
                if out is None:
                    chunk_path = os.path.join(
                        chunk_dir, "chunk_b{:06d}_f{:06d}.fastq.gz".format(bin_index, file_index)
                    )
                    out = gzip.open(chunk_path, "ab", compresslevel=compresslevel)
                    handles[bin_index] = out
                    if len(handles) > MAX_OPEN_CHUNKS:
                        _old_bin, old_handle = handles.popitem(last=False)
                        old_handle.close()
                else:
                    handles.move_to_end(bin_index)
                out.write(header + b"\n" + seq + b"\n+\n" + qual + b"\n")
                reads_per_bin[bin_index] += 1
                bp_per_bin[bin_index] += len(seq)
    except (EOFError, zlib.error, gzip.BadGzipFile) as err:
        raise ValueError("Corrupt or truncated input file {}: {}".format(path, err)) from err
    finally:
        for out in handles.values():
            out.close()

    return reads_per_bin, bp_per_bin


def format_number(value: float) -> str:
    """Render 2.0 as '2' and 0.5 as '0.5' for use in file names."""
    return str(int(value)) if float(value).is_integer() else str(value)


class NanoTimeSort:
    """Bin Nanopore reads into cumulative sequencing-time interval files."""

    def __init__(
        self,
        input_path: str,
        output_folder: str,
        interval: str,
        prefix: str = "interval",
        threads: int = 1,
        compresslevel: int = 4,
        max_time: Optional[str] = None,
    ):
        if os.sep in prefix or (os.altsep and os.altsep in prefix):
            raise ValueError("Prefix must not contain path separators: '{}'".format(prefix))
        self.input_path = input_path
        self.output_folder = output_folder
        self.prefix = prefix
        self.threads = max(1, threads)
        self.compresslevel = compresslevel
        self.bin_size, self.units = self._parse_interval(interval)
        self.bin_seconds = self.bin_size * UNIT_SECONDS[self.units]
        if max_time is None:
            self.max_seconds = None
            self.max_label = None
        else:
            max_value, max_units = self._parse_interval(max_time, flag="--max-time")
            self.max_seconds = max_value * UNIT_SECONDS[max_units]
            # Label the truncated last bin the way the user wrote the cutoff.
            self.max_label = format_number(max_value) + max_units

    @staticmethod
    def _parse_interval(interval: str, flag: str = "--interval") -> Tuple[float, str]:
        units = interval[-1].lower()
        if units not in UNIT_SECONDS:
            raise ValueError(
                "Invalid unit '{}' in {} value '{}'. Use one of: {}".format(
                    units, flag, interval, ", ".join(sorted(UNIT_SECONDS))
                )
            )
        try:
            bin_size = float(interval[:-1])
        except ValueError:
            raise ValueError("Invalid {} value: '{}'".format(flag, interval)) from None
        if bin_size <= 0:
            raise ValueError("{} must be greater than zero.".format(flag))
        return bin_size, units

    def run(self) -> List[str]:
        """Execute the full pipeline. Returns the list of output file paths."""
        overall_start = time()
        os.makedirs(self.output_folder, exist_ok=True)

        fastq_list = find_fastq_files(self.input_path)
        if not fastq_list:
            raise FileNotFoundError("No FASTQ files found in {}".format(self.input_path))
        workers = min(self.threads, len(fastq_list))

        # Pass 1: find the time range of the run.
        print("Scanning {} FASTQ file(s)...".format(len(fastq_list)), end="", flush=True)
        start = time()
        scans = self._map(scan_file, fastq_list, workers)
        total_reads = sum(s.reads for s in scans)
        total_missing = sum(s.missing for s in scans)
        total_malformed = sum(s.malformed for s in scans)
        t_mins = [s.t_min for s in scans if s.t_min is not None]
        t_maxs = [s.t_max for s in scans if s.t_max is not None]
        print(" {} reads in {}".format(total_reads, elapsed_time(time() - start)))
        if total_missing:
            print(
                "Warning: {} read(s) without a 'start_time=' or 'st:Z:' header field "
                "were skipped.".format(total_missing),
                file=sys.stderr,
            )
        if total_malformed:
            print(
                "Warning: {} read(s) with an unparseable start time value were skipped "
                "(field present but not a valid timestamp).".format(total_malformed),
                file=sys.stderr,
            )
        if not t_mins:
            raise ValueError(
                "No read start times found. Headers must contain a Guppy/MinKNOW "
                "'start_time=' field or a Dorado 'st:Z:' tag."
            )

        t_min, t_max = min(t_mins), max(t_maxs)
        run_seconds = (t_max - t_min).total_seconds()
        cutoff = None
        # Strictly '<': a cutoff equal to the run length must not silently
        # drop the read(s) acquired at exactly t_max.
        if self.max_seconds is not None and self.max_seconds < run_seconds:
            cutoff = self.max_seconds
            # The small tolerance guards against float imprecision creating a
            # spurious extra bin (e.g. (1.1*3600)/(0.1*3600) -> 11.000000000000002,
            # whose ceil would be 12); the extra bin's name would collide with
            # the real last bin's cutoff label.
            num_bins = max(1, ceil(cutoff / self.bin_seconds - 1e-9))
            print(
                "Run spans {}; binning only the first {} -> {} interval(s) of {}{}".format(
                    elapsed_time(run_seconds), self.max_label, num_bins,
                    format_number(self.bin_size), self.units,
                )
            )
        else:
            num_bins = int(run_seconds // self.bin_seconds) + 1
            print(
                "Run spans {} -> {} intervals of {}{}".format(
                    elapsed_time(run_seconds), num_bins, format_number(self.bin_size), self.units
                )
            )

        # Pass 2: compress each read once into its interval chunk.
        print("Binning reads...", end="", flush=True)
        start = time()
        with tempfile.TemporaryDirectory(prefix="nanotimesort_", dir=self.output_folder) as chunk_dir:
            # With a cutoff, whole files whose earliest read is already past
            # it never need to be decompressed again.
            jobs = []
            skipped_files = 0
            for i, path in enumerate(fastq_list):
                if scans[i].t_min is None:
                    continue  # no binnable reads in this file at any time
                if cutoff is not None and (scans[i].t_min - t_min).total_seconds() >= cutoff:
                    skipped_files += 1
                    continue
                jobs.append(
                    (path, i, t_min, self.bin_seconds, num_bins, cutoff, chunk_dir, self.compresslevel)
                )
            results = self._starmap(chunk_file, jobs, workers)
            reads_per_bin = [0] * num_bins
            bp_per_bin = [0] * num_bins
            for file_reads, file_bp in results:
                for b in range(num_bins):
                    reads_per_bin[b] += file_reads[b]
                    bp_per_bin[b] += file_bp[b]
            skipped_note = (
                " ({} file(s) past the cutoff skipped)".format(skipped_files) if skipped_files else ""
            )
            print(" done in {}{}".format(elapsed_time(time() - start), skipped_note))

            # Pass 3: assemble cumulative outputs by gzip member concatenation.
            print("Writing cumulative interval files...", end="", flush=True)
            start = time()
            last_label = self.max_label if cutoff is not None else None
            outputs = self._assemble(chunk_dir, num_bins, reads_per_bin, bp_per_bin, last_label)
            print(" done in {}".format(elapsed_time(time() - start)))

        print("Total run time: {}".format(elapsed_time(time() - overall_start)))
        return outputs

    def _assemble(
        self,
        chunk_dir: str,
        num_bins: int,
        reads_per_bin: List[int],
        bp_per_bin: List[int],
        last_label: Optional[str] = None,
    ) -> List[str]:
        chunk_names = sorted(os.listdir(chunk_dir))
        outputs: List[str] = []
        seen_names = set()
        previous_path: Optional[str] = None
        cumulative_reads = 0
        cumulative_bp = 0

        for b in range(num_bins):
            cumulative_reads += reads_per_bin[b]
            cumulative_bp += bp_per_bin[b]
            if last_label is not None and b == num_bins - 1:
                # A --max-time cutoff truncates the last bin: name it after
                # the actual cutoff so the file never claims more time than
                # it contains.
                span = last_label
            else:
                span = format_number((b + 1) * self.bin_size) + self.units
            out_name = "{}_0-{}_{}reads_{}bp.fastq.gz".format(
                self.prefix, span, cumulative_reads, cumulative_bp
            )
            out_path = os.path.join(self.output_folder, out_name)
            # Opening a colliding name 'wb' would silently truncate the
            # previous cumulative file before it is read back; fail loudly
            # instead if bin labels ever collide.
            if out_name in seen_names:
                raise RuntimeError(
                    "Internal error: output name collision for '{}'".format(out_name)
                )
            seen_names.add(out_name)
            prefix = "chunk_b{:06d}_".format(b)
            with open(out_path, "wb") as out:
                if previous_path is not None:
                    with open(previous_path, "rb") as prev:
                        shutil.copyfileobj(prev, out, _COPY_BUFFER)
                for name in chunk_names:
                    if name.startswith(prefix):
                        with open(os.path.join(chunk_dir, name), "rb") as chunk:
                            shutil.copyfileobj(chunk, out, _COPY_BUFFER)
            outputs.append(out_path)
            previous_path = out_path

        return outputs

    def _map(self, func, items, workers):
        if workers <= 1 or len(items) <= 1:
            return [func(item) for item in items]
        with ProcessPoolExecutor(max_workers=workers) as pool:
            return list(pool.map(func, items))

    def _starmap(self, func, jobs, workers):
        if workers <= 1 or len(jobs) <= 1:
            return [func(*job) for job in jobs]
        with ProcessPoolExecutor(max_workers=workers) as pool:
            futures = [pool.submit(func, *job) for job in jobs]
            return [f.result() for f in futures]


def elapsed_time(seconds: float) -> str:
    """Format a duration in seconds as a compact '1d2h3m4s' style string."""
    if seconds < 1:
        return "{:.2f}s".format(seconds)
    minutes, secs = divmod(int(round(seconds)), 60)
    hours, minutes = divmod(minutes, 60)
    days, hours = divmod(hours, 24)
    parts = [("d", days), ("h", hours), ("m", minutes), ("s", secs)]
    return "".join("{}{}".format(value, name) for name, value in parts if value) or "0s"
