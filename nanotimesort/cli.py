"""Command-line interface for nanoTimeSort."""

from __future__ import annotations

import sys
from argparse import ArgumentParser, RawTextHelpFormatter
from multiprocessing import cpu_count

from . import __version__
from .binner import NanoTimeSort


def build_parser() -> ArgumentParser:
    cpu = cpu_count()
    parser = ArgumentParser(
        prog="nanotimesort",
        description="Bin Oxford Nanopore reads by cumulative sequencing time intervals.\n"
                    "Supports Guppy/MinKNOW ('start_time=') and Dorado ('st:Z:') FASTQ headers.",
        formatter_class=RawTextHelpFormatter,
    )
    parser.add_argument(
        "-f", "--fastq", metavar="/basecalled/folder/", required=True,
        help="Input folder (searched recursively) or single FASTQ file.\n"
             "Accepts .fastq, .fq, .fastq.gz and .fq.gz.",
    )
    parser.add_argument(
        "-o", "--output", metavar="/output/folder/", required=True,
        help="Output folder. Created if it does not exist.",
    )
    parser.add_argument(
        "-i", "--interval", metavar="1h", required=True,
        help="Time interval for the bins, e.g. '1h', '30m' or '90s'.\n"
             "Bins are cumulative: with '-i 1h', the second file also\n"
             "contains the reads of the first hour.",
    )
    parser.add_argument(
        "-m", "--max-time", metavar="2h", default=None,
        help="Only bin reads acquired up to this elapsed time; later reads\n"
             "are discarded. E.g. '-i 2h -m 2h' produces a single file with\n"
             "the first two hours of the run. Same format as --interval.\n"
             "Default: bin the whole run.",
    )
    parser.add_argument(
        "-p", "--prefix", metavar="my_sample", default="interval",
        help="Output file prefix. Files are named like\n"
             "'my_sample_0-1h_123reads_456789bp.fastq.gz'.\n"
             "Default: interval",
    )
    parser.add_argument(
        "-t", "--threads", metavar=str(cpu), type=int, default=cpu,
        help="Number of FASTQ files to process in parallel.\n"
             "Default: {}".format(cpu),
    )
    parser.add_argument(
        "-c", "--compression-level", metavar="4", type=int, default=4,
        choices=range(1, 10),
        help="Gzip compression level for output files (1=fastest, 9=smallest).\n"
             "Default: 4",
    )
    parser.add_argument(
        "-v", "--version", action="version", version="nanoTimeSort v{}".format(__version__),
    )
    return parser


def main(argv=None) -> int:
    args = build_parser().parse_args(argv)
    try:
        binner = NanoTimeSort(
            input_path=args.fastq,
            output_folder=args.output,
            interval=args.interval,
            prefix=args.prefix,
            threads=args.threads,
            compresslevel=args.compression_level,
            max_time=args.max_time,
        )
        binner.run()
    except (ValueError, FileNotFoundError) as err:
        print("Error: {}".format(err), file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
