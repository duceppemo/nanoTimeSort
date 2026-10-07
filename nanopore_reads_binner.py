#!/usr/bin/env python3
"""Deprecated entry point kept for backward compatibility.

Use the 'nanotimesort' command instead (pip install .). This wrapper maps the
old script invocation onto the new CLI.
"""

import sys

from nanotimesort.cli import main

if __name__ == "__main__":
    print(
        "Note: nanopore_reads_binner.py is deprecated; use the 'nanotimesort' "
        "command instead.",
        file=sys.stderr,
    )
    sys.exit(main())
