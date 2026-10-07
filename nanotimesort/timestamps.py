"""Extraction and parsing of read start times from Nanopore FASTQ headers.

Two header dialects are supported:

* Guppy / MinKNOW (key=value fields)::

    @<read_id> runid=... read=... ch=... start_time=2019-07-16T19:51:22Z
    @<read_id> ... start_time=2025-01-13T10:45:28.681306+00:00 ...

* Dorado (SAM-style tags)::

    @<read_id> qs:f:21.3 du:f:12.44 ch:i:942 st:Z:2023-09-01T11:13:45.731+00:00 ...
"""

from __future__ import annotations

import re
from datetime import datetime, timezone
from typing import Optional, Union

# Guppy/MinKNOW style field.
_GUPPY_PREFIX = b"start_time="
# Dorado SAM-tag style field (st:Z:<ISO 8601>).
_DORADO_PREFIX = b"st:Z:"

_FRACTION_RE = re.compile(r"\.(\d+)")


def extract_start_time(header: bytes) -> Optional[datetime]:
    """Return the read start time from a FASTQ header line, or None if absent.

    :param header: raw FASTQ header line (bytes, with or without trailing newline)
    :return: timezone-aware datetime (UTC assumed when the timestamp is naive)
    """
    for item in header.split():
        if item.startswith(_GUPPY_PREFIX):
            return parse_timestamp(item[len(_GUPPY_PREFIX):])
        if item.startswith(_DORADO_PREFIX):
            return parse_timestamp(item[len(_DORADO_PREFIX):])
    return None


def parse_timestamp(raw: Union[bytes, str]) -> datetime:
    """Parse an ISO 8601 / RFC 3339 timestamp into a timezone-aware datetime.

    Handles 'Z' suffixes and unusual fractional-second precision on Python
    versions where datetime.fromisoformat() is strict (< 3.11).
    """
    if isinstance(raw, bytes):
        raw = raw.decode("ascii")
    try:
        dt = datetime.fromisoformat(raw)
    except ValueError:
        dt = datetime.fromisoformat(_normalize(raw))
    if dt.tzinfo is None:
        dt = dt.replace(tzinfo=timezone.utc)
    return dt


def _normalize(s: str) -> str:
    """Rewrite a timestamp so strict fromisoformat() implementations accept it."""
    if s.endswith(("Z", "z")):
        s = s[:-1] + "+00:00"

    def _pad(match: "re.Match[str]") -> str:
        digits = match.group(1)[:6]
        return "." + digits.ljust(6, "0")

    return _FRACTION_RE.sub(_pad, s, count=1)
