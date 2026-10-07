from datetime import datetime, timezone

from nanotimesort.timestamps import extract_start_time, parse_timestamp

GUPPY_HEADER = (
    b"@c041234f-1234-4a81-9a5a-1234567890ab runid=abc123 read=42 ch=133 "
    b"start_time=2019-07-16T19:51:22Z flow_cell_id=FAK12345"
)
MINKNOW_HEADER = (
    b"@bd8655fb-383c-45cc-bff3-eb1dc86533e0 runid=abc parent_read_id=bd8655fb "
    b"start_time=2025-01-13T10:45:28.681306+00:00 protocol_group_id=test"
)
DORADO_HEADER = (
    b"@0000813e-1111-4c35-8f57-222233334444 qs:f:21.5 du:f:12.44 ns:i:62205 "
    b"ts:i:10 mx:i:3 ch:i:942 st:Z:2023-09-01T11:13:45.731+00:00 rn:i:9566 "
    b"fn:Z:PAO12345_pass_0.pod5 sm:f:421.3 sd:f:93.0 sv:Z:quantile dx:i:0 "
    b"RG:Z:abc_dna_r10.4.1_e8.2_400bps_hac@v4.2.0"
)


def test_guppy_header():
    dt = extract_start_time(GUPPY_HEADER)
    assert dt == datetime(2019, 7, 16, 19, 51, 22, tzinfo=timezone.utc)


def test_minknow_header_with_microseconds():
    dt = extract_start_time(MINKNOW_HEADER)
    assert dt == datetime(2025, 1, 13, 10, 45, 28, 681306, tzinfo=timezone.utc)


def test_dorado_header():
    dt = extract_start_time(DORADO_HEADER)
    assert dt == datetime(2023, 9, 1, 11, 13, 45, 731000, tzinfo=timezone.utc)


def test_header_without_time_returns_none():
    assert extract_start_time(b"@read1 length=100") is None


def test_naive_timestamp_assumed_utc():
    dt = parse_timestamp("2023-09-01T11:13:45")
    assert dt.tzinfo is not None
    assert dt.utcoffset().total_seconds() == 0


def test_z_suffix_with_milliseconds():
    dt = parse_timestamp("2023-09-01T11:13:45.7Z")
    assert dt == datetime(2023, 9, 1, 11, 13, 45, 700000, tzinfo=timezone.utc)


def test_normalize_z_suffix():
    from nanotimesort.timestamps import _normalize

    assert _normalize("2019-07-16T19:51:22Z") == "2019-07-16T19:51:22+00:00"
    assert _normalize("2019-07-16T19:51:22z") == "2019-07-16T19:51:22+00:00"


def test_normalize_fractional_digits():
    from nanotimesort.timestamps import _normalize

    # Padded to 6 digits for strict parsers (Python < 3.11)...
    assert _normalize("2023-09-01T11:13:45.731+00:00") == "2023-09-01T11:13:45.731000+00:00"
    # ...and excess digits truncated to 6.
    assert _normalize("2023-09-01T11:13:45.1234567890Z") == "2023-09-01T11:13:45.123456+00:00"


def test_non_utc_offset_converted_consistently():
    # 12:00+02:00 is 10:00 UTC, i.e. EARLIER than 11:00Z.
    a = parse_timestamp("2024-05-01T12:00:00+02:00")
    b = parse_timestamp("2024-05-01T11:00:00Z")
    assert a < b
    assert (b - a).total_seconds() == 3600


def test_first_time_field_wins():
    # Both dialects present: the first field encountered is used.
    header = b"@r1 start_time=2024-05-01T10:00:00Z st:Z:2024-05-01T12:00:00Z"
    dt = extract_start_time(header)
    assert dt == datetime(2024, 5, 1, 10, 0, 0, tzinfo=timezone.utc)


def test_similar_field_names_do_not_match():
    # Only an exact 'start_time=' prefix counts, not e.g. parent fields.
    header = b"@r1 parent_start_time=2024-05-01T10:00:00Z other=1"
    assert extract_start_time(header) is None


def test_malformed_timestamp_returns_none():
    assert extract_start_time(b"@r1 start_time=notadate ch=1") is None
    assert extract_start_time(b"@r1 st:Z:2024-99-99T99:99:99Z") is None
    assert extract_start_time(b"@r1 start_time=") is None


def test_malformed_timestamp_raises_in_direct_parse():
    import pytest

    with pytest.raises(ValueError):
        parse_timestamp("notadate")
