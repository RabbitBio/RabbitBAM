#!/usr/bin/env python3
"""Compare complete BAM record bytes and order against a validated output."""
import argparse
import contextlib
import gzip
import itertools
import struct


def exact(stream, size):
    data = stream.read(size)
    if len(data) != size:
        raise ValueError(f"Truncated BAM: wanted {size} bytes, got {len(data)}")
    return data


def integer(stream):
    return struct.unpack("<i", exact(stream, 4))[0]


def header(stream):
    if exact(stream, 4) != b"BAM\x01":
        raise ValueError("Invalid BAM magic")
    text_size = integer(stream)
    if text_size < 0:
        raise ValueError("Negative header size")
    exact(stream, text_size)
    references = []
    count = integer(stream)
    if count < 0:
        raise ValueError("Negative reference count")
    for _ in range(count):
        name_size = integer(stream)
        if name_size <= 0:
            raise ValueError("Invalid reference name size")
        name = exact(stream, name_size)
        length = integer(stream)
        if not name.endswith(b"\0") or length < 0:
            raise ValueError("Invalid reference entry")
        references.append((name, length))
    return references


def records(stream):
    while True:
        length_data = stream.read(4)
        if not length_data:
            return
        if len(length_data) != 4:
            raise ValueError("Truncated BAM record size")
        length = struct.unpack("<i", length_data)[0]
        if length < 32:
            raise ValueError(f"Invalid record size {length}")
        yield exact(stream, length)


def verify(paths, expected_records, expected_duplicates):
    with contextlib.ExitStack() as stack:
        streams = [stack.enter_context(gzip.open(path, "rb")) for path in paths]
        references = [header(stream) for stream in streams]
        if any(value != references[0] for value in references[1:]):
            raise ValueError("Reference dictionaries differ")
        count = duplicates = 0
        previous = (-1, -1)
        for count, group in enumerate(itertools.zip_longest(*(records(stream) for stream in streams)), 1):
            if any(record is None for record in group):
                raise ValueError(f"Record counts differ at record {count}")
            if any(record != group[0] for record in group[1:]):
                raise ValueError(f"Record bytes or order differ at record {count}")
            tid, pos = struct.unpack_from("<ii", group[0])
            flag = struct.unpack_from("<H", group[0], 14)[0]
            key = (tid if tid >= 0 else 1 << 31, pos)
            if key < previous:
                raise ValueError(f"Coordinate order decreases at record {count}")
            previous = key
            duplicates += bool(flag & 0x400)
        if count != expected_records or duplicates != expected_duplicates:
            raise ValueError(f"Unexpected totals: records={count}, duplicates={duplicates}")
    print(f"PASS files={len(paths)} records={count} duplicates={duplicates}; "
          "all BAM record bytes and order match; coordinate sorted; gzip CRC checked")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("paths", nargs="+", help="Validated baseline followed by new outputs")
    parser.add_argument("--records", type=int, default=2668351)
    parser.add_argument("--duplicates", type=int, default=15695)
    args = parser.parse_args()
    if len(args.paths) < 2:
        parser.error("Provide a baseline and at least one new output")
    verify(args.paths, args.records, args.duplicates)


if __name__ == "__main__":
    main()
