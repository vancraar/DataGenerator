#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Convert MultiAPEX DataGenerator CSV output into PEPC gravity binary input.

PEPC's gravity frontend expects a raw binary stream with one particle after
another. Each particle stores 3D position, 3D velocity, and mass as double
precision values in that order.
"""

from __future__ import annotations

import argparse
import csv
import struct
from pathlib import Path

REQUIRED_COLUMNS = ("mass", "pos_x", "pos_y", "pos_z", "vel_x", "vel_y", "vel_z")
PACK_FORMAT = "<7d"
RECORD_SIZE = struct.calcsize(PACK_FORMAT)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        prog="convert_to_pepc_binary",
        description="Convert DataGenerator CSV output to PEPC gravity binary input.",
    )
    parser.add_argument("--input", "-i", required=True, help="Input CSV from DataGenerator")
    parser.add_argument("--output", "-o", required=True, help="Output PEPC binary file")
    return parser.parse_args()


def convert_csv_to_pepc_binary(input_csv: str | Path, output_bin: str | Path) -> int:
    """Convert a DataGenerator CSV file into PEPC's binary particle stream.

    Returns the number of particles written.
    """
    input_csv = Path(input_csv)
    output_bin = Path(output_bin)

    with input_csv.open("r", newline="") as csv_file, output_bin.open("wb") as bin_file:
        reader = csv.DictReader(csv_file)
        if reader.fieldnames is None:
            raise ValueError("CSV file is missing a header row")

        missing_columns = [column for column in REQUIRED_COLUMNS if column not in reader.fieldnames]
        if missing_columns:
            raise ValueError(f"CSV file is missing required columns: {', '.join(missing_columns)}")

        count = 0
        for row in reader:
            record = struct.pack(
                PACK_FORMAT,
                float(row["pos_x"]),
                float(row["pos_y"]),
                float(row["pos_z"]),
                float(row["vel_x"]),
                float(row["vel_y"]),
                float(row["vel_z"]),
                float(row["mass"]),
            )
            if len(record) != RECORD_SIZE:
                raise AssertionError("internal packing error")
            bin_file.write(record)
            count += 1

    return count


def main() -> int:
    args = parse_args()
    count = convert_csv_to_pepc_binary(args.input, args.output)
    print(f"Wrote {count} particles to {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
