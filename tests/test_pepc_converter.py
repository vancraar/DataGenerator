import csv
import struct
import tempfile
import unittest
from pathlib import Path

from convert_to_pepc_binary import convert_csv_to_pepc_binary


class PepcConverterTests(unittest.TestCase):
    def test_converts_csv_rows_to_pepc_binary_stream(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            csv_path = tmpdir / "input.csv"
            bin_path = tmpdir / "output.bin"

            with csv_path.open("w", newline="") as f:
                writer = csv.writer(f)
                writer.writerow(["id", "mass", "pos_x", "pos_y", "pos_z", "vel_x", "vel_y", "vel_z"])
                writer.writerow([0, 1.5, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0])
                writer.writerow([1, 8.5, -1.0, -2.0, -3.0, -4.0, -5.0, -6.0])

            count = convert_csv_to_pepc_binary(csv_path, bin_path)

            self.assertEqual(count, 2)
            self.assertEqual(bin_path.stat().st_size, 2 * 56)

            payload = bin_path.read_bytes()
            self.assertEqual(
                payload,
                struct.pack(
                    "<7d",
                    2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 1.5,
                ) + struct.pack(
                    "<7d",
                    -1.0, -2.0, -3.0, -4.0, -5.0, -6.0, 8.5,
                ),
            )


if __name__ == "__main__":
    unittest.main()
